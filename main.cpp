#include "engine/easy.h"
#include <iostream>
#include <map>
#include <ranges>
#include <set>
#include <variant>

#include "metrics.h"
#include "tuples.h"

using namespace arctic;  //NOLINT


/// Simulation time step
double g_dt = 0.0001;
/// Current simulation time
double g_t = 0.0;
/// Time multiplier for simulation speed control
double g_t_mult = 1.0 / 128.0;

std::mt19937 g_random(548);

TInstant Now() {
    return g_t;
}

// GUI elements
GuiFactory gf;
std::shared_ptr<GuiTheme> g_theme;
std::shared_ptr<Panel> g_gui;
std::shared_ptr<Scrollbar> g_tablet_metrics_scrollbar;
std::shared_ptr<Scrollbar> g_usage_scrollbar;
std::shared_ptr<Checkbox> g_no_moves_checkbox;
std::shared_ptr<Checkbox> g_guess_usage_checkbox;
std::shared_ptr<Checkbox> g_fill_nodes_checkbox;
std::shared_ptr<Checkbox> g_stick_checkbox;
std::shared_ptr<Checkbox> g_priority_checkbox;
Si32 g_max_tablet_metrics = 100'000'000; 


// Global simulation state

// Visual elements
Font g_font;

/// @brief Increases simulation speed
void SpeedUp() {
    g_t_mult = std::min(32.0, g_t_mult * 2.0);
}

/// @brief Decreases simulation speed
void SlowDown() {
    g_t_mult = std::max(1.0 / 128.0, g_t_mult / 2.0);
}

using TNodeId = Ui32;
using TResourceRawValues = std::tuple<Si64, Si64, Si64, Si64>;
using TResourceNormalizedValues = std::tuple<double, double, double, double>;
using TFullTabletId = Ui32;
using TDataCenterId = Ui32;
using TTabletCategoryId = Ui32;

const TResourceRawValues MAXIMUM_VALUES = {1'000'000'000, 64'000'000'000ull, 1'000'000'000, 1'000'000};

enum EResource : size_t {
    CPU,
    Memory,
    Network,
    Counter,
};

enum class EResourceToBalance {
    ComputeResources,
    Counter,
    CPU,
    Memory,
    Network,
};

enum ETabletType {
    Dummy,
    DataShard,
    ColumnShard,
    Coordinator,
};

Rgba GetColor(ETabletType type) {
    switch (type) {
        case ETabletType::Dummy: return Rgba(255, 200, 100);
        case ETabletType::DataShard: return Rgba(0, 64, 255);
        case ETabletType::ColumnShard: return Rgba(200, 200, 200);
        case ETabletType::Coordinator: return Rgba(200, 200, 0);
    }
}

TResourceNormalizedValues NormalizeRawValues(const TResourceRawValues& values, const TResourceRawValues& maximum = MAXIMUM_VALUES) {
    return safe_div(values, maximum);
}

template <typename ResourceTypes>
static auto ExtractResourceUsage(const ResourceTypes& values, EResourceToBalance resource = EResourceToBalance::ComputeResources) {
    switch (resource) {
            case EResourceToBalance::CPU: return std::get<EResource::CPU>(values);
    case EResourceToBalance::Memory: return std::get<EResource::Memory>(values);
    case EResourceToBalance::Network: return std::get<EResource::Network>(values);
    case EResourceToBalance::Counter: return std::get<EResource::Counter>(values);
    case EResourceToBalance::ComputeResources: return max(values);
    }
}

template<std::ranges::range TRange>
std::ranges::range_value_t<TRange> StableSum(const TRange& values) {
    using TValue = std::ranges::range_value_t<TRange>;
    TValue sum{};
    TValue correction{};
    for (const auto& x : values) {
        TValue y = x - correction;
        TValue tmp = sum + y;
        correction = (tmp - sum) - y;
        sum = tmp;
    }
    return sum;
}

template <typename... ResourceTypes>
inline std::tuple<ResourceTypes...> GetStDev(std::vector<std::tuple<ResourceTypes...>> values) {
    std::tuple<ResourceTypes...> sum;
    if (values.empty())
        return sum;
    sum = StableSum(values);
    auto mean = sum / values.size();
    auto quadraticDev = [&] (const std::tuple<ResourceTypes...>& value) {
        auto diff = value - mean;
        return diff * diff;
    };
    for (auto& value : values) {
        value = quadraticDev(value);
    }        
    sum = StableSum(values);
    auto div = sum / values.size();
    auto st_dev = sqrt(div);
    return tuple_cast<ResourceTypes...>::cast(st_dev);
}

TResourceNormalizedValues GetMinNodeUsageToBalance() {
    return {0.3, 0.3, 0.3, 0};
}

Ui64 GetMaxResourceCounter() {
    return 1'000'000;
}

TResourceNormalizedValues GetMinScatterToBalance() {
    return {0.5, 0.5, 0.5, 0.02};
}

double GetMaxNodeUsageToKick() {
    return 0.9;
}

double GetNodeUsageRangeToKick() {
    return 0.2;
}

double GetResourceOvercommitment() {
    return 3;
}

Ui64 GetMaxMovementsOnAutoBalancer() {
    return 1;
}

Ui64 GetMaxMovementsOnEmergencyBalancer() {
    return 2;
}

TDuration GetTabletKickCooldownPeriod() {
    return 60;
}

TDuration GetMinPeriodBetweenBalance() {
    return 0.2;
}

TInstant g_last_balancer_trigger;
TFullTabletId g_last_tablet_moved;
TNodeId g_last_node_moved;

void UpdateResourceValues(TNodeId nodeId, const TResourceRawValues& delta);

struct TTabletInfo;

struct TTabletCategoryInfo {
    TTabletCategoryId Id;
    std::unordered_set<TFullTabletId> Tablets;
    Ui64 MaxDisconnectTimeout = 0;
    bool StickTogetherInDC = false;

    TTabletCategoryInfo(TTabletCategoryId id, bool stickTogether)
        : Id(id)
        , StickTogetherInDC(stickTogether)
    {}
};

TTabletCategoryInfo g_system_category(0, true);

struct TTabletInfo {
    TResourceRawValues ResourceValues;
    TTabletMetricsAggregates ResourceMetricsAggregates;
    TFullTabletId Id;
    ETabletType Type;
    TNodeId NodeId;
    TInstant LastBalancerDecisionTime = 0;
    double Weight = 0;
    TTabletCategoryInfo* Category = nullptr; 
    bool Alive = true;
    TNodeId FailedNodeId = 0;

    bool IsGoodForBalancer(TInstant now) {
        return (now - LastBalancerDecisionTime > GetTabletKickCooldownPeriod()) || LastBalancerDecisionTime == 0;
    }

    bool HasAllowedMetric(EResourceToBalance resource) const {
        if (Type == ETabletType::ColumnShard) {
            return resource == EResourceToBalance::Counter;
        }
        return true;
    }
    
    bool HasMetric(EResourceToBalance resource) const {
        if (!HasAllowedMetric(resource)) {
            return false;
        }
        return ExtractResourceUsage(ResourceValues, resource) > 0;
    }

    TFullTabletId GetFullTabletId() const {
        return Id;
    }

    double GetWeight(EResourceToBalance resourceToBalance = EResourceToBalance::ComputeResources) const {
        auto result = ExtractResourceUsage(NormalizeRawValues(ResourceValues), resourceToBalance);
        //*Log() << "Tablet " << Id << " has weight " << result << " in resource " << (int)resourceToBalance << std::endl;
        return std::max(result, Weight);
    }

    bool IsAliveOnNode(TNodeId nodeId) const {
        return NodeId == nodeId;
    }

    template <typename ResourcesType>
    static double GetUsage(const ResourcesType& current, const ResourcesType& maximum, EResourceToBalance resource = EResourceToBalance::ComputeResources) {
        auto normValues = NormalizeRawValues(current, maximum);
        return ExtractResourceUsage(normValues, resource);
    }

    void MakeBalancerDecision(TInstant now) {
        LastBalancerDecisionTime = now;
    }

    Si64 GetCounterValue() {
        if (std::get<EResource::CPU>(ResourceValues) > 0) {
            return 0;
        }
        if (std::get<EResource::Memory>(ResourceValues) > 0) {
            return 0;
        }
        if (std::get<EResource::Network>(ResourceValues) > 0) {
            return 0;
        }
        return 1;
    }

    void UpdateResourceUsage(const TMetrics& metrics) {
        TInstant now = Now();
        auto before = ResourceValues;
        if (metrics.HasCPU()) {
            ResourceMetricsAggregates.MaximumCPU.SetValue(metrics.GetCPU(), now);
        } else {
            ResourceMetricsAggregates.MaximumCPU.AdvanceTime(now);
        }
        std::get<EResource::CPU>(ResourceValues) = ResourceMetricsAggregates.MaximumCPU.GetValue();
        if (metrics.HasMemory()) {
                ResourceMetricsAggregates.MaximumMemory.SetValue(metrics.GetMemory(), now);
        } else {
            ResourceMetricsAggregates.MaximumMemory.AdvanceTime(now);
        }
        std::get<EResource::Memory>(ResourceValues) = ResourceMetricsAggregates.MaximumMemory.GetValue();
        if (metrics.HasNetwork()) {
            ResourceMetricsAggregates.MaximumNetwork.SetValue(metrics.GetNetwork(), now);
        } else {
            ResourceMetricsAggregates.MaximumNetwork.AdvanceTime(now);
        }
        std::get<EResource::Network>(ResourceValues) = ResourceMetricsAggregates.MaximumNetwork.GetValue();
        std::get<EResource::Counter>(ResourceValues) = GetCounterValue();
        const auto& after = ResourceValues;
        UpdateResourceValues(NodeId, after - before);
    }

    template <typename ResourcesType>
    void FilterRawValues(ResourcesType& values) const {
        if (std::get<EResource::Counter>(ResourceValues) == 0) {
            std::get<EResource::Counter>(values) = 0;
        }
        if (std::get<EResource::CPU>(ResourceValues) == 0) {
            std::get<EResource::CPU>(values) = 0;
        }
        if (std::get<EResource::Memory>(ResourceValues) == 0) {
            std::get<EResource::Memory>(values) = 0;
        }
        if (std::get<EResource::Network>(ResourceValues) == 0) {
            std::get<EResource::Network>(values) = 0;
        }
    }

    bool IsAlive() const {
        return Alive && NodeId != 0;
    }

    double GetPriority() const {
        return (Type == ETabletType::Coordinator) + GetWeight();
    }

    TTabletInfo(TFullTabletId id, const TResourceRawValues& metrics, ETabletType type = ETabletType::Dummy)
        : Id(id)
        , ResourceValues(metrics)
        , Type(type)
    {
        if (Type == ETabletType::Coordinator) {
            Category = &g_system_category;
            Category->Tablets.insert(Id);
        }
    }

    Rgba GetColor() const {
        Ui32 scale = 255 * GetWeight() * 7;
        if (Type == ETabletType::ColumnShard) {
            scale = 255;
        }
        return Scale(::GetColor(Type), std::min(scale, 255u));
    }
};

struct TNodeInfo {
    enum ENodeKind : Ui32 {
        ModernCPU,
        OldCPU,
    };
    TResourceRawValues RawValues;
    NMetrics::TFastRiseAverageValue<double, 20> AveragedNodeTotalUsage;
    std::set<TTabletInfo*> Tablets;
    TNodeId Id;
    Ui32 NodeKind = ModernCPU;
    mutable bool Drawn = false;
    bool Alive = true;
    TDataCenterId DataCenterId;
    TInstant StartTime = Now();
    std::vector<TInstant> RestartTimestamps; 

    bool IsAlive() const {
        return Alive;
    }

    double GetNodeUsage(const TResourceNormalizedValues& normValues, EResourceToBalance resource = EResourceToBalance::ComputeResources) const {
        double usage = ExtractResourceUsage(normValues, resource);
        if (resource == EResourceToBalance::ComputeResources && AveragedNodeTotalUsage.IsValueStable()) {
            usage = std::max(usage, AveragedNodeTotalUsage.GetValue());
        }
        return usage;
    }

    double GetNodeUsage(EResourceToBalance resource = EResourceToBalance::ComputeResources) const {
        auto normValues = NormalizeRawValues(RawValues, GetResourceMaximumValues());
        return GetNodeUsage(normValues, resource);
    }

    TResourceRawValues GetResourceMaximumValues() const {
        auto maxValues = MAXIMUM_VALUES;
        if (NodeKind == OldCPU) {
            std::get<EResource::CPU>(maxValues) /= 2;
        }
        return maxValues;
    }

    bool IsOverloaded() const {
        auto maxValues = GetResourceMaximumValues();
        auto normValues = NormalizeRawValues(RawValues, maxValues);
        return GetNodeUsage(normValues) >= GetMaxNodeUsageToKick();
    }

    double GetNodeUsageForTablet(const TTabletInfo& tablet) const {
        // what it would like when tablet will run on this node?
        TResourceRawValues nodeValues = RawValues;
        TResourceRawValues tabletValues = tablet.ResourceValues;
        tablet.FilterRawValues(nodeValues);
        tablet.FilterRawValues(tabletValues);
        auto current = tablet.IsAliveOnNode(Id) ? nodeValues : nodeValues + tabletValues;
        auto maximum = GetResourceMaximumValues();
        // basically, this is: return max(a / b);
        double usage = TTabletInfo::GetUsage(current, maximum);
        usage = std::max(usage, GetNodeUsage() + (tablet.IsAliveOnNode(Id) ? 0 : tablet.Weight));
        return usage;
    }

    void UpdateResourceTotalUsage(double usage) {
        if (NodeKind == OldCPU) {
            usage *= 1.2;
        }
        usage *= g_usage_scrollbar->GetValue() / 100.0;
        AveragedNodeTotalUsage.Push(usage);
        if (g_guess_usage_checkbox->IsChecked() && Tablets.size() == 1) {
            (*Tablets.begin())->Weight = usage;
        }
    }

    void AddTablet(TTabletInfo& tablet) {
        Tablets.insert(&tablet);
        RawValues += tablet.ResourceValues;
        tablet.NodeId = Id;
    }

    void RemoveTablet(TTabletInfo& tablet) {
        Tablets.erase(&tablet);
        RawValues -= tablet.ResourceValues;
        tablet.NodeId = 0;
    }

    TDataCenterId GetDataCenter() const {
        return DataCenterId;
    }

    Si32 GetPriorityForTablet(const TTabletInfo& tablet, std::unordered_map<TDataCenterId, Si32>& dcPriority) const {
        Si32 priority = 0;

        if (tablet.FailedNodeId == Id) {
            --priority;
        }
        
        priority += dcPriority[GetDataCenter()];
        if (g_priority_checkbox->IsChecked()) {
            priority -= GetRestarts() / 3;
        }
        return priority;
    }

    Ui64 GetRestarts() const {
        return RestartTimestamps.size();
    }

    TNodeInfo(TNodeId id, TDataCenterId dc = 0) : Id(id), DataCenterId(dc) {}
};


enum EBalancerType {
    Emergency,
    ScatterCPU,
    ScatterMemory,
    ScatterNetwork,
    ScatterCounter,
    Scatter,
};

struct TBalancerSettings {
    EBalancerType Type;
    int MaxMovements = 0;
    const std::vector<TNodeId> FilterNodeIds = {};
    EResourceToBalance ResourceToBalance = EResourceToBalance::ComputeResources;
};

void StartHiveBalancer(TBalancerSettings&& settings);

void BalanceNodes(std::vector<TNodeInfo*>& nodes, EResourceToBalance resourceToBalance) {
    std::sort(nodes.begin(), nodes.end(), [resourceToBalance](const TNodeInfo* a, const TNodeInfo* b) -> bool {
        return a->GetNodeUsage(resourceToBalance) > b->GetNodeUsage(resourceToBalance);
    });
}

void BalanceTablets(std::vector<TTabletInfo*>::iterator first, std::vector<TTabletInfo*>::iterator last, EResourceToBalance resourceToBalance) {
    std::uniform_real_distribution randGen;
    std::vector<std::pair<double, TTabletInfo*>> weights;
    weights.reserve(last - first);
    for (auto it = first; it != last; ++it) {
        double weight = (*it)->GetWeight(resourceToBalance);
        weights.emplace_back(weight * randGen(g_random), *it);
    }
    std::sort(weights.begin(), weights.end(), [](const auto& a, const auto& b) -> bool {
        return a.first > b.first;
    });
    for (size_t n = 0; n < weights.size(); ++n) {
        first[n] = weights[n].second;
    }
}

struct THive {
    struct THiveStats {
        struct TNodeStat {
            TNodeId NodeId;
            double Usage;
            TResourceNormalizedValues ResourceNormValues;

            TNodeStat(TNodeId node, double usage, TResourceNormalizedValues values)
                : NodeId(node)
                , Usage(usage)
                , ResourceNormValues(values)
            {
            }
        };

        double MinUsage;
        TNodeId MinUsageNodeId;
        double MaxUsage;
        TNodeId MaxUsageNodeId;
        double Scatter;
        TResourceNormalizedValues ScatterByResource;
        std::vector<TNodeStat> Values;
    };

    struct TNoNodeFound {};
    using TBestNodeResult = std::variant<TNodeInfo*, TNoNodeFound>;

    struct TSelectedNode {
        double Usage;
        TNodeInfo* Node;

        TSelectedNode(double usage, TNodeInfo* node)
            : Usage(usage)
            , Node(node)
        {}

        bool operator <(const TSelectedNode& b) const {
            return Usage < b.Usage;
        }
    };

    std::map<TNodeId, TNodeInfo> Nodes;
    std::map<TFullTabletId, TTabletInfo> Tablets;
    EBalancerType LastTrigger;
    Ui64 LastMovements;

    TNodeInfo* FindNode(TNodeId nodeId) {
        auto it = Nodes.find(nodeId);
        if (it == Nodes.end()) {
            return nullptr;
        }
        return &it->second;
    }

    TTabletInfo* FindTablet(TFullTabletId nodeId) {
        auto it = Tablets.find(nodeId);
        if (it == Tablets.end()) {
            return nullptr;
        }
        return &it->second;
    }

    THiveStats GetStats() const {
        THiveStats stats = {};
        stats.Values.reserve(Nodes.size());
        for (const auto& ni : Nodes) {
            if (ni.second.IsAlive()) {
                auto nodeValues = NormalizeRawValues(ni.second.RawValues, ni.second.GetResourceMaximumValues());
                stats.Values.emplace_back(ni.first, ni.second.GetNodeUsage(nodeValues), nodeValues);
            }
        }
        if (stats.Values.empty()) {
            return stats;
        }
        auto it = std::minmax_element(stats.Values.begin(), stats.Values.end(), [](const THiveStats::TNodeStat& a, const THiveStats::TNodeStat& b) -> bool {
            return a.Usage < b.Usage;
        });
        stats.MaxUsage = it.second->Usage;
        stats.MaxUsageNodeId = it.second->NodeId;
        stats.MinUsage = it.first->Usage;
        stats.MinUsageNodeId = it.first->NodeId;

        TResourceNormalizedValues minValues = stats.Values.front().ResourceNormValues;
        TResourceNormalizedValues maxValues = stats.Values.front().ResourceNormValues;
        for (size_t i = 1; i < stats.Values.size(); ++i) {
            minValues = piecewise_min(minValues, stats.Values[i].ResourceNormValues);
            maxValues = piecewise_max(maxValues, stats.Values[i].ResourceNormValues);
        }

        auto minValuesToBalance = GetMinNodeUsageToBalance();
        maxValues = piecewise_max(maxValues, minValuesToBalance);
        minValues = piecewise_max(minValues, minValuesToBalance);
        auto discrepancy = maxValues - minValues;
        auto& counterDiscrepancy = std::get<EResource::Counter>(discrepancy);
        if (counterDiscrepancy * GetMaxResourceCounter() <= 1.5) {
            // We should ignore counter discrepancy of one - it cannot be fixed by balancer
            // Value 1.5 is used to avoid rounding errors
            counterDiscrepancy = 0;
        }
        stats.ScatterByResource = safe_div(discrepancy, maxValues);
        stats.Scatter = max(stats.ScatterByResource);

        return stats;
    }

    double GetScatter() const {
        THiveStats stats = GetStats();
        return stats.Scatter;
    }

    double GetUsage() const {
        THiveStats stats = GetStats();
        return stats.MaxUsage;
    }

    std::optional<EResourceToBalance> CheckScatter(const TResourceNormalizedValues& scatterByResource) const {
        auto minScatterToBalance = GetMinScatterToBalance();
        auto cmp = piecewise_compare(scatterByResource, minScatterToBalance);
        if (std::get<EResource::Counter>(cmp) == std::partial_ordering::greater) {
            return EResourceToBalance::Counter;
        }
        if (std::get<EResource::CPU>(cmp) == std::partial_ordering::greater) {
            return EResourceToBalance::CPU;
        }
        if (std::get<EResource::Memory>(cmp) == std::partial_ordering::greater) {
            return EResourceToBalance::Memory;
        }
        if (std::get<EResource::Network>(cmp) == std::partial_ordering::greater) {
            return EResourceToBalance::Network;
        }
        return std::nullopt;
    }

    void ProcessTabletBalancer() {
        TInstant now = Now();
        if (now < g_last_balancer_trigger + GetMinPeriodBetweenBalance()) {
            return;
        }
        THiveStats stats = GetStats();
        // *Log() << g_t << " ProcessTabletBalancer: max usage = " << stats.MaxUsage << ", min usage = " << stats.MinUsage << ", scatter = ";
        Out(*Log(),  stats.ScatterByResource);
        *Log() << std::endl;

        double minUsageToKick = GetMaxNodeUsageToKick() - GetNodeUsageRangeToKick();
        auto scatteredResource = CheckScatter(stats.ScatterByResource);
        EBalancerType balancerType = EBalancerType::Emergency;
        if (scatteredResource) {
            switch (*scatteredResource) {
                case EResourceToBalance::Counter:
                    balancerType = EBalancerType::ScatterCounter;
                    break;
                case EResourceToBalance::CPU:
                    balancerType = EBalancerType::ScatterCPU;
                    break;
                case EResourceToBalance::Memory:
                    balancerType = EBalancerType::ScatterMemory;
                    break;
                case EResourceToBalance::Network:
                    balancerType = EBalancerType::ScatterNetwork;
                    break;
                case EResourceToBalance::ComputeResources:
                    balancerType = EBalancerType::Scatter;
                    break;
            }
        }
        if (stats.MaxUsage >= GetMaxNodeUsageToKick() && stats.MinUsage < minUsageToKick) {
            std::vector<TNodeId> overloadedNodes;
            for (const auto& [nodeId, nodeInfo] : Nodes) {
                if (nodeInfo.IsAlive() && nodeInfo.IsOverloaded()) {
                    overloadedNodes.emplace_back(nodeId);
                }
            }

            if (!overloadedNodes.empty()) {
                if (!g_no_moves_checkbox->IsChecked() || LastTrigger == balancerType || LastMovements != 0) {
                    StartHiveBalancer({
                        .Type = EBalancerType::Emergency,
                        .MaxMovements = (int)GetMaxMovementsOnEmergencyBalancer(),
                        .FilterNodeIds = std::move(overloadedNodes),
                    });
                    return;
                }
            }
        }

        if (scatteredResource) {
            StartHiveBalancer({
                .Type = balancerType,
                .MaxMovements = (int)GetMaxMovementsOnAutoBalancer(),
                .ResourceToBalance = *scatteredResource,
            });
            return;
        }
    }

    TNodeInfo* SelectNode(const std::vector<THive::TSelectedNode>& selectedNodes) {
        if (selectedNodes.empty()) {
            return nullptr;
        }
        double sumUsage = 0;
        double maxUsage = 0;
        for (const TSelectedNode& selectedNode : selectedNodes) {
            double usage = selectedNode.Usage;
            sumUsage += usage;
            maxUsage = std::max(maxUsage, usage);
        }
        double sumAvail = maxUsage * selectedNodes.size() - sumUsage;
        if (sumAvail > 0) {
            std::uniform_real_distribution rnd;
            double pos = rnd(g_random) * sumAvail;
            for (const TSelectedNode& selectedNode : selectedNodes) {
                double avail = maxUsage - selectedNode.Usage;
                if (pos < avail) {
                    return selectedNode.Node;
                } else {
                    pos -= avail;
                }
            }
        }
        std::uniform_int_distribution<size_t> rnd(0, selectedNodes.size() - 1);
        return selectedNodes[rnd(g_random)].Node;
    }

    std::vector<THive::TSelectedNode> SelectMaxPriorityNodes(std::vector<TSelectedNode> selectedNodes, const TTabletInfo& tablet, std::unordered_map<TDataCenterId, Si32>& dcPriority) const
    {
        Si32 priority = std::numeric_limits<Si32>::min();
        for (const TSelectedNode& selectedNode : selectedNodes) {
            priority = std::max(priority, selectedNode.Node->GetPriorityForTablet(tablet, dcPriority));
        }

        auto it = std::partition(selectedNodes.begin(), selectedNodes.end(), [&] (const TSelectedNode& selectedNode) {
            return selectedNode.Node->GetPriorityForTablet(tablet, dcPriority) == priority;
        });

        selectedNodes.erase(it, selectedNodes.end());

        return selectedNodes;
    }

    TBestNodeResult FindBestNode(const TTabletInfo& tablet) {
        std::vector<TDataCenterId> dcs;
        if (tablet.Category && tablet.Category->StickTogetherInDC && g_stick_checkbox->IsChecked()) {
            std::unordered_map<TDataCenterId, Ui32> dcTablets;
            for (TFullTabletId tabletId : tablet.Category->Tablets) {
                TTabletInfo* tab = FindTablet(tabletId);
                if (tab->IsAlive()) {
                    TDataCenterId dc = FindNode(tab->NodeId)->GetDataCenter();
                    dcTablets[dc]++;
                    *Log() << "Increment for dc " << dc << std::endl;
                }
            }
            if (!dcTablets.empty()) {
                for (const auto& [dc, count] : dcTablets) {
                    dcs.push_back(dc);
                }
                std::sort(dcs.begin(), dcs.end(), [&](TDataCenterId a, TDataCenterId b) -> bool {
                    return dcTablets[a] > dcTablets[b];
                });
            }
        }

        std::vector<std::vector<TNodeInfo*>> candidateGroups;
        candidateGroups.resize(dcs.size() + 1);
        std::unordered_map<TDataCenterId, std::vector<TNodeInfo*>*> indexDC2Group;
        std::unordered_map<TDataCenterId, Si32> dcPriority;
        for (size_t numGroup = 0; numGroup < dcs.size(); ++numGroup) {
            indexDC2Group[dcs[numGroup]] = candidateGroups.data() + numGroup;
            if (g_priority_checkbox->IsChecked()) {
                dcPriority[dcs[numGroup]] = dcs.size() - numGroup;
            }
        }
        for (auto it = Nodes.begin(); it != Nodes.end(); ++it) {
            TNodeInfo* nodeInfo = &it->second;
            if (nodeInfo->IsAlive()) {
                TDataCenterId dataCenterId = nodeInfo->GetDataCenter();
                auto itDataCenter = indexDC2Group.find(dataCenterId);
                if (!g_priority_checkbox->IsChecked() && itDataCenter != indexDC2Group.end()) {
                    itDataCenter->second->push_back(nodeInfo);
                    *Log() << "Node " << nodeInfo->Id << " in dc " << dataCenterId << " in group " << (Ui64)itDataCenter->second << std::endl;
                } else {
                    candidateGroups.back().push_back(nodeInfo);
                    *Log() << "Node " << nodeInfo->Id << " in dc " << dataCenterId << " in group [default]" << std::endl; 
                }
            }
        }

        std::vector<TSelectedNode> selectedNodes;

        for (auto itCandidateNodes = candidateGroups.begin(); itCandidateNodes != candidateGroups.end(); ++itCandidateNodes) {
            const std::vector<TNodeInfo*>& candidateNodes(*itCandidateNodes);

            selectedNodes.reserve(candidateNodes.size());

            for (auto it = candidateNodes.begin(); it != candidateNodes.end(); ++it) {
                TNodeInfo& nodeInfo = *(*it);
                double usage = nodeInfo.GetNodeUsageForTablet(tablet);
                selectedNodes.emplace_back(usage, &nodeInfo);
            }

            if (!selectedNodes.empty()) {
                break;
            }
        }

        TNodeInfo* selectedNode = nullptr;
        if (!selectedNodes.empty()) {
            selectedNodes = SelectMaxPriorityNodes(std::move(selectedNodes), tablet, dcPriority);
            selectedNode = SelectNode(selectedNodes);
        }
        if (selectedNode != nullptr) {
            return selectedNode;
        } else {
            return TNoNodeFound();
        }
    }
    
    bool IsTabletMoveExpedient(const TTabletInfo& tablet, const TNodeInfo& node) const {
        if (g_guess_usage_checkbox->IsChecked() && Nodes.at(tablet.NodeId).GetNodeUsageForTablet(tablet) <= node.GetNodeUsageForTablet(tablet)) {
            return false;
        }
        if (Nodes.at(tablet.NodeId).IsOverloaded() && !node.IsOverloaded()) {
            return true;
        }
        std::vector<TResourceNormalizedValues> values;
        std::size_t oldNode = std::numeric_limits<std::size_t>::max();
        std::size_t newNode = std::numeric_limits<std::size_t>::max();
        values.reserve(Nodes.size());
        for (const auto& ni : Nodes) {
            if (ni.second.IsAlive()) {
                if (ni.first == node.Id)
                    newNode = values.size();
                if (ni.first == tablet.NodeId)
                    oldNode = values.size();
                values.push_back(NormalizeRawValues(ni.second.RawValues, ni.second.GetResourceMaximumValues()));
            }
        }

        if (oldNode == std::numeric_limits<std::size_t>::max()
                || newNode == std::numeric_limits<std::size_t>::max()) {
            return false;
        }

        auto tabletResources = tablet.ResourceValues;
    //    NMetrics::TResourceMetrics::TResourceNormalizedValues oldValues = values[oldNode];
    //    NMetrics::TResourceMetrics::TResourceNormalizedValues newValues = values[newNode] + NMetrics::TResourceMetrics::Normalize(tabletResources, node.GetResourceMaximumValues());
    //    return sum(newValues) < sum(oldValues);

        TResourceNormalizedValues beforeStDev = GetStDev(values);
        values[oldNode] -= NormalizeRawValues(tabletResources);
        values[newNode] += NormalizeRawValues(tabletResources);
        TResourceNormalizedValues afterStDev = GetStDev(values);
        tablet.FilterRawValues(beforeStDev);
        tablet.FilterRawValues(afterStDev);
        double before = max(beforeStDev);
        double after = max(afterStDev);
        bool result = after < before;
        if (!result) {
            *Log() << "Move " << tablet.Id << " from " << tablet.NodeId << " to " << node.Id << " is not expedient" << std::endl;
        }
        return result;
    }

    void MoveTablet(TTabletInfo& tablet, TNodeId from, TNodeId to) {
        *Log() << g_t << ": Balancer moving tablet " << tablet.Id << " from node " << from << " to node " << to << std::endl; 
        g_last_tablet_moved = tablet.Id;
        g_last_node_moved = from;
        TNodeInfo* fromNode = FindNode(from);
        TNodeInfo* toNode = FindNode(to);
        if (fromNode) {
            fromNode->RemoveTablet(tablet);
        }
        toNode->AddTablet(tablet);
        tablet.Alive = true;
    }

    void FillNode(TNodeId nodeId) {
        for (auto& [tabletId, tablet] : Tablets) {
            auto result = FindBestNode(tablet);
            if (std::holds_alternative<TNodeInfo*>(result)) {
                TNodeInfo* node = std::get<TNodeInfo*>(result);
                if (node->Id == nodeId) {
                    MoveTablet(tablet, tablet.NodeId, nodeId);
                }
            }
        }
    }

    void KillNode(TNodeId nodeId) {
        auto it = Nodes.find(nodeId);
        if (it == Nodes.end()) {
            return;
        }
        auto& node = it->second;
        node.Alive = false;
        std::vector<TTabletInfo*> tablets(node.Tablets.begin(), node.Tablets.end());
        for (auto* tablet : tablets) {
            tablet->Alive = false;
            tablet->FailedNodeId = nodeId;
        }
        std::sort(tablets.begin(), tablets.end(), [](const auto& lhs, const auto& rhs) {
            return lhs->GetPriority() < rhs->GetPriority();     
        });
        for (auto* tablet : tablets) {
            auto result = FindBestNode(*tablet);
            if (std::holds_alternative<TNodeInfo*>(result)) {
                TNodeInfo* node = std::get<TNodeInfo*>(result);
                MoveTablet(*tablet, nodeId, node->Id);
            }
        }
        node.Drawn = false;
    }

    void RessurectNode(TNodeId nodeId) {
        auto it = Nodes.find(nodeId);
        if (it == Nodes.end()) {
            return;
        }
        it->second.Alive = true;
        it->second.Drawn = false;
        it->second.StartTime = Now();
        it->second.RestartTimestamps.push_back(Now());
    }

};

THive g_hive;
std::unordered_set<TFullTabletId> g_burst_tablets;

void ChangeTabletMetrics() {
    g_max_tablet_metrics = g_tablet_metrics_scrollbar->GetValue();
    const auto max_cpu = g_max_tablet_metrics;
    std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
    for (auto& [tabletId, tablet] : g_hive.Tablets) {
        if (tablet.Type != ETabletType::ColumnShard) {
            tablet.UpdateResourceUsage({max_cpu, std::nullopt, std::nullopt});
        }
    }
}

struct THiveBalancer {
    THive* Hive = &g_hive;
    int Movements = 0;
    TBalancerSettings Settings;
    std::vector<TNodeId> Nodes;
    std::vector<TNodeId>::iterator NextNode;
    std::vector<TFullTabletId> Tablets;
    std::vector<TFullTabletId>::iterator NextTablet;

    bool CanKickNextTablet() const {
        return (Settings.MaxMovements == 0 || Movements < Settings.MaxMovements);
    }

    void BalanceNodes() {
        std::vector<TNodeInfo*> nodes;
        if (!Settings.FilterNodeIds.empty()) {
            nodes.reserve(Settings.FilterNodeIds.size());
            for (TNodeId nodeId : Settings.FilterNodeIds) {
                TNodeInfo* node = Hive->FindNode(nodeId);
                if (node != nullptr && node->IsAlive()) {
                    nodes.emplace_back(node);
                }
            }
        } else {
            nodes.reserve(Hive->Nodes.size());
            for (auto& [nodeId, nodeInfo] : Hive->Nodes) {
                if (nodeInfo.IsAlive()) {
                    nodes.emplace_back(&nodeInfo);
                }
            }
        }

        ::BalanceNodes(nodes, Settings.ResourceToBalance);

        Nodes.reserve(nodes.size());
        for (auto node : nodes) {
            Nodes.push_back(node->Id);
        }


        NextNode = Nodes.begin();
        Tablets.clear();
    }

    std::optional<TFullTabletId> GetNextTablet(TInstant now) {
        for (; Tablets.empty() || NextTablet == Tablets.end(); ++NextNode) {
            if (NextNode == Nodes.end()) {
                return std::nullopt;
            }
            TNodeInfo* node = Hive->FindNode(*NextNode);
            if (node == nullptr) {
                continue;
            }
            const std::set<TTabletInfo*>& nodeTablets = node->Tablets;
            std::vector<TTabletInfo*> tablets;
            tablets.reserve(nodeTablets.size());
            for (TTabletInfo* tablet : nodeTablets) {
                if (tablet->IsGoodForBalancer(now) && 
                    tablet->HasMetric(Settings.ResourceToBalance)) {
                    tablets.emplace_back(tablet);
                }
            }
            if (!tablets.empty()) {
                BalanceTablets(tablets.begin(), tablets.end(), Settings.ResourceToBalance); 
                Tablets.clear();
                Tablets.reserve(tablets.size());
                for (auto tablet : tablets) {
                    Tablets.push_back(tablet->GetFullTabletId());
                }
            }
            NextTablet = Tablets.begin();
        }
        return *(NextTablet++);
    }

    void KickNextTablet() {
        if (Settings.MaxMovements != 0 && Movements >= Settings.MaxMovements) {
            return;
        }

        TInstant now = Now();

        while (CanKickNextTablet()) {
            std::optional<TFullTabletId> tabletId = GetNextTablet(now);
            if (!tabletId) {
                break;
            }
            TTabletInfo* tablet = Hive->FindTablet(*tabletId);
            if (tablet == nullptr) {
                continue;
            }
            THive::TBestNodeResult result = Hive->FindBestNode(*tablet);
            if (std::holds_alternative<TNodeInfo*>(result)) {
                TNodeInfo* node = std::get<TNodeInfo*>(result);
                if (node->Id != tablet->NodeId && Hive->IsTabletMoveExpedient(*tablet, *node)) {
                    tablet->MakeBalancerDecision(now);
                    ++Movements;
                    Hive->MoveTablet(*tablet, tablet->NodeId, node->Id);
                }
            }
        }
    }

    THiveBalancer(TBalancerSettings&& settings) : Settings(std::move(settings)) {}

    void Run() {
        BalanceNodes();
        KickNextTablet();
        g_hive.LastTrigger = Settings.Type;
        g_hive.LastMovements = Movements;
    }
};

void StartHiveBalancer(TBalancerSettings&& settings) {
    *Log() << Now() << "Start Hive balancer, resource = " << (int)settings.ResourceToBalance << std::endl;
    g_last_balancer_trigger = Now();
    THiveBalancer balancer(std::move(settings));
    balancer.Run();
}

void UpdateResourceValues(TNodeId nodeId, const TResourceRawValues& delta) {
    TNodeInfo* node = g_hive.FindNode(nodeId);
    if (node == nullptr) {
        return;
    }
    node->RawValues += delta;
}

TTabletInfo& MakeDS() {
    const auto max_cpu = g_max_tablet_metrics;
    std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
    const Ui64 max_memory = (double)g_max_tablet_metrics * std::get<EResource::Memory>(MAXIMUM_VALUES) / std::get<EResource::CPU>(MAXIMUM_VALUES);
    std::uniform_int_distribution<Ui64> get_memory(0, max_memory);
    TFullTabletId id = g_hive.Tablets.size() + 1;
    g_hive.Tablets.emplace(id, TTabletInfo(id, {get_cpu(g_random), get_memory(g_random), 0, 0}, ETabletType::DataShard));
    return g_hive.Tablets.at(id);
}

TTabletInfo& MakeCS() {
    TFullTabletId id = g_hive.Tablets.size() + 1;
    g_hive.Tablets.emplace(id, TTabletInfo(id, {0, 0, 0, 1}, ETabletType::ColumnShard));
    return g_hive.Tablets.at(id);
}

TTabletInfo& MakeC() {
    const auto max_cpu = g_max_tablet_metrics;
    std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
    const Ui64 max_memory = (double)g_max_tablet_metrics * std::get<EResource::Memory>(MAXIMUM_VALUES) / std::get<EResource::CPU>(MAXIMUM_VALUES);
    std::uniform_int_distribution<Ui64> get_memory(0, max_memory);
    TFullTabletId id = g_hive.Tablets.size() + 1;
    g_hive.Tablets.emplace(id, TTabletInfo(id, {get_cpu(g_random), get_memory(g_random), 0, 0}, ETabletType::Coordinator));
    return g_hive.Tablets.at(id);
}

void CreateDS() {
    TTabletInfo& tablet = MakeDS();
    auto result = g_hive.FindBestNode(tablet);
    if (std::holds_alternative<TNodeInfo*>(result)) {
        g_hive.MoveTablet(tablet, 0, std::get<TNodeInfo*>(result)->Id);
    }
}

void CreateCS() {
    TTabletInfo& tablet = MakeCS();
    auto result = g_hive.FindBestNode(tablet);
    if (std::holds_alternative<TNodeInfo*>(result)) {
        g_hive.MoveTablet(tablet, 0, std::get<TNodeInfo*>(result)->Id);
    }
}

void CreateC() {
    TTabletInfo& tablet = MakeC();
    auto result = g_hive.FindBestNode(tablet);
    if (std::holds_alternative<TNodeInfo*>(result)) {
        g_hive.MoveTablet(tablet, 0, std::get<TNodeInfo*>(result)->Id);
    }
}

struct IScenario {
    virtual void CreateInitialState() {
    }

    virtual void UpdateModel() {
        g_t += g_dt;
        g_hive.ProcessTabletBalancer(); 

        std::uniform_real_distribution<double> coordinator_usage(0.8, 1);
        //std::uniform_real_distribution<double> normal_usage(0.3, 0.7);
        for (auto& [nodeId, node] : g_hive.Nodes) {
            unsigned coordinators = 0;
            for (const auto* tablet : node.Tablets) {
                coordinators += (tablet->Type == ETabletType::Coordinator);
            } 
            if (coordinators) {
                node.UpdateResourceTotalUsage(coordinator_usage(g_random) * coordinators);
            } else {
                node.UpdateResourceTotalUsage(0);
            }
        }
    }
};

struct BurstScenario : IScenario {
    void CreateInitialState() override {
        const size_t num_nodes = 8;
        const size_t num_tablets = 42;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i));
        }

        for (TFullTabletId i = 1; i <= num_tablets; ++i) {
            MakeDS();
        }

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() override {
        IScenario::UpdateModel();

        const auto max_cpu = std::get<EResource::CPU>(MAXIMUM_VALUES) * g_hive.Nodes.size() * 6 / g_hive.Tablets.size() / 5;
        std::poisson_distribution is_burst(g_dt / g_hive.Tablets.size());
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            if (is_burst(g_random)) {
                //uint64_t cpu = std::get<EResource::CPU>(MAXIMUM_VALUES) * .84;
                tablet.UpdateResourceUsage({max_cpu, std::nullopt, std::nullopt});
                g_burst_tablets.insert(tabletId);
                *Log() << g_t << ": burst on tablet " << tabletId << " on node " << tablet.NodeId << std::endl;
            } else if (g_burst_tablets.contains(tabletId)) {
                std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
                tablet.UpdateResourceUsage({get_cpu(g_random), std::nullopt, std::nullopt});
                g_burst_tablets.erase(tabletId);
            }
        }
    }
};

struct HtapScenario : IScenario {
    void CreateInitialState() override {
        const size_t num_nodes = 8;
        const size_t num_tablets = 42;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            TNodeInfo node(i);
            g_hive.Nodes.emplace(i, node);
        }
        for (TFullTabletId i = 1; i <= num_tablets; ++i) {
            if (i % 2 == 0) {
                MakeDS();
            } else {
                MakeCS();
            }
        }

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() override {
        IScenario::UpdateModel();

        const auto max_cpu = g_max_tablet_metrics;
        std::poisson_distribution is_burst(g_dt / g_hive.Tablets.size());
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            if (tablet.Type == ETabletType::DataShard && is_burst(g_random)) {
                //uint64_t cpu = std::get<EResource::CPU>(MAXIMUM_VALUES) * .84;
                tablet.UpdateResourceUsage({max_cpu, std::nullopt, std::nullopt});
                g_burst_tablets.insert(tabletId);
                *Log() << g_t << ": burst on tablet " << tabletId << " on node " << tablet.NodeId << std::endl;
            } else if (g_burst_tablets.contains(tabletId)) {
                std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
                tablet.UpdateResourceUsage({get_cpu(g_random), std::nullopt, std::nullopt});
                g_burst_tablets.erase(tabletId);
            }
        }

        for (auto& [nodeId, node] : g_hive.Nodes) {
            node.UpdateResourceTotalUsage(node.Tablets.size() / 6.0);
        }
    }
};

struct OneFatTabletScenario : IScenario {
    void CreateInitialState() override {
        const size_t num_nodes = 8;
        const size_t num_tablets = 42;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i));
        }

        for (TFullTabletId i = 1; i < num_tablets; ++i) {
            MakeDS();
        }
        MakeC();

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() override {
        IScenario::UpdateModel();
    }
};

struct SmallScenario : IScenario {
    void CreateInitialState() override {
        const size_t num_nodes = 3;
        const size_t num_tablets = 25;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i));
        }
        for (TFullTabletId i = 1; i <= num_tablets; ++i) {
            MakeDS();
        }

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() {
        IScenario::UpdateModel(); 
    }
};

struct ColumnShardScenario : IScenario {
    void CreateInitialState() {
        const size_t num_nodes = 8;
        const size_t num_tablets = 42;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i));
        }

        for (TFullTabletId i = 1; i <= num_tablets; ++i) {
            g_hive.Tablets.emplace(i, TTabletInfo(i, {0, 0, 0, 1}, ETabletType::ColumnShard));
        }

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() {
        IScenario::UpdateModel();

        for (auto& [nodeId, node] : g_hive.Nodes) {
            node.UpdateResourceTotalUsage(node.Tablets.size() / 10.0);
        }
    }
};

struct BadDcScenario : IScenario {
    void CreateInitialState() override {
        const size_t num_nodes = 9;
        const size_t num_tablets = 70;
        g_usage_scrollbar->SetValue(0);
        g_tablet_metrics_scrollbar->SetValue(0.025 * std::get<EResource::CPU>(MAXIMUM_VALUES));
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i, (i - 1) / 3 + 1));
        }

        auto& firstC = MakeC(); // to put it in the first dc
        g_hive.FindNode(0)->AddTablet(firstC);
        for (int i = 0; i < num_tablets; ++i) {
            if (i % 14 == 0) {
                CreateC();
            } else {
                CreateDS();
            }
        }
    }

    void UpdateModel() override {
        IScenario::UpdateModel();
        std::uniform_int_distribution<TNodeId> pickNode(1, 3);
        std::poisson_distribution verify(5 * g_dt);
        if (verify(g_random)) {
            TNodeId nodeId = pickNode(g_random);
            g_hive.KillNode(nodeId);
            g_hive.RessurectNode(nodeId);
        }
    }
};

void AddNode() {
    TNodeId id = g_hive.Nodes.size() + 1;
    g_hive.Nodes.emplace(id, TNodeInfo(id));
    if (g_fill_nodes_checkbox->IsChecked()) {
        g_hive.FillNode(id);
    }
}

void DrawModel() {
    Si32 x = 50;
    for (const auto& [nodeId, node] : g_hive.Nodes) {
        Rgba color;
        if (!node.IsAlive()) {
            color = Rgba{180, 0, 64}; 
        } else if (node.NodeKind) {
            color = Rgba{32, 255, 64};
        } else {
            color = Rgba{0, 255, 100};
        }
        const std::vector<Rgba> dcColors = {color, Rgba(200, 0, 200), Rgba(0, 180, 120), Rgba(255, 128, 0)};
        DrawRectangle(Vec2Si32{x, 400}, Vec2Si32{x + 150, 550}, dcColors[node.GetDataCenter()]);
        DrawRectangle(Vec2Si32{x + 10, 400 + 10}, Vec2Si32{x + 150 - 10, 550 - 10}, color);
        Si32 x_tablet = x + 10;
        Si32 y_tablet = 540;
        for (const auto* tablet: node.Tablets) {
            Rgba color = tablet->GetColor();
            if (tablet->Id == g_last_tablet_moved) {
                DrawCircle(Vec2Si32{x_tablet, y_tablet}, 10, Rgba(255, 0, 0));
                DrawCircle(Vec2Si32{x_tablet, y_tablet}, 5, color);
            } else {
                DrawCircle(Vec2Si32{x_tablet, y_tablet}, 10, color);
            }
            x_tablet += 21;
            if (x_tablet > x + 140) {
                x_tablet = x + 10;
                y_tablet -= 21;
            }
        }
        std::stringstream text;
        text << static_cast<int>(node.GetNodeUsage() * 100) << "%";
        text << " " << Now() - node.StartTime;
        if (g_last_node_moved == node.Id) {
            text << " (*)";
        }
        g_font.Draw(text.str().c_str(), x, 580, kTextOriginTop);

        if (!node.Drawn) {
            std::shared_ptr<Button> button;
            button = gf.MakeButton();
            if (node.NodeKind == TNodeInfo::ENodeKind::ModernCPU) {
                button->SetText("Modern CPU");
            } else {
                button->SetText("Old CPU");
            }
            button->SetPos(Vec2Si32(x, 320));
            button->SetWidth(175);
            button->OnButtonClick = [id = nodeId]() {
                auto& node = g_hive.Nodes.at(id);
                node.NodeKind ^= 1;
                node.Drawn = false;
            };
            g_gui->AddChild(button);
            button = gf.MakeButton();
            button->SetPos(Vec2Si32(x, 250));
            button->SetWidth(175);
            if (node.IsAlive()) {
                button->SetText("Kill");
                button->OnButtonClick = [id = nodeId]() {
                    g_hive.KillNode(id);
                };
            } else {
                button->SetText("Start");
                button->OnButtonClick = [id = nodeId]() {
                    g_hive.RessurectNode(id);
                };
            }
            g_gui->AddChild(button);
            node.Drawn = true;
        }
        x += 200; 
    }
}

/// @brief Main entry point for the simulation
void EasyMain() {
  // Load GUI theme and font
  g_theme = std::make_shared<GuiTheme>();
  g_theme->Load("data/gui_theme.xml");
  g_font.Load("data/arctic_one_bmf.fnt");

  // Resize screen to match the desired resolution
  ResizeScreen(1920, 1080);

  // Create GUI elements
  gf.theme_ = g_theme;
  g_gui = gf.MakeTransparentPanel();
  std::shared_ptr<Button> button;

  button = gf.MakeButton();
  button->SetText("Slow Down");
  button->SetPos(Vec2Si32(16, 935));
  button->SetWidth(200);
  button->OnButtonClick = SlowDown;
  g_gui->AddChild(button);

  button = gf.MakeButton();
  button->SetText("Speed Up");
  button->SetPos(Vec2Si32(16+208, 935));
  button->SetWidth(200);
  button->OnButtonClick = SpeedUp;
  g_gui->AddChild(button);

  button = gf.MakeButton();
  button->SetText("Add node");
  button->SetPos(Vec2Si32(16+2*208, 935));
  button->SetWidth(200);
  button->OnButtonClick = AddNode;
  g_gui->AddChild(button);

  g_tablet_metrics_scrollbar = gf.MakeHorizontalScrollbar();
  g_tablet_metrics_scrollbar->SetMinValue(0);
  g_tablet_metrics_scrollbar->SetMaxValue(1'000'000'000);
  g_tablet_metrics_scrollbar->SetPos(Vec2Si32(16+3*208, 935));
  g_tablet_metrics_scrollbar->SetValue(g_max_tablet_metrics);
  g_tablet_metrics_scrollbar->OnScrollChange = ChangeTabletMetrics;
  g_gui->AddChild(g_tablet_metrics_scrollbar);

  g_usage_scrollbar = gf.MakeHorizontalScrollbar();
  g_usage_scrollbar->SetMinValue(0);
  g_usage_scrollbar->SetMaxValue(200);
  g_usage_scrollbar->SetPos(Vec2Si32(16+4*208, 935));
  g_usage_scrollbar->SetValue(100);
  g_gui->AddChild(g_usage_scrollbar);

  button = gf.MakeButton();
  button->SetText("Create DS");
  button->SetPos(Vec2Si32(16+5*208, 935));
  button->SetWidth(200);
  button->OnButtonClick = CreateDS;
  g_gui->AddChild(button);

  button = gf.MakeButton();
  button->SetText("Create CS");
  button->SetPos(Vec2Si32(16+6*208, 935));
  button->SetWidth(200);
  button->OnButtonClick = CreateCS;
  g_gui->AddChild(button);

  button = gf.MakeButton();
  button->SetText("Create C");
  button->SetPos(Vec2Si32(16+7*208, 935));
  button->SetWidth(200);
  button->OnButtonClick = CreateC;
  g_gui->AddChild(button);

  g_no_moves_checkbox = gf.MakeCheckbox();
  g_no_moves_checkbox->SetText("check no moves");
  g_no_moves_checkbox->SetPos(Vec2Si32(16+8*208, 935));
  g_gui->AddChild(g_no_moves_checkbox);

  g_guess_usage_checkbox = gf.MakeCheckbox();
  g_guess_usage_checkbox->SetText("guess usage");
  g_guess_usage_checkbox->SetPos(Vec2Si32(16+8*208, 985));
  g_gui->AddChild(g_guess_usage_checkbox);

  g_fill_nodes_checkbox = gf.MakeCheckbox();
  g_fill_nodes_checkbox->SetText("fill nodes");
  g_fill_nodes_checkbox->SetPos(Vec2Si32(16+8*208, 885));
  g_gui->AddChild(g_fill_nodes_checkbox);

  g_stick_checkbox = gf.MakeCheckbox();
  g_stick_checkbox->SetText("stick together");
  g_stick_checkbox->SetPos(Vec2Si32(16+8*208, 835));
  g_stick_checkbox->SetChecked(true);
  g_gui->AddChild(g_stick_checkbox);

  g_priority_checkbox = gf.MakeCheckbox();
  g_priority_checkbox->SetText("updated priority");
  g_priority_checkbox->SetPos(Vec2Si32(16+8*208, 785));
  g_priority_checkbox->SetChecked(false);
  g_gui->AddChild(g_priority_checkbox);



  // Create initial simulation state
  std::unique_ptr<IScenario> scenario;
  auto txt = ReadFile("scenario.txt");
  switch (txt[0]) {
      default:
      case '0':
          scenario = std::make_unique<BurstScenario>();
          break;
      case '1':
          scenario = std::make_unique<HtapScenario>();
          break;
      case '2':
          scenario = std::make_unique<OneFatTabletScenario>();
          break;
      case '3':
          scenario = std::make_unique<ColumnShardScenario>();
          break;
      case '4':
          scenario = std::make_unique<SmallScenario>();
          break;
      case '5':
          scenario = std::make_unique<BadDcScenario>();
          break;
  }

  scenario->CreateInitialState();

  // Main simulation loop
  double rt0 = Time();
  double rt1 = Time();
  double t_target = 0.0;

  while (!IsKeyDownward(kKeyEscape)) {
    // Update simulation time
    rt0 = rt1;
    rt1 = Time();
    double rdt = rt1 - rt0;
    t_target = g_t + rdt * g_t_mult;

    // Clear screen
    Clear();

    // Update simulation state until target time is reached or timeout
    while (g_t < t_target && Time() - rt1 < 0.1) {
      scenario->UpdateModel();
    }

    // Apply GUI input
    for (Si32 i = 0; i < InputMessageCount(); ++i) {
      g_gui->ApplyInput(GetInputMessage(i), nullptr);
    }

    // Draw simulation state
    DrawModel();
    g_gui->Draw(Vec2Si32(0, 0));

    // Draw time info
    char text[128];
    snprintf(text, sizeof(text), "Model time: %f s\nTarget multiplier: %f", g_t, g_t_mult);
    g_font.Draw(text, 20, ScreenSize().y - 20, kTextOriginTop);
    g_font.Draw("tablet metrics", 16+3*208, 1000, kTextOriginTop);
    g_font.Draw("total usage", 16+4*208, 1000, kTextOriginTop);

    // Show frame
    ShowFrame();
  }
}
