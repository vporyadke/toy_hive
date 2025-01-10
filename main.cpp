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
std::shared_ptr<GuiTheme> g_theme;
std::shared_ptr<Panel> g_gui;


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
    return {0.1, 0.1, 0.1, 0};
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

struct TTabletInfo {
    TResourceRawValues ResourceValues;
    TTabletMetricsAggregates ResourceMetricsAggregates;
    TFullTabletId Id;
    ETabletType Type;
    TNodeId NodeId;
    TInstant LastBalancerDecisionTime = 0;

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
        return result;
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

    TTabletInfo(TFullTabletId id, const TResourceRawValues& metrics, ETabletType type = ETabletType::Dummy)
        : Id(id)
        , ResourceValues(metrics)
        , Type(type)
    {
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
    TResourceRawValues RawValues;
    NMetrics::TFastRiseAverageValue<double, 20> AveragedNodeTotalUsage;
    std::set<TTabletInfo*> Tablets;
    TNodeId Id;

    bool IsAlive() const {
        return true;
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
        return MAXIMUM_VALUES; 
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
        return usage;
    }

    void UpdateResourceTotalUsage(double usage) {
        AveragedNodeTotalUsage.Push(usage);
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

    TNodeInfo(TNodeId id) : Id(id) {}
};

struct TBalancerSettings {
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
        *Log() << "GetStats: ";
        Out(*Log(), minValues);
        *Log() << " vs ";
        Out(*Log(), maxValues);
        *Log() << std::endl;
        auto discrepancy = maxValues - minValues;
        auto& counterDiscrepancy = std::get<EResource::Counter>(discrepancy);
        *Log() << counterDiscrepancy * GetMaxResourceCounter() << std::endl;
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
        *Log() << g_t << " ProcessTabletBalancer: max usage = " << stats.MaxUsage << ", min usage = " << stats.MinUsage << ", scatter = ";
        Out(*Log(),  stats.ScatterByResource);
        *Log() << std::endl;

        double minUsageToKick = GetMaxNodeUsageToKick() - GetNodeUsageRangeToKick();
        if (stats.MaxUsage >= GetMaxNodeUsageToKick() && stats.MinUsage < minUsageToKick) {
            std::vector<TNodeId> overloadedNodes;
            for (const auto& [nodeId, nodeInfo] : Nodes) {
                if (nodeInfo.IsAlive() && nodeInfo.IsOverloaded()) {
                    overloadedNodes.emplace_back(nodeId);
                }
            }

            if (!overloadedNodes.empty()) {
                StartHiveBalancer({
                    .MaxMovements = (int)GetMaxMovementsOnEmergencyBalancer(),
                    .FilterNodeIds = std::move(overloadedNodes),
                });
                return;
            }
        }

        auto scatteredResource = CheckScatter(stats.ScatterByResource);
        if (scatteredResource) {
            StartHiveBalancer({
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

    TBestNodeResult FindBestNode(const TTabletInfo& tablet) {
        std::vector<TNodeInfo*> candidateNodes;
        for (auto it = Nodes.begin(); it != Nodes.end(); ++it) {
            TNodeInfo* nodeInfo = &it->second;
            if (nodeInfo->IsAlive()) {
                candidateNodes.push_back(nodeInfo);
            }
        }

        std::vector<TSelectedNode> selectedNodes;
        selectedNodes.reserve(candidateNodes.size());

        for (auto it = candidateNodes.begin(); it != candidateNodes.end(); ++it) {
            TNodeInfo& nodeInfo = *(*it);
            double usage = nodeInfo.GetNodeUsageForTablet(tablet);
            selectedNodes.emplace_back(usage, &nodeInfo);
        }

        TNodeInfo* selectedNode = nullptr;
        if (!selectedNodes.empty()) {
            selectedNode = SelectNode(selectedNodes);
        }
        if (selectedNode != nullptr) {
            return selectedNode;
        } else {
            return TNoNodeFound();
        }
    }
    
    bool IsTabletMoveExpedient(const TTabletInfo& tablet, const TNodeInfo& node) const {
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
        fromNode->RemoveTablet(tablet);
        toNode->AddTablet(tablet);
    }

};

THive g_hive;
std::unordered_set<TFullTabletId> g_burst_tablets;

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

struct BurstScenario {
    void CreateInitialState() {
        const size_t num_nodes = 8;
        const size_t num_tablets = 42;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i));
        }

        const auto max_cpu = std::get<EResource::CPU>(MAXIMUM_VALUES) * num_nodes * 6 / num_tablets / 5;
        std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
        for (TFullTabletId i = 1; i <= num_tablets; ++i) {
            g_hive.Tablets.emplace(i, TTabletInfo(i, {get_cpu(g_random), 0, 0, 0}));
        }

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() {
        g_t += g_dt;
        g_hive.ProcessTabletBalancer(); 

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

struct HtapScenario {
    void CreateInitialState() {
        const size_t num_nodes = 8;
        const size_t num_tablets = 42;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i));
        }

        const auto max_cpu = std::get<EResource::CPU>(MAXIMUM_VALUES) * num_nodes * 4 / num_tablets;
        std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
        const auto max_memory = std::get<EResource::Memory>(MAXIMUM_VALUES) * num_nodes * 12 / num_tablets / 5;
        std::uniform_int_distribution<Ui64> get_memory(0, max_memory);
        for (TFullTabletId i = 1; i <= num_tablets; ++i) {
            if (i % 2 == 0) {
                g_hive.Tablets.emplace(i, TTabletInfo(i, {get_cpu(g_random), get_memory(g_random), 0, 0}, ETabletType::DataShard));
            } else {
                g_hive.Tablets.emplace(i, TTabletInfo(i, {0, 0, 0, 1}, ETabletType::ColumnShard));
            }
        }

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() {
        g_t += g_dt;
        g_hive.ProcessTabletBalancer(); 

        const auto max_cpu = std::get<EResource::CPU>(MAXIMUM_VALUES) * g_hive.Nodes.size() * 12 / g_hive.Tablets.size() / 5;
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
            double weight = 0;
            for (const auto* tablet : node.Tablets) {
                weight += tablet->GetWeight();
            }
            node.UpdateResourceTotalUsage(weight * 1.1);
        }
    }
};

struct OneFatTabletScenario {
    void CreateInitialState() {
        const size_t num_nodes = 8;
        const size_t num_tablets = 42;
        for (TNodeId i = 1; i <= num_nodes; ++i) {
            g_hive.Nodes.emplace(i, TNodeInfo(i));
        }

        const auto max_cpu = std::get<EResource::CPU>(MAXIMUM_VALUES) * num_nodes / num_tablets;
        std::uniform_int_distribution<Ui64> get_cpu(0, max_cpu);
        const auto max_memory = std::get<EResource::Memory>(MAXIMUM_VALUES) * num_nodes / num_tablets;
        std::uniform_int_distribution<Ui64> get_memory(0, max_memory);
        for (TFullTabletId i = 1; i < num_tablets; ++i) {
            g_hive.Tablets.emplace(i, TTabletInfo(i, {get_cpu(g_random), get_memory(g_random), 0, 0}, ETabletType::DataShard));
        }
        g_hive.Tablets.emplace(42, TTabletInfo(42, {get_cpu(g_random), get_memory(g_random), 0, 0}, ETabletType::Coordinator));

        std::uniform_int_distribution<TNodeId> pick_node(1, num_nodes);
        for (auto& [tabletId, tablet] : g_hive.Tablets) {
            g_hive.FindNode(pick_node(g_random))->AddTablet(tablet);
        }
    }

    void UpdateModel() {
        g_t += g_dt;
        g_hive.ProcessTabletBalancer(); 

        std::uniform_real_distribution<double> coordinator_usage(0.8, 1);
        std::uniform_real_distribution<double> normal_usage(0.3, 0.7);
        auto coordinator_node = g_hive.Tablets.at(42).NodeId;
        for (auto& [nodeId, node] : g_hive.Nodes) {
            if (nodeId == coordinator_node) {
                node.UpdateResourceTotalUsage(coordinator_usage(g_random));
            } else {
                node.UpdateResourceTotalUsage(normal_usage(g_random));
            }
        }
    }
};

void DrawModel() {
    char text[128];

    Si32 x = 50;
    for (const auto& [nodeId, node] : g_hive.Nodes) {
        DrawRectangle(Vec2Si32{x, 400}, Vec2Si32{x + 150, 550}, Rgba{32, 255, 64});
        Si32 x_tablet = x + 10;
        Si32 y_tablet = 540;
        for (const auto* tablet: node.Tablets) {
            Rgba color = tablet->GetColor();
            if (tablet->Id == g_last_tablet_moved) {
                DrawCircle(Vec2Si32{x_tablet, y_tablet}, 10, Rgba(255, 0, 0));
                DrawCircle(Vec2Si32{x_tablet, y_tablet}, 8, color);
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
        if (g_last_node_moved == node.Id) {
            text << " (*)";
        }
        g_font.Draw(text.str().c_str(), x, 580, kTextOriginTop);
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
  GuiFactory gf;
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


  // Create initial simulation state
  OneFatTabletScenario scenario; // Change the type here for different scenarios
  scenario.CreateInitialState();

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
      scenario.UpdateModel();
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

    // Show frame
    ShowFrame();
  }
}
