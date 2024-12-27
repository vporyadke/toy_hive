#pragma once
#include <algorithm>
#include <array>
#include <cmath>

namespace std {

template <>
struct make_signed<double> {
    using type = double;
};

}

using TInstant = double;
using TDuration = double;

TInstant Now();

struct TMetrics {
    std::optional<uint64_t> CPU;
    std::optional<uint64_t> Memory;
    std::optional<uint64_t> Network;

    bool HasCPU() const {
        return CPU.has_value();
    }

    bool HasMemory() const {
        return Memory.has_value();
    }

    bool HasNetwork() const {
        return Network.has_value();
    }

    uint64_t GetCPU() const {
        return CPU.value_or(0);
    }

    uint64_t GetMemory() const {
        return Memory.value_or(0);
    }

    uint64_t GetNetwork() const {
        return Network.value_or(0);
    }
};

namespace NMetrics {

// gauge resource type, supports absolute and delta population
template <typename ValueType>
class TGaugeValue {
public:
    using TType = ValueType;

    TGaugeValue(ValueType value = ValueType())
        : Value(value)
    {}

    std::make_signed_t<ValueType> Set(ValueType value) {
        std::make_signed_t<ValueType> diff =
                static_cast<std::make_signed_t<ValueType>>(value) -
                static_cast<std::make_signed_t<ValueType>>(Value);
        Value = value;
        return diff;
    }

    void Increment(std::make_signed_t<ValueType> value) { Value += value; }
    ValueType GetValue() const { return Value; }
    TGaugeValue operator +(const TGaugeValue& o) const { return TGaugeValue(Value + o.GetValue()); }
    TGaugeValue operator -(const TGaugeValue& o) const { return TGaugeValue(Value - o.GetValue()); }
    constexpr bool IsValueReady() const { return Value != ValueType(); }
    constexpr bool IsValueObsolete(TInstant) const { return false; }

protected:
    ValueType Value;
};

template <typename ValueType, size_t MaxCount = 20>
class TFastRiseAverageValue {
public:
    using TType = ValueType;

    TFastRiseAverageValue()
        : AccumulatorValue()
        , AccumulatorCount()
    {}

    void Push(ValueType value) {
        if (IsValueReady()) {
            ValueType currentValue = GetValue();
            if (value > currentValue) {
                ValueType newValue = (currentValue + value) / 2;
                AccumulatorCount++;
                AccumulatorValue = newValue * AccumulatorCount;
                return;
            }
        }
        if (AccumulatorCount >= MaxCount) {
            AccumulatorValue = AccumulatorValue / 2;
            AccumulatorCount /= 2;
        }
        AccumulatorValue += value;
        AccumulatorCount++;
    }

    ValueType GetValue() const {
        if (AccumulatorCount == 0) {
            return ValueType();
        }
        return AccumulatorValue / AccumulatorCount;
    }

    bool IsValueReady() const {
        return AccumulatorCount > 0;
    }

    bool IsValueStable() const {
        return AccumulatorCount >= 2;
    }

protected:
    ValueType AccumulatorValue;
    size_t AccumulatorCount;
};

struct TProto {
    TInstant LastBucketStartTime;
    std::vector<uint64_t> Values;
    uint64_t AllTimeMaximum;
};

class TMaximumValueVariableWindowUI64 : public TProto {
public:
    using TType = uint64_t;

    void SetValue(TType value, TInstant now = Now()) {
        if (TProto::AllTimeMaximum > 0 || MaximumValue > 0) { // ignoring initial value
            TProto::AllTimeMaximum = std::max(value, TProto::AllTimeMaximum);
        }
        TDuration elapsedCurrentBucket = now - TProto::LastBucketStartTime;
        if (TProto::Values.size() == 0 || elapsedCurrentBucket >= BucketDuration) {
            size_t bucketsPassed = 0;
            // new bucket
            if (TProto::Values.size() == 0) {
                TProto::LastBucketStartTime = now;
                bucketsPassed = 1;
            } else {
                bucketsPassed = elapsedCurrentBucket / BucketDuration;
                TProto::LastBucketStartTime = TProto::LastBucketStartTime + bucketsPassed * BucketDuration;
            }

            if (bucketsPassed >= BucketCount) {
                TProto::Values.clear();
                TProto::Values.push_back(value);
                MaximumValue = value;
            } else if (TProto::Values.size() + bucketsPassed > BucketCount) {
                auto shift = TProto::Values.size() + bucketsPassed - BucketCount;
                for (size_t pass = TProto::Values.size(); pass < BucketCount; ++pass) {
                    TProto::Values.push_back(value);
                }
                auto newBegin = TProto::Values.begin();
                std::advance(newBegin, shift);
                std::rotate(TProto::Values.begin(), newBegin, TProto::Values.end());
                auto last = TProto::Values.end();
                std::advance(last, -shift);
                std::fill(last, TProto::Values.end(), value);
                MaximumValue = *std::max_element(TProto::Values.begin(), TProto::Values.end());
            } else {
                for (size_t pass = 0; pass < bucketsPassed; ++pass) {
                    TProto::Values.push_back(value);
                }
                MaximumValue = std::max(MaximumValue, value);
            }
        } else {
            auto& lastBucketValue(*std::prev(TProto::Values.end()));
            lastBucketValue = std::max(lastBucketValue, value);
            MaximumValue = std::max(MaximumValue, value);
        }
    }

    void AdvanceTime(TInstant now) {
        // Nothing changed, last value is still relevant
        TType lastValue = {};
        if (!TProto::Values.empty()) {
            lastValue = *std::prev(TProto::Values.end());
        }
        SetValue(lastValue, now);
    }

    TType GetValue() const {
        return MaximumValue;
    }

protected:
    TType MaximumValue = {};
    size_t BucketCount = 24;
    TDuration BucketDuration = 60 / BucketCount;
};

} // NMetrics

using TMetricsMaximum = NMetrics::TMaximumValueVariableWindowUI64;

struct TTabletMetricsAggregates {
    TMetricsMaximum MaximumCPU;
    TMetricsMaximum MaximumMemory;
    TMetricsMaximum MaximumNetwork;
};
