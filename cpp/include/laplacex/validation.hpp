#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace laplacex {
inline void validate_samples(const std::vector<double>& t,
                             const std::vector<double>& f, std::size_t minimum = 4) {
    if (t.size() != f.size() || t.size() < minimum)
        throw std::invalid_argument("Need equal-length time and signal arrays with enough samples");
    for (std::size_t i = 0; i < t.size(); ++i)
        if (!std::isfinite(t[i]) || !std::isfinite(f[i]) || t[i] <= 0 ||
            (i && t[i] <= t[i-1]))
            throw std::invalid_argument("Need finite samples and strictly increasing t > 0");
}
} // namespace laplacex
