#include "laplacex/critical.hpp"
#include "laplacex/validation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace laplacex {
namespace {
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>; // contiguous columns

void validate_decomposition(const Vector& t, const Vector& f, const DecomposeConfig& c) {
    validate_samples(t, f);
    if (c.n_scales < 2 || c.n_log_points < 2 || c.min_ridge_len < 1 ||
        !std::isfinite(c.scale_min) || !std::isfinite(c.scale_max) ||
        c.scale_min <= 0 || c.scale_max < c.scale_min ||
        !std::isfinite(c.ridge_threshold) || c.ridge_threshold < 0 ||
        c.ridge_threshold > 1 || !std::isfinite(c.noise_sigma) ||
        !std::isfinite(c.snr_threshold) || c.snr_threshold < 0)
        throw std::invalid_argument("Invalid decomposition configuration");
    if (c.wavelet != "mexh")
        throw std::invalid_argument("Only the mexh wavelet is implemented");
    if (c.model_selection != "bic" && c.model_selection != "aicc")
        throw std::invalid_argument("model_selection must be 'bic' or 'aicc'");
}

double dot(const Vector& a, const Vector& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

// Householder reduction preserves the least-squares problem without forming
// normal equations. LM reuses this reduction when trying different damping.
void triangularize(Matrix& a, Vector& rhs) {
    const int n = static_cast<int>(rhs.size()), k = static_cast<int>(a.size());
    for (int j = 0; j < std::min(n, k); ++j) {
        double ss = 0;
        for (int i = j; i < n; ++i) ss += a[j][i]*a[j][i];
        double norm = std::sqrt(ss);
        if (norm < 1e-14) continue;
        double alpha = -std::copysign(norm, a[j][j]);
        Vector v(n-j);
        for (int i = j; i < n; ++i) v[i-j] = a[j][i]/norm;
        v[0] -= alpha/norm;
        double beta = 2/dot(v, v);
        for (int q = j+1; q < k; ++q) {
            double product = 0;
            for (int i = j; i < n; ++i) product += v[i-j]*a[q][i];
            for (int i = j; i < n; ++i) a[q][i] -= beta*product*v[i-j];
        }
        double product = 0;
        for (int i = j; i < n; ++i) product += v[i-j]*rhs[i];
        for (int i = j; i < n; ++i) rhs[i] -= beta*product*v[i-j];
        a[j][j] = alpha;
    }
    // The reflector's tail is workspace, not part of R.
    for (int j = 0; j < k; ++j)
        std::fill(a[j].begin()+std::min(j+1, n), a[j].end(), 0);
}

// Scaled Householder QR: no normal equations, including the damped LM solve.
Vector least_squares(Matrix a, Vector rhs) {
    const int n = static_cast<int>(rhs.size()), k = static_cast<int>(a.size());
    Vector scaling(k), result(k, 0);
    for (int j = 0; j < k; ++j) {
        scaling[j] = std::sqrt(dot(a[j], a[j]));
        if (scaling[j] > 0)
            for (auto& value : a[j]) value /= scaling[j];
    }
    triangularize(a, rhs);
    for (int j = std::min(n, k)-1; j >= 0; --j) {
        if (std::abs(a[j][j]) < 1e-13) continue;
        double value = rhs[j];
        for (int q = j+1; q < k; ++q) value -= a[q][j]*result[q];
        result[j] = value/a[j][j];
    }
    for (int j = 0; j < k; ++j)
        result[j] = scaling[j] > 0 ? result[j]/scaling[j] : 0;
    return result;
}

Matrix basis(const Vector& t, const Vector& rates) {
    Matrix e(rates.size(), Vector(t.size()));
    for (std::size_t j = 0; j < rates.size(); ++j)
        for (std::size_t i = 0; i < t.size(); ++i) e[j][i] = std::exp(-rates[j]*t[i]);
    return e;
}

struct Fit {
    Vector rates, amplitudes, residual;
    Matrix e;
    double rss = std::numeric_limits<double>::infinity();
};

void evaluate(Fit& fit, const Vector& t, const Vector& y) {
    fit.e = basis(t, fit.rates);
    fit.residual = y;
    long double rss = 0;
    for (std::size_t i = 0; i < t.size(); ++i) {
        long double prediction = 0;
        for (std::size_t j = 0; j < fit.rates.size(); ++j)
            prediction += static_cast<long double>(fit.amplitudes[j])*fit.e[j][i];
        fit.residual[i] = y[i]-static_cast<double>(prediction);
        rss += static_cast<long double>(fit.residual[i])*fit.residual[i];
    }
    fit.rss = static_cast<double>(rss);
}

Fit refine(const Vector& t, const Vector& y, Vector rates, bool optimize) {
    Fit fit;
    fit.rates = std::move(rates);
    fit.e = basis(t, fit.rates);
    fit.amplitudes = least_squares(fit.e, y);
    evaluate(fit, t, y);
    if (!optimize) return fit;
    const int n = static_cast<int>(t.size()), k = static_cast<int>(fit.rates.size());
    const double low = .005/t.back(), high = 8/t.front();
    double damping = 1e-3;
    const double floor = 1e-30*dot(y, y);
    for (int iteration = 0; iteration < 100; ++iteration) {
        Matrix jacobian(2*k, Vector(n));
        Vector norms(2*k);
        for (int j = 0; j < k; ++j)
            for (int i = 0; i < n; ++i) {
                jacobian[j][i] = fit.e[j][i];
                jacobian[j+k][i] = -fit.amplitudes[j]*fit.rates[j]*t[i]*fit.e[j][i];
            }
        for (int j = 0; j < 2*k; ++j) {
            norms[j] = std::max(1e-100, std::sqrt(dot(jacobian[j], jacobian[j])));
            for (auto& value : jacobian[j]) value /= norms[j];
        }
        Vector projected = fit.residual;
        int rows = n;
        if (k > 1) {
            triangularize(jacobian, projected);
            rows = 2*k;
            for (auto& column : jacobian) column.resize(rows);
            projected.resize(rows);
        }
        bool accepted = false;
        const double previous = fit.rss;
        for (int trial = 0; trial < 12; ++trial) {
            Matrix augmented = jacobian;
            Vector rhs = projected;
            rhs.resize(rows+2*k, 0);
            for (int j = 0; j < 2*k; ++j) {
                augmented[j].resize(rows+2*k, 0);
                augmented[j][rows+j] = std::sqrt(damping);
            }
            auto step = least_squares(std::move(augmented), std::move(rhs));
            Fit candidate = fit;
            double log_step = 0;
            for (int j = 0; j < k; ++j)
                log_step = std::max(log_step, std::abs(step[j+k]/norms[j+k]));
            double factor = log_step > 1.5 ? 1.5/log_step : 1;
            for (int j = 0; j < k; ++j) {
                candidate.amplitudes[j] += factor*step[j]/norms[j];
                candidate.rates[j] = std::clamp(
                    fit.rates[j]*std::exp(factor*step[j+k]/norms[j+k]), low, high);
            }
            evaluate(candidate, t, y);
            if (std::isfinite(candidate.rss) && candidate.rss < fit.rss) {
                fit = std::move(candidate);
                damping = std::max(1e-16, damping*.15);
                accepted = true;
                break;
            }
            damping = std::min(1e12, damping*8);
        }
        if (!accepted || fit.rss <= floor) break;
        if (previous-fit.rss <= 2e-15*fit.rss) break;
    }
    // Refit amplitudes at the final rates without a damping penalty.
    auto amplitudes = least_squares(fit.e, y);
    Fit polished = fit;
    polished.amplitudes = std::move(amplitudes);
    evaluate(polished, t, y);
    if (polished.rss < fit.rss) fit = std::move(polished);
    return fit;
}

double median(Vector v) {
    if (v.empty()) return 0;
    auto middle = v.begin()+v.size()/2;
    std::nth_element(v.begin(), middle, v.end());
    double value = *middle;
    if (v.size()%2 == 0) value = .5*(value+*std::max_element(v.begin(), middle));
    return value;
}

// Local linear interpolation removes trends even on irregular time grids.
double estimate_noise(const Vector& t, const Vector& y) {
    Vector differences;
    for (std::size_t i = 1; i+1 < t.size(); ++i) {
        double weight = (t[i]-t[i-1])/(t[i+1]-t[i-1]);
        differences.push_back((y[i]-(1-weight)*y[i-1]-weight*y[i+1])
                             /std::sqrt(1+(1-weight)*(1-weight)+weight*weight));
    }
    double center = median(differences);
    for (auto& v : differences) v = std::abs(v-center);
    return median(std::move(differences))/.6744897501960817;
}

struct Candidate { double rate, strength; int votes = 1; };

// Mexican-Hat zero crossings locate extrema of the Gaussian-smoothed log
// derivative. Constant edge extension and exact DC cancellation remove the
// artificial boundary modes that arose from transforming t*f(t).
struct LogWavelets {
    struct Kernel { Vector derivative, mexican, noise_norm; double scale; };
    Vector x, alpha;
    std::vector<std::size_t> left;
    std::vector<Kernel> kernels;
    const DecomposeConfig& cfg;
    double dx;
    Vector noise_weights;

    LogWavelets(const Vector& t, const DecomposeConfig& c) : cfg(c) {
        int n = c.n_log_points;
        noise_weights.resize(t.size());
        dx = (std::log(t.back())-std::log(t.front()))/(n-1);
        x.resize(n); alpha.resize(n); left.resize(n);
        std::size_t j = 0;
        for (int i = 0; i < n; ++i) {
            x[i] = std::log(t.front())+i*dx;
            double ti = std::exp(x[i]);
            while (j+2 < t.size() && t[j+1] < ti) ++j;
            left[i] = j;
            alpha[i] = std::clamp((ti-t[j])/(t[j+1]-t[j]), 0.0, 1.0);
        }
        for (int jscale = 0; jscale < c.n_scales; ++jscale) {
            double a = std::exp(std::log(c.scale_min)+(std::log(c.scale_max)
                               -std::log(c.scale_min))*jscale/(c.n_scales-1));
            a = std::max(a, 1.25*dx);
            int radius = static_cast<int>(std::ceil(std::min(double(n-1), 4.5*(a/dx))));
            Kernel kernel; kernel.scale = a;
            kernel.derivative.resize(radius+1); kernel.mexican.resize(radius+1);
            for (int d = 1; d <= radius; ++d) {
                double u = d*dx/a;
                double gaussian = std::exp(-.5*u*u)*dx/(a*std::sqrt(2*std::acos(-1.0)));
                kernel.derivative[d] = d*dx/(a*a)*gaussian;
                kernel.mexican[d] = (1-u*u)*gaussian;
            }
            // Propagate white sample noise through interpolation AND filtering.
            // Oversampled early times share observations; treating interpolated
            // grid values as independent gives incorrect detection thresholds.
            // Only crossings surviving the relative threshold need a noise
            // norm. Cache those exact variances across residual searches.
            kernel.noise_norm.assign(n, -1);
            kernels.push_back(std::move(kernel));
        }
    }

    double noise_norm(Kernel& kernel, int i) {
        if (kernel.noise_norm[i] < 0) {
            int n = static_cast<int>(x.size());
            int radius = static_cast<int>(kernel.derivative.size())-1;
            auto& weights = noise_weights;
            std::size_t first = left[std::max(0, i-radius)];
            std::size_t last = left[std::min(n-1, i+radius)]+1;
            std::fill(weights.begin()+first, weights.begin()+last+1, 0);
            for (int d = 1; d <= radius; ++d) {
                int l = std::max(0, i-d), r = std::min(n-1, i+d);
                double w = kernel.derivative[d];
                weights[left[l]] += w*(1-alpha[l]);
                weights[left[l]+1] += w*alpha[l];
                weights[left[r]] -= w*(1-alpha[r]);
                weights[left[r]+1] -= w*alpha[r];
            }
            double variance = 0;
            for (std::size_t j = first; j <= last; ++j) variance += weights[j]*weights[j];
            kernel.noise_norm[i] = std::sqrt(variance);
        }
        return kernel.noise_norm[i];
    }

    std::vector<Candidate> candidates(const Vector& residual, double sigma) {
        int n = static_cast<int>(x.size());
        Vector values(n);
        for (int i = 0; i < n; ++i) {
            auto j = left[i];
            values[i] = residual[j]*(1-alpha[i])+residual[j+1]*alpha[i];
        }
        const int padding = static_cast<int>(kernels.back().derivative.size())-1;
        Vector extended(n+2*padding, values.back());
        std::fill(extended.begin(), extended.begin()+padding, values.front());
        std::copy(values.begin(), values.end(), extended.begin()+padding);
        std::vector<Candidate> ridges;
        for (auto& kernel : kernels) {
            Vector derivative(n, 0), mexican(n, 0);
            int radius = static_cast<int>(kernel.derivative.size())-1;
            for (int i = 0; i < n; ++i) {
                for (int d = 1; d <= radius; ++d) {
                    double l = extended[padding+i-d], r = extended[padding+i+d];
                    derivative[i] += kernel.derivative[d]*(l-r);
                    mexican[i] += kernel.mexican[d]*(l+r-2*values[i]);
                }
            }
            double maximum = 0;
            for (double value : derivative) maximum = std::max(maximum, std::abs(value));
            if (maximum == 0) continue;
            std::vector<int> used;
            for (int i = 1; i+1 < n; ++i) {
                double sign = derivative[i] >= 0 ? 1 : -1;
                if (std::abs(derivative[i]) < cfg.ridge_threshold*maximum ||
                    sign*mexican[i] < 0 || sign*mexican[i+1] > 0 ||
                    mexican[i] == mexican[i+1]) continue;
                if (std::abs(derivative[i]) < cfg.snr_threshold*sigma*noise_norm(kernel, i)) continue;
                double fraction = mexican[i]/(mexican[i]-mexican[i+1]);
                double log_rate = -(x[i]+fraction*dx)-.5*kernel.scale*kernel.scale;
                if (!std::isfinite(log_rate)) continue;
                // Link from fine to coarse, retaining branches that merge later.
                int best = -1; double distance = .3;
                for (int k = 0; k < static_cast<int>(ridges.size()); ++k) {
                    double delta = std::abs(std::log(ridges[k].rate)-log_rate);
                    if (delta < distance &&
                        std::find(used.begin(), used.end(), k) == used.end()) {
                        best = k; distance = delta;
                    }
                }
                if (best < 0) {
                    ridges.push_back({std::exp(log_rate), std::abs(derivative[i]), 1});
                    used.push_back(static_cast<int>(ridges.size()-1));
                } else {
                    auto& ridge = ridges[best];
                    ++ridge.votes;
                    // Fine-scale localization has less smoothing bias.
                    ridge.strength = std::max(ridge.strength, std::abs(derivative[i]));
                    used.push_back(best);
                }
            }
        }
        ridges.erase(std::remove_if(ridges.begin(), ridges.end(), [&](const Candidate& r) {
            return r.votes < cfg.min_ridge_len;
        }), ridges.end());
        return ridges;
    }
};

Matrix orthonormal(const Matrix& e) {
    Matrix q;
    for (auto v : e) {
        for (int pass = 0; pass < 2; ++pass)
            for (const auto& column : q) {
                double projection = dot(v, column);
                for (std::size_t i = 0; i < v.size(); ++i) v[i] -= projection*column[i];
            }
        double norm = std::sqrt(dot(v, v));
        if (norm < 1e-12) continue;
        for (auto& value : v) value /= norm;
        q.push_back(std::move(v));
    }
    return q;
}

double score(const Fit& fit, int n, double floor, const DecomposeConfig& cfg) {
    double parameters = 2*fit.rates.size();
    if (cfg.model_selection == "bic")
        return n*std::log(std::max(fit.rss, floor)/n)+parameters*std::log(double(n));
    // Compare independently refitted nested models, not truncated amplitudes.
    return n*std::log(std::max(fit.rss, floor)/n)+2*parameters
           +2*parameters*(parameters+1)/std::max(1.0, n-parameters-1);
}

Fit prune(Fit current, const Vector& t, const Vector& y, double sigma,
          const DecomposeConfig& cfg, double floor) {
    if (!cfg.use_aic) return current;
    bool removed = true;
    while (removed && current.rates.size() > 1) {
        removed = false;
        for (std::size_t j = 0; j < current.rates.size(); ++j) {
            bool weak = std::abs(current.amplitudes[j]) <=
                        std::max(1e-10, cfg.snr_threshold*sigma);
            // Conditional amplitude uncertainty flags overlapping/cancelling
            // modes, including large amplitudes with nearly dependent columns.
            double fit_sigma = std::min(sigma, std::sqrt(current.rss/
                std::max(1.0, double(t.size())-2*current.rates.size())));
            if (fit_sigma > 1e-12) {
                Matrix others;
                for (std::size_t k = 0; k < current.rates.size(); ++k) {
                    if (k != j) others.push_back(current.e[k]);
                    Vector derivative(t.size());
                    for (std::size_t i = 0; i < t.size(); ++i)
                        derivative[i] = current.rates[k]*t[i]*current.e[k][i];
                    others.push_back(std::move(derivative));
                }
                auto q = orthonormal(others);
                auto innovation = current.e[j];
                for (int pass = 0; pass < 2; ++pass)
                    for (const auto& column : q) {
                        double product = dot(innovation, column);
                        for (std::size_t i = 0; i < t.size(); ++i)
                            innovation[i] -= product*column[i];
                    }
                double precision = std::sqrt(dot(innovation, innovation));
                if (std::abs(current.amplitudes[j])*precision <
                    std::max(cfg.snr_threshold, std::sqrt(2*std::log(double(t.size()))))*fit_sigma)
                    weak = true;
            }
            for (std::size_t k = 0; k < current.rates.size(); ++k)
                if (k != j && std::abs(std::log(current.rates[j]/current.rates[k])) < .1)
                    weak = true;
            if (!weak) continue;
            auto rates = current.rates; rates.erase(rates.begin()+j);
            auto sub = refine(t, y, std::move(rates), !cfg.disable_als);
            if (score(sub, static_cast<int>(t.size()), floor, cfg) <
                score(current, static_cast<int>(t.size()), floor, cfg)) {
                current = std::move(sub); removed = true; break;
            }
        }
    }
    return current;
}

std::vector<ExponentialComponent> fit_signal(const Vector& t, const Vector& f,
        const DecomposeConfig& cfg, int budget, double energy_tol, bool wavelets) {
    validate_decomposition(t, f, cfg);
    if (budget < 1 || budget > 32 || !std::isfinite(energy_tol) || energy_tol <= 0 || energy_tol >= 1)
        throw std::invalid_argument("Need 1 <= max_components <= 32 and 0 < energy_tol < 1");
    const int n = static_cast<int>(t.size());
    budget = std::min(budget, (n-1)/2);
    double amplitude_scale = 0;
    for (double value : f) amplitude_scale = std::max(amplitude_scale, std::abs(value));
    if (amplitude_scale == 0) return {};
    Vector y = f;
    for (auto& value : y) value /= amplitude_scale;
    const double energy = dot(y, y), floor = 1e-28*energy;
    const double sigma = cfg.noise_sigma >= 0 ? cfg.noise_sigma/amplitude_scale : estimate_noise(t, y);
    const double mean = std::accumulate(y.begin(), y.end(), 0.0)/n;
    double centered = 0, time_mean = std::accumulate(t.begin(), t.end(), 0.0)/n;
    double covariance = 0, time_variance = 0;
    for (int i = 0; i < n; ++i) {
        centered += (y[i]-mean)*(y[i]-mean);
        covariance += (t[i]-time_mean)*(y[i]-mean);
        time_variance += (t[i]-time_mean)*(t[i]-time_mean);
    }
    if (centered <= std::max(floor, (n+3*std::sqrt(2.0*n))*sigma*sigma) &&
        std::abs(covariance) <= std::max(1e-14, cfg.snr_threshold*sigma*std::sqrt(time_variance))) {
        if (std::abs(mean) <= std::max(cfg.snr_threshold, std::sqrt(2*std::log(double(n))))
                             *sigma/std::sqrt(double(n))) return {};
        ExponentialComponent constant;
        constant.amplitude = mean*amplitude_scale;
        constant.decay_time = constant.log_time = std::numeric_limits<double>::infinity();
        constant.significance = 1;
        return {constant};
    }
    // Exact single-mode shortcut, shared by wavelet and rate-grid routes.
    // A cheap log-linearity screen is followed by a fit on EVERY observation;
    // accept only at the common numerical floor, never at a noise tolerance.
    if (!cfg.disable_als && y.front() != 0) {
        int last = n-1;
        while (last > 0 && std::abs(y[last]) < 1e-4) --last;
        if (last > 0 && y.front()*y[last] > 0) {
            double rate = -std::log(std::abs(y[last]/y.front()))/(t[last]-t.front());
            bool exponential = rate > 0 && rate*t.front() <= 8;
            for (int j = 1; j < 16 && exponential; ++j) {
                int i = j*last/16;
                double prediction = y.front()*std::exp(-rate*(t[i]-t.front()));
                if (std::abs(prediction-y[i]) > 1e-10) exponential = false;
            }
            if (exponential) {
                auto single = refine(t, y, {rate}, true);
                if (single.rss <= floor) {
                    ExponentialComponent component;
                    component.rate = single.rates[0];
                    component.amplitude = single.amplitudes[0]*amplitude_scale;
                    component.decay_time = 1/component.rate;
                    component.log_time = std::log(component.decay_time);
                    component.significance = 1;
                    return {component};
                }
            }
        }
    }
    // Kernels and interpolation weights are built once, reused for every residual.
    std::unique_ptr<LogWavelets> detector;
    if (wavelets) detector = std::make_unique<LogWavelets>(t, cfg);
    Fit current; current.rss = energy; current.residual = y;
    if (wavelets) {
        // Initialize ALL persistent branches together. Sequentially fitting
        // one or two exponentials to a three-mode signal wastes nonlinear solves.
        auto ridges = detector->candidates(y, sigma);
        std::sort(ridges.begin(), ridges.end(), [](const auto& a, const auto& b) {
            return a.votes == b.votes ? a.strength > b.strength : a.votes > b.votes;
        });
        Vector initial;
        for (const auto& ridge : ridges) {
            if (ridge.votes < std::max(cfg.min_ridge_len, cfg.n_scales/2)) continue;
            bool distinct = true;
            for (double rate : initial)
                if (std::abs(std::log(rate/ridge.rate)) < .3) distinct = false;
            if (distinct) initial.push_back(ridge.rate);
            if (static_cast<int>(initial.size()) == budget) break;
        }
        if (!initial.empty())
            current = prune(refine(t, y, std::move(initial), !cfg.disable_als),
                            t, y, sigma, cfg, floor);
    }
    for (int mode = static_cast<int>(current.rates.size()); mode < budget; ++mode) {
        if (current.rss <= floor) break;
        auto proposals = wavelets ? detector->candidates(current.residual, sigma)
                                 : std::vector<Candidate>{};
        if (!wavelets) {
            for (int j = 0; j < 64; ++j)
                proposals.push_back({std::exp(std::log(.005/t.back())
                    +j*std::log(1600*t.back()/t.front())/63), 0, 1});
        }
        Matrix tangent = current.e;
        for (std::size_t j = 0; j < current.rates.size(); ++j) {
            Vector derivative(n);
            for (int i = 0; i < n; ++i)
                derivative[i] = t[i]*current.rates[j]*current.e[j][i];
            tangent.push_back(std::move(derivative));
        }
        // Condition the new mode on BOTH existing amplitudes and rates.
        // This estimates joint-fit improvement instead of a fixed-rate addition.
        const auto q = orthonormal(tangent);
        auto conditional_strength = [&](double rate) {
            Vector v(n);
            for (int i = 0; i < n; ++i) v[i] = std::exp(-rate*t[i]);
            for (int pass = 0; pass < 2; ++pass)
                for (const auto& column : q) {
                    double projection = dot(v, column);
                    for (int i = 0; i < n; ++i) v[i] -= projection*column[i];
                }
            double norm = dot(v, v), product = dot(v, current.residual);
            return norm > 1e-20 ? product*product/norm : 0;
        };
        for (auto& proposal : proposals) proposal.strength = conditional_strength(proposal.rate);
        // A merged ridge can leave a broad residual lobe on the wrong side of
        // the fitted rate. Search both log-rate sides of that ridge's span,
        // using the existing Jacobian projection rather than nonlinear fits.
        // Keep the original proposals as alternatives; no count is prescribed.
        if (wavelets && current.rates.size() == 1 && !cfg.disable_als && !proposals.empty()) {
            const auto strongest = std::max_element(proposals.begin(), proposals.end(),
                [](const Candidate& a, const Candidate& b) { return a.strength < b.strength; });
            double center = std::log(current.rates[0]);
            double span = std::min(2.0, std::abs(std::log(strongest->rate)-center)+.3);
            for (int side : {-1, 1}) {
                double low = std::max(std::log(.005/t.back()), center+(side < 0 ? -span : .03));
                double high = std::min(std::log(8/t.front()), center+(side < 0 ? -.03 : span));
                if (low >= high) continue;
                constexpr double fraction = .6180339887498948482;
                double a = high-fraction*(high-low), b = low+fraction*(high-low);
                double fa = conditional_strength(std::exp(a)), fb = conditional_strength(std::exp(b));
                for (int iteration = 0; iteration < 18; ++iteration) {
                    if (fa < fb) {
                        low = a; a = b; fa = fb;
                        b = low+fraction*(high-low); fb = conditional_strength(std::exp(b));
                    } else {
                        high = b; b = a; fb = fa;
                        a = high-fraction*(high-low); fa = conditional_strength(std::exp(a));
                    }
                }
                proposals.push_back({std::exp(fa > fb ? a : b), std::max(fa, fb), 1});
            }
        }
        std::sort(proposals.begin(), proposals.end(), [](const Candidate& a, const Candidate& b) {
            return a.strength > b.strength;
        });
        Fit best = current;
        Vector tried;
        for (const auto& proposal : proposals) {
            bool duplicate = false;
            for (double rate : current.rates)
                if (std::abs(std::log(proposal.rate/rate)) < .015) duplicate = true;
            for (double rate : tried)
                if (std::abs(std::log(proposal.rate/rate)) < .15) duplicate = true;
            if (duplicate) continue;
            if (proposal.strength < std::max(cfg.snr_threshold*cfg.snr_threshold,
                                            2*std::log(double(n)))*sigma*sigma) break;
            Vector rates = current.rates; rates.push_back(proposal.rate);
            auto candidate = refine(t, y, std::move(rates), !cfg.disable_als);
            bool identifiable = true;
            for (std::size_t j = 0; j < candidate.rates.size(); ++j) {
                if (std::abs(candidate.amplitudes[j]) > 1e4) identifiable = false;
                for (std::size_t k = j+1; k < candidate.rates.size(); ++k)
                    if (std::abs(std::log(candidate.rates[j]/candidate.rates[k])) < 1e-5)
                        identifiable = false;
            }
            if (identifiable && candidate.rss < best.rss) best = std::move(candidate);
            tried.push_back(proposal.rate);
            if (best.rss <= std::max(floor, 1.3*n*sigma*sigma) || tried.size() >= 3) break;
        }
        if (best.rates.size() == current.rates.size()) break;
        if (mode > 0 && cfg.use_aic && score(best, n, floor, cfg) >= score(current, n, floor, cfg)) break;
        current = std::move(best);
        if (!cfg.use_aic && current.rss < energy_tol*energy) break;
    }
    // Backward refits remove redundant modes jointly; never truncate a fitted sum.
    current = prune(std::move(current), t, y, sigma, cfg, floor);
    std::vector<ExponentialComponent> out;
    double total = 0;
    for (double amplitude : current.amplitudes) total += std::abs(amplitude);
    for (std::size_t j = 0; j < current.rates.size(); ++j) {
        ExponentialComponent component;
        component.rate = current.rates[j]; component.decay_time = 1/component.rate;
        component.log_time = std::log(component.decay_time);
        component.amplitude = current.amplitudes[j]*amplitude_scale;
        component.significance = total > 0 ? std::abs(current.amplitudes[j])/total : 0;
        out.push_back(component);
    }
    std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) {
        return a.significance > b.significance;
    });
    return out;
}
} // namespace

std::vector<ExponentialComponent> decompose(const Vector& t, const Vector& f,
        const DecomposeConfig& cfg, int max_components) {
    return fit_signal(t, f, cfg, max_components, 1e-6, true);
}

std::vector<ExponentialComponent> decompose_deflate(const Vector& t, const Vector& f,
        const DecomposeConfig& cfg, int max_components, double energy_tol) {
    return fit_signal(t, f, cfg, max_components, energy_tol, false);
}

std::vector<double> reconstruct(
    const std::vector<ExponentialComponent>& components,
    const std::vector<double>& t_eval)
{
    std::vector<double> out(t_eval.size(), 0.0);
    for (std::size_t i = 0; i < t_eval.size(); ++i) {
        double tt = t_eval[i];
        for (const auto& c : components) {
            double val = c.amplitude * std::exp(-c.rate * tt);
            if (c.frequency != 0.0)
                val *= std::cos(c.frequency * tt + c.phase);
            out[i] += val;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// reconstruct_laplace
// ---------------------------------------------------------------------------
std::vector<std::complex<double>> reconstruct_laplace(
    const std::vector<ExponentialComponent>& components,
    const std::vector<std::complex<double>>& s_eval)
{
    std::vector<std::complex<double>> out(s_eval.size(), 0.0);
    for (std::size_t i = 0; i < s_eval.size(); ++i) {
        std::complex<double> s = s_eval[i];
        for (const auto& c : components) {
            auto shifted = s + c.rate;
            if (c.frequency == 0.0) {
                if (std::abs(shifted) == 0)
                    throw std::domain_error("Laplace image evaluated at a component pole");
                out[i] += c.amplitude / shifted;
            } else {
                // L{A exp(-rate*t) cos(omega*t+phase)} includes BOTH poles.
                auto denom = shifted * shifted + c.frequency * c.frequency;
                if (std::abs(denom) == 0)
                    throw std::domain_error("Laplace image evaluated at a component pole");
                out[i] += c.amplitude * (shifted * std::cos(c.phase)
                            - c.frequency * std::sin(c.phase)) / denom;
            }
        }
    }
    return out;
}

} // namespace laplacex
