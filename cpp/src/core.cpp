#include "laplacex/core.hpp"

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <thread>
#include <future>
#include <numeric>
#include <limits>
#include <vector>    // already included via core.hpp, but for clarity

namespace laplacex {

// ---------------------------------------------------------------------------
// Gauss-Kronrod G7-K15 nodes and weights (standard interval [-1, 1])
// ---------------------------------------------------------------------------
static const double GK_NODES[15] = {
     0.0000000000000000,
     0.2077849550078985, -0.2077849550078985,
     0.4058451513773972, -0.4058451513773972,
     0.5860872354676911, -0.5860872354676911,
     0.7415311855993945, -0.7415311855993945,
     0.8648644233597691, -0.8648644233597691,
     0.9491079123427585, -0.9491079123427585,
     0.9914553711208126, -0.9914553711208126
};

static const double GK_W[15] = {
     0.2094821410847278,
     0.2044329400752989,  0.2044329400752989,
     0.1903505780647854,  0.1903505780647854,
     0.1690047266392679,  0.1690047266392679,
     0.1406532597155259,  0.1406532597155259,
     0.1047900103222502,  0.1047900103222502,
     0.0630920926299786,  0.0630920926299786,
     0.0229353220105292,  0.0229353220105292
};

static const double G7_W[15] = {
     0.4179591836734694,
     0.0,                 0.0,
     0.3818300505051189,  0.3818300505051189,
     0.0,                 0.0,
     0.2797053914892767,  0.2797053914892767,
     0.0,                 0.0,
     0.1294849661688697,  0.1294849661688697,
     0.0,                 0.0
};

static void validate_config(const TransformConfig& cfg) {
    if (!std::isfinite(cfg.tol) || cfg.tol <= 0 ||
        !std::isfinite(cfg.rel_tol) || cfg.rel_tol < 0 ||
        !std::isfinite(cfg.t_max) || cfg.t_max <= 0 ||
        cfg.max_levels < 0 || cfg.max_levels > 50)
        throw std::invalid_argument("Invalid quadrature tolerance, depth or t_max");
}

static void validate_frequency(std::complex<double> s) {
    if (!std::isfinite(s.real()) || !std::isfinite(s.imag()) || s.real() < 0)
        throw std::domain_error("Forward transform requires finite s with Re(s) >= 0");
}

// Evaluate integrand in [a, b] using G7-K15; return (kronrod, gauss) pair
static std::pair<std::complex<double>, std::complex<double>>
gk15(const std::function<std::complex<double>(double)>& g, double a, double b)
{
    double mid  = 0.5 * (a + b);
    double half = 0.5 * (b - a);

    std::complex<double> kron{}, gauss{};
    for (int i = 0; i < 15; ++i) {
        double x  = mid + half * GK_NODES[i];
        auto   fx = g(x);
        if (!std::isfinite(fx.real()) || !std::isfinite(fx.imag()))
            throw std::domain_error("Non-finite quadrature integrand");
        kron  += GK_W[i] * fx;
        gauss += G7_W[i] * fx;
    }
    kron  *= half;
    gauss *= half;
    return {kron, gauss};
}

// Adaptive recursive integration of g over [a, b]
static std::complex<double> adaptive_gk(
    const std::function<std::complex<double>(double)>& g,
    double a, double b,
    double tol, double rel_tol, int depth, int max_depth)
{
    auto [kron, gauss] = gk15(g, a, b);
    double err = std::abs(kron - gauss);
    double ref = std::abs(kron);

    const double target = std::max(tol, rel_tol * ref);
    // Do not bisect further when the embedded estimate reaches roundoff.
    if (err <= target || err <= 32 * std::numeric_limits<double>::epsilon() * ref) {
        return kron;
    }
    if (depth >= max_depth)
        throw std::runtime_error("Quadrature did not converge; increase max_levels or relax tol");

    double mid = 0.5 * (a + b);
    double t2  = tol * 0.5;
    return adaptive_gk(g, a, mid, t2, rel_tol, depth + 1, max_depth)
         + adaptive_gk(g, mid, b, t2, rel_tol, depth + 1, max_depth);
}

// ---------------------------------------------------------------------------
// Forward transform
// ---------------------------------------------------------------------------
std::complex<double> forward(
    const std::function<double(double)>& f,
    std::complex<double> s,
    const TransformConfig& cfg)
{
    validate_config(cfg);
    validate_frequency(s);
    double s_abs  = std::abs(s);
    double center = std::min((s_abs > 0) ? 1.0 / s_abs : 1.0, cfg.t_max);

    std::vector<double> breaks;
    breaks.push_back(0.0);
    for (double x = std::min(1e-6, center * 0.25); x > 0 && x < center; x *= 4.0)
        breaks.push_back(x);
    breaks.push_back(center);
    for (double x = center * 4.0; x < cfg.t_max; x *= 4.0) breaks.push_back(x);
    breaks.push_back(cfg.t_max);
    std::sort(breaks.begin(), breaks.end());
    breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());

    std::complex<double> result{}, correction{};
    const double local_tol = cfg.tol / (breaks.size() - 1);
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        double a = breaks[i], b = breaks[i + 1];
        auto g = [&](double t) -> std::complex<double> {
            return f(t) * std::exp(-s * t);
        };
        std::complex<double> value;
        if (a == 0.0) {
            // t=u^2 removes the common integrable 1/sqrt(t) endpoint singularity.
            auto endpoint = [&](double u) { return g(u * u) * (2 * u); };
            value = adaptive_gk(endpoint, 0, std::sqrt(b), local_tol,
                                cfg.rel_tol, 0, cfg.max_levels);
        } else {
            value = adaptive_gk(g, a, b, local_tol, cfg.rel_tol, 0, cfg.max_levels);
        }
        auto adjusted = value - correction;
        auto next = result + adjusted;
        correction = (next - result) - adjusted;
        result = next;
    }
    return result;
}

// ---------------------------------------------------------------------------
// Batch forward (parallelised over s values) – used only from native C++,
// Python interface uses sequential loop to avoid GIL issues.
// ---------------------------------------------------------------------------
std::vector<std::complex<double>> forward_batch(
    const std::function<double(double)>& f,
    const std::vector<std::complex<double>>& s_vals,
    const TransformConfig& cfg)
{
    std::size_t n = s_vals.size();
    std::vector<std::complex<double>> out(n);

    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::future<void>> futs;

    auto worker = [&](std::size_t from, std::size_t to) {
        for (std::size_t i = from; i < to; ++i)
            out[i] = forward(f, s_vals[i], cfg);
    };

    std::size_t chunk = (n + nthreads - 1) / nthreads;
    for (unsigned t = 0; t < nthreads; ++t) {
        std::size_t from = t * chunk;
        std::size_t to   = std::min(from + chunk, n);
        if (from < to)
            futs.push_back(std::async(std::launch::async, worker, from, to));
    }
    for (auto& fut : futs) fut.get();
    return out;
}

// ---------------------------------------------------------------------------
// Cubic-spline interpolation (natural BC)
// ---------------------------------------------------------------------------
struct CubicSpline {
    std::vector<double> t, a, b, c, d;

    CubicSpline(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size())
            throw std::invalid_argument("t and f must have equal lengths");
        int n = static_cast<int>(x.size()) - 1;
        if (n < 1) throw std::invalid_argument("Need at least 2 points");
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || x[i] < 0 ||
                (i && x[i] <= x[i-1]))
                throw std::invalid_argument("Need finite samples with strictly increasing t >= 0");
        }

        t = x; a = y;
        std::vector<double> h(n), alpha(n), l(n+1), mu(n+1), z(n+1);
        b.resize(n); c.resize(n+1); d.resize(n);

        for (int i = 0; i < n; ++i) h[i] = x[i+1] - x[i];
        for (int i = 1; i < n; ++i)
            alpha[i] = 3*(y[i+1]-y[i])/h[i] - 3*(y[i]-y[i-1])/h[i-1];

        l[0] = 1; mu[0] = z[0] = 0;
        for (int i = 1; i < n; ++i) {
            l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i];
        }
        l[n] = 1; z[n] = c[n] = 0;
        for (int j = n-1; j >= 0; --j) {
            c[j] = z[j] - mu[j]*c[j+1];
            b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3;
            d[j] = (c[j+1]-c[j]) / (3*h[j]);
        }
    }

    double eval(double x) const {
        int n = static_cast<int>(t.size()) - 1;
        if (x > t.back()) return 0.0;
        int lo = 0, hi = n - 1;
        while (lo < hi) { int mid = (lo+hi)/2; if (t[mid+1] < x) lo=mid+1; else hi=mid; }
        double dx = x - t[lo];
        return a[lo] + b[lo]*dx + c[lo]*dx*dx + d[lo]*dx*dx*dx;
    }
};

// ---------------------------------------------------------------------------
// Forward from discrete samples
// ---------------------------------------------------------------------------
std::vector<std::complex<double>> forward_discrete(
    const std::vector<double>& t,
    const std::vector<double>& f,
    const std::vector<std::complex<double>>& s_vals,
    const TransformConfig& cfg)
{
    validate_config(cfg);
    CubicSpline spline(t, f);
    std::vector<std::complex<double>> out;
    out.reserve(s_vals.size());
    // Respect every spline knot: global adaptive sampling can miss narrow features.
    for (auto s : s_vals) {
        validate_frequency(s);
        auto g = [&](double x) { return spline.eval(x) * std::exp(-s * x); };
        std::complex<double> sum{}, correction{};
        double left = 0.0;
        for (double knot : t) {
            const double right = std::min(knot, cfg.t_max);
            if (right > left) {
                auto value = adaptive_gk(g, left, right, cfg.tol / t.size(),
                                         cfg.rel_tol, 0, cfg.max_levels);
                auto adjusted = value - correction;
                auto next = sum + adjusted;
                correction = (next - sum) - adjusted;
                sum = next;
            }
            left = right;
            if (knot >= cfg.t_max) break;
        }
        out.push_back(sum);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Inverse: Weeks-Laguerre expansion, with midpoint Fourier coefficients.
// See Abate, Choudhury & Whitt (1996), DOI 10.1287/ijoc.8.4.413.
// ---------------------------------------------------------------------------
std::vector<double> inverse_weeks(
    const std::function<std::complex<double>(std::complex<double>)>& F,
    const std::vector<double>& t_eval,
    double sigma, double b, int N,
    const TransformConfig& cfg)
{
    validate_config(cfg);
    if (N < 2 || N > 4096 || !std::isfinite(sigma) ||
        !std::isfinite(b) || b < 0)
        throw std::invalid_argument("Weeks requires 2 <= N <= 4096, finite sigma, and b >= 0");
    for (double t : t_eval)
        if (!std::isfinite(t) || t < 0)
            throw std::invalid_argument("Inverse evaluation times must be finite and nonnegative");
    if (t_eval.empty()) return {};
    const double tmax = std::max(1e-12, *std::max_element(t_eval.begin(), t_eval.end()));
    if (sigma == 0) sigma = 4.0 / tmax;
    const int samples = 4 * N;
    const double pi = std::acos(-1.0);
    std::vector<double> coefficients;
    double best_score = std::numeric_limits<double>::infinity();
    double best_b = b;
    const bool auto_b = b == 0;
    // Tail minimisation chooses a useful scale for both decays and oscillations.
    // Coefficients are computed once per candidate and shared by all times.
    for (int candidate = 0; candidate < (auto_b ? 17 : 1); ++candidate) {
        double scale = auto_b ? std::pow(2.0, -1.0 + candidate * 0.5) / tmax : b;
        std::vector<long double> sums(N, 0.0L);
        for (int j = 0; j < samples; ++j) {
            double theta = 2 * pi * (j + 0.5) / samples;
            std::complex<double> z = std::polar(1.0, theta);
            auto one_minus_z = 1.0 - z;
            // Cotangent form keeps the contour's real part exactly sigma.
            auto value = F({sigma, scale / std::tan(theta * 0.5)});
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                throw std::domain_error("Non-finite Laplace image on the Weeks contour");
            auto g = 2.0 * scale / one_minus_z * value;
            std::complex<double> phase(1.0, 0.0);
            auto step = std::conj(z);
            for (int k = 0; k < N; ++k) {
                sums[k] += (g * phase).real();
                phase *= step;
            }
        }
        std::vector<double> a(N);
        double tail = 0, norm = 0;
        for (int k = 0; k < N; ++k) {
            a[k] = static_cast<double>(sums[k] / samples);
            norm += std::abs(a[k]);
            if (k >= 3 * N / 4) tail += std::abs(a[k]);
        }
        // Relative tail size prevents an under-resolved, low-amplitude expansion
        // from winning solely because it misses most of an oscillatory signal.
        const double roundoff = 32 * std::numeric_limits<double>::epsilon();
        double score = norm > 0 ? tail / norm + roundoff : 0;
        if (score < best_score) {
            best_score = score;
            coefficients = std::move(a);
            best_b = scale;
        }
        if (tail + roundoff * norm <
            cfg.tol * std::exp(-std::max(0.0, sigma * tmax)) * 0.25) break;
    }
    std::vector<double> out(t_eval.size());
    for (std::size_t i = 0; i < t_eval.size(); ++i) {
        // Scaled Laguerre recurrence avoids separately overflowing L_n(2bt).
        const long double x = 2.0L * best_b * t_eval[i];
        long double prev = std::exp((sigma - best_b) * t_eval[i]);
        long double curr = (1.0L - x) * prev;
        long double sum = coefficients[0] * prev + coefficients[1] * curr;
        for (int k = 1; k < N - 1; ++k) {
            long double next = ((2 * k + 1 - x) * curr - k * prev) / (k + 1);
            sum += coefficients[k+1] * next;
            prev = curr;
            curr = next;
        }
        out[i] = static_cast<double>(sum);
        if (!std::isfinite(out[i]))
            throw std::overflow_error("Weeks reconstruction overflow; adjust sigma and b");
    }
    return out;
}

} // namespace laplacex
