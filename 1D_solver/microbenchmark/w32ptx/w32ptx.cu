// w32ptx.cu -- driver for the extracted w32_poly64 modes.
//
// CLI contract is identical to ../weno_micro.f90's driver so the same sweep
// plumbing parses it:   ./w32ptx <mode> [nx] [nrepeat] [nlaunch]
// Extra long options after the positionals: --dump FILE, --list-kernels,
// --accuracy.
//
// Defaults: nx = 4194304, nrepeat = 1, nlaunch = 10.
//
// nrepeat stays 1 for every comparable row. In the Fortran harness nvfortran
// hoisted the inner loop for register-only modes but not for shared-memory ones,
// which silently made rows incomparable; that is why run_weno_nsys.sh rejects
// --nrepeat outright and why the loop here carries `#pragma unroll 1`.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <string>

#include "kernels.cuh"

// ---------------------------------------------------------------------------
// Emit one kernel per mode, then one launch wrapper per mode.
// ---------------------------------------------------------------------------
#define MODE(name, R, Rc, Body)                                                       \
    extern "C" __global__ __launch_bounds__(128) void k_##name(                       \
        int n, int nrepeat, const double* __restrict__ x, double* __restrict__ out)    \
    {                                                                                 \
        w32::kernel_body<R, w32::Rc, w32::Body>(n, nrepeat, x, out);                   \
    }
#include "modes.inc"
#undef MODE

#define MODE(name, R, Rc, Body)                                                       \
    static void launch_##name(int n, int nrepeat, const double* x, double* out,        \
                              dim3 grid, dim3 block)                                  \
    {                                                                                 \
        k_##name<<<grid, block>>>(n, nrepeat, x, out);                                 \
    }
#include "modes.inc"
#undef MODE

struct ModeEntry {
    const char* name;
    const char* kernel;
    int R;
    void (*launch)(int, int, const double*, double*, dim3, dim3);
};

static const ModeEntry kModes[] = {
#define MODE(name, R, Rc, Body) {#name, "k_" #name, R, &launch_##name},
#include "modes.inc"
#undef MODE
};
static const int kNumModes = (int)(sizeof(kModes) / sizeof(kModes[0]));

// ---------------------------------------------------------------------------
// Error handling. cudaGetLastError() FIRST, then cudaDeviceSynchronize():
// sync alone returns 0 when the launch itself was rejected, which is how a
// wrong-architecture build once reported an untouched buffer as checksum=0.0.
// ---------------------------------------------------------------------------
static void check_launch(const char* what)
{
    cudaError_t e = cudaGetLastError();
    if (e != cudaSuccess) {
        std::fprintf(stderr, "FAIL: launch error after %s: %s (%d)\n", what,
                     cudaGetErrorString(e), (int)e);
        std::exit(3);
    }
    e = cudaDeviceSynchronize();
    if (e != cudaSuccess) {
        std::fprintf(stderr, "FAIL: sync error after %s: %s (%d)\n", what,
                     cudaGetErrorString(e), (int)e);
        std::exit(3);
    }
}

#define CUDA_OK(call)                                                                 \
    do {                                                                              \
        cudaError_t _e = (call);                                                       \
        if (_e != cudaSuccess) {                                                       \
            std::fprintf(stderr, "FAIL: %s -> %s (%d)\n", #call,                       \
                         cudaGetErrorString(_e), (int)_e);                             \
            std::exit(3);                                                              \
        }                                                                             \
    } while (0)

/// The architecture this binary actually contains device code for. A stale
/// CASE_GPU_CC/ARCH builds cleanly and then fails at every launch with
/// cudaErrorInvalidPtx (218) because LTO device code cannot be JIT'd across
/// architectures. Refusing up front beats diagnosing it at the checksum.
static int compiled_arch()
{
#if defined(W32PTX_ARCH)
    return W32PTX_ARCH;
#else
    return 0;
#endif
}

static void arch_guard()
{
    int dev = 0;
    CUDA_OK(cudaGetDevice(&dev));
    cudaDeviceProp prop{};
    CUDA_OK(cudaGetDeviceProperties(&prop, dev));
    const int runtime_cc = prop.major * 10 + prop.minor;
    const int built = compiled_arch();
    std::printf("gpu_cc_compiled=%d gpu_cc_runtime=%d\n", built, runtime_cc);
    if (built != 0 && built != runtime_cc) {
        std::fprintf(stderr,
                     "FAIL: binary built for sm_%d but device is sm_%d (%s).\n"
                     "      Rebuild with: make ARCH=%d\n",
                     built, runtime_cc, prop.name, runtime_cc);
        std::exit(3);
    }
}

static const ModeEntry* find_mode(const char* name)
{
    for (int i = 0; i < kNumModes; ++i)
        if (std::strcmp(kModes[i].name, name) == 0) return &kModes[i];
    return nullptr;
}

// ---------------------------------------------------------------------------
// Accuracy gate.
//
// checksum_all provably cannot see a weight error: init_input is piecewise
// linear and every WENO candidate is exact on linear data, so the documented
// order-destroying `dd` mirror bug moved the whole-array sum by 3e-11 out of
// 1.27e6 and individual values by <= 1 ulp. And ../../report/check_weno_order.py
// --w32-split, while it gates the algebra, models the division with numpy
// float32 IEEE arithmetic and therefore cannot see rcp.approx at all.
//
// So the gate has to run the REAL device code. Two tests, both biases scored
// separately -- the historical bug was right-bias-only.
// ---------------------------------------------------------------------------

static void run_mode_on(const ModeEntry* m, int n, const std::vector<double>& hx,
                        std::vector<double>& hout)
{
    double *x = nullptr, *out = nullptr;
    CUDA_OK(cudaMalloc(&x, sizeof(double) * (size_t)n * 3));
    CUDA_OK(cudaMalloc(&out, sizeof(double) * (size_t)n * 6));
    CUDA_OK(cudaMemcpy(x, hx.data(), sizeof(double) * (size_t)n * 3, cudaMemcpyHostToDevice));

    dim3 block(128), grid((unsigned)((n + 127) / 128));
    w32::poison<<<dim3((unsigned)(((size_t)n * 6 + 127) / 128)), block>>>(out, (size_t)n * 6, -1.0);
    check_launch("poison");
    m->launch(n, 1, x, out, grid, block);
    check_launch(m->name);

    hout.resize((size_t)n * 6);
    CUDA_OK(cudaMemcpy(hout.data(), out, sizeof(double) * (size_t)n * 6, cudaMemcpyDeviceToHost));
    CUDA_OK(cudaFree(x));
    CUDA_OK(cudaFree(out));
}

/// (a) Max relative deviation on a smooth non-linear field. Same field as
/// ../check_df_accuracy.sh so the two gates are commensurable.
///
/// Gated RELATIVE TO _ref, not against FP64. Measuring an FP32-weights split
/// against FP64 does not answer the question we have: on a smooth field the
/// dense beta quadratic forms cancel catastrophically (terms of order 1e6
/// producing a beta of order 1e-8, a relative cancellation far beyond FP32's
/// 1.2e-7 precision), so the FP32 weights carry ~1e-7 relative error by
/// construction. That is why ../../report/check_weno_order.py carries
/// W32_FLOOR = 3e-7. Measured here: _ref itself sits at 9e-08 .. 3e-07.
///
/// So the FP32 split's own error is the baseline, and the question for a rung is
/// only whether its reciprocal made things materially worse. Threshold 1.5x.
static bool accuracy_deviation()
{
    const int n = 4096;
    std::vector<double> hx((size_t)n * 3);
    for (int i = 0; i < n; ++i) {
        const double v = std::exp(0.3 * std::sin(0.017 * i))
                       + 0.5 * std::cos(0.0031 * std::pow((double)i, 1.3));
        hx[(size_t)0 * n + i] = v;
        hx[(size_t)1 * n + i] = v;
        hx[(size_t)2 * n + i] = v;
    }

    std::printf("\n-- (a) max relative deviation vs weno64, smooth field, n=%d --\n", n);
    std::printf("   _ref is the FP32-split baseline (informational; ~1e-7 is the\n");
    std::printf("   FP32 weight-noise floor). Other rungs gated at <= 1.5x _ref.\n");
    std::printf("%-30s %14s %14s %9s %9s   %s\n", "mode", "left", "right",
                "vs_ref_L", "vs_ref_R", "verdict");

    bool ok = true;
    for (int R = 3; R <= 5; ++R) {
        char refname[64];
        std::snprintf(refname, sizeof(refname), "weno64_seq%d", 2 * R - 1);
        const ModeEntry* ref64 = find_mode(refname);
        if (!ref64) continue;
        std::vector<double> href;
        run_mode_on(ref64, n, hx, href);

        const int valid = n - (2 * R - 1);
        double refmax = 0.0;
        for (int c = 0; c < 6; ++c)
            for (int i = 0; i < valid; ++i) refmax = fmax(refmax, fabs(href[(size_t)c * n + i]));

        // Baseline first: the _ref rung at this width.
        double base_l = 0.0, base_r = 0.0;
        bool have_base = false;
        for (int pass_i = 0; pass_i < 2; ++pass_i) {
            for (int mi = 0; mi < kNumModes; ++mi) {
                const ModeEntry* m = &kModes[mi];
                if (m->R != R) continue;
                if (std::strncmp(m->name, "w32_poly64_", 11) != 0) continue;
                const bool is_ref = (std::strstr(m->name, "_ref") != nullptr);
                if ((pass_i == 0) != is_ref) continue;

                std::vector<double> hv;
                run_mode_on(m, n, hx, hv);
                double dl = 0.0, dr = 0.0;
                for (int c = 0; c < 6; ++c) {
                    double& d = (c % 2 == 0) ? dl : dr;
                    for (int i = 0; i < valid; ++i)
                        d = fmax(d, fabs(hv[(size_t)c * n + i] - href[(size_t)c * n + i]));
                }
                const double rl = dl / refmax, rr = dr / refmax;

                if (is_ref) {
                    base_l = rl; base_r = rr; have_base = true;
                    std::printf("%-30s %14.3e %14.3e %9s %9s   %s\n", m->name, rl, rr,
                                "-", "-", "baseline");
                } else {
                    const double xl = have_base && base_l > 0 ? rl / base_l : 0.0;
                    const double xr = have_base && base_r > 0 ? rr / base_r : 0.0;
                    const bool good = (xl <= 1.5) && (xr <= 1.5);
                    ok = ok && good;
                    std::printf("%-30s %14.3e %14.3e %8.2fx %8.2fx   %s\n", m->name, rl, rr,
                                xl, xr, good ? "PASS" : "FAIL");
                }
            }
        }
    }
    return ok;
}

/// (b) Order test on the real device code. This is the gate that catches a
/// too-coarse reciprocal AND the mirrored-dd class of bug, neither of which any
/// checksum can see. Replicates micro_split_w32_test in
/// ../../report/check_weno_order.py:465.
///
/// Face geometry: the left bias uses cells i..i+2R-2 (centre i+R-1) and gives
/// v- at the RIGHT edge of cell i+R-1; the right bias uses cells i+1..i+2R-1
/// (centre i+R) and gives v+ at the LEFT edge of cell i+R. Same face,
/// x = (i+R)*h. Ghost cells are periodic (sin has period 2pi), which both
/// reproduces numpy's np.roll and keeps all N faces inside the kernel's
/// `i >= n-(2R-1)` guard.
static bool accuracy_order()
{
    const int Ns[] = {16, 32, 64, 128, 256};
    const int nN = 5;

    std::printf("\n-- (b) order on device code, u(x)=sin(x), both biases --\n");
    std::printf("%-30s %5s %12s %12s %7s %7s   %s\n", "mode", "bias", "L1(min)", "L1(max N)",
                "rate", "need", "verdict");

    bool ok = true;
    for (int mi = 0; mi < kNumModes; ++mi) {
        const ModeEntry* m = &kModes[mi];
        if (std::strncmp(m->name, "w32_poly64_", 11) != 0) continue;
        const int R = m->R;

        double l1[2][nN];
        for (int k = 0; k < nN; ++k) {
            const int N = Ns[k];
            const double h = 2.0 * M_PI / N;
            const int n = N + 2 * R;
            std::vector<double> hx((size_t)n * 3);
            for (int j = 0; j < n; ++j) {
                const int jp = j % N;  // periodic ghost cells
                const double vb = (std::cos(jp * h) - std::cos((jp + 1) * h)) / h;
                hx[(size_t)0 * n + j] = vb;
                hx[(size_t)1 * n + j] = vb;
                hx[(size_t)2 * n + j] = vb;
            }
            std::vector<double> hv;
            run_mode_on(m, n, hx, hv);

            double sl = 0.0, sr = 0.0;
            for (int i = 0; i < N; ++i) {
                const double exact = std::sin((i + R) * h);
                sl += fabs(hv[(size_t)0 * n + i] - exact);
                sr += fabs(hv[(size_t)1 * n + i] - exact);
            }
            l1[0][k] = sl / N;
            l1[1][k] = sr / N;
        }

        for (int b = 0; b < 2; ++b) {
            double best_rate = 0.0, lmin = l1[b][0];
            for (int k = 0; k < nN; ++k) lmin = fmin(lmin, l1[b][k]);
            for (int k = 1; k < nN; ++k) {
                if (l1[b][k] > 3e-7 && l1[b][k - 1] > 3e-7) {
                    const double r = std::log2(l1[b][k - 1] / l1[b][k]);
                    best_rate = fmax(best_rate, r);
                }
            }
            const double need = 2.0 * R - 1.5;
            const bool pass = (best_rate > need) || (lmin <= 5e-7);
            ok = ok && pass;
            std::printf("%-30s %5s %12.4e %12.4e %7.2f %7.2f   %s\n", m->name,
                        b == 0 ? "left" : "right", lmin, l1[b][nN - 1], best_rate, need,
                        pass ? "PASS" : "FAIL");
        }
    }
    return ok;
}

static int run_accuracy()
{
    const bool a = accuracy_deviation();
    const bool b = accuracy_order();
    std::printf("\n-> %s\n", (a && b) ? "PASS" : "FAIL");
    return (a && b) ? 0 : 1;
}

// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
    std::string mode_name = "w32_poly64_seq9_ref";
    int n = 4194304, nrepeat = 1, nlaunch = 10;
    const char* dump_path = nullptr;

    // Positionals, then long options.
    int pos = 0;
    for (int a = 1; a < argc; ++a) {
        if (std::strcmp(argv[a], "--list-kernels") == 0) {
            for (int i = 0; i < kNumModes; ++i)
                std::printf("%s\t%s\n", kModes[i].name, kModes[i].kernel);
            return 0;
        }
        if (std::strcmp(argv[a], "--dump") == 0 && a + 1 < argc) { dump_path = argv[++a]; continue; }
        if (std::strcmp(argv[a], "--accuracy") == 0 || std::strcmp(argv[a], "accuracy") == 0) {
            arch_guard();
            return run_accuracy();
        }
        if (argv[a][0] == '-') {
            std::fprintf(stderr, "FAIL: unknown option %s\n", argv[a]);
            return 2;
        }
        switch (pos++) {
            case 0: mode_name = argv[a]; break;
            case 1: n = std::atoi(argv[a]); break;
            case 2: nrepeat = std::atoi(argv[a]); break;
            case 3: nlaunch = std::atoi(argv[a]); break;
            default: break;
        }
    }

    arch_guard();

    const ModeEntry* m = find_mode(mode_name.c_str());
    if (!m) {
        std::fprintf(stderr, "FAIL: unknown mode '%s' (try --list-kernels)\n", mode_name.c_str());
        return 2;
    }
    if (n < 64 || nlaunch < 1 || nrepeat < 1) {
        std::fprintf(stderr, "FAIL: bad nx/nrepeat/nlaunch\n");
        return 2;
    }

    double *x = nullptr, *out = nullptr;
    CUDA_OK(cudaMalloc(&x, sizeof(double) * (size_t)n * 3));
    CUDA_OK(cudaMalloc(&out, sizeof(double) * (size_t)n * 6));

    dim3 block(128), grid((unsigned)((n + 127) / 128));

    // Poison once, before the warm-up -- not per launch.
    w32::poison<<<dim3((unsigned)(((size_t)n * 6 + 127) / 128)), block>>>(out, (size_t)n * 6, -1.0);
    check_launch("poison");
    w32::init_input<<<grid, block>>>(n, x);
    check_launch("init_input");

    // Untimed warm-up, or launch 1 of the timed loop pays JIT//first-touch.
    m->launch(n, nrepeat, x, out, grid, block);
    check_launch(m->name);

    cudaEvent_t ev0, ev1;
    CUDA_OK(cudaEventCreate(&ev0));
    CUDA_OK(cudaEventCreate(&ev1));
    float ms_min = FLT_MAX, ms_sum = 0.0f;
    for (int il = 0; il < nlaunch; ++il) {
        CUDA_OK(cudaEventRecord(ev0, 0));
        m->launch(n, nrepeat, x, out, grid, block);
        CUDA_OK(cudaEventRecord(ev1, 0));
        CUDA_OK(cudaEventSynchronize(ev1));
        float ms = 0.0f;
        CUDA_OK(cudaEventElapsedTime(&ms, ev0, ev1));
        if (ms < ms_min) ms_min = ms;
        ms_sum += ms;
    }
    CUDA_OK(cudaEventDestroy(ev0));
    CUDA_OK(cudaEventDestroy(ev1));
    check_launch(m->name);

    // Clocks cannot be locked without root here, so min is the figure to
    // compare across interleaved rounds; a mean far above it means throttling.
    const float ms_avg = ms_sum / (float)nlaunch;
    std::printf("time_min_us=%14.3f\n", ms_min * 1.0e3f);
    std::printf("time_avg_us=%14.3f\n", ms_avg * 1.0e3f);
    std::printf("nlaunch=%d\n", nlaunch);
    std::printf("time_ratio=%.4f\n", ms_min > 0.0f ? ms_avg / ms_min : 0.0f);

    std::vector<double> hout((size_t)n * 6);
    CUDA_OK(cudaMemcpy(hout.data(), out, sizeof(double) * (size_t)n * 6, cudaMemcpyDeviceToHost));

    if (dump_path) {
        FILE* f = std::fopen(dump_path, "wb");
        if (!f) { std::fprintf(stderr, "FAIL: cannot open %s\n", dump_path); return 5; }
        std::fwrite(hout.data(), sizeof(double), (size_t)n * 6, f);
        std::fclose(f);
        std::printf("dump=%s\n", dump_path);
    }

    // Checksum window. The Fortran driver sums hall(1:n-5,:) for EVERY mode,
    // so the recorded WENO7/9 numbers include 2 and 4 rows of -1.0 poison
    // (-12 and -24 baked in). Reproduced here on purpose, and printed so a
    // regression is visible instead of silent.
    const int wend = n - 5;
    double csum = 0.0, small = 0.0;
    long long nbad = 0;
    unsigned long long hash = 1469598103934665603ULL;  // FNV-1a-64
    for (int c = 0; c < 6; ++c) {
        for (int i = 0; i < wend; ++i) {
            const double v = hout[(size_t)c * n + i];
            csum += v;
            if (!(v == v && fabs(v) <= DBL_MAX)) ++nbad;
            unsigned long long bits;
            std::memcpy(&bits, &v, sizeof(bits));
            for (int b = 0; b < 8; ++b) {
                hash ^= (bits >> (8 * b)) & 0xFFULL;
                hash *= 1099511628211ULL;
            }
        }
        for (int i = 0; i < 8 && i < wend; ++i) small += hout[(size_t)c * n + i];
    }

    std::printf("mode=%s\n", m->name);
    std::printf("kernel=%s\n", m->kernel);
    std::printf("n=%d nrepeat=%d\n", n, nrepeat);
    std::printf("checksum=%16.8E\n", small);
    std::printf("checksum_all=%24.16E\n", csum);
    std::printf("checksum_bits=%llu\n", hash);
    std::printf("checksum_window_end=%d\n", wend);
    std::printf("nonfinite=%lld\n", nbad);

    CUDA_OK(cudaFree(x));
    CUDA_OK(cudaFree(out));

    if (nbad > 0) {
        std::fprintf(stderr, "FAIL: non-finite values in output\n");
        return 4;
    }
    return 0;
}
