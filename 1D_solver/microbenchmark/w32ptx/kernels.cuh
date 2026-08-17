// kernels.cuh -- one shared face body, one kernel template, and the ablations.
//
// EVERY variant goes through the SAME face body and the SAME kernel template.
// Only the reciprocal policy (weno_w32.cuh) changes between rungs. That is the
// experimental control: a timing difference cannot be an artifact of two
// different kernel bodies, an accidentally different load pattern, or a
// different epilogue.
#pragma once

#include "weno_w32.cuh"

namespace w32 {

// ---------------------------------------------------------------------------
// Stencil selectors. R == r == the number of candidate substencils.
//   R = 3 -> WENO5-Z, R = 4 -> WENO7-Z, R = 5 -> WENO9-Z
// Width per bias is 2R-1; the left and right biases together span 2R cells,
// which is why exactly 2R doubles are loaded per variable (LDG = 6R).
// ---------------------------------------------------------------------------
template <int R, class Rc> struct WSel;
template <int R> struct PSel;
template <int R> struct W64Sel;

template <class Rc> struct WSel<3, Rc> {
    __device__ __forceinline__ static void left(const double* v, float* w)
    { weights5_32_left<Rc>(v[0], v[1], v[2], v[3], v[4], w[0], w[1], w[2]); }
    __device__ __forceinline__ static void right(const double* v, float* w)
    { weights5_32_right<Rc>(v[0], v[1], v[2], v[3], v[4], w[0], w[1], w[2]); }
};
template <class Rc> struct WSel<4, Rc> {
    __device__ __forceinline__ static void left(const double* v, float* w)
    { weights7_32_left<Rc>(v[0], v[1], v[2], v[3], v[4], v[5], v[6], w[0], w[1], w[2], w[3]); }
    __device__ __forceinline__ static void right(const double* v, float* w)
    { weights7_32_right<Rc>(v[0], v[1], v[2], v[3], v[4], v[5], v[6], w[0], w[1], w[2], w[3]); }
};
template <class Rc> struct WSel<5, Rc> {
    __device__ __forceinline__ static void left(const double* v, float* w)
    { weights9_32_left<Rc>(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                           w[0], w[1], w[2], w[3], w[4]); }
    __device__ __forceinline__ static void right(const double* v, float* w)
    { weights9_32_right<Rc>(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                            w[0], w[1], w[2], w[3], w[4]); }
};

template <> struct PSel<3> {
    __device__ __forceinline__ static void left(const double* v, double* p)
    { poly5_64_left(v[0], v[1], v[2], v[3], v[4], p[0], p[1], p[2]); }
    __device__ __forceinline__ static void right(const double* v, double* p)
    { poly5_64_right(v[0], v[1], v[2], v[3], v[4], p[0], p[1], p[2]); }
};
template <> struct PSel<4> {
    __device__ __forceinline__ static void left(const double* v, double* p)
    { poly7_64_left(v[0], v[1], v[2], v[3], v[4], v[5], v[6], p[0], p[1], p[2], p[3]); }
    __device__ __forceinline__ static void right(const double* v, double* p)
    { poly7_64_right(v[0], v[1], v[2], v[3], v[4], v[5], v[6], p[0], p[1], p[2], p[3]); }
};
template <> struct PSel<5> {
    __device__ __forceinline__ static void left(const double* v, double* p)
    { poly9_64_left(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                    p[0], p[1], p[2], p[3], p[4]); }
    __device__ __forceinline__ static void right(const double* v, double* p)
    { poly9_64_right(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                     p[0], p[1], p[2], p[3], p[4]); }
};

template <> struct W64Sel<3> {
    __device__ __forceinline__ static void left(const double* v, double* w)
    { weights5_64_left(v[0], v[1], v[2], v[3], v[4], w[0], w[1], w[2]); }
    __device__ __forceinline__ static void right(const double* v, double* w)
    { weights5_64_right(v[0], v[1], v[2], v[3], v[4], w[0], w[1], w[2]); }
};
template <> struct W64Sel<4> {
    __device__ __forceinline__ static void left(const double* v, double* w)
    { weights7_64_left(v[0], v[1], v[2], v[3], v[4], v[5], v[6], w[0], w[1], w[2], w[3]); }
    __device__ __forceinline__ static void right(const double* v, double* w)
    { weights7_64_right(v[0], v[1], v[2], v[3], v[4], v[5], v[6], w[0], w[1], w[2], w[3]); }
};
template <> struct W64Sel<5> {
    __device__ __forceinline__ static void left(const double* v, double* w)
    { weights9_64_left(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                       w[0], w[1], w[2], w[3], w[4]); }
    __device__ __forceinline__ static void right(const double* v, double* w)
    { weights9_64_right(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                        w[0], w[1], w[2], w[3], w[4]); }
};

// ---------------------------------------------------------------------------
// Load helper. Fortran x(n,3) is column-major, so x(i,f) -> x[f*n + i].
// Loading all 2R cells that the two biases span together is what reproduces
// nvfortran's CSE of the 2R-2 shared narrowing conversions: ptxas sees the same
// value feeding both the left and the right weights call.
// ---------------------------------------------------------------------------
template <int R>
__device__ __forceinline__ void load_stencil(const double* __restrict__ x, int n, int i, int f,
                                             double* v)
{
#pragma unroll
    for (int k = 0; k < 2 * R; ++k) v[k] = x[(size_t)f * n + i + k];
}

/// Left-to-right sum of products, matching the Fortran
/// `q = real(w0,8)*p0 + real(w1,8)*p1 + ...` association exactly. Seeding with
/// the first term rather than 0.0 avoids an extra leading DADD.
template <int R>
__device__ __forceinline__ double combine(const float* w, const double* p)
{
    double s = (double)w[0] * p[0];
#pragma unroll
    for (int k = 1; k < R; ++k) s += (double)w[k] * p[k];
    return s;
}

template <int R>
__device__ __forceinline__ double combine64(const double* w, const double* p)
{
    double s = w[0] * p[0];
#pragma unroll
    for (int k = 1; k < R; ++k) s += w[k] * p[k];
    return s;
}

// ---------------------------------------------------------------------------
// The face bodies. q[] column map matches the Fortran: q[2f] = left bias
// (v-), q[2f+1] = right bias (v+), f = 0,1,2 = rho,u,p.
// ---------------------------------------------------------------------------

/// FP32 weights x FP64 polynomials -- the kernel under study.
template <int R, class Rc>
__device__ __forceinline__ void face_w32_poly64(const double* __restrict__ x, int n, int i,
                                                double* q)
{
#pragma unroll
    for (int f = 0; f < 3; ++f) {
        double v[2 * R];
        load_stencil<R>(x, n, i, f, v);
        float w[R];
        double p[R];
        WSel<R, Rc>::left(v, w);
        PSel<R>::left(v, p);
        q[2 * f] = combine<R>(w, p);
        WSel<R, Rc>::right(v + 1, w);
        PSel<R>::right(v + 1, p);
        q[2 * f + 1] = combine<R>(w, p);
    }
}

/// Ablation: the FP32 weight stream ALONE. Combined in FP32 so the FP64 pipe
/// stays near zero (only 6 widening F2F for the stores survive). This is the
/// denominator for the co-issue question -- t(w32_poly64) - t(w32_only) is the
/// marginal cost of adding the FP64 half, and both terms sit far above the
/// ~200 us DRAM floor, unlike poly64_only which is buried in it.
template <int R, class Rc>
__device__ __forceinline__ void face_w32_only(const double* __restrict__ x, int n, int i,
                                              double* q)
{
#pragma unroll
    for (int f = 0; f < 3; ++f) {
        double v[2 * R];
        load_stencil<R>(x, n, i, f, v);
        float w[R];
        WSel<R, Rc>::left(v, w);
        float s = w[0];
#pragma unroll
        for (int k = 1; k < R; ++k) s += (float)(k + 1) * w[k];
        q[2 * f] = (double)s;
        WSel<R, Rc>::right(v + 1, w);
        s = w[0];
#pragma unroll
        for (int k = 1; k < R; ++k) s += (float)(k + 1) * w[k];
        q[2 * f + 1] = (double)s;
    }
}

/// Full FP64 -- the in-binary accuracy reference. Never a timing baseline.
template <int R>
__device__ __forceinline__ void face_weno64(const double* __restrict__ x, int n, int i, double* q)
{
#pragma unroll
    for (int f = 0; f < 3; ++f) {
        double v[2 * R];
        load_stencil<R>(x, n, i, f, v);
        double w[R], p[R];
        W64Sel<R>::left(v, w);
        PSel<R>::left(v, p);
        q[2 * f] = combine64<R>(w, p);
        W64Sel<R>::right(v + 1, w);
        PSel<R>::right(v + 1, p);
        q[2 * f + 1] = combine64<R>(w, p);
    }
}

// ---------------------------------------------------------------------------
// The kernel template. `Body` is a functor tag selecting one of the face
// bodies above; everything else -- guard, repeat loop, DCE guard, epilogue --
// is shared verbatim so no variant can differ in scaffolding.
//
// `#pragma unroll 1` on the repeat loop: nvfortran unrolls `do k=1,nrepeat` by
// 4 with a predicated remainder, so ~13 of WENO9's 37 DADDs in the reference
// are unroll copies never executed at nrepeat=1. Pinning the unroll here keeps
// the C++ count clean and makes the reference delta a single documented number
// rather than a moving target.
// ---------------------------------------------------------------------------
struct BodyW32Poly64 {
    template <int R, class Rc>
    __device__ __forceinline__ static void run(const double* __restrict__ x, int n, int i, double* q)
    { face_w32_poly64<R, Rc>(x, n, i, q); }
};
struct BodyW32Only {
    template <int R, class Rc>
    __device__ __forceinline__ static void run(const double* __restrict__ x, int n, int i, double* q)
    { face_w32_only<R, Rc>(x, n, i, q); }
};
struct BodyWeno64 {
    template <int R, class Rc>
    __device__ __forceinline__ static void run(const double* __restrict__ x, int n, int i, double* q)
    { face_weno64<R>(x, n, i, q); }
};

template <int R, class Rc, class Body>
__device__ __forceinline__ void kernel_body(int n, int nrepeat,
                                            const double* __restrict__ x,
                                            double* __restrict__ out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n - (2 * R - 1)) return;

    double q[6];
    double acc = 0.0;

    // do/while rather than for(k=0;k<nrepeat;++k): the caller guarantees
    // nrepeat >= 1, so there is no zero-trip path and q needs no initialiser.
    // Worth 8 instructions at WENO9 _rcp (1160 vs 1168) and 0 at _ref -- i.e.
    // marginal, kept for tidiness rather than for speed.
    //
    // The repeat loop does NOT amplify arithmetic: NVVM hoists the whole
    // loop-invariant face body out of it and leaves only the ~6 DADD accumulator,
    // measured at ~3.5 us/iteration for every mode regardless of its arithmetic.
    // So nrepeat is a launch-shape knob only, run_w32ptx.sh rejects it, and
    // amplification would need a runtime-opaque index offset instead.
    int k = 0;
#pragma unroll 1
    do {
        Body::template run<R, Rc>(x, n, i, q);
        acc += q[0] + q[1] + q[2] + q[3] + q[4] + q[5];
    } while (++k < nrepeat);

    // The 1e-30*acc term is what stops NVVM deleting the repeat loop. The
    // factor keeps it numerically invisible in the checksum. All six columns
    // must be stored or half the kernel disappears.
    out[i] = q[0] + 1.0e-30 * acc;
    out[(size_t)1 * n + i] = q[1];
    out[(size_t)2 * n + i] = q[2];
    out[(size_t)3 * n + i] = q[3];
    out[(size_t)4 * n + i] = q[4];
    out[(size_t)5 * n + i] = q[5];
}

// ---------------------------------------------------------------------------
// Setup kernels
// ---------------------------------------------------------------------------

/// Exact transcription of ../weno_micro.f90:1171-1187. NOTE the 1-based index:
/// Fortran's `i` runs 1..n and `z = dble(i)/dble(n)`. Using a 0-based i here
/// would shift the whole field and change every checksum.
__global__ void init_input(int n, double* __restrict__ x)
{
    const int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    if (i0 >= n) return;
    const int i = i0 + 1;
    const double z = (double)i / (double)(n > 1 ? n : 1);
    const double stepv = (i > n / 2) ? 1.0 : 0.0;
    x[(size_t)0 * n + i0] = 1.0 + 0.03 * z + 0.15 * stepv;
    x[(size_t)1 * n + i0] = 0.2 + 0.02 * z - 0.04 * stepv;
    x[(size_t)2 * n + i0] = 1.0 - 0.01 * z + 0.30 * stepv;
}

/// Poison the output so a kernel that never launched cannot masquerade as a
/// valid result. This is what exposed cudaErrorInvalidPtx (218) in the Fortran
/// harness -- cudaDeviceSynchronize() alone returns 0 when the launch itself
/// was rejected, and an untouched zero buffer was reported as checksum=0.0.
__global__ void poison(double* __restrict__ p, size_t cnt, double val)
{
    const size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < cnt) p[i] = val;
}

}  // namespace w32
