// weno_w32.cuh -- WENO-Z arithmetic for the w32_poly64 split, extracted from
// ../weno_micro.f90 (weights{5,7,9}_32_{left,right} at lines 167-254 / 420-522 /
// 738-868, poly{5,7,9}_64_{left,right} at 236-254 / 524-544 / 870-892).
//
// FP32 smoothness-indicator weights x FP64 candidate polynomials. The weights
// carry ~88% of the arithmetic and every division; the polynomials and the final
// combine stay FP64. The split is preserved on purpose -- folding the weight
// normalisation into the polynomial scaling the way the solver does would move
// FP64 work across the very boundary this benchmark measures.
//
// WHY THE RECIPROCAL IS A POLICY
// nvfortran -fast does not relax FP32 division, so every `min(tau/(b+eps), cap)`
// becomes the full IEEE div.rn.f32 sequence: MUFU.RCP + FCHK + 5 FFMA + a
// predicated BRA + CALL.REL.NOINC into a ~101-instruction out-of-line stub +
// BSYNC. Measured per face: 24 / 30 / 36 CALL sites at WENO5 / 7 / 9. Those are
// basic-block boundaries, and ptxas cannot co-schedule across one -- the same
// mechanism that made deleting the solver's order-degrading boundary ladder the
// thing that unblocked ILP above ORDER=2.
//
// So the reciprocal is the only degree of freedom between the rungs. Everything
// else -- betas, tau, dd, eps, the clamp, the polynomials, the combine -- is
// written once and shared, which is what makes a timing difference attributable.
//
// TRAPS ENCODED HERE, each of which has silently cost a debugging session:
//  * `dd` is MIRRORED between biases. Until 2026-08-17 both sides used the
//    left-biased d, making v+ THIRD order instead of fifth. Instruction counts
//    are identical, so no timing showed it and checksum_all could not see it
//    (piecewise-linear input; every candidate is exact on linear data, so the
//    whole-array sum moved by 3e-11 out of 1.27e6). ../report/check_weno_order.py
//    is what catches it.
//  * The `_right` betas are written in MIRRORED index order with the b-index
//    assignment reversed. Do not tidy them into forward order.
//  * `+ eps` is load-bearing. The beta coefficient rows sum to exactly zero, so
//    on constant data beta == 0 AND tau == 0; with a bare reciprocal that is
//    0 * rcp(0) = 0 * Inf = NaN. init_input is a ramp, so no checksum sees this.
//  * `fminf(., ratio_cap)` is load-bearing. If beta_k == 0 while tau != 0 then
//    r_k = Inf, a_k = Inf, invs = 0, and every weight becomes NaN.
//  * Ratios are UNBATCHED, normalisation IS batched. Batching the ratios needs
//    1/(c0*c1*...) and eps^r underflows to zero in FP32 (eps^5 = 1e-100), giving
//    1/0 = Inf then 0*Inf = NaN on a smooth plateau. The FP64 twins may batch.
//  * WENO5 keeps the textbook sum-of-two-squares beta form. That is a WENO5
//    special case; the published WENO7/9 betas are dense quadratic forms.
//    Rewriting r=3 would change WENO5's instruction count and invalidate every
//    prior WENO5 timing.
#pragma once

#include <cuda_runtime.h>

namespace w32 {

// ---------------------------------------------------------------------------
// Inline-PTX reciprocal primitives
//
// All are NON-volatile on purpose. Each is a pure function of its input whose
// result is consumed, so NVVM will not delete it, and leaving off `volatile`
// keeps CSE and loop-invariant hoisting -- both of which we want. No memory is
// touched, so no "memory" clobber. Marking these volatile would block the CSE
// of identical reciprocals and inflate the MUFU count.
//
// The braces around multi-instruction bodies are mandatory: without them a
// second instantiation in the same function redeclares the .reg names and ptxas
// fails with "Redefinition of variable".
//
// %0 is written only by the final instruction in every body, so these are safe
// even if the register allocator aliases an output with an input. Do NOT rely on
// "=&f" -- nvcc's inline-PTX guide does not define the earlyclobber modifier.
//
// Float literals are written in the 0f........ hex form throughout. PTX parses
// decimal literals as double and then converts, so `0.1` in a .f32 instruction
// is not reliably the float you meant. 0f3F800000 == 1.0f.
// ---------------------------------------------------------------------------

/// One MUFU.RCP on the SFU pipe. Max relative error 2^-23.
///
/// THE `.ftz` ON THE INSTRUCTION IS LOAD-BEARING, and it is not the same thing
/// as the -ftz compile flag. Measured in isolation on sm_89:
///   rcp.approx.f32      -> FSETP.GEU + FSETP.GT + FSEL + FSEL + FMUL
///                          + MUFU.RCP + FMUL   = 7 instructions
///   rcp.approx.ftz.f32  -> MUFU.RCP            = 1 instruction
/// Without `.ftz`, ptxas must keep denormal inputs and overflowing results
/// correct, so it emits an inline predicated range-scaling sequence around the
/// MUFU. At WENO9 that is 6 extra instructions x 36 reciprocals per face = ~216
/// instructions, which swallowed most of the saving this rung exists to produce
/// (measured: the FP32 stream fell by only 29 instead of ~175). Compiling with
/// -ftz=true does NOT fix it -- tried, the guards remain and the total grows.
///
/// Numerically inert on this data: every reciprocal argument here is either
/// beta+eps >= 1e-20 or sum(alpha) >= 1, both far above the 1.18e-38 min normal,
/// and the results (<= 1e20 and <= 1) are normal too. No denormal can arise, so
/// there is nothing for `.ftz` to flush. The `+eps` and the ratio_cap clamp
/// still cover the pathological beta == 0 cases.
__device__ __forceinline__ float rcp_approx(float x)
{
    float y;
    asm("rcp.approx.ftz.f32 %0, %1;" : "=f"(y) : "f"(x));
    return y;
}

/// a/b to ~2 ulp as MUFU.RCP + FMUL, spelled out rather than delegated to
/// `div.approx.ftz.f32`.
///
/// Equivalent either way -- verified in isolation on sm_89, both forms are
/// exactly two instructions (`div.approx.ftz.f32` gives MUFU.RCP + FMUL.FTZ,
/// this gives MUFU.RCP + FMUL). The explicit form is kept only because it puts
/// the `.ftz` decision at the same site as the rcp, where the comment above
/// explains why it matters. No performance claim attaches to the choice.
__device__ __forceinline__ float div_approx(float a, float b)
{
    float q;
    asm("{\n\t"
        ".reg .f32 y;\n\t"
        "rcp.approx.ftz.f32 y,  %2;\n\t"
        "mul.rn.f32     %0, %1, y;\n\t"
        "}"
        : "=f"(q) : "f"(a), "f"(b));
    return q;
}

/// div.full.f32 -- full-range approximate divide, ~2 ulp. More instructions
/// than div.approx but still branch-free and still stub-free. This is the
/// attribution control: it removes the basic-block boundaries while removing
/// almost none of the FP32 arithmetic.
__device__ __forceinline__ float div_full(float a, float b)
{
    float q;
    asm("div.full.ftz.f32 %0, %1, %2;" : "=f"(q) : "f"(a), "f"(b));
    return q;
}

/// rcp.approx + one Newton-Raphson step: ~0.5-1 ulp, branch-free.
///   e  = 1 - x*r0      (fma, so the residual costs one rounding)
///   r1 = r0 + r0*e
/// Needed for the NORMALISATION specifically. A raw rcp.approx there puts its
/// 2^-23 error on all r weights as a common multiplicative factor, so sum(w)
/// becomes 1 + 1.2e-7 instead of 1. Nothing damps that -- unlike a ratio error,
/// which is multiplied by an O(h^r) candidate difference -- so it would floor
/// the reconstruction above the accuracy gate. No CUDA intrinsic exposes this:
/// every __frcp_* is the branchy IEEE one, which is what we are removing.
///
/// NOTE the explicit `neg.f32`. PTX rejects an operand negation on `fma`
/// ("Operand negation not allowed for instruction 'fma'") even though SASS
/// supports it -- the reference division sequence itself contains
/// `FFMA R0, -R63, R64, 1`. Writing the negation as its own instruction lets
/// ptxas fold it back into the FFMA source modifier, so it costs nothing in
/// SASS while keeping the PTX legal.
__device__ __forceinline__ float rcp_refined(float x)
{
    float y;
    asm("{\n\t"
        ".reg .f32 r0, e, nx;\n\t"
        "rcp.approx.ftz.f32 r0, %1;\n\t"
        "neg.f32        nx, %1;\n\t"
        "fma.rn.f32     e,  nx, r0, 0f3F800000;\n\t"
        "fma.rn.f32     %0, r0, e,  r0;\n\t"
        "}"
        : "=f"(y) : "f"(x));
    return y;
}

// ---------------------------------------------------------------------------
// Reciprocal policies -- the ONLY difference between the rungs.
//   ratio(tau, bpe) -> tau / (beta + eps)
//   norm(s)         -> 1 / sum(alpha)
// ---------------------------------------------------------------------------

/// Baseline: what nvfortran -fast emits. IEEE div.rn.f32 at every site, with
/// its FCHK range test and out-of-line slow-path CALL.
struct RecipRef {
    static constexpr const char* tag = "ref";
    __device__ __forceinline__ static float ratio(float tau, float bpe) { return tau / bpe; }
    __device__ __forceinline__ static float norm(float s) { return 1.0f / s; }
};

/// PRIMARY rung. Approximate ratios (2 ulp is far below the FP32 cancellation
/// already present in the dense beta forms), Newton-refined normalisation so the
/// partition of unity survives. Removes 5r FP32 FMA-pipe instructions per call
/// and all 24/30/36 basic-block boundaries.
struct RecipRcp {
    static constexpr const char* tag = "rcp";
    __device__ __forceinline__ static float ratio(float tau, float bpe) { return div_approx(tau, bpe); }
    __device__ __forceinline__ static float norm(float s) { return rcp_refined(s); }
};

/// Fallback if `rcp` fails the order gate: near-IEEE everywhere, still
/// branch-free. ~1 ulp vs div.rn's 0.5 ulp.
struct RecipRcpN {
    static constexpr const char* tag = "rcpn";
    __device__ __forceinline__ static float ratio(float tau, float bpe) { return tau * rcp_refined(bpe); }
    __device__ __forceinline__ static float norm(float s) { return rcp_refined(s); }
};

/// Attribution control: branch-free but ~as many instructions as the baseline.
/// t(divfull) - t(ref) isolates the removal of the basic-block boundaries;
/// t(rcp) - t(divfull) isolates the removal of the arithmetic.
struct RecipFull {
    static constexpr const char* tag = "divfull";
    __device__ __forceinline__ static float ratio(float tau, float bpe) { return div_full(tau, bpe); }
    __device__ __forceinline__ static float norm(float s) { return div_full(1.0f, s); }
};

// ---------------------------------------------------------------------------
// Shared constants
// ---------------------------------------------------------------------------
__device__ __constant__ constexpr float EPS32 = 1.0e-20f;
__device__ __constant__ constexpr float RATIO_CAP = 1.0e9f;

// ===========================================================================
// WENO5-Z  (r = 3)
// ===========================================================================

/// Textbook sum-of-two-squares betas. Optimal weights 1/10, 6/10, 3/10.
template <class Rc>
__device__ __forceinline__ void weights5_32_left(
    double v1, double v2, double v3, double v4, double v5,
    float& w0, float& w1, float& w2)
{
    constexpr float c13_12 = 13.0f / 12.0f;
    constexpr float dd0 = 0.1f, dd1 = 0.6f, dd2 = 0.3f;

    const float x1 = (float)v1, x2 = (float)v2, x3 = (float)v3;
    const float x4 = (float)v4, x5 = (float)v5;

    const float e0 = x1 - 2.0f * x2 + x3;
    const float f0 = x1 - 4.0f * x2 + 3.0f * x3;
    const float e1 = x2 - 2.0f * x3 + x4;
    const float f1 = x2 - x4;
    const float e2 = x3 - 2.0f * x4 + x5;
    const float f2 = 3.0f * x3 - 4.0f * x4 + x5;

    const float b0 = c13_12 * (e0 * e0) + 0.25f * (f0 * f0);
    const float b1 = c13_12 * (e1 * e1) + 0.25f * (f1 * f1);
    const float b2 = c13_12 * (e2 * e2) + 0.25f * (f2 * f2);

    const float tau = fabsf(b0 - b2);

    const float r0 = fminf(Rc::ratio(tau, b0 + EPS32), RATIO_CAP);
    const float r1 = fminf(Rc::ratio(tau, b1 + EPS32), RATIO_CAP);
    const float r2 = fminf(Rc::ratio(tau, b2 + EPS32), RATIO_CAP);

    const float a0 = dd0 * (1.0f + r0 * r0);
    const float a1 = dd1 * (1.0f + r1 * r1);
    const float a2 = dd2 * (1.0f + r2 * r2);

    const float invs = Rc::norm(a0 + a1 + a2);
    w0 = a0 * invs;
    w1 = a1 * invs;
    w2 = a2 * invs;
}

/// Right bias: identical betas, MIRRORED optimal weights (3/10, 6/10, 1/10).
template <class Rc>
__device__ __forceinline__ void weights5_32_right(
    double v1, double v2, double v3, double v4, double v5,
    float& w0, float& w1, float& w2)
{
    constexpr float c13_12 = 13.0f / 12.0f;
    constexpr float dd0 = 0.3f, dd1 = 0.6f, dd2 = 0.1f;

    const float x1 = (float)v1, x2 = (float)v2, x3 = (float)v3;
    const float x4 = (float)v4, x5 = (float)v5;

    const float e0 = x1 - 2.0f * x2 + x3;
    const float f0 = x1 - 4.0f * x2 + 3.0f * x3;
    const float e1 = x2 - 2.0f * x3 + x4;
    const float f1 = x2 - x4;
    const float e2 = x3 - 2.0f * x4 + x5;
    const float f2 = 3.0f * x3 - 4.0f * x4 + x5;

    const float b0 = c13_12 * (e0 * e0) + 0.25f * (f0 * f0);
    const float b1 = c13_12 * (e1 * e1) + 0.25f * (f1 * f1);
    const float b2 = c13_12 * (e2 * e2) + 0.25f * (f2 * f2);

    const float tau = fabsf(b0 - b2);

    const float r0 = fminf(Rc::ratio(tau, b0 + EPS32), RATIO_CAP);
    const float r1 = fminf(Rc::ratio(tau, b1 + EPS32), RATIO_CAP);
    const float r2 = fminf(Rc::ratio(tau, b2 + EPS32), RATIO_CAP);

    const float a0 = dd0 * (1.0f + r0 * r0);
    const float a1 = dd1 * (1.0f + r1 * r1);
    const float a2 = dd2 * (1.0f + r2 * r2);

    const float invs = Rc::norm(a0 + a1 + a2);
    w0 = a0 * invs;
    w1 = a1 * invs;
    w2 = a2 * invs;
}

/// 1/6 folded into the coefficients, so each candidate is a division-free chain.
__device__ __forceinline__ void poly5_64_left(
    double v1, double v2, double v3, double v4, double v5,
    double& p0, double& p1, double& p2)
{
    constexpr double c1 = 1.0 / 6.0, c2 = 2.0 / 6.0, c5 = 5.0 / 6.0;
    constexpr double c7 = 7.0 / 6.0, c11 = 11.0 / 6.0;
    p0 =  c2 * v1 - c7 * v2 + c11 * v3;
    p1 = -c1 * v2 + c5 * v3 +  c2 * v4;
    p2 =  c2 * v3 + c5 * v4 -  c1 * v5;
}

__device__ __forceinline__ void poly5_64_right(
    double v1, double v2, double v3, double v4, double v5,
    double& p0, double& p1, double& p2)
{
    constexpr double c1 = 1.0 / 6.0, c2 = 2.0 / 6.0, c5 = 5.0 / 6.0;
    constexpr double c7 = 7.0 / 6.0, c11 = 11.0 / 6.0;
    p0 = -c1 * v1 + c5 * v2 +  c2 * v3;
    p1 =  c2 * v2 + c5 * v3 -  c1 * v4;
    p2 = c11 * v3 - c7 * v4 +  c2 * v5;
}

// ===========================================================================
// WENO7-Z  (r = 4).  Betas scaled by 240; the common factor cancels in
// tau/(beta+eps).  Optimal weights {1,12,18,4}/35, mirrored for the right bias.
// ===========================================================================

template <class Rc>
__device__ __forceinline__ void weights7_32_left(
    double v1, double v2, double v3, double v4, double v5, double v6, double v7,
    float& w0, float& w1, float& w2, float& w3)
{
    constexpr float dd0 = 1.0f / 35.0f, dd1 = 12.0f / 35.0f;
    constexpr float dd2 = 18.0f / 35.0f, dd3 = 4.0f / 35.0f;

    const float x1 = (float)v1, x2 = (float)v2, x3 = (float)v3, x4 = (float)v4;
    const float x5 = (float)v5, x6 = (float)v6, x7 = (float)v7;

    const float b0 = x1 * (547.0f * x1 - 3882.0f * x2 + 4642.0f * x3 - 1854.0f * x4)
                   + x2 * (7043.0f * x2 - 17246.0f * x3 + 7042.0f * x4)
                   + x3 * (11003.0f * x3 - 9402.0f * x4)
                   + x4 * (2107.0f * x4);
    const float b1 = x2 * (267.0f * x2 - 1642.0f * x3 + 1602.0f * x4 - 494.0f * x5)
                   + x3 * (2843.0f * x3 - 5966.0f * x4 + 1922.0f * x5)
                   + x4 * (3443.0f * x4 - 2522.0f * x5)
                   + x5 * (547.0f * x5);
    const float b2 = x3 * (547.0f * x3 - 2522.0f * x4 + 1922.0f * x5 - 494.0f * x6)
                   + x4 * (3443.0f * x4 - 5966.0f * x5 + 1602.0f * x6)
                   + x5 * (2843.0f * x5 - 1642.0f * x6)
                   + x6 * (267.0f * x6);
    const float b3 = x4 * (2107.0f * x4 - 9402.0f * x5 + 7042.0f * x6 - 1854.0f * x7)
                   + x5 * (11003.0f * x5 - 17246.0f * x6 + 4642.0f * x7)
                   + x6 * (7043.0f * x6 - 3882.0f * x7)
                   + x7 * (547.0f * x7);

    const float tau = fabsf(1.0f * b0 + 3.0f * b1 - 3.0f * b2 - 1.0f * b3);

    const float r0 = fminf(Rc::ratio(tau, b0 + EPS32), RATIO_CAP);
    const float r1 = fminf(Rc::ratio(tau, b1 + EPS32), RATIO_CAP);
    const float r2 = fminf(Rc::ratio(tau, b2 + EPS32), RATIO_CAP);
    const float r3 = fminf(Rc::ratio(tau, b3 + EPS32), RATIO_CAP);

    const float a0 = dd0 * (1.0f + r0 * r0);
    const float a1 = dd1 * (1.0f + r1 * r1);
    const float a2 = dd2 * (1.0f + r2 * r2);
    const float a3 = dd3 * (1.0f + r3 * r3);

    const float invs = Rc::norm(a0 + a1 + a2 + a3);
    w0 = a0 * invs;
    w1 = a1 * invs;
    w2 = a2 * invs;
    w3 = a3 * invs;
}

/// Mirrored index order with the b-index assignment reversed -- transcribed
/// verbatim from ../weno_micro.f90:489-504. Note tau's sign pattern is also
/// flipped here (bit-identical under fabsf, but kept for source fidelity).
template <class Rc>
__device__ __forceinline__ void weights7_32_right(
    double v1, double v2, double v3, double v4, double v5, double v6, double v7,
    float& w0, float& w1, float& w2, float& w3)
{
    constexpr float dd0 = 4.0f / 35.0f, dd1 = 18.0f / 35.0f;
    constexpr float dd2 = 12.0f / 35.0f, dd3 = 1.0f / 35.0f;

    const float x1 = (float)v1, x2 = (float)v2, x3 = (float)v3, x4 = (float)v4;
    const float x5 = (float)v5, x6 = (float)v6, x7 = (float)v7;

    const float b3 = x7 * (547.0f * x7 - 3882.0f * x6 + 4642.0f * x5 - 1854.0f * x4)
                   + x6 * (7043.0f * x6 - 17246.0f * x5 + 7042.0f * x4)
                   + x5 * (11003.0f * x5 - 9402.0f * x4)
                   + x4 * (2107.0f * x4);
    const float b2 = x6 * (267.0f * x6 - 1642.0f * x5 + 1602.0f * x4 - 494.0f * x3)
                   + x5 * (2843.0f * x5 - 5966.0f * x4 + 1922.0f * x3)
                   + x4 * (3443.0f * x4 - 2522.0f * x3)
                   + x3 * (547.0f * x3);
    const float b1 = x5 * (547.0f * x5 - 2522.0f * x4 + 1922.0f * x3 - 494.0f * x2)
                   + x4 * (3443.0f * x4 - 5966.0f * x3 + 1602.0f * x2)
                   + x3 * (2843.0f * x3 - 1642.0f * x2)
                   + x2 * (267.0f * x2);
    const float b0 = x4 * (2107.0f * x4 - 9402.0f * x3 + 7042.0f * x2 - 1854.0f * x1)
                   + x3 * (11003.0f * x3 - 17246.0f * x2 + 4642.0f * x1)
                   + x2 * (7043.0f * x2 - 3882.0f * x1)
                   + x1 * (547.0f * x1);

    const float tau = fabsf(-1.0f * b0 - 3.0f * b1 + 3.0f * b2 + 1.0f * b3);

    const float r0 = fminf(Rc::ratio(tau, b0 + EPS32), RATIO_CAP);
    const float r1 = fminf(Rc::ratio(tau, b1 + EPS32), RATIO_CAP);
    const float r2 = fminf(Rc::ratio(tau, b2 + EPS32), RATIO_CAP);
    const float r3 = fminf(Rc::ratio(tau, b3 + EPS32), RATIO_CAP);

    const float a0 = dd0 * (1.0f + r0 * r0);
    const float a1 = dd1 * (1.0f + r1 * r1);
    const float a2 = dd2 * (1.0f + r2 * r2);
    const float a3 = dd3 * (1.0f + r3 * r3);

    const float invs = Rc::norm(a0 + a1 + a2 + a3);
    w0 = a0 * invs;
    w1 = a1 * invs;
    w2 = a2 * invs;
    w3 = a3 * invs;
}

__device__ __forceinline__ void poly7_64_left(
    double v1, double v2, double v3, double v4, double v5, double v6, double v7,
    double& p0, double& p1, double& p2, double& p3)
{
    p0 = -(1.0 / 4.0) * v1 + (13.0 / 12.0) * v2 - (23.0 / 12.0) * v3 + (25.0 / 12.0) * v4;
    p1 =  (1.0 / 12.0) * v2 - (5.0 / 12.0) * v3 + (13.0 / 12.0) * v4 + (1.0 / 4.0) * v5;
    p2 = -(1.0 / 12.0) * v3 + (7.0 / 12.0) * v4 + (7.0 / 12.0) * v5 - (1.0 / 12.0) * v6;
    p3 =  (1.0 / 4.0) * v4 + (13.0 / 12.0) * v5 - (5.0 / 12.0) * v6 + (1.0 / 12.0) * v7;
}

__device__ __forceinline__ void poly7_64_right(
    double v1, double v2, double v3, double v4, double v5, double v6, double v7,
    double& p0, double& p1, double& p2, double& p3)
{
    p3 = -(1.0 / 4.0) * v7 + (13.0 / 12.0) * v6 - (23.0 / 12.0) * v5 + (25.0 / 12.0) * v4;
    p2 =  (1.0 / 12.0) * v6 - (5.0 / 12.0) * v5 + (13.0 / 12.0) * v4 + (1.0 / 4.0) * v3;
    p1 = -(1.0 / 12.0) * v5 + (7.0 / 12.0) * v4 + (7.0 / 12.0) * v3 - (1.0 / 12.0) * v2;
    p0 =  (1.0 / 4.0) * v4 + (13.0 / 12.0) * v3 - (5.0 / 12.0) * v2 + (1.0 / 12.0) * v1;
}

// ===========================================================================
// WENO9-Z  (r = 5).  Betas scaled by 10080.  Optimal weights
// {1/126, 10/63, 10/21, 20/63, 5/126}, mirrored for the right bias.
// tau's sign pattern is the same on both biases here.
// ===========================================================================

template <class Rc>
__device__ __forceinline__ void weights9_32_left(
    double v1, double v2, double v3, double v4, double v5,
    double v6, double v7, double v8, double v9,
    float& w0, float& w1, float& w2, float& w3, float& w4)
{
    constexpr float dd0 = 1.0f / 126.0f, dd1 = 10.0f / 63.0f, dd2 = 10.0f / 21.0f;
    constexpr float dd3 = 20.0f / 63.0f, dd4 = 5.0f / 126.0f;

    const float x1 = (float)v1, x2 = (float)v2, x3 = (float)v3, x4 = (float)v4;
    const float x5 = (float)v5, x6 = (float)v6, x7 = (float)v7, x8 = (float)v8;
    const float x9 = (float)v9;

    const float b0 = x1 * (45316.0f * x1 - 417002.0f * x2 + 729726.0f * x3 - 576014.0f * x4 + 172658.0f * x5)
                   + x2 * (965926.0f * x2 - 3408792.0f * x3 + 2716916.0f * x4 - 822974.0f * x5)
                   + x3 * (3042786.0f * x3 - 4924152.0f * x4 + 1517646.0f * x5)
                   + x4 * (2041126.0f * x4 - 1299002.0f * x5)
                   + x5 * (215836.0f * x5);
    const float b1 = x2 * (13816.0f * x2 - 121742.0f * x3 + 198426.0f * x4 - 140474.0f * x5 + 36158.0f * x6)
                   + x3 * (277126.0f * x3 - 929952.0f * x4 + 674036.0f * x5 - 176594.0f * x6)
                   + x4 * (812586.0f * x4 - 1223952.0f * x5 + 330306.0f * x6)
                   + x5 * (485446.0f * x5 - 280502.0f * x6)
                   + x6 * (45316.0f * x6);
    const float b2 = x3 * (13816.0f * x3 - 102002.0f * x4 + 135846.0f * x5 - 77894.0f * x6 + 16418.0f * x7)
                   + x4 * (209926.0f * x4 - 598152.0f * x5 + 358196.0f * x6 - 77894.0f * x7)
                   + x5 * (462306.0f * x5 - 598152.0f * x6 + 135846.0f * x7)
                   + x6 * (209926.0f * x6 - 102002.0f * x7)
                   + x7 * (13816.0f * x7);
    const float b3 = x4 * (45316.0f * x4 - 280502.0f * x5 + 330306.0f * x6 - 176594.0f * x7 + 36158.0f * x8)
                   + x5 * (485446.0f * x5 - 1223952.0f * x6 + 674036.0f * x7 - 140474.0f * x8)
                   + x6 * (812586.0f * x6 - 929952.0f * x7 + 198426.0f * x8)
                   + x7 * (277126.0f * x7 - 121742.0f * x8)
                   + x8 * (13816.0f * x8);
    const float b4 = x5 * (215836.0f * x5 - 1299002.0f * x6 + 1517646.0f * x7 - 822974.0f * x8 + 172658.0f * x9)
                   + x6 * (2041126.0f * x6 - 4924152.0f * x7 + 2716916.0f * x8 - 576014.0f * x9)
                   + x7 * (3042786.0f * x7 - 3408792.0f * x8 + 729726.0f * x9)
                   + x8 * (965926.0f * x8 - 417002.0f * x9)
                   + x9 * (45316.0f * x9);

    const float tau = fabsf(1.0f * b0 + 2.0f * b1 - 6.0f * b2 + 2.0f * b3 + 1.0f * b4);

    const float r0 = fminf(Rc::ratio(tau, b0 + EPS32), RATIO_CAP);
    const float r1 = fminf(Rc::ratio(tau, b1 + EPS32), RATIO_CAP);
    const float r2 = fminf(Rc::ratio(tau, b2 + EPS32), RATIO_CAP);
    const float r3 = fminf(Rc::ratio(tau, b3 + EPS32), RATIO_CAP);
    const float r4 = fminf(Rc::ratio(tau, b4 + EPS32), RATIO_CAP);

    const float a0 = dd0 * (1.0f + r0 * r0);
    const float a1 = dd1 * (1.0f + r1 * r1);
    const float a2 = dd2 * (1.0f + r2 * r2);
    const float a3 = dd3 * (1.0f + r3 * r3);
    const float a4 = dd4 * (1.0f + r4 * r4);

    const float invs = Rc::norm(a0 + a1 + a2 + a3 + a4);
    w0 = a0 * invs;
    w1 = a1 * invs;
    w2 = a2 * invs;
    w3 = a3 * invs;
    w4 = a4 * invs;
}

/// Mirrored index order with the b-index assignment reversed -- verbatim from
/// ../weno_micro.f90:823-847.
template <class Rc>
__device__ __forceinline__ void weights9_32_right(
    double v1, double v2, double v3, double v4, double v5,
    double v6, double v7, double v8, double v9,
    float& w0, float& w1, float& w2, float& w3, float& w4)
{
    constexpr float dd0 = 5.0f / 126.0f, dd1 = 20.0f / 63.0f, dd2 = 10.0f / 21.0f;
    constexpr float dd3 = 10.0f / 63.0f, dd4 = 1.0f / 126.0f;

    const float x1 = (float)v1, x2 = (float)v2, x3 = (float)v3, x4 = (float)v4;
    const float x5 = (float)v5, x6 = (float)v6, x7 = (float)v7, x8 = (float)v8;
    const float x9 = (float)v9;

    const float b4 = x9 * (45316.0f * x9 - 417002.0f * x8 + 729726.0f * x7 - 576014.0f * x6 + 172658.0f * x5)
                   + x8 * (965926.0f * x8 - 3408792.0f * x7 + 2716916.0f * x6 - 822974.0f * x5)
                   + x7 * (3042786.0f * x7 - 4924152.0f * x6 + 1517646.0f * x5)
                   + x6 * (2041126.0f * x6 - 1299002.0f * x5)
                   + x5 * (215836.0f * x5);
    const float b3 = x8 * (13816.0f * x8 - 121742.0f * x7 + 198426.0f * x6 - 140474.0f * x5 + 36158.0f * x4)
                   + x7 * (277126.0f * x7 - 929952.0f * x6 + 674036.0f * x5 - 176594.0f * x4)
                   + x6 * (812586.0f * x6 - 1223952.0f * x5 + 330306.0f * x4)
                   + x5 * (485446.0f * x5 - 280502.0f * x4)
                   + x4 * (45316.0f * x4);
    const float b2 = x7 * (13816.0f * x7 - 102002.0f * x6 + 135846.0f * x5 - 77894.0f * x4 + 16418.0f * x3)
                   + x6 * (209926.0f * x6 - 598152.0f * x5 + 358196.0f * x4 - 77894.0f * x3)
                   + x5 * (462306.0f * x5 - 598152.0f * x4 + 135846.0f * x3)
                   + x4 * (209926.0f * x4 - 102002.0f * x3)
                   + x3 * (13816.0f * x3);
    const float b1 = x6 * (45316.0f * x6 - 280502.0f * x5 + 330306.0f * x4 - 176594.0f * x3 + 36158.0f * x2)
                   + x5 * (485446.0f * x5 - 1223952.0f * x4 + 674036.0f * x3 - 140474.0f * x2)
                   + x4 * (812586.0f * x4 - 929952.0f * x3 + 198426.0f * x2)
                   + x3 * (277126.0f * x3 - 121742.0f * x2)
                   + x2 * (13816.0f * x2);
    const float b0 = x5 * (215836.0f * x5 - 1299002.0f * x4 + 1517646.0f * x3 - 822974.0f * x2 + 172658.0f * x1)
                   + x4 * (2041126.0f * x4 - 4924152.0f * x3 + 2716916.0f * x2 - 576014.0f * x1)
                   + x3 * (3042786.0f * x3 - 3408792.0f * x2 + 729726.0f * x1)
                   + x2 * (965926.0f * x2 - 417002.0f * x1)
                   + x1 * (45316.0f * x1);

    const float tau = fabsf(1.0f * b0 + 2.0f * b1 - 6.0f * b2 + 2.0f * b3 + 1.0f * b4);

    const float r0 = fminf(Rc::ratio(tau, b0 + EPS32), RATIO_CAP);
    const float r1 = fminf(Rc::ratio(tau, b1 + EPS32), RATIO_CAP);
    const float r2 = fminf(Rc::ratio(tau, b2 + EPS32), RATIO_CAP);
    const float r3 = fminf(Rc::ratio(tau, b3 + EPS32), RATIO_CAP);
    const float r4 = fminf(Rc::ratio(tau, b4 + EPS32), RATIO_CAP);

    const float a0 = dd0 * (1.0f + r0 * r0);
    const float a1 = dd1 * (1.0f + r1 * r1);
    const float a2 = dd2 * (1.0f + r2 * r2);
    const float a3 = dd3 * (1.0f + r3 * r3);
    const float a4 = dd4 * (1.0f + r4 * r4);

    const float invs = Rc::norm(a0 + a1 + a2 + a3 + a4);
    w0 = a0 * invs;
    w1 = a1 * invs;
    w2 = a2 * invs;
    w3 = a3 * invs;
    w4 = a4 * invs;
}

__device__ __forceinline__ void poly9_64_left(
    double v1, double v2, double v3, double v4, double v5,
    double v6, double v7, double v8, double v9,
    double& p0, double& p1, double& p2, double& p3, double& p4)
{
    p0 =  (1.0 / 5.0) * v1 - (21.0 / 20.0) * v2 + (137.0 / 60.0) * v3 - (163.0 / 60.0) * v4 + (137.0 / 60.0) * v5;
    p1 = -(1.0 / 20.0) * v2 + (17.0 / 60.0) * v3 - (43.0 / 60.0) * v4 + (77.0 / 60.0) * v5 + (1.0 / 5.0) * v6;
    p2 =  (1.0 / 30.0) * v3 - (13.0 / 60.0) * v4 + (47.0 / 60.0) * v5 + (9.0 / 20.0) * v6 - (1.0 / 20.0) * v7;
    p3 = -(1.0 / 20.0) * v4 + (9.0 / 20.0) * v5 + (47.0 / 60.0) * v6 - (13.0 / 60.0) * v7 + (1.0 / 30.0) * v8;
    p4 =  (1.0 / 5.0) * v5 + (77.0 / 60.0) * v6 - (43.0 / 60.0) * v7 + (17.0 / 60.0) * v8 - (1.0 / 20.0) * v9;
}

__device__ __forceinline__ void poly9_64_right(
    double v1, double v2, double v3, double v4, double v5,
    double v6, double v7, double v8, double v9,
    double& p0, double& p1, double& p2, double& p3, double& p4)
{
    p4 =  (1.0 / 5.0) * v9 - (21.0 / 20.0) * v8 + (137.0 / 60.0) * v7 - (163.0 / 60.0) * v6 + (137.0 / 60.0) * v5;
    p3 = -(1.0 / 20.0) * v8 + (17.0 / 60.0) * v7 - (43.0 / 60.0) * v6 + (77.0 / 60.0) * v5 + (1.0 / 5.0) * v4;
    p2 =  (1.0 / 30.0) * v7 - (13.0 / 60.0) * v6 + (47.0 / 60.0) * v5 + (9.0 / 20.0) * v4 - (1.0 / 20.0) * v3;
    p1 = -(1.0 / 20.0) * v6 + (9.0 / 20.0) * v5 + (47.0 / 60.0) * v4 - (13.0 / 60.0) * v3 + (1.0 / 30.0) * v2;
    p0 =  (1.0 / 5.0) * v5 + (77.0 / 60.0) * v4 - (43.0 / 60.0) * v3 + (17.0 / 60.0) * v2 - (1.0 / 20.0) * v1;
}

// ===========================================================================
// FP64 weights -- the ACCURACY REFERENCE only, never a timing baseline.
//
// Mathematically the same scheme as the FP32 twins above. Deliberately NOT a
// transcription of ../weno_micro.f90's weights*_64_*: those batch the ratio
// inversion (1/(c0*c1*c2)) as a division-count optimisation, which is safe in
// FP64 but is not what we need here. Unbatched is simpler and agrees to ~1e-16,
// which is 8 decades below the 5e-8 accuracy gate.
// ===========================================================================
constexpr double EPS64 = 1.0e-20;

__device__ __forceinline__ void weights5_64_left(
    double v1, double v2, double v3, double v4, double v5,
    double& w0, double& w1, double& w2)
{
    const double e0 = v1 - 2.0 * v2 + v3, f0 = v1 - 4.0 * v2 + 3.0 * v3;
    const double e1 = v2 - 2.0 * v3 + v4, f1 = v2 - v4;
    const double e2 = v3 - 2.0 * v4 + v5, f2 = 3.0 * v3 - 4.0 * v4 + v5;
    const double b0 = (13.0 / 12.0) * e0 * e0 + 0.25 * f0 * f0;
    const double b1 = (13.0 / 12.0) * e1 * e1 + 0.25 * f1 * f1;
    const double b2 = (13.0 / 12.0) * e2 * e2 + 0.25 * f2 * f2;
    const double tau = fabs(b0 - b2);
    const double r0 = tau / (b0 + EPS64), r1 = tau / (b1 + EPS64), r2 = tau / (b2 + EPS64);
    const double a0 = 0.1 * (1.0 + r0 * r0), a1 = 0.6 * (1.0 + r1 * r1), a2 = 0.3 * (1.0 + r2 * r2);
    const double invs = 1.0 / (a0 + a1 + a2);
    w0 = a0 * invs; w1 = a1 * invs; w2 = a2 * invs;
}

__device__ __forceinline__ void weights5_64_right(
    double v1, double v2, double v3, double v4, double v5,
    double& w0, double& w1, double& w2)
{
    const double e0 = v1 - 2.0 * v2 + v3, f0 = v1 - 4.0 * v2 + 3.0 * v3;
    const double e1 = v2 - 2.0 * v3 + v4, f1 = v2 - v4;
    const double e2 = v3 - 2.0 * v4 + v5, f2 = 3.0 * v3 - 4.0 * v4 + v5;
    const double b0 = (13.0 / 12.0) * e0 * e0 + 0.25 * f0 * f0;
    const double b1 = (13.0 / 12.0) * e1 * e1 + 0.25 * f1 * f1;
    const double b2 = (13.0 / 12.0) * e2 * e2 + 0.25 * f2 * f2;
    const double tau = fabs(b0 - b2);
    const double r0 = tau / (b0 + EPS64), r1 = tau / (b1 + EPS64), r2 = tau / (b2 + EPS64);
    const double a0 = 0.3 * (1.0 + r0 * r0), a1 = 0.6 * (1.0 + r1 * r1), a2 = 0.1 * (1.0 + r2 * r2);
    const double invs = 1.0 / (a0 + a1 + a2);
    w0 = a0 * invs; w1 = a1 * invs; w2 = a2 * invs;
}

__device__ __forceinline__ void weights7_64_left(
    double v1, double v2, double v3, double v4, double v5, double v6, double v7,
    double& w0, double& w1, double& w2, double& w3)
{
    const double b0 = v1 * (547.0 * v1 - 3882.0 * v2 + 4642.0 * v3 - 1854.0 * v4)
                    + v2 * (7043.0 * v2 - 17246.0 * v3 + 7042.0 * v4)
                    + v3 * (11003.0 * v3 - 9402.0 * v4) + v4 * (2107.0 * v4);
    const double b1 = v2 * (267.0 * v2 - 1642.0 * v3 + 1602.0 * v4 - 494.0 * v5)
                    + v3 * (2843.0 * v3 - 5966.0 * v4 + 1922.0 * v5)
                    + v4 * (3443.0 * v4 - 2522.0 * v5) + v5 * (547.0 * v5);
    const double b2 = v3 * (547.0 * v3 - 2522.0 * v4 + 1922.0 * v5 - 494.0 * v6)
                    + v4 * (3443.0 * v4 - 5966.0 * v5 + 1602.0 * v6)
                    + v5 * (2843.0 * v5 - 1642.0 * v6) + v6 * (267.0 * v6);
    const double b3 = v4 * (2107.0 * v4 - 9402.0 * v5 + 7042.0 * v6 - 1854.0 * v7)
                    + v5 * (11003.0 * v5 - 17246.0 * v6 + 4642.0 * v7)
                    + v6 * (7043.0 * v6 - 3882.0 * v7) + v7 * (547.0 * v7);
    const double tau = fabs(b0 + 3.0 * b1 - 3.0 * b2 - b3);
    const double r0 = tau / (b0 + EPS64), r1 = tau / (b1 + EPS64);
    const double r2 = tau / (b2 + EPS64), r3 = tau / (b3 + EPS64);
    const double a0 = (1.0 / 35.0) * (1.0 + r0 * r0), a1 = (12.0 / 35.0) * (1.0 + r1 * r1);
    const double a2 = (18.0 / 35.0) * (1.0 + r2 * r2), a3 = (4.0 / 35.0) * (1.0 + r3 * r3);
    const double invs = 1.0 / (a0 + a1 + a2 + a3);
    w0 = a0 * invs; w1 = a1 * invs; w2 = a2 * invs; w3 = a3 * invs;
}

__device__ __forceinline__ void weights7_64_right(
    double v1, double v2, double v3, double v4, double v5, double v6, double v7,
    double& w0, double& w1, double& w2, double& w3)
{
    const double b3 = v7 * (547.0 * v7 - 3882.0 * v6 + 4642.0 * v5 - 1854.0 * v4)
                    + v6 * (7043.0 * v6 - 17246.0 * v5 + 7042.0 * v4)
                    + v5 * (11003.0 * v5 - 9402.0 * v4) + v4 * (2107.0 * v4);
    const double b2 = v6 * (267.0 * v6 - 1642.0 * v5 + 1602.0 * v4 - 494.0 * v3)
                    + v5 * (2843.0 * v5 - 5966.0 * v4 + 1922.0 * v3)
                    + v4 * (3443.0 * v4 - 2522.0 * v3) + v3 * (547.0 * v3);
    const double b1 = v5 * (547.0 * v5 - 2522.0 * v4 + 1922.0 * v3 - 494.0 * v2)
                    + v4 * (3443.0 * v4 - 5966.0 * v3 + 1602.0 * v2)
                    + v3 * (2843.0 * v3 - 1642.0 * v2) + v2 * (267.0 * v2);
    const double b0 = v4 * (2107.0 * v4 - 9402.0 * v3 + 7042.0 * v2 - 1854.0 * v1)
                    + v3 * (11003.0 * v3 - 17246.0 * v2 + 4642.0 * v1)
                    + v2 * (7043.0 * v2 - 3882.0 * v1) + v1 * (547.0 * v1);
    const double tau = fabs(-b0 - 3.0 * b1 + 3.0 * b2 + b3);
    const double r0 = tau / (b0 + EPS64), r1 = tau / (b1 + EPS64);
    const double r2 = tau / (b2 + EPS64), r3 = tau / (b3 + EPS64);
    const double a0 = (4.0 / 35.0) * (1.0 + r0 * r0), a1 = (18.0 / 35.0) * (1.0 + r1 * r1);
    const double a2 = (12.0 / 35.0) * (1.0 + r2 * r2), a3 = (1.0 / 35.0) * (1.0 + r3 * r3);
    const double invs = 1.0 / (a0 + a1 + a2 + a3);
    w0 = a0 * invs; w1 = a1 * invs; w2 = a2 * invs; w3 = a3 * invs;
}

__device__ __forceinline__ void weights9_64_left(
    double v1, double v2, double v3, double v4, double v5,
    double v6, double v7, double v8, double v9,
    double& w0, double& w1, double& w2, double& w3, double& w4)
{
    const double b0 = v1 * (45316.0 * v1 - 417002.0 * v2 + 729726.0 * v3 - 576014.0 * v4 + 172658.0 * v5)
                    + v2 * (965926.0 * v2 - 3408792.0 * v3 + 2716916.0 * v4 - 822974.0 * v5)
                    + v3 * (3042786.0 * v3 - 4924152.0 * v4 + 1517646.0 * v5)
                    + v4 * (2041126.0 * v4 - 1299002.0 * v5) + v5 * (215836.0 * v5);
    const double b1 = v2 * (13816.0 * v2 - 121742.0 * v3 + 198426.0 * v4 - 140474.0 * v5 + 36158.0 * v6)
                    + v3 * (277126.0 * v3 - 929952.0 * v4 + 674036.0 * v5 - 176594.0 * v6)
                    + v4 * (812586.0 * v4 - 1223952.0 * v5 + 330306.0 * v6)
                    + v5 * (485446.0 * v5 - 280502.0 * v6) + v6 * (45316.0 * v6);
    const double b2 = v3 * (13816.0 * v3 - 102002.0 * v4 + 135846.0 * v5 - 77894.0 * v6 + 16418.0 * v7)
                    + v4 * (209926.0 * v4 - 598152.0 * v5 + 358196.0 * v6 - 77894.0 * v7)
                    + v5 * (462306.0 * v5 - 598152.0 * v6 + 135846.0 * v7)
                    + v6 * (209926.0 * v6 - 102002.0 * v7) + v7 * (13816.0 * v7);
    const double b3 = v4 * (45316.0 * v4 - 280502.0 * v5 + 330306.0 * v6 - 176594.0 * v7 + 36158.0 * v8)
                    + v5 * (485446.0 * v5 - 1223952.0 * v6 + 674036.0 * v7 - 140474.0 * v8)
                    + v6 * (812586.0 * v6 - 929952.0 * v7 + 198426.0 * v8)
                    + v7 * (277126.0 * v7 - 121742.0 * v8) + v8 * (13816.0 * v8);
    const double b4 = v5 * (215836.0 * v5 - 1299002.0 * v6 + 1517646.0 * v7 - 822974.0 * v8 + 172658.0 * v9)
                    + v6 * (2041126.0 * v6 - 4924152.0 * v7 + 2716916.0 * v8 - 576014.0 * v9)
                    + v7 * (3042786.0 * v7 - 3408792.0 * v8 + 729726.0 * v9)
                    + v8 * (965926.0 * v8 - 417002.0 * v9) + v9 * (45316.0 * v9);
    const double tau = fabs(b0 + 2.0 * b1 - 6.0 * b2 + 2.0 * b3 + b4);
    const double r0 = tau / (b0 + EPS64), r1 = tau / (b1 + EPS64), r2 = tau / (b2 + EPS64);
    const double r3 = tau / (b3 + EPS64), r4 = tau / (b4 + EPS64);
    const double a0 = (1.0 / 126.0) * (1.0 + r0 * r0), a1 = (10.0 / 63.0) * (1.0 + r1 * r1);
    const double a2 = (10.0 / 21.0) * (1.0 + r2 * r2), a3 = (20.0 / 63.0) * (1.0 + r3 * r3);
    const double a4 = (5.0 / 126.0) * (1.0 + r4 * r4);
    const double invs = 1.0 / (a0 + a1 + a2 + a3 + a4);
    w0 = a0 * invs; w1 = a1 * invs; w2 = a2 * invs; w3 = a3 * invs; w4 = a4 * invs;
}

__device__ __forceinline__ void weights9_64_right(
    double v1, double v2, double v3, double v4, double v5,
    double v6, double v7, double v8, double v9,
    double& w0, double& w1, double& w2, double& w3, double& w4)
{
    const double b4 = v9 * (45316.0 * v9 - 417002.0 * v8 + 729726.0 * v7 - 576014.0 * v6 + 172658.0 * v5)
                    + v8 * (965926.0 * v8 - 3408792.0 * v7 + 2716916.0 * v6 - 822974.0 * v5)
                    + v7 * (3042786.0 * v7 - 4924152.0 * v6 + 1517646.0 * v5)
                    + v6 * (2041126.0 * v6 - 1299002.0 * v5) + v5 * (215836.0 * v5);
    const double b3 = v8 * (13816.0 * v8 - 121742.0 * v7 + 198426.0 * v6 - 140474.0 * v5 + 36158.0 * v4)
                    + v7 * (277126.0 * v7 - 929952.0 * v6 + 674036.0 * v5 - 176594.0 * v4)
                    + v6 * (812586.0 * v6 - 1223952.0 * v5 + 330306.0 * v4)
                    + v5 * (485446.0 * v5 - 280502.0 * v4) + v4 * (45316.0 * v4);
    const double b2 = v7 * (13816.0 * v7 - 102002.0 * v6 + 135846.0 * v5 - 77894.0 * v4 + 16418.0 * v3)
                    + v6 * (209926.0 * v6 - 598152.0 * v5 + 358196.0 * v4 - 77894.0 * v3)
                    + v5 * (462306.0 * v5 - 598152.0 * v4 + 135846.0 * v3)
                    + v4 * (209926.0 * v4 - 102002.0 * v3) + v3 * (13816.0 * v3);
    const double b1 = v6 * (45316.0 * v6 - 280502.0 * v5 + 330306.0 * v4 - 176594.0 * v3 + 36158.0 * v2)
                    + v5 * (485446.0 * v5 - 1223952.0 * v4 + 674036.0 * v3 - 140474.0 * v2)
                    + v4 * (812586.0 * v4 - 929952.0 * v3 + 198426.0 * v2)
                    + v3 * (277126.0 * v3 - 121742.0 * v2) + v2 * (13816.0 * v2);
    const double b0 = v5 * (215836.0 * v5 - 1299002.0 * v4 + 1517646.0 * v3 - 822974.0 * v2 + 172658.0 * v1)
                    + v4 * (2041126.0 * v4 - 4924152.0 * v3 + 2716916.0 * v2 - 576014.0 * v1)
                    + v3 * (3042786.0 * v3 - 3408792.0 * v2 + 729726.0 * v1)
                    + v2 * (965926.0 * v2 - 417002.0 * v1) + v1 * (45316.0 * v1);
    const double tau = fabs(b0 + 2.0 * b1 - 6.0 * b2 + 2.0 * b3 + b4);
    const double r0 = tau / (b0 + EPS64), r1 = tau / (b1 + EPS64), r2 = tau / (b2 + EPS64);
    const double r3 = tau / (b3 + EPS64), r4 = tau / (b4 + EPS64);
    const double a0 = (5.0 / 126.0) * (1.0 + r0 * r0), a1 = (20.0 / 63.0) * (1.0 + r1 * r1);
    const double a2 = (10.0 / 21.0) * (1.0 + r2 * r2), a3 = (10.0 / 63.0) * (1.0 + r3 * r3);
    const double a4 = (1.0 / 126.0) * (1.0 + r4 * r4);
    const double invs = 1.0 / (a0 + a1 + a2 + a3 + a4);
    w0 = a0 * invs; w1 = a1 * invs; w2 = a2 * invs; w3 = a3 * invs; w4 = a4 * invs;
}

}  // namespace w32
