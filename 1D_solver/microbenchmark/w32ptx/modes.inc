// modes.inc -- the single source of truth for every mode.
//
// Included multiple times by w32ptx.cu with MODE() redefined each time, so the
// kernel definitions, the launch wrappers, the dispatch table and the
// --list-kernels output can never drift apart. (In the Fortran harness the
// mode->kernel map lives separately in analyze_weno_static.py and has to be
// maintained by hand.)
//
//   MODE(name, R, RecipPolicy, Body)
//     R = 3 -> WENO5-Z,  4 -> WENO7-Z,  5 -> WENO9-Z
//
// Naming: <family>_seq<width>_<rung>. The Fortran equivalent of
// w32_poly64_seq5_ref is `w32_poly64_seq` -- no 5 suffix. Keep that alias in
// mind when joining the two summary.csv files.

// ---- the kernel under study: FP32 weights x FP64 polynomials ----
// _ref     : IEEE div.rn.f32, what nvfortran -fast emits. The baseline.
// _rcp     : div.approx.f32 ratios + Newton-refined normalisation. PRIMARY.
// _divfull : div.full.f32 everywhere -- branch-free, ~same instruction count.
//            The attribution control between _ref and _rcp.
// _rcpn    : Newton-refined everywhere, near-IEEE. Fallback rung.
MODE(w32_poly64_seq5_ref,     3, RecipRef,   BodyW32Poly64)
MODE(w32_poly64_seq5_rcp,     3, RecipRcp,   BodyW32Poly64)
MODE(w32_poly64_seq5_divfull, 3, RecipFull,  BodyW32Poly64)
MODE(w32_poly64_seq5_rcpn,    3, RecipRcpN,  BodyW32Poly64)

MODE(w32_poly64_seq7_ref,     4, RecipRef,   BodyW32Poly64)
MODE(w32_poly64_seq7_rcp,     4, RecipRcp,   BodyW32Poly64)
MODE(w32_poly64_seq7_divfull, 4, RecipFull,  BodyW32Poly64)
MODE(w32_poly64_seq7_rcpn,    4, RecipRcpN,  BodyW32Poly64)

MODE(w32_poly64_seq9_ref,     5, RecipRef,   BodyW32Poly64)
MODE(w32_poly64_seq9_rcp,     5, RecipRcp,   BodyW32Poly64)
MODE(w32_poly64_seq9_divfull, 5, RecipFull,  BodyW32Poly64)
MODE(w32_poly64_seq9_rcpn,    5, RecipRcpN,  BodyW32Poly64)

// ---- ablation: the FP32 weight stream alone (combined in FP32) ----
// t(w32_poly64_seqX_Y) - t(w32_only_seqX_Y) is the marginal cost of the FP64
// half at rung Y. Both terms are far above the DRAM floor, which is what makes
// this usable where poly64_only is not.
MODE(w32_only_seq5_ref,       3, RecipRef,   BodyW32Only)
MODE(w32_only_seq5_rcp,       3, RecipRcp,   BodyW32Only)
MODE(w32_only_seq7_ref,       4, RecipRef,   BodyW32Only)
MODE(w32_only_seq7_rcp,       4, RecipRcp,   BodyW32Only)
MODE(w32_only_seq9_ref,       5, RecipRef,   BodyW32Only)
MODE(w32_only_seq9_rcp,       5, RecipRcp,   BodyW32Only)

// ---- accuracy reference: full FP64 weights x FP64 polynomials ----
// Never quote as a timing baseline; it exists to be the truth the FP32 splits
// are measured against in --accuracy.
MODE(weno64_seq5,             3, RecipRef,   BodyWeno64)
MODE(weno64_seq7,             4, RecipRef,   BodyWeno64)
MODE(weno64_seq9,             5, RecipRef,   BodyWeno64)
