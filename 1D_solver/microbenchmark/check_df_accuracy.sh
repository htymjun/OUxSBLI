#!/usr/bin/env bash
# Accuracy gate for the family G double-float candidate polynomials.
#
# WHY THIS EXISTS, and why summary.csv cannot do its job: init_input builds a
# piecewise-linear field, and every WENO candidate polynomial is exact on linear
# data. A DF implementation whose lo-word arithmetic has been optimised away
# still reproduces the FP64 checksum_all BIT-FOR-BIT. That is not a theoretical
# concern -- it is exactly what happened to poly9_df_* (the fltflt-operator arm),
# which passes the checksum and is only 2.1e-08 accurate.
#
# So this probe evaluates the candidate polynomials on a SMOOTH NON-LINEAR field
# and compares each DF arm against poly9_64_*. It builds a single translation
# unit out of the real weno_micro.f90 module text -- not a copy of the routines
# -- so it cannot drift from what the benchmark actually measures. That matters:
# the failure being guarded against is a whole-translation-unit optimisation
# effect, and it does NOT reproduce when the same routines are lifted into a
# small standalone module.
#
# Run after ANY change to dfr_dot5, poly9_dfr_*, poly9_df_*, the fltflt library,
# or the compiler flags.
#
#   bash 1D_solver/microbenchmark/check_df_accuracy.sh [--gpu-cc 89]
#
# Expected:
#   poly9_dfr_* (relaxed, hand-written)  ~1e-14   PASS
#   poly9_df_*  (fltflt operators)       ~2e-08   known-bad, reported, not a gate
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FC="${FC:-nvfortran}"
GPU_CC="${CASE_GPU_CC:-89}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

while [ $# -gt 0 ]; do
  case "$1" in
    --gpu-cc) GPU_CC="${2:?missing --gpu-cc value}"; shift 2;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done

# The fltflt module object. Reuse the benchmark's build if present, else compile.
FLTFLT_OBJ=""
for cand in "$SCRIPT_DIR"/build*/CMakeFiles/weno_micro.dir/**/fltflt.f90.o; do
  [ -f "$cand" ] && FLTFLT_OBJ="$cand" && break
done
if [ -z "$FLTFLT_OBJ" ]; then
  ( cd "$WORK" && $FC -cuda -fast -Mfma -Mpreprocess \
      -gpu=rdc,lto,cc${GPU_CC},maxregcount:128 \
      -I"$SCRIPT_DIR/../src" -c "$SCRIPT_DIR/../src/fltflt.f90" -o fltflt.o )
  FLTFLT_OBJ="$WORK/fltflt.o"
fi

# Assemble one translation unit: the real module text, plus one probe kernel per
# polynomial arm, plus a driver. One arm per kernel keeps register pressure low
# so that spilling cannot be confused with a precision effect.
awk '/^end module weno_micro_kernels/{
  split("64 df dfr", a, " ")
  for (j=1;j<=3;j++) {
    for (s=1;s<=2;s++) {
      v=a[j]; sd=(s==1)?"left":"right"; sn=(s==1)?"l":"r"
      print "  attributes(global) subroutine kp_" v "_" sn "(n, x, o)"
      print "    integer, value :: n"
      print "    real(8), device :: x(n), o(n,5)"
      print "    integer :: i"
      print "    real(8) :: a0,a1,a2,a3,a4"
      print "    i = (blockIdx%x-1)*blockDim%x + threadIdx%x"
      print "    if (i > n-9) return"
      print "    call poly9_" v "_" sd "(x(i),x(i+1),x(i+2),x(i+3),x(i+4),x(i+5),x(i+6),x(i+7),x(i+8), a0,a1,a2,a3,a4)"
      print "    o(i,1)=a0; o(i,2)=a1; o(i,3)=a2; o(i,4)=a3; o(i,5)=a4"
      print "  end subroutine kp_" v "_" sn
    }
  }
}{print} /^end module weno_micro_kernels/{exit}' "$SCRIPT_DIR/weno_micro.f90" > "$WORK/probe.f90"

cat >> "$WORK/probe.f90" <<'EOF'
program main
  use cudafor
  use weno_micro_kernels
  implicit none
  integer, parameter :: n = 4096
  real(8), allocatable, device :: dx(:), a(:,:), b(:,:), c(:,:)
  real(8) :: hx(n), h64(n,5), hdf(n,5), hdfr(n,5), den, edf, edfr
  integer :: i, side, ierr
  type(dim3) :: g, blk
  logical :: ok
  allocate(dx(n), a(n,5), b(n,5), c(n,5))
  ! smooth but genuinely non-linear: every candidate differs here, unlike the
  ! piecewise-linear init_input the checksum uses
  do i = 1, n
    hx(i) = exp(0.3d0*sin(0.017d0*dble(i))) + 0.5d0*cos(0.0031d0*dble(i)**1.3d0)
  enddo
  dx = hx
  blk = dim3(128,1,1); g = dim3((n+127)/128,1,1)
  ok = .true.
  do side = 1, 2
    if (side == 1) then
      call kp_64_l <<<g,blk>>>(n, dx, a)
      call kp_df_l <<<g,blk>>>(n, dx, b)
      call kp_dfr_l<<<g,blk>>>(n, dx, c)
    else
      call kp_64_r <<<g,blk>>>(n, dx, a)
      call kp_df_r <<<g,blk>>>(n, dx, b)
      call kp_dfr_r<<<g,blk>>>(n, dx, c)
    endif
    ierr = cudaGetLastError()
    if (ierr /= 0) then
      print *, ' launch failed: ', cudaGetErrorString(ierr)
      error stop 3
    endif
    ierr = cudaDeviceSynchronize()
    h64 = a; hdf = b; hdfr = c
    den  = maxval(abs(h64(1:n-9,:)))
    edf  = maxval(abs(hdf (1:n-9,:) - h64(1:n-9,:))) / den
    edfr = maxval(abs(hdfr(1:n-9,:) - h64(1:n-9,:))) / den
    if (side == 1) then
      print "(a)", '  left-biased candidates (v-):'
    else
      print "(a)", '  right-biased candidates (v+):'
    endif
    print "(a,es11.4,a)", '    poly9_dfr_* relaxed, hand-written : ', edfr, '   (gate: < 1e-12)'
    print "(a,es11.4,a)", '    poly9_df_*  fltflt operators      : ', edf,  '   (known-bad, not gated)'
    if (edfr >= 1.d-12) ok = .false.
  enddo
  print "(a)", ''
  if (ok) then
    print "(a)", '  -> PASS: the relaxed DF candidates are genuine double-float on both biases'
  else
    print "(a)", '  -> FAIL: relaxed DF lost its lo word. The error-free transforms were'
    print "(a)", '           compiled away -- check that dfr_dot5 still uses __fadd_rn for'
    print "(a)", '           every add/sub, and that the input split is written inline'
    print "(a)", '           rather than via fltflt_init (see poly9_df_left comments).'
  endif
  deallocate(dx, a, b, c)
  if (.not. ok) error stop 1
end program main
EOF

echo "target architecture: cc${GPU_CC}"
( cd "$WORK" && $FC -cuda -fast -Mfma -Mpreprocess \
    -gpu=rdc,lto,cc${GPU_CC},maxregcount:128 \
    -I"$SCRIPT_DIR/../src" probe.f90 "$FLTFLT_OBJ" -o probe 2>&1 \
  | grep -vE 'NVFORTRAN-W-0473|^probe\.f90:$|^\s*[0-9]+ inform|GNU-stack|deprecated' || true )

"$WORK/probe"
