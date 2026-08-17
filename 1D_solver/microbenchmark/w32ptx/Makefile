# Makefile for the extracted w32_poly64 PTX benchmark.
#
# Deliberately NOT CMake. The Fortran side's CASE_GPU_CC is a cached variable
# that defaults to `native`, and a stale cached value produces a clean build
# that then fails at EVERY launch with cudaErrorInvalidPtx (218). Here ARCH is
# a required argument with no default, so that failure mode cannot occur.

NVCC ?= nvcc
ARCH ?=

ifeq ($(strip $(ARCH)),)
$(error ARCH is required. Use `make ARCH=89` (RTX 4060) or `make ARCH=80` (A100). \
There is no default on purpose: a wrong architecture builds cleanly and then fails \
at every launch with cudaErrorInvalidPtx)
endif

# -prec-div=true (the default, stated explicitly) is what reproduces
# nvfortran -fast's IEEE div.rn.f32 with its FCHK range test and out-of-line
# slow-path CALL. That sequence IS the baseline being optimised away, so
# -use_fast_math must never appear here: it would silently turn the _ref rung
# into the _rcp rung and the whole ladder would measure nothing. `make regs`
# asserts the FCHK count for exactly this reason.
#
# -ftz=false matches the reference, where FADD/FMUL/FFMA carry no .FTZ
# (measured: 156 plain FADD vs 10 FADD.FTZ, plus 30 FMNMX.FTZ -- the residual
# is a modifier-only difference, recorded rather than chased).
#
# -maxrregcount=128 matches the Fortran side's -gpu=maxregcount:128, raised
# from 96 because at 96 the WENO9 kernels pinned to the cap and spilled 64-84 B,
# injecting LDL/STL into the very instruction stream being measured.
NVCCFLAGS := -O3 -std=c++17 -arch=sm_$(ARCH) -maxrregcount=128 \
             --fmad=true -prec-div=true -prec-sqrt=true -ftz=false \
             -lineinfo -DW32PTX_ARCH=$(ARCH) \
             -Xptxas -v -Xptxas -warn-spills

SRC  := w32ptx.cu
DEPS := weno_w32.cuh kernels.cuh modes.inc
BIN  := w32ptx

.PHONY: all ptx sass regs accuracy clean

all: $(BIN)

$(BIN): $(SRC) $(DEPS)
	@echo "=== building for sm_$(ARCH) ==="
	$(NVCC) $(NVCCFLAGS) -o $@ $(SRC) 2>&1 | tee build.log
	@echo "=== register / spill summary ==="
	@grep -E 'ptxas info|registers|spill' build.log || true
	@# Match only a NONZERO spill: ptxas always prints "0 bytes spill stores".
	@if grep -qE '[1-9][0-9]* bytes spill' build.log; then \
	  echo "ERROR: register spill detected -- LDL/STL in the measured stream invalidates every timing."; \
	  grep -E '[1-9][0-9]* bytes spill' build.log; \
	  exit 1; \
	fi
	@echo "no spill"

# NVVM's output, before ptxas touches it. The only way to confirm what an asm
# block actually became at PTX level -- and to catch a decimal float literal
# that should have been written 0f........
ptx: $(SRC) $(DEPS)
	$(NVCC) $(NVCCFLAGS) -ptx -o w32ptx.ptx $(SRC)
	@echo "wrote w32ptx.ptx"
	@echo "--- decimal .f32 literals (should be empty) ---"
	@grep -nE '\.f32[^;]*[0-9]\.[0-9]' w32ptx.ptx || echo "  (none)"

sass: $(BIN)
	cuobjdump -sass $(BIN) > w32ptx.sass
	@echo "wrote w32ptx.sass"

# The gate that catches -use_fast_math contamination and non-inlined helpers.
# _ref must retain its division CALLs; _rcp must have none.
regs: sass
	@echo "mode                            FCHK  CALL  MUFU   LDG   STG  LDL  STL"
	@# EXACT match on the symbol, never a substring: `_rcp` is a prefix of
	@# `_rcpn`, so a regex match silently concatenates the two kernels and every
	@# count comes out as their sum (LDG 36 instead of 18, STG 12 instead of 6).
	@for k in $$(./$(BIN) --list-kernels | cut -f2); do \
	  awk -v K="$$k" '/Function :/{f=($$NF==K)} f' w32ptx.sass > .k.sass; \
	  printf '%-30s' "$$k"; \
	  for m in FCHK CALL MUFU LDG STG LDL STL; do \
	    printf '%6s' "$$(grep -c "\\b$$m" .k.sass)"; \
	  done; echo; \
	done; rm -f .k.sass

accuracy: $(BIN)
	./$(BIN) accuracy

clean:
	rm -f $(BIN) build.log w32ptx.ptx w32ptx.sass .k.sass
