FC = mpif90

OUTPUT_OPTIN = -cuda -Xptxas -warn-spills -acc -fast -Minfo=accel -lnvhpcwrapnvtx

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = SBI8
vpath %f90 ../src:../../src

OBJ =	mod_globals.o \
      set_coordinate.o \
      set.o \
      calc_common_dim.o \
      calc_mat.o \
      calc_physical_quantities.o \
      calc_sutherland.o \
      calc_visc.o \
      calc_common.o \
      calc_hybrid.o \
      calc_slau.o \
      calc_term.o \
      calc_keep.o \
      calc_muscl.o \
      calc_flux.o \
			calc_les.o \
      calc_steps.o \
      calc_rescale.o \
			calc_para.o \
      print.o \
      calc_time_dev.o \
      main.o

a.out: $(OBJ)
		$(FC) -cuda -Xptxas -warn-spills -acc -fast -Minfo -lnvhpcwrapnvtx -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set.o: mod_globals.mod
calc_physical_quantities.o: mod_globals.mod
calc_common.o: mod_globals.mod
calc_term.o: calc_common_dim.mod
calc_keep.o: mod_globals.mod calc_common_dim.mod calc_term.mod calc_mat.mod
calc_visc.o: mod_globals.mod calc_sutherland.mod
calc_mat.o: mod_globals.mod calc_common_dim.mod
calc_slau.o: mod_globals.mod calc_common.mod calc_hybrid.mod
calc_roe.o: mod_globals.mod calc_common.mod calc_common_dim.mod calc_physical_quantities.mod calc_mat.mod
calc_muscl.o: mod_globals.mod
calc_flux.o: mod_globals.mod calc_keep.mod calc_slau.mod calc_roe.mod calc_muscl.mod
calc_hybrid.o: mod_globals.mod
calc_steps.o: mod_globals.mod
calc_rescale.o: mod_globals.mod
calc_time_dev.o: mod_globals.mod calc_physical_quantities.mod set.mod calc_flux.mod calc_hybrid.mod calc_visc.mod calc_steps.mod calc_les.mod calc_rescale.mod calc_para.mod print.mod
print.o: mod_globals.mod
main.o: main.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

