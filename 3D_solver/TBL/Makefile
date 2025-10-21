FC = mpif90

OUTPUT_OPTIN = -cpp -cuda -acc -fast -Minfo=accel -Mchkptr -gpu=ptxinfo -lnvhpcwrapnvtx

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = shock_tube
vpath %f90 ../src:../../src

OBJ =	mod_globals.o \
      mod_constant.o \
      set_coordinate.o \
      calc_physical_quantities.o \
      calc_rand.o \
      store_shared.o \
      calc_visc.o \
      calc_les.o \
      calc_hybrid.o \
      calc_slau.o \
      calc_keep.o \
      calc_roe.o \
      calc_muscl.o \
      calc_flux.o \
      calc_flux_base.o \
      calc_steps.o \
      calc_rescale.o \
      calc_para.o \
      set_bc_common.o \
      set_bc_tbl_sbli.o \
      set_init_common.o \
      set.o \
      print.o \
      calc_time_dev.o \
      main.o

a.out: $(OBJ)
		$(FC) -cuda -acc -fast -Minfo -lnvhpcwrapnvtx -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

mod_constant.o: mod_globals.mod
set_bc_common.o: mod_globals.mod
set_bc_tbl_sbli.o: mod_globals.mod mod_constant.mod
set_init_common.o: mod_globals.mod mod_constant.mod
set.o: mod_globals.mod mod_constant.mod set_bc_common.mod set_init_common.mod set_coordinate.mod calc_para.mod
calc_physical_quantities.o: mod_globals.mod mod_constant.mod
calc_keep.o: mod_globals.mod mod_constant.mod
calc_roe.o: mod_constant.mod 
calc_visc_common.o: mod_constant.mod
store_shared.o: mod_globals.mod
calc_visc.o: mod_globals.mod mod_constant.mod store_shared.mod
calc_slau.o: mod_globals.mod mod_constant.mod calc_hybrid.mod
calc_muscl.o: mod_globals.mod
calc_flux.o: mod_globals.mod mod_constant.mod calc_keep.mod calc_slau.mod calc_roe.mod calc_hybrid.mod calc_muscl.mod
calc_flux_base.o: mod_globals.mod calc_physical_quantities.mod calc_hybrid.mod calc_flux.mod calc_visc.mod calc_les.mod set.mod
calc_hybrid.o: mod_globals.mod
calc_steps.o: mod_globals.mod
calc_rescale.o: mod_globals.mod mod_constant.mod
calc_time_dev.o: mod_globals.mod set.mod calc_flux_base.mod calc_steps.mod calc_rescale.mod calc_para.mod print.mod
print.o: mod_globals.mod
main.o: main.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

clean:
		rm -f *.mod *.o *.out


