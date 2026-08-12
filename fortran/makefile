# Build the Hornbostel DLCQ QCD(1+1) solver.
#
# -mcmodel=medium is MANDATORY, not an optimization: the COMMON blocks total
# ~2.87 GB of static BSS (HAM0(3,6902,6902), HAM, HNU, HNORM, plus the main
# program's local Z(6902,6902)).  Without it the link fails with R_X86_64_PC32
# relocation overflows.
#
# The original makefile used `f77`, which on this system is a symlink to
# gfortran-13.  Naming gfortran directly makes that explicit and portable.
#
# The original also ran `strip` on the binary.  That is dropped so the
# executable stays debuggable; it costs nothing at runtime.

FC      = gfortran
FFLAGS  = -mcmodel=medium

.PHONY: all clean wide

all: qcdf

qcdf: qcdf.f
	$(FC) $(FFLAGS) -o $@ $<

# Variant with the colour-array caps raised so high Fock sectors are not
# corrupted; see docs/fortran-color-overflow.md. Generated, never edited.
wide: qcdf_wide

qcdf_wide.f: qcdf.f widen.py
	python3 widen.py

qcdf_wide: qcdf_wide.f
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f qcdf qcdf_wide qcdf_wide.f
