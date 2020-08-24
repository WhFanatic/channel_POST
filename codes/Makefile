
FC = mpiifort
OPT = -O3

EXEC = CHANNEL_WHN

FILELIST = global.F90 \
CHANNEL_000_MAIN.F90 \
CHANNEL_010_INIT.F90 \
CHANNEL_011_PREMODIFYP.F90 \
CHANNEL_011_PREMODIFYV.F90 \
CHANNEL_012_INITFIELD.F90 \
CHANNEL_020_VOR.F90 \
CHANNEL_021_FP.F90 \
CHANNEL_022_BOUNDARYP.F90 \
CHANNEL_023_UTP.F90 \
CHANNEL_024_FR.F90 \
CHANNEL_030_NONLINEAR90.F90 \
CHANNEL_040_FIXFR.F90 \
CHANNEL_041_DISPLACE.F90 \
CHANNEL_042_WRITEFILE.F90 \
CHANNEL_050_CHEBTRANS.F90 \
CHANNEL_056_MODIFYDIV.F90 \
CHANNEL_071_TRANSPOSE_2.F90 \
CHANNEL_072_DIAGONAL.F90 \
CHANNEL_073_PRENEW.F90 \
CHANNEL_074_VISNEW.F90 \
CHANNEL_075_DIV.F90

FFTLIB = -DMKL_SINGLE \
-I"/apps/intel/Compiler/11.1/069/f/mkl"/include \
-I"/apps/intel/Compiler/11.1/069/f/mkl"/include/fftw

MKLLIB = -I"/apps/intel/Compiler/11.1/069/f/mkl"/include \
-I"/apps/intel/Compiler/11.1/069/f/mkl"/include/fftw\
-L/work1/cuigx1_work/zhouzisong/fftw/fftw2xf/lib/em64t -lfftw2xf_intel \
-L"/apps/intel/Compiler/11.1/069/f/mkl/lib/em64t" \
"/apps/intel/Compiler/11.1/069/f/mkl/lib/em64t"/libmkl_intel_lp64.so \
-Wl,--start-group "/apps/intel/Compiler/11.1/069/f/mkl/lib/em64t"/libmkl_intel_thread.so \
"/apps/intel/Compiler/11.1/069/f/mkl/lib/em64t"/libmkl_core.so \
-Wl,--end-group -L"/apps/intel/Compiler/11.1/069/f/mkl/lib/em64t" \
-liomp5 -lpthread -lm 

OBJS = $(FILELIST:.f90=.o)

%.o:%.F90
	$(FC) $(OPT) $(FFTLIB) -c  $<

default: $(OBJS)
	$(FC) -O3 $(FFTLIB) $(FILELIST) $(MKLLIB) -o $(EXEC)

clean:
	rm -f *.out *.o work.* *.pc *.mod

