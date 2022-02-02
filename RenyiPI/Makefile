# GNU
#F77  = /usr/local/pgi/linux86/bin/pgf77 -O3 -r8
#F77 = g77 -O3 -w
#LIB = -L/sw/lib/ -llapack
#F77 = gfortran -g  -w

# XLF   
#F77 = /opt/ibmcmp/xlf/8.1/bin/f77 -O2 -WF,-DOSX -qmaxmem=-1 -qdpc -qextname 
#F77 = /opt/ibmcmp/xlf/8.1/bin/f77 -g -C -qdpc -qmaxmem=-1 #-qextname 
#F77 = /opt/ibmcmp/xlf/8.1/bin/f77 -O3 -WF,-DOSX -qmaxmem=-1 -qdpc #-qextname 
#XTRALIB = -Wl,-framework -Wl,Accelerate 

# Linux INTEL
#F77 = /opt/mpich-1.2.6/intel/bin/mpif77 -g -C -w -r8
#F77 = /opt/mpich-1.2.6/intel/bin/mpif77 -O3 -w -r8
#F77 = /opt/mpich-1.2.6/intel/bin/mpif77 -O3 -r8
F90 = gfortran
F77 = gfortran 

# DEBUGGING OPTIONS
#FLAGS = -O3 -w -ffixed-form
#FLAGS90 =  -cpp -g -C -ffree-form -fbounds-check -fopenmp -finit-local-zero -Wall -Wextra -fbacktrace -ffree-line-length-none
# OPTIMIZED OPTIONS
FLAGS = -O3 -ffixed-form -ffast-math -ftree-vectorize -funroll-loops -fbounds-check
FLAGS90 = -O3 -cpp -ffree-form -ffree-line-length-none -ffast-math -ftree-vectorize -funroll-loops -fbounds-check \
          -fopenmp -finit-local-zero


# blas/lapack linking
#ATLAS
#LINK = -framework Accelerate -bind_at_load
#LINK =  ~/lib/lapack-3.5.0/liblapack.a \
#this one above was commented before
#~/lib/lapack-3.5.0/librefblas.a \
#LINK=-L/home/tommaso/Desktop/libraries -lrefblas -llapack
#LINK = -lifport -lifcore -limf -lm -lirc -lirc_s -ldl \
#-L/home/casula/intel/mkl8.0/lib/32 -lmkl_lapack -lguide -lmkl_ia32

# fftw3 linking
INCLUDE=-I/usr/include/
#INCLUDE=-I/usr/share/doc/packages/fftw3-devel
#INCLUDE=-I/usr/local/include/
FFTFLAGS=-L/usr/lib/x86_64-linux-gnu -lfftw3
FFTFLAGS=-L/usr/share/doc/packages/fftw3-devel - lfftw3
FFTFLAGS=-L/usr/lib64/ruby/2.1.0/x86_64-linux-gnu -lfftw3
#LINK=-L/home/tommaso/Desktop/libraries/lapack-master -llapack -lrefblas
LINK=-L/srv/nfs4/opt/matlab/R2011a/bin/glnxa6i4/ -llapack -lblas
#this one above vas uncommented before
#LINK=-L/home/tommaso/Desktop/libraries -llapack -lrefblas
#FFTFLAGS=-L/usr/local/lib -lfftw3

# general linking
LDFLAGS=

# include flags
INC = $(INCLUDE)

# main program
MD= scnd.o md_variables.o input_tom.o cumul.o cumul1.o potential.o checkpoint.o \
    sprint.o zeroav.o readc2.o mdnvt_tom.o init_ceriotti_tom.o dynmatrix.o pot.o \
    setvel.o force0.o propNH.o propP.o propX.o prop_ceriotti_tom.o blockMD_tom.o init_pioud_tom.o \
    reweight_pioud_tom.o set_turboq.o root2mat.o dsyev2.o dscalzero.o  \
    dealloc.o vel_verlet_tom.o prop_pioud_tom.o 

# tools
PHONON_FD = md_variables.o potential.o minmz_and_phonons_fd.o input_tools_tom.o 

AUTOCORR = md_variables.o kinds.o fft_methods.o correlations.o autocorr.o

DIFFUSION = md_variables.o diffusion.o

GR = md_variables.o estimator.o readgr_qe.o input_tools_tom.o

TEMP = md_variables.o estimator.o readtemp.o input_tools_tom.o

.SUFFIXES: .f90 .F .x .o

.F.o:   
	$(F77) $(INC) -c $?

.f.o:   
	$(F77) $(FLAGS) -c $?

.f90.o: 
	$(F90) $(FLAGS90) $(INC) -c $?

all : pilmdt.x phfd.x post_autocorr.x post_diffusion.x post_gr.x post_temp.x

pilmdt.x:$(MD)
	$(F90) $(FLAGS90) $(MD) -o $@ $(LDFLAGS) \
	$(LINK)

phfd.x:$(PHONON_FD)
	$(F90) $(FLAGS90) $(PHONON_FD) -o $@ $(LDFLAGS) \
	$(LINK)

post_autocorr.x:$(AUTOCORR)
	$(F90) $(FLAGS90) $(AUTOCORR) -o $@ $(FFTFLAGS) $(LDFLAGS) \
	$(LINK)

post_diffusion.x:$(DIFFUSION)
	$(F90) $(FLAGS90) $(DIFFUSION) -o $@ $(LDFLAGS) \
	$(LINK)

post_gr.x:$(GR)
	$(F90) $(FLAGS90) $(GR) -o $@ $(FFTFLAGS) $(LDFLAGS) \
	$(LINK)

post_temp.x:$(TEMP)
	$(F90) $(FLAGS90) $(TEMP) -o $@ $(FFTFLAGS) $(LDFLAGS) \
	$(LINK)

clean:
	rm -f *.o *.mod *~

cleanall:
	rm -f *.o *.mod *.x
