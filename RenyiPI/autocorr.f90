program read_autocorr

  use md_variables, only: filen, lfnl, nbeadMD
  use correlations
  use fft_methods

  implicit none 

  logical :: ifex
  integer :: i,j,k,icount,ii,nw
  real(8)  :: c0, c_infinity, delta_t
  integer  :: iskipr,istepinit,istepmax, ndim, nav, ntf, total_time_steps,countj
  real(8), dimension(:), ALLOCATABLE :: forces,tbuf
  real(8), dimension(:,:), allocatable:: pbuf
  real(8), dimension(:), allocatable :: vec_khintchine, corr_function
  complex(8), dimension(:), allocatable :: spectrum
  logical :: yesquantum

  istepmax=0
  write(6,*) '  istepinit ?'  
  read(5,*) istepinit 
  if(istepinit.lt.0) then
     write(6,*) ' Max generation from command line '
     read(5,*) istepmax
     istepinit=-istepinit
  endif


  yesquantum=.true.

  open(222,file='prefix.input',status='old',form='formatted')

  read(222,*) filen
  lfnl=index(filen,' ')-1
  write(6,*) 'name of the run: ',filen(1:lfnl)
  close(222)

! does this file exist?
  inquire(file=filen(1:lfnl)//'.rs',exist=ifex)
  if(ifex) then
     open(unit=100,file=filen(1:lfnl)//'.rs',status='old',form ='formatted')
     write(6,*) 'reading .rs file'
     read (100,*) nbeadMD
     if(nbeadMD.eq.1) yesquantum=.false.
     close(100)

  else
     write(6,*) 'checkpoint file not found!'
     stop

  endif

  
  open(23,file='sigma.dat',status='old')


  if(yesquantum) then 
     ndim=5
     nav=4
     write(6,*) 'Quantum molecular dynamics'
  else
     ndim=3
     nav=1
     write(6,*) 'Classical molecular dynamics'
  endif
  allocate(forces(ndim))

  if(istepmax.eq.0) then 
       
     i=0
     icount = 0
     do while(i.ge.0)
        read(23,*,end=100)   
        icount=icount+1
     enddo
100  continue
     istepmax=icount
  else
     icount=istepmax
  endif
      
  write(*,*) 'num lines files, tot step statistic =',icount , icount-istepinit+1
  if(icount-istepinit+1.le.0) then 
     write(*,*) 'ERROR istepinit > icount'
     stop
  endif

  total_time_steps=icount-istepinit+1
  ntf=2*total_time_steps

  write(*,*) total_time_steps,'total_time_steps'


  allocate(pbuf(nav,total_time_steps),tbuf(nav))

  pbuf=0.d0
  tbuf=0.d0

  rewind(23)

  countj=0
  icount=0

  do while(icount.lt.istepmax)

     read(23,*,end=200) (forces(ii),ii=1,ndim)
	
     icount=icount+1
     
     if(icount.ge.istepinit) then 

        countj=countj+1

        pbuf(1,countj) = forces(1)  ! total energy

        if(yesquantum) then 
           pbuf(2,countj)=forces(1)-forces(4) ! potential energy 
           pbuf(3,countj)=forces(4)   ! Kin energy virial 
           pbuf(4,countj)=forces(5)  ! Kin energy primitive
        endif

     endif
     
  enddo
200 continue     

  close(23)
  

! compute autocorrelation time

  allocate(vec_khintchine(ntf),spectrum(ntf/2+1),corr_function(ntf))
  vec_khintchine=0.d0
  spectrum=cmplx(0.d0,0.d0)

  do ii=1,nav

     do i=1,total_time_steps
        vec_khintchine(2*i-1)=pbuf(ii,i)
     enddo

     call khintchine_array(vec_khintchine,spectrum,ntf,1)

     call fft_bw_vector(spectrum,corr_function,ntf)

!normalization
! c0 = c(1)
     c0=corr_function(1)

! c_infinity taken as average value over the second-half interval
     c_infinity=0.d0
     delta_t=total_time_steps/2-total_time_steps/4+1
     do i=total_time_steps/4,total_time_steps/2
        c_infinity=c_infinity+corr_function(2*i-1)/delta_t
     enddo

     do i=1,total_time_steps
        pbuf(ii,i)=(corr_function(2*i-1)-c_infinity)/(c0-c_infinity)
     enddo

! time integration (in time step units!!)
! to avoid nasty fluctuations stop the integration when integrand becomes negative
     tbuf(ii)=0.d0
     i=1
     do while (i.le.total_time_steps.and.pbuf(ii,i).ge.0.d0)
        tbuf(ii)=tbuf(ii)+pbuf(ii,i)
        i=i+1
     enddo

  enddo

  open(100,file='corr_function.dat',form='formatted',status='unknown')

  do i=1,total_time_steps
     write(100,*) i-1,(pbuf(ii,i),ii=1,nav)
  enddo
  
  close(100)


  if(yesquantum) then
     write(6,*) 'autocorrelation time (in time step units)'
     write(6,*) 'total energy autocorrelation',tbuf(1)
     write(6,*) 'potential energy autocorrelation',tbuf(2)
     write(6,*) 'virial kinetic energy autocorrelation',tbuf(3)
     write(6,*) 'primitive kinetic energy autocorrelation',tbuf(4)
  else
     write(6,*) 'autocorrelation time (in time step units)'
     write(6,*) 'potential energy autocorrelation',tbuf(1)
  endif
     

  deallocate(pbuf,forces,tbuf)
  deallocate(vec_khintchine,spectrum,corr_function)
   
end program read_autocorr
