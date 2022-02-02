program readtemp
  
  use estimator
  use md_variables, only : natoms,ndimMD,tempMD,output_directory,lonl
  implicit none

  integer i, j, b, ix, ig, skip, du,&
       flag_select, done, icut, i1, &
       n_atom1, n_atom2,  repetition, ibinit, iend, &
       lbin, iflagerr, nbin, ngen, nion, ieskinr, d

  real(8) tmp(5), temp_max, Delta, norm, rrr, den, rr, &
       atom1, atom2, r_cut, fant, diff,diffs, &
       cellscale(3), omega,dum
  
  logical :: ifex
  
  character*7 startpoint
  real(8), dimension(:), allocatable :: temperatura
  real(8), dimension(:), allocatable :: temp_hist_one
  type(avg_array) :: temp_hist,temp_histJN
  type(esti_array) :: temp_histbin
  
   
  call input_tools_tom 
  nion=natoms
  ieskinr=ndimMD*natoms
  
  repetition=0
  
  temp_max=5.d0*tempMD
  
! get nion, atomic_weight
! allocate rion
! allocate my_iond(nion,nion),iond_cart(3,nion,nion)
  
  open(unit=22,file=output_directory(1:lonl)//'/local_temp.dat',form='formatted',status='old',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) "Opening local_temp.dat error!"
    stop
  endif
 
  ngen=0
  do while(.true.)
    read(22,*,iostat=iflagerr)
    if(iflagerr/=0) exit
    ngen=ngen+1
  enddo
  rewind(22)


  allocate(temperatura(ngen))
  
  do i1=1,ngen
      read(22,*) temperatura(i1)
  end do
  

  open(unit=25,file=output_directory(1:lonl)//'/temp_hist.out',form='formatted',status='unknown')

  !count the number of total available snapshots

  write(*,*)
  write(*,*) ' Number of generation = ',ngen

  ! all default options
  iflagerr=0

  omega=1.d0

  
  write(*,*) ' Histogram points: (if < 0 set default = 1000;'
  write(*,*) '                    if > 0 normal) '
  read(5,*) nbin
  write(*,*) 'bin lenght, ibinit, iend (if <= 0 ngen),  skip:'
  read(5,*) lbin, ibinit,iend,skip

  if(nbin .lt. 0) then
    nbin=1000
  endif


  if(iend .gt. 0 .and. iend .lt. ngen) ngen=iend

  if(mod(ngen-ibinit+1,lbin).ne.0) then
    ibinit=ibinit+mod(ngen-ibinit+1,lbin)
    write(6,*) ' Warning changing ibinit to match bin length &
    &  of last iterations, new ibinit',ibinit
  endif


  Delta=temp_max/dble(nbin)

  icut = int(r_cut/Delta)

  !all the gr allocation here
  allocate(temp_hist_one(nbin))
  call alloc(temp_hist,nbin)

  !skip the initial ibinit-1 records
  do i=1,ibinit-1
    read(22,*,iostat=iflagerr)    
  enddo
  
  done=0

  !main loop for all records
  do ig=ibinit,ngen,skip

    if(mod(ig,100).eq.0) write(*,*) ' Generation Number = ',ig
    done=done+1

    temp_hist_one(:)=0.d0
    ix=int(temperatura(ig)/Delta)+1
    temp_hist_one(ix)=temp_hist_one(ix)+1.
   
    call push(temp_hist, temp_hist_one)
  enddo !main loop ends

  call calc(temp_hist)


  write(*,*) done, 'generations done, norm = ', norm

  write(*,*) ' Begin Jacknife!'
  rewind(22)

  !skip the initial ibinit-1 records
  do i=1,ibinit-1
    read(22,*,iostat=iflagerr)
    if(iflagerr/=0) exit
  enddo


  call alloc(temp_histJN,nbin)
  call alloc(temp_histbin,nbin)

  done=0

  do ig=ibinit,ngen,skip

    write(*,*) ' Generation Number = ',ig
    done=done+1


    temp_hist_one(:)=0.d0
    ix=int(temperatura(ig)/Delta)+1
    !write(*,*) ix,temperatura(ig)
    temp_hist_one(ix)=temp_hist_one(ix)+1.

    call push(temp_histJN, temp_hist_one)
    if(temp_histJN%num==lbin) then
      call push(temp_histbin,(temp_hist%summation-temp_histJN%summation)/(temp_hist%num-temp_histJN%num))
      call reset(temp_histJN)
    endif
  enddo !main loop ends

  call calc(temp_histbin)
  write(*,*) temp_histbin%num,' Jacknife bins!'

  write(25,*) '#  Temperature, temp. distribution, error'
!  write(25,*) ' Normalization =',sum(temp_hist%average(1:nbin))
  do b=1,nbin
    rrr=(dble(b)-0.5d0)*Delta
    write(25,'(F18.10,2E18.8)') rrr, temp_hist%average(b), &
      & sqrt(temp_histbin%vari(b)*max(temp_histbin%num-1,1))
   
  enddo

  close(25)
  close(22)

  deallocate(temp_hist_one)
  deallocate(temperatura)
  call free(temp_hist)
  call free(temp_histJN)
  call free(temp_histbin)

  stop

end program
