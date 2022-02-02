program readgr_qe
  
  use estimator
  use md_variables, only : nbeadMD,natoms,ndimMD,tempMD,output_directory,lonl,avec,nunitcells, &
                           avecsp,ncellsx,ncellsy,ncellsz
  implicit none

  integer i, i2, b, ix, ig, skip, du,&
       flag_select, done,  i1, ii1,ij1,ik1,iddmax1, ii2,ij2,ik2,iddmax2, &
       n_atom1, n_atom2,  repetition, ibinit, iend, &
       lbin, iflagerr, nbin, ngen, nion, ieskinr

  real(8)  pos_max, Delta, norm, rrr, den, rr, &
       atom1, atom2, r_cut, fant, diff, diffs, &
       cellscale(3), omega,dum, vec(3), vec1(3), vec2(3), &
       avecspd(3,3)
  
  logical :: ifex,repimage
  
  character*7 startpoint,dumat
  character*1 rptf
  real(8), dimension(:,:), allocatable :: pos
  real(8), dimension(:), allocatable :: pos_hist_one, poseq, poseq_hist_one
  type(avg_array) :: pos_hist,pos_histJN
  type(esti_array) :: pos_histbin
  
   
  call input_tools_tom 
  nion=natoms
  ieskinr=ndimMD*natoms
     
  repetition=0
  
  open(unit=12,file=output_directory(1:lonl)//'/configinit.xyz',form='formatted',status='old',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) "Opening configinit.xyz error!"
    stop
  endif
  
  IF (nbeadMD .gt.1 ) then
    open(unit=22,file=output_directory(1:lonl)//'/positions_cen.dat',form='formatted',status='old',iostat=iflagerr)
  ELSE
    open(unit=22,file=output_directory(1:lonl)//'/positions.dat',form='formatted',status='old',iostat=iflagerr)
  END IF
  if(iflagerr.ne.0) then
    write(*,*) "Opening positions_cen.dat error!"
    stop
  endif
 
  ngen=0
  do while(.true.)
    read(22,*,iostat=iflagerr)
    if(iflagerr/=0) exit
    ngen=ngen+1
  enddo
  rewind(22)


  allocate(pos(ngen,nion*ndimMD),poseq(nion*ndimMD))
  
  read(12,*)
  read(12,*) 
  do i1=1,nion
    read(12,*) dumat,poseq(1+(i1-1)*ndimMD:ndimMD+(i1-1)*ndimMD) 
  end do
  poseq = poseq / 0.5291772d0

  do i1=1,ngen
      read(22,*) pos(i1,:)
  end do
  

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
  write(*,*) 'do you want to replicate the supercell?? digit t or f'
  read(5,*) rptf
  if(rptf.eq.'t') repimage=.true.
  if(rptf.eq.'f') repimage=.false.

  if(nbin .lt. 0 .and. repimage .eqv. .false.) then
    nbin=1000
  else if(nbin .lt. 0 .and. repimage .eqv. .true.) then
    nbin=2000
  endif
  
  iddmax1=0; iddmax2=0
  pos_max=25.0
  avecspd = avecsp
  if(repimage .eqv. .true.) then
    iddmax1=1; iddmax2=1
    pos_max=50.0
    avecspd = avecsp * 2.d0
    open(unit=25,file=output_directory(1:lonl)//'/gr_hist_repl.out',form='formatted',status='unknown')
  else
    open(unit=25,file=output_directory(1:lonl)//'/gr_hist_sing.out',form='formatted',status='unknown')
  endif


  if(iend .gt. 0 .and. iend .lt. ngen) ngen=iend

  if(mod(ngen-ibinit+1,lbin).ne.0) then
    ibinit=ibinit+mod(ngen-ibinit+1,lbin)
    write(6,*) ' Warning changing ibinit to match bin length &
    &  of last iterations, new ibinit',ibinit
  endif


  Delta=pos_max/dble(nbin)


  allocate(pos_hist_one(nbin),poseq_hist_one(nbin))

 !!! equilibrium configuration
  poseq_hist_one = 0.0
  
  do ii1=0,iddmax1
  do ij1=0,iddmax1
  do ik1=0,iddmax1
    do i1=1,natoms
       vec1(:) = poseq(1+(i1-1)*ndimMD : ndimMD+(i1-1)*ndimMD ) + &
                 dble(ii1)*avecsp(:,1) + dble(ij1)*avecsp(:,2) + dble(ik1)*avecsp(:,3)

      do ii2=0,iddmax2
      do ij2=0,iddmax2
      do ik2=0,iddmax2
        do i2=1,natoms
           vec2(:) = poseq(1+(i2-1)*ndimMD : ndimMD+(i2-1)*ndimMD ) + &
                     dble(ii2)*avecsp(:,1) + dble(ij2)*avecsp(:,2) + dble(ik2)*avecsp(:,3)
     
          if (ii2.eq.ii1 .and. ij2.eq.ij1 .and. ik2.eq.ik1) then
            if (i2.ne.i1) then
              vec=0.0

              call mimage(vec2(:) - vec1(:) , avecspd, vec )

              ix=int(norm2(vec)/Delta)+1
              poseq_hist_one(ix)=poseq_hist_one(ix)+1.
            end if

          else
              vec=0.0

              call mimage(vec2(:) - vec1(:) , avecspd, vec )

              ix=int(norm2(vec)/Delta)+1
              poseq_hist_one(ix)=poseq_hist_one(ix)+1.

          end if
        end do
      end do
      end do
      end do
    end do
  end do
  end do
  end do
 !!! MD configurations
  
  done=0
  call alloc(pos_hist,nbin)
  call reset(pos_hist)
  
  do ig=ibinit,ngen,skip

    if(mod(ig,100).eq.0) write(*,*) ' Generation Number = ',ig
    done=done+1

    pos_hist_one(:)=0.d0


    do ii1=0,iddmax1
    do ij1=0,iddmax1
    do ik1=0,iddmax1
      do i1=1,natoms
        vec1(:) = pos(ig,1+(i1-1)*ndimMD : ndimMD+(i1-1)*ndimMD ) + &
                  dble(ii1)*avecsp(:,1) + dble(ij1)*avecsp(:,2) + dble(ik1)*avecsp(:,3)

        do ii2=0,iddmax2
        do ij2=0,iddmax2
        do ik2=0,iddmax2
          do i2=1,natoms
            vec2(:) = pos(ig,1+(i2-1)*ndimMD : ndimMD+(i2-1)*ndimMD ) + &
                      dble(ii2)*avecsp(:,1) + dble(ij2)*avecsp(:,2) + dble(ik2)*avecsp(:,3)

            if (ii2.eq.ii1 .and. ij2.eq.ij1 .and. ik2.eq.ik1) then
              if (i2.ne.i1) then
                vec=0.0

                call mimage(vec2(:) - vec1(:) , avecspd, vec )

                ix=int(norm2(vec)/Delta)+1
                pos_hist_one(ix)=pos_hist_one(ix)+1.
              end if

            else
              vec=0.0

              call mimage(vec2(:) - vec1(:) , avecspd, vec )

              ix=int(norm2(vec)/Delta)+1
              pos_hist_one(ix)=pos_hist_one(ix)+1.

            end if
          end do
        end do
        end do
        end do
      end do
    end do
    end do
    end do

    call push(pos_hist, pos_hist_one)
  enddo !main loop ends

  call calc(pos_hist)


  write(*,*) done, 'generations done, norm = ', norm

  write(*,*) ' Begin Jacknife!'
  rewind(22)


  call alloc(pos_histJN,nbin)
  call alloc(pos_histbin,nbin)
  call reset(pos_histJN)
  call reset(pos_histbin)

  done=0

  do ig=ibinit,ngen,skip

    if(mod(ig,100).eq.0) write(*,*) ' Generation Number = ',ig
    done=done+1

    pos_hist_one(:)=0.d0

    do ii1=0,iddmax1
    do ij1=0,iddmax1
    do ik1=0,iddmax1
      do i1=1,natoms
        vec1(:) = pos(ig,1+(i1-1)*ndimMD : ndimMD+(i1-1)*ndimMD ) + &
                  dble(ii1)*avecsp(:,1) + dble(ij1)*avecsp(:,2) + dble(ik1)*avecsp(:,3)

        do ii2=0,iddmax2
        do ij2=0,iddmax2
        do ik2=0,iddmax2
          do i2=1,natoms
            vec2(:) = pos(ig,1+(i2-1)*ndimMD : ndimMD+(i2-1)*ndimMD ) + &
                      dble(ii2)*avecsp(:,1) + dble(ij2)*avecsp(:,2) + dble(ik2)*avecsp(:,3)

            if (ii2.eq.ii1 .and. ij2.eq.ij1 .and. ik2.eq.ik1) then
              if (i2.ne.i1) then
                vec=0.0

                call mimage(vec2(:) - vec1(:) , avecspd, vec )

                ix=int(norm2(vec)/Delta)+1
                pos_hist_one(ix)=pos_hist_one(ix)+1.
              end if

            else
              vec=0.0

              call mimage(vec2(:) - vec1(:) , avecspd, vec )

              ix=int(norm2(vec)/Delta)+1
              pos_hist_one(ix)=pos_hist_one(ix)+1.

            end if
          end do
        end do
        end do
        end do
      end do
    end do
    end do
    end do

    call push(pos_histJN, pos_hist_one)
    if(pos_histJN%num==lbin) then
      call push(pos_histbin,(pos_hist%summation - pos_histJN%summation)/(pos_hist%num - pos_histJN%num))
      call reset(pos_histJN)
    endif
  enddo !main loop ends

  call calc(pos_histbin)
  write(*,*) pos_histbin%num,' Jacknife bins!'

  write(25,*) '#  r, g(r), error'
!  write(25,*) ' Normalization =',sum(temp_hist%average(1:nbin))
  do b=1,nbin
    rrr=(dble(b)-0.5d0)*Delta
    write(25,'(F18.10,3E18.8)') rrr, pos_hist%average(b), &
      & sqrt(pos_histbin%vari(b)*max(pos_histbin%num-1,1)), poseq_hist_one(b)
   
  enddo

  close(25)
  close(22)

  deallocate(pos_hist_one,poseq_hist_one)
  deallocate(pos,poseq)
  call free(pos_hist)
  call free(pos_histJN)
  call free(pos_histbin)

  stop

end program


subroutine mimage(vector,vecsp,my_minimum_image)
    double precision,dimension(3), intent(in) :: vector
    double precision, dimension(3,3),intent(in) :: vecsp(3,3)
    double precision, dimension(3) :: my_minimum_image
    double precision dist_min
    integer i1,i2,i3

    dist_min = norm2(vector)
    my_minimum_image = vector
    do i1 = -2, 2
      do i2 = -2, 2
        do i3 = -2, 2
          if (norm2(i1*vecsp(:,1) + i2*vecsp(:,2) + i3*vecsp(:,3) + vector(:)) .le. dist_min) then
            dist_min = norm2(i1*vecsp(:,1) + i2*vecsp(:,2) + i3*vecsp(:,3) + vector(:))
            my_minimum_image(:) = vector(:) + i1*vecsp(:,1) + i2*vecsp(:,2) + i3*vecsp(:,3)
          end if
        end do
      end do
    end do

    !!write(*,*) 'i_min',i1_min, i2_min, i3_min

end subroutine mimage
