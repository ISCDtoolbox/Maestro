subroutine readc2 
!    *******************************************************************
!    ** subroutine to read in the configuration from unit 30          **
!    *******************************************************************
  use md_variables
  implicit none
  integer :: i,l,k
  open(unit=30,file='input_files_dir/configinit')
  write(*,*) 'reading the initial configuration file'
  write(unit_dot_out,*) 'reading the initial configuration file'
      
  do k=1,nbeadMD 
    do i=1,n  
      read (30,*) (rpos(l,i,k),l=1,ndimMD)
    enddo
    rewind(30)
  enddo

  close (30)
  
  rpos_init(:,:,:)=rpos(:,:,:)
  
  return

end subroutine
