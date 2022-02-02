
module fft_methods

  use kinds
  
  implicit none

contains

  subroutine fft_fw_array(input,output,n,m)
    
    include 'fftw3.f'
    
    integer                        :: n,m
    real(kind=dp), dimension(n,m) :: input
    complex(kind=dp), dimension(n/2+1,m) :: output
    
    real(kind=dp), dimension(n) :: in
    complex(kind=dp), dimension(n/2+1) :: out
    
    integer(kind=dp) :: plan
    
    integer :: i
    
    call dfftw_plan_dft_r2c_1d(plan,N,in,out,FFTW_ESTIMATE)
    
    do i=1,m
       in=input(:,i)
       call dfftw_execute(plan)
       output(:,i)=out
    end do
    call dfftw_destroy_plan(plan)
    
  end subroutine fft_fw_array


  subroutine fft_bw_array(input,output,n,m)
    
    include 'fftw3.f'
    
    integer                        :: n,m
    complex(kind=dp), dimension(n/2+1,m) :: input
    real(kind=dp), dimension(n,m) :: output
    
    complex(kind=dp), dimension(n/2+1) :: in
    real(kind=dp), dimension(n) :: out
    
    integer(kind=dp) :: plan
    
    integer :: i
    
    call dfftw_plan_dft_c2r_1d(plan,N,in,out,FFTW_ESTIMATE)
    
    do i=1,m
       in=input(:,i)
       call dfftw_execute(plan)
       output(:,i)=out
    end do
    call dfftw_destroy_plan(plan)
    
  end subroutine fft_bw_array

  subroutine fft_fw_vector(input,output,n)
    
    include 'fftw3.f'
    
    integer                        :: n
    real(kind=dp), dimension(n) :: input
    complex(kind=dp), dimension(n/2+1) :: output
    
    real(kind=dp), dimension(n) :: in
    complex(kind=dp), dimension(n/2+1) :: out
    
    integer(kind=dp) :: plan
    
    call dfftw_plan_dft_r2c_1d(plan,n,in,out,FFTW_ESTIMATE)
    
    in=input
    call dfftw_execute(plan)
    output=out

    call dfftw_destroy_plan(plan)
    
  end subroutine fft_fw_vector

  subroutine fft_bw_vector(input,output,n)
    
    include 'fftw3.f'
    
    integer                        :: n
    complex(kind=dp), dimension(n/2+1) :: input
    real(kind=dp), dimension(n) :: output
    
    complex(kind=dp), dimension(n/2+1) :: in
    real(kind=dp), dimension(n) :: out
    
    integer(kind=dp) :: plan
    
    call dfftw_plan_dft_c2r_1d(plan,N,in,out,FFTW_ESTIMATE)
    
    in=input
    call dfftw_execute(plan)
    output=out

    call dfftw_destroy_plan(plan)
    
  end subroutine fft_bw_vector

  subroutine fft_fw_3d(array,n,m,l)
    
    include 'fftw3.f'

    integer :: n,m,l
    complex(kind=dp), dimension(n,m,l), intent(inout) :: array

    integer(kind=8) :: plan

    call dfftw_plan_dft_3d(plan, n,m,l, array,array, &
         FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, array, array)
    call dfftw_destroy_plan(plan)

  end subroutine fft_fw_3d

  subroutine fft_bw_3d(array,n,m,l)
    
    include 'fftw3.f'

    integer :: n,m,l
    complex(kind=dp), dimension(n,m,l), intent(inout) :: array

    integer(kind=8) :: plan

    call dfftw_plan_dft_3d(plan, n,m,l, array,array,&
         FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, array, array)
    call dfftw_destroy_plan(plan)

  end subroutine fft_bw_3d

end module fft_methods
