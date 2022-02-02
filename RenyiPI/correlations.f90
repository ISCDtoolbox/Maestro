module correlations

  use kinds
  use fft_methods

  implicit none

  interface write_spectrum
     module procedure write_spectrum_array,write_spectrum_vector
  end interface

contains

  subroutine khintchine_array(array,spectrum,ntf,ndim)
    integer, intent(in) :: ntf,ndim
    real(kind=dp), dimension(ntf,ndim), intent(in) :: array
    complex(kind=dp), dimension(ntf/2+1,ndim,ndim), intent(out) :: spectrum

    complex(kind=dp), dimension(ntf/2+1,ndim) :: tf_array
    integer id1,id2

    call fft_fw_array(array,tf_array,ntf,ndim)

    do id1=1,ndim
       do id2=1,ndim
          spectrum(:,id1,id2)=&
               conjg(tf_array(:,id1))*tf_array(:,id2)
       end do
    end do

  end subroutine khintchine_array

  subroutine khintchine_array_multi(array,spectrum,ntf,ndim,naver)
    integer, intent(in) :: ntf,ndim,naver
    real(kind=dp), dimension(ntf,ndim,naver), intent(in) :: array
    complex(kind=dp), dimension(ntf/2+1,ndim,ndim), intent(out) :: spectrum

    complex(kind=dp), dimension(ntf/2+1,ndim) :: tf_array
    integer iw,id1,id2,iav

    spectrum=0.0_dp

    do iav=1,naver

       call fft_fw_array(array(:,:,iav),tf_array,ntf,ndim)

       do iw=1,ntf/2+1
          do id1=1,ndim
             do id2=1,ndim
                spectrum(iw,id1,id2)=spectrum(iw,id1,id2)+&
                     conjg(tf_array(iw,id1))*tf_array(iw,id2)
             end do
          end do
       end do

    end do

    spectrum=spectrum/dble(naver)

  end subroutine khintchine_array_multi

  subroutine khintchine_self(array,spectrum,ntf,ndim)
    integer, intent(in) :: ntf,ndim
    real(kind=dp), dimension(ntf,ndim), intent(in) :: array
    complex(kind=dp), dimension(ntf/2+1,ndim), intent(out) :: spectrum

    complex(kind=dp), dimension(:,:), allocatable :: tf_array
    integer id

    allocate(tf_array(ntf/2+1,ndim))

    call fft_fw_array(array,tf_array,ntf,ndim)

    do id=1,ndim
       spectrum(:,id)=&
            conjg(tf_array(:,id))*tf_array(:,id)
    end do

  end subroutine khintchine_self

  subroutine write_spectrum_array(unit,spectrum,ntf,ndim,domega)
    integer, intent(in) :: unit
    integer, intent(in) :: ntf,ndim
    complex(kind=dp), dimension(ntf/2+1,ndim,ndim), intent(in) :: spectrum
    real(kind=dp), intent(in),optional :: domega

    real(kind=dp) :: omega,dw
    integer       :: iw

    dw=1.0_dp
    if (present(domega)) dw=domega

    do iw=1,ntf/2+1
       omega=dw*dble(iw-1)
       write(unit,*) omega,real(spectrum(iw,:,:))
    end do

  end subroutine write_spectrum_array

  subroutine write_spectrum_diag(unit,spectrum,ntf,ndim,domega,intensity)
    integer, intent(in) :: unit
    integer, intent(in) :: ntf,ndim
    complex(kind=dp), dimension(ntf/2+1,ndim,ndim), intent(in) :: spectrum
    real(kind=dp), intent(in),optional :: domega
    real(kind=dp), dimension(ndim), intent(in), optional :: intensity
    real(kind=dp) :: omega,dw
    integer       :: iw,i

    real(kind=dp), dimension(ndim) :: intens
    character(len=20) :: myformat
    write(myformat,'(a,I4,a)')'(',(ndim+1),'F20.10)'

    dw=1.0_dp
    if (present(domega)) dw=domega

    intens=1.0_dp
    if (present(intensity)) intens=intensity

    do iw=1,ntf/2+1
       omega=dw*dble(iw-1)
       write(unit,myformat) omega,(real(spectrum(iw,i,i))*intens(i),i=1,ndim)
    end do

  end subroutine write_spectrum_diag

  subroutine write_spectrum_vector(unit,spectrum,ntf,ndim,domega)
    integer, intent(in) :: unit
    integer, intent(in) :: ntf,ndim
    complex(kind=dp), dimension(ntf/2+1,ndim), intent(in) :: spectrum
    real(kind=dp), intent(in),optional :: domega

    real(kind=dp) :: omega,dw
    integer       :: iw,i
    character(len=20) :: myformat
    write(myformat,'(a,I4,a)')'(',(ndim+1),'F20.10)'

    dw=1.0_dp
    if (present(domega)) dw=domega

    do iw=1,ntf/2+1
       omega=dw*dble(iw-1)
       write(unit,myformat) omega,(real(spectrum(iw,i)),i=1,ndim)
    end do

  end subroutine write_spectrum_vector

end module correlations
