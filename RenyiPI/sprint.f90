subroutine sprint(unit,jndex,dim,string)
      
  use md_variables

  integer :: jndex,dim,k, unit, l
  real(8) :: cnorm,a1
  character(*) :: string

  cnorm=max(anorm(jndex)-1.d0,1.d0)
  k=jndex-1
 !do l=1,dim
    k=k+2
    if (av(k) == 0.d0) then
      a1=0.d0
    else
      a1=sqrt(av(k)/cnorm)
    endif
    write(unit,'(3g20.10)') avp(k-1),av(k-1),a1
  !enddo
  call flush(unit)

  return
end

