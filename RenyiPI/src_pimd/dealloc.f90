subroutine dealloc

  use md_variables
  
  if (irun .eq. 3) then
     deallocate(cost1,cost2)
     deallocate(tmes_bead)
     deallocate(mass_ion)
     deallocate(rpos_old,forceMD_old)
     if (sigmacov .ne. 0) then
       deallocate(fk)
       deallocate(cov)
       deallocate(alpha_qmc)
       deallocate(alphaqmc_eig)
       deallocate(gamma_eigen)
     endif     
     if (yesquantum) then  
       deallocate(omega_mode,friction_mode) 
       deallocate(cmatrix)  
       deallocate(ptilde)
       deallocate(rtilde_mode) 
       deallocate(fbead,rcentroid,rcentroidtilde)
     endif
  endif
  
  
  if (irun .eq. 4) then
     deallocate(psip)
     deallocate(scalpar)
     deallocate(cov_pimd)
     deallocate(fk)
     deallocate(tmes_bead)
     deallocate(rcentroid,rcentroidtilde)
     deallocate(mass_ion)
     if (yesquantum) then 
       deallocate(fbead)
       deallocate(kdyn,kdyn_eig)
     endif
     if (averaged_cov.and.sigmacov.ne.0.d0) deallocate(cov)
  end if
  
  return
  
end subroutine
