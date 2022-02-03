module md_variables

      implicit none

       INTEGER  :: &
       unit_dot_ep,  &
       unit_dot_ep2, &
       unit_dot_epsr, &
       unit_dot_eplr, &
       unit_dot_ek,   &
       unit_dot_positions, &
       unit_dot_positions_cen, &
       unit_dot_forces,    &
       unit_dot_forces_cen,    &
       unit_dot_velocities, &
       unit_dot_velocities_cen, &
       unit_dot_xyz,        &
       unit_dot_localtemp,  &
       unit_dot_sigma,      &
       unit_dot_out

       LOGICAL :: restart_pimd
! Parameters
       integer,parameter :: nmax=7,mtmax=3001,ngofrtb=200
       integer,parameter :: maxnab=500,maxm=14,mav=2000+nmax,mtypes=2

       real(8), parameter :: pi=3.14159265358979d0, kbm1=315775.024d0  ! in K/Ha units
       real(8), parameter :: eps_cov = 1.d-10   ! cut off of lower eigenvalues for the inversion of the covariance matrix
       double precision, parameter :: hbar2pi=4.135667516d0,a0=0.529177d0
       double precision, parameter :: ev2cm1=8065.543937d0
       double precision, parameter :: ha2cm1 = 219474.63137d0
       double precision, parameter :: timeau2fs = 0.0241888432658d0

! Variables

       integer :: natoms,ndimMD,nblocks,nstep,iprint,iratio,n,nspecies,gMD,nbeadMD
       integer :: nstep_block, numqsymmpoints
       integer :: itarget,irun,lfnl,lonl,lonll,ntmax,ifl,nav,icl,ikin,nunitcells,ibrav
       integer, dimension (:), allocatable :: igr,ipot,partner,indx
       integer :: ncellsx,ncellsy,ncellsz

       real(8) :: avec(3,3),avecsp(3,3)
       real(8) :: dens,rcut,rg,rlist,rmin,tempMD,drmax,delta_harm,delta_force
       real(8) :: vol,S,PS,Q,delt,sint,mtot,gammaMD,sig1,sigmacov,sigmavar
       real(8) :: ekinq, ekinqp, epot_centroid,kspring,c0_k,c1_k,c2_k,tfakeMD,deltahtilde,epot_old
       real(8), dimension (:), allocatable :: amas,av,avp,anorm,anormp,avpre,anormpre
       real(8), dimension (:), allocatable :: sigma,sig2,sig3,sig4,dynmat_eig,tmes_bead,dynmatforce_eig
       real(8), dimension (:), allocatable :: friction_mode,omega_mode,cost1,cost2,gamma_eigen
       real(8), dimension (:), allocatable :: alpha_qmc,alphaqmc_eig
       real(8), dimension (:,:), allocatable :: rcm,vcm,dynmat,covMD,cmatrix,dynmat_force0
       real(8), dimension (:,:), allocatable :: el,rtilde,fnoiseMD,rcentroid,rcentroidtilde,forcedyn
       real(8), dimension (:,:,:), allocatable :: rpos,rpos_init,forceMD,vel,pimp,etaMD,nuMD,velocity
       real(8), dimension (:,:,:), allocatable :: rtilde_mode, ptilde,forceMD_old,rpos_old

       logical :: nh,fixcm,verbose,yesglobal
       character(len=50) :: filen,potential  !!filen is the name of the system
       character(len=150) :: output_dir, output_directory,pot_type
       character(len=2), dimension(:), allocatable :: ion_name

! Variables added by Miha Srdinsek
       integer :: renyi_order
       real(8) :: renyi_lam_s, renyi_lam_j

! Myreweight and former TurboRVB variables

       integer :: ieskin,iflagerr,info,nbead,nion
       integer :: iscramax,rank
       real(8) :: dt,cost,temp,delta0,delta0k,delta0q,ener_true,friction,epstion,dynmat_cutoff
       real(8) :: normcorr,scalecov,sigma_true,tmes,friction_qmc
       real(8), dimension(:,:), allocatable :: fk,mass_ion,fbead,kdyn,cov_pimd
       real(8), dimension(:), allocatable :: psip,scalpar,kdyn_eig,cov,cov_old
       logical :: yesquantum,yesturboq,yessecond,velocity_langevin,averaged_cov

end module md_variables
