module qc_system
  implicit none

  !--- (JOB INFO) ---
  integer, parameter :: job_id_ener = 0
  integer, parameter :: job_id_engrad = 1
  integer, parameter :: job_id_opt = 2
  integer, parameter :: job_id_md   = 3
  integer, parameter :: job_id_freq = 4
  integer, parameter :: job_id_rescale = 5

  integer, parameter :: theory_id_scf = 0
  integer, parameter :: theory_id_mp2 = 1
  integer, parameter :: theory_id_ccsd = 2
  integer, parameter :: theory_id_ccsd_pt = 3
  integer, parameter :: embed_id_esp = 0
  integer, parameter :: embed_id_dipole = 1
  integer, parameter :: embed_id_atomdipole = 2
  character(len=200) :: code, embedding_field, scratch_top, my_scratch
  
  ! optimizer options
  integer, parameter :: opt_id_grad = 0
  integer, parameter :: opt_id_lbfgs = 1
  logical :: l_freezecell
  integer :: optimizer_id
  integer :: opt_maxiter
  real*8 :: atom_gmax, lat_gmax
  
  real*8, parameter :: delta_bohr = 1.0d-3
  real*8, parameter :: delta_ang  = delta_bohr*0.529177249d0 
  integer :: n_hess(3)
  logical :: l_bq, l_grd, l_frq, l_opt, l_md, l_bsse
  character(len=10) :: theory
  integer :: theory_id, embed_id

  ! rescale info
  real*8 :: rescale_param(6)
  
  !--- (FRAGMENT) ---
  integer          :: nfrg
  !--- (GEOM)---
  integer          :: natom, nmol
  character(len=2), dimension(:), allocatable :: at_atnm
  integer, dimension(:),  allocatable :: mol_iatom, mol_nsite, mol_netcharge, mol_ichg, mol_nchg
  real*8, dimension(:),   allocatable :: at_mass
  real*8, dimension(:,:), allocatable :: chg_pos, chg_pos_new, chg_coef
  real*8, dimension(:,:), allocatable :: gm_pos, gm_pos_old  ! 
  real*8, dimension(:,:), allocatable :: d_rr, corr_d_rr !, d_rq  ! gradient; AK : corr_d_rr -  MP2 correction to gradient for MTS-MD
  !real*8, dimension(:,:,:), allocatable :: d_mbq            ! monomer d_bq
  !real*8, dimension(:,:,:), allocatable :: d_mrr            ! monomer d_rr
  real*8, dimension(:,:),   allocatable :: d_mr, corr_d_mr   ! monomer d_rr
  real*8, dimension(:,:),   allocatable :: d_bq0
  real*8, dimension(:,:),   allocatable :: d_mbqi
  real*8, dimension(:,:),   allocatable :: d_mbqj
  real*8 :: m_virt(3,3)
  
  !--- (MD)----
  integer :: nequil, nprod, nsave, nstep
  real*8, dimension(:,:), allocatable :: vel
  integer :: pstart, pend, pstep
  real*8 :: temp0, pext0
  real*8 :: dt, virt(3,3)
  real*8, dimension(3,3) :: Pr, pmom2
  
  integer :: maxnum_pair = 500
  integer :: maxnum_bq = 8000
  
  real*8 :: lat_a, lat_b, lat_c !M!
  real*8 :: lat_alpha, lat_beta, lat_gamma !M!
  integer :: lat_axis !M!
  real*8, dimension(3,3) :: lat, lati, lat_old !M!
  logical :: a_dim, b_dim, c_dim
  
  real*8  :: rcut_qq2 ! rcut for QM-QM
  real*8  :: rcut_bq2 ! rcut for QM-BQ
  real*8  :: rcut_lr2 ! rcut for QM-LR !M!

  !--- (BIM )--
  !real*8, dimension(:),   allocatable :: en_list
  real*8 :: en_kin

  ! logical 
  logical :: ensemble_nvt
  logical :: is_temp_tol
  logical :: use_respa, xirespa, xorespa
  character(len=120) :: fconfig
  character(len=30)  :: basis_type
  character(len=30)  :: jname
  character(len=240) :: jcmd

  ! --- (Nose Hoover Method)---
  ! real*8, dimension(2) :: gnh, vnh, xnh, qmass
  ! real*8 :: ekin0, gfree0, glogv, vlogv, xlogv, pmass
  ! 
  real*8 :: time_system
  real*8 :: ekin0, gfree0, glogv, vlogv, xlogv, pmass
  integer :: nnos, nresn, nresp                             ! NRESP - number of multiple timesteps for Verlet
  real*8, dimension(:), allocatable :: vnh, xnh, gnh, qmass ! thermostat coordinates, velocities and masses
  integer, parameter :: nyosh = 3                           ! order of Yoshida integrator
  ! w_i for Yoshida scheme
  real*8, dimension(nyosh) :: wyosh = (/1.0/(2.0 - 2.0**(1.0/3.0)), 1.0 - 2.0 /(2.0 - 2.0**(1.0/3.0)), 1.0/(2.0 - 2.0**(1.0/3.0))/) 
  real*8, dimension(nyosh) :: wdti2, wdti4, wdti8
  
  

end module qc_system
