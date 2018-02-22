subroutine qc_md_run (irst)
  use qc_system
  use qc_mpi
  use qc_bim
  implicit none
  integer, intent(in) :: irst
  real*8 :: temp, umon, udim, ulr
  integer :: istep, istart0


  if (sys_master) then
     open (45, file='trajectory.xyz', access='append', status='unknown')
  end if

  call md_nose_init

  istart0 = 0
  istep   = 0
  if (irst .eq. 0) then
     call read_system 
     call md_vel_init
  else
     call md_restart (istart0)
  end if

  !-- Parrallel JOB
  call qc_mpi_barrier 
  call qc_mpi_real_bcast (gm_pos, 3*natom)

  if (irst .eq. 0) then
     ! Obtain Gradient
     call qc_bim_pot (umon, udim, ulr)
     istep = 0
     if (sys_master) then
        call md_kinen
        call md_report (istep+istart0, umon, udim, ulr, temp)
     end if

  end if
  
  !--- begin MD procedure

  do istep = 1, nprod
     
     call md_integrate (umon, udim, ulr)
     
     if (sys_master) then
        call md_report (istep+istart0, umon, udim, ulr, temp)
        if (mod(istep+istart0, nsave) .eq. 0) &
             call md_save (istep+istart0)
        if (mod(istep+istart0,5) .eq. 0) then
           call md_trajectory (istep+istart0, umon+udim+ulr, temp)
        end if
     end if
     
  end do

  if (sys_master) close (45)

end subroutine qc_md_run


subroutine md_integrate (Umon, Udim, Ulr)
  use qc_system
  use consts
  use qc_mpi
  use qc_bim
  implicit none
  real*8, intent(out)  :: Umon, Udim, Ulr
  integer :: iat, ii, im, is, ia
  real*8  :: dt2
  real*8, dimension(3) :: arr, vrr
  real*8 :: p(6)
  real*8, dimension(:,:), allocatable :: gm_pos0

  allocate(gm_pos0(3,natom))
  dt2 = 0.5d0*dt

  if (sys_master) then

     call md_nose_chain 
     
     do iat = 1, natom
        arr(1:3)    = -d_rr(1:3,iat) / at_mass(iat)
        vrr(1:3)    = vel(1:3,iat) + arr(1:3)*dt2
        gm_pos(1:3,iat) = gm_pos(1:3,iat)  + dt*vrr(1:3)*bohr2ang
        vel(1:3,iat)= vrr(1:3)
     end do

     
     call md_pbc
  end if
  call qc_mpi_real_bcast (gm_pos, 3*natom)
    
  ! Cell shape/size update: use lattice gradient
  call MDlat_grad(p)
  p(1:3) = p(1:3) / bohr2ang

  ! store scaled/internal coordinates
  ii = 1
  do im = 1,nmol
      is = mol_iatom(im)
      gm_pos0(1:3, ii) = matmul(lati,gm_pos(1:3, is+1))
      ii = ii + 1
      do ia = is+2,is+mol_nsite(im)
          gm_pos0(1:3,ii) = (gm_pos(1:3, ia) - gm_pos(1:3, is+1))*ang2bohr
          ii = ii + 1
      end do
  end do
    
  ! Update cell
  if (a_dim) lat_a      = lat_a*(1.0d0-0.01d0*p(1))
  if (b_dim) lat_b      = lat_b*(1.0d0-0.01d0*p(2))
  if (c_dim) lat_c      = lat_c*(1.0d0-0.01d0*p(3))
  
  call MDupdate_geometry(gm_pos0)

  call qc_bim_pot  (Umon, Udim, Ulr)
  
  if (sys_master) then

     do iat = 1, natom
        arr(1:3) = -d_rr(1:3,iat)/at_mass(iat)
        vel(1:3,iat) = vel(1:3,iat) + dt2*arr(1:3)
     end do
     
     call md_nose_chain

  end if
  deallocate(gm_pos0)

end subroutine md_integrate


subroutine md_nose_init
  use qc_system
  use consts  ! 
  implicit none
  real(8) :: time_system, omega_system, omega2

  Ekin0 = temp0*boltz   ! K-> A.U.

  gfree0 = 3.0d0 * natom
  
  xnh(1) = 0.0d0
  gnh(1) = 0.0d0
  vnh(1) = 0.0d0
  !time_system  = 8.881d0/au_time ! fs --> A.U.
  time_system  = 20.000d0/au_time ! fs --> A.U.
  omega_system = 2.0d0*pi/time_system
  omega2       = omega_system**2

  qmass(1) = gfree0 * Ekin0 / omega2
  xnh(2) = 0.0d0
  gnh(2) = 0.0d0
  qmass(2) = Ekin0 / omega2
  
  glogv = 0.0d0
  vlogv = 0.0d0
  xlogv = 0.0d0

  pmass = (gfree0 + 3.0d0) * Ekin0 * 0.1d0 !* 100.0d0

  return

end subroutine md_nose_init



subroutine md_nose_chain 
  use qc_system
  use consts
  implicit none
  integer :: iat
  real*8 :: dt2, dt4, dt8
  real*8 :: scale

  call md_kinen

  dt2 = dt / 2.0d0
  dt4 = dt / 4.0d0
  dt8 = dt / 8.0d0
  scale = 1.0D0

  if (ensemble_nvt) then ! Nose-Hoover Method for constant temperature
     gnh(2) = (qmass(1)*vnh(1)*vnh(1) - Ekin0) / qmass(2)
     vnh(2) = vnh(2) + gnh(2)*dt4
     vnh(1) = vnh(1) * exp(-vnh(2)*dt8)

     gnh(1) = (2.0d0*en_kin - gfree0*Ekin0)/qmass(1)

     vnh(1) = vnh(1) + gnh(1)*dt4
     vnh(1) = vnh(1) * exp(-vnh(2)*dt8)

     xnh(1) = xnh(1) + vnh(1)*dt2
     xnh(2) = xnh(2) + vnh(2)*dt2

     scale = exp(-vnh(1)*dt2)

     en_kin = en_kin * scale*scale

     vnh(1) = vnh(1)*exp(-vnh(2)*dt8)
     gnh(1) = (2.0d0*en_kin - gfree0*Ekin0)/qmass(1)
     vnh(1) = vnh(1) + gnh(1)*dt4
     vnh(1) = vnh(1)*exp(-vnh(2)*dt8)

     gnh(2) = (qmass(1)*vnh(1)*vnh(1) - Ekin0)/qmass(2)
     vnh(2) = vnh(2) + gnh(2)*dt4
  end if

  ! Andersen-Hoover method combined with Nose-Hoover chain

  do iat = 1, natom
     vel(1:3,iat) = vel(1:3, iat)*scale
  end do
  pmom2 = pmom2 * scale*scale

  return

end subroutine md_nose_chain



subroutine md_kinen
    use qc_system
    use consts
    implicit none
    integer :: iat
    integer :: i, j

    pmom2 = 0.0d0

    do iat = 1, natom
        do j = 1,3
        do i = 1,3
            pmom2(i,j) = pmom2(i,j) + at_mass(iat)*vel(i,iat)*vel(j,iat)
        end do
        end do
    end do

    en_kin = 0.5d0 * (pmom2(1,1) + pmom2(2,2) + pmom2(3,3))

    return

end subroutine md_kinen


subroutine md_save (istep)
  use qc_system
  use qc_mpi
  implicit none
  integer, intent(in) :: istep
  integer :: iat
  
  if (sys_master) then
     open (10, file='mdrr.sav', status='unknown')
     open (11, file='mdrr.xyz', status='unknown')
     write (10, *) natom
     write (10, '(i9)') istep 
     write (11, *) natom
     write (11, '(i9)') istep 
     
     ! Coordinate
     do iat = 1, natom
        write (10, *) at_atnm(iat), gm_pos(1:3,iat)
        write (11, *) at_atnm(iat), gm_pos(1:3,iat)
     end do
     ! Velocity
     do iat = 1, natom
        write (10, *) vel(1:3,iat)
     end do
     ! Gradient
     do iat = 1, natom
        write (10, *) d_rr(1:3,iat)
     end do

     !do iat = 1, natom
     !   write (10, *) at_chg(iat)
     !end do

     ! Nose-Hoover Temperature
     
     write (10, *) xnh(1), gnh(1), vnh(1), qmass(1)
     write (10, *) xnh(2), gnh(2), vnh(2), qmass(2)
     ! Pressure
     write (10, *) xlogv, glogv, vlogv, pmass

     close(10)
     close(11)
  end if

  
end subroutine md_save


subroutine md_trajectory (istep,Upot,temp)
  use qc_system
  use consts
  use qc_mpi
  implicit none
  integer, intent(in) :: istep
  real*8,  intent(in) :: Upot, temp
  real*8  :: time_ps
  integer :: iat
  
  if (sys_master) then
     time_ps = istep*dt*au_time/1.0e3 ! ps
    write (45, '(i0)') natom
    write (45, '(f0.8, 5f13.8, i3, A6, f16.8, A6, f10.2, A9, f9.3)') lat_a, lat_b, lat_c, lat_alpha, &
    lat_beta, lat_gamma, lat_axis, 'Upot=', Upot, 'Temp=', temp, 'Time(ps) ', time_ps
    
     ! Coordinate
     do iat = 1, natom
        write (45, '(A,3f16.8)') at_atnm(iat), gm_pos(1:3,iat)
     end do
     
     call flush(45)
  end if

end subroutine md_trajectory



subroutine md_restart (istart)
  use qc_system
  use qc_mpi
  implicit none
  integer, intent(out) :: istart
  integer :: iat

  if (sys_master) then
     open (10, file='./mdrr.sav', status='old')
     read (10, *) !natom
     read (10, *) istart
     
     ! Coordinate
     do iat = 1, natom
        read (10, *) at_atnm(iat), gm_pos(1:3,iat)
     end do
     
     do iat = 1, natom
        read (10, *) vel(1:3, iat) 
     end do
     
     do iat = 1, natom
        read (10, *) d_rr(1:3, iat) 
     end do
     
     !do iat = 1, natom
     !   read (10, *) at_chg(iat)
     !end do

     read (10, *) xnh(1), gnh(1), vnh(1), qmass(1)
     read (10, *) xnh(2), gnh(2), vnh(2), qmass(2)
     ! Pressure
     read (10, *) xlogv, glogv, vlogv, pmass
     
     close(10)

     call md_kinen

  end if

  call qc_mpi_barrier

  call qc_mpi_int_bcast1 (istart)
  !call qc_mpi_real_bcast (at_chg, natom)
  call qc_mpi_real_bcast (d_rr, 3*natom)

end subroutine md_restart


subroutine md_vel_init
  use qc_system
  use consts
  use qc_mpi
  implicit none
  integer :: iat
  real*8 :: vxx, vyy, vzz
  real*8 :: rand_f
  real*8 :: tmp, scale
  real*8, dimension(3) :: pmom

  if (sys_master) then
     rand_f = 3.0d0**15
     
     scale = sqrt(Ekin0)
     pmom = 0.0d0
     do iat = 1, natom
        call random (vxx, vyy, vzz,rand_f)
        vel(1,iat) = vxx*scale
        vel(2,iat) = vyy*scale
        vel(3,iat) = vzz*scale
        pmom = pmom + at_mass(iat)*vel(1:3,iat)
     end do
     
     pmom = pmom / natom

     ! correct velocities
     do iat = 1, natom
        vel(1:3,iat) = vel(1:3,iat) - pmom(1:3)/at_mass(iat)
     end do
     
     call md_kinen
     
     tmp  = 2.0D0*en_kin/(boltz*gfree0)
     scale = sqrt (temp0/tmp)
     
     do iat = 1, natom
        vel(1:3,iat) = vel(1:3,iat)*scale
     end do
  end if

end subroutine md_vel_init


subroutine random (vxx, vyy, vzz, rand_f)
  implicit none
  real*8 :: vxx, vyy, vzz
  real*8 :: rand_b, rand_f, rand_r, rand_ri
  real*8 :: ra, re, rg, pid

  rand_b  = 5.0d0**3
  rand_r  = 2.0d0**24
  rand_ri = 1.0d0/rand_r
  pid = 6.28318530717958647692528677d0

  ra = rand_b * rand_f
  rand_f=ra-dint(ra*rand_ri)*rand_r
  re=rand_f*rand_ri
  ra=rand_b*rand_f
  rand_f=ra-dint(ra*rand_ri)*rand_r
  rg=rand_f*rand_ri
  vxx=sqrt(-2.0d0*log(re))*sin(pid*rg)
  ra=rand_b*rand_f
  rand_f=ra-dint(ra*rand_ri)*rand_r
  re=rand_f*rand_ri
  ra=rand_b*rand_f
  rand_f=ra-dint(ra*rand_ri)*rand_r
  rg=rand_f*rand_ri
  vyy=sqrt(-2.0d0*log(re))*sin(pid*rg)
  ra=rand_b*rand_f
  rand_f=ra-dint(ra*rand_ri)*rand_r
  re=rand_f*rand_ri
  ra=rand_b*rand_f
  rand_f=ra-dint(ra*rand_ri)*rand_r
  rg=rand_f*rand_ri
  vzz=sqrt(-2.0d0*log(re))*sin(pid*rg)
  
end subroutine random



subroutine md_report (istep, Umon, Udim, Ulr, temp)
  use qc_system
  use qc_lattice
  use consts
  use qc_mpi
  implicit none
  integer, intent(in) :: istep
  real*8,  intent(in) :: Umon, Udim, Ulr
  real*8,  intent(out):: temp
  real*8 :: Upot, Ten, scale
  real*8 :: time_ps, press
  real*8 :: qkin1, qkin2, qpot1, qpot2, vol
  integer :: i, j

  if (sys_master) then
     Upot = Umon + Udim + Ulr
     Ten = Upot + en_kin
     temp = 2.0d0*en_kin/(boltz*gfree0)

     ! Ang --> Bohr
     vol  = volume(lat)*ang2bohr**3 !M!
     
     do i = 1,3
     do j = 1,3
        Pr(i,j) = (pmom2(i,j) + virt(i,j))/vol
    end do
    end do
     
     press = (Pr(1,1)+Pr(2,2)+Pr(3,3))/3.0d0 * au_bar ! E_h/bohr**3 --> Bar

     time_ps = istep*dt*au_time/1.0e3 ! ps
     
     if (ensemble_nvt) then
        qkin1 = 0.5d0 * qmass(1)*vnh(1)*vnh(1)
        qkin2 = 0.5d0 * qmass(2)*vnh(2)*vnh(2)
        qpot1 = gfree0* Ekin0 * xnh(1)
        qpot2 = Ekin0 * xnh(2)
        Ten = Ten + qkin1 + qkin2 + qpot1 + qpot2
     end if
     
     !write (6,'(I8,F10.3,5F14.6,F10.2)') &
     write (6,'(A, I8,F10.3,5F14.6,F10.2,F12.2)') &
          '@ ', istep, time_ps,Upot/Nmol, Umon/Nmol, Udim/Nmol, Ulr/Nmol, &
          Ten/Nmol, temp, press
     !  end if
     
     if (is_temp_tol) then
     !if (.false.) then
        if (temp .lt. temp0 - 50.0d0 .or. temp .gt. temp0 + 50.0d0) then
           scale = sqrt(temp0/temp)
           vel = vel * scale
        end if
     end if
     
     call flush(6)
  end if

end subroutine md_report



subroutine md_pbc
    use qc_system
    implicit none
    integer :: iat, im, nsite, is

    gm_pos = matmul(lati, gm_pos) ! scale coordinates
    
    do im = 1, nmol ! fold molecules back into central cell
        iat   = mol_iatom(im)
        nsite = mol_nsite(im)
        ! X
        if (gm_pos(1,iat+1) .lt. 0.0d0) then
            do is = 1, nsite
                gm_pos(1,iat+is) = gm_pos(1,iat+is) + 1.0d0
            end do
        else if (gm_pos(1,iat+1) .ge. 1.0d0) then
            do is = 1, nsite
                gm_pos(1,iat+is) = gm_pos(1,iat+is) - 1.0d0
            end do
        end if
        ! Y
        if (gm_pos(2,iat+1) .lt. 0.0d0) then
            do is = 1, nsite
                gm_pos(2,iat+is) = gm_pos(2,iat+is) + 1.0d0
            end do
        else if (gm_pos(2,iat+1) .ge. 1.0d0) then
            do is = 1, nsite
                gm_pos(2,iat+is) = gm_pos(2,iat+is) - 1.0d0
            end do
        end if
        ! Z
        if (gm_pos(3,iat+1) .lt. 0.0d0) then
            do is = 1, nsite
                gm_pos(3,iat+is) = gm_pos(3,iat+is) + 1.0d0
            end do
        else if (gm_pos(3,iat+1) .ge. 1.0d0) then
            do is = 1, nsite
                gm_pos(3,iat+is) = gm_pos(3,iat+is) - 1.0d0
            end do
        end if
    end do

    gm_pos = matmul(lat, gm_pos) ! unscale coordinates
    

end subroutine md_pbc

subroutine MDlat_grad (lat_constants_grad)
    use qc_system
    use qc_lattice
    use qc_mpi
    use consts
    implicit none
    real*8 :: stress_bar(3,3), vol
    real*8 :: g_lat(3,3), d_al(3,3), d_be(3,3), d_ga(3,3)
    real*8, intent(out) :: lat_constants_grad(6)

    ! g_lat: enthalpy derivatives wrt all 9 lattice vector components
    vol = volume(lat) * ang2bohr**3
    stress_bar = virt/vol*au_bar ! bar
    stress_bar = 0.5d0*( stress_bar + transpose(stress_bar) ) ! symmetrized (no rotation)
    stress_bar(1,1) = stress_bar(1,1) - pext0
    stress_bar(2,2) = stress_bar(2,2) - pext0
    stress_bar(3,3) = stress_bar(3,3) - pext0
    g_lat = (-1.0d0*vol/au_bar)*matmul(stress_bar, transpose(lati)/ang2bohr) ! au/bohr

    ! get derivatives wrt a,b,c,alpha,beta,gamma
    call lat_angle_grad(d_al, d_be, d_ga)
    lat_constants_grad(1) = dot_product(g_lat(:,1), lat(:,1)/lat_a) ! dH/da (au/bohr)
    lat_constants_grad(2) = dot_product(g_lat(:,2), lat(:,2)/lat_b) ! dH/db
    lat_constants_grad(3) = dot_product(g_lat(:,3), lat(:,3)/lat_c) ! dH/dc
    lat_constants_grad(4) = sum( (g_lat/bohr2ang)*d_al)
    lat_constants_grad(5) = sum( (g_lat/bohr2ang)*d_be)
    lat_constants_grad(6) = sum( (g_lat/bohr2ang)*d_ga)

    ! avoid need NaN if in less than 3 dimensions
    if (.not. a_dim) lat_constants_grad(1) = 0.0d0
    if (.not. b_dim) lat_constants_grad(2) = 0.0d0
    if (.not. c_dim) lat_constants_grad(3) = 0.0d0
    if (.not. (a_dim.and.b_dim)) lat_constants_grad(4) = 0.0d0
    if (.not. (a_dim.and.b_dim.and.c_dim)) lat_constants_grad(5) = 0.0d0
    if (.not. (a_dim.and.b_dim.and.c_dim)) lat_constants_grad(6) = 0.0d0
end subroutine MDlat_grad

subroutine MDupdate_geometry (gm_pos0)
    use qc_system
    use qc_mpi
    use qc_lattice
    use consts
    implicit none
    integer :: ii, im, is, ia
    real*8, intent(in) :: gm_pos0(3,natom)

    call create_lattice_vectors
    call update_scaling_matrix
    ii = 1
    do im = 1,nmol
        is = mol_iatom(im)
        gm_pos(1:3,is+1) = matmul(lat,gm_pos0(1:3,ii))
        ii = ii + 1
        do ia = is+2, is+mol_nsite(im)
            gm_pos(1:3,ia) = gm_pos(1:3,is+1) + gm_pos0(1:3,ii)*bohr2ang
            ii = ii + 1
        end do
    end do
    call qc_mpi_barrier
    call qc_mpi_real_bcast1(lat_a)
    call qc_mpi_real_bcast1(lat_b)
    call qc_mpi_real_bcast1(lat_c)
    call qc_mpi_real_bcast1(lat_alpha)
    call qc_mpi_real_bcast1(lat_beta)
    call qc_mpi_real_bcast1(lat_gamma)
    call qc_mpi_real_bcast (lat, 9)
    call qc_mpi_real_bcast (lati, 9)
    call qc_mpi_real_bcast (gm_pos, 3*natom)
end subroutine MDupdate_geometry
