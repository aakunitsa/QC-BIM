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
  call qc_mpi_real_bcast (gm_pos, 3*natom) ! makes sense if we're doing restart

  if (irst .eq. 0) then
     ! Obtain Gradient
     call qc_bim_pot (umon, udim, ulr) ! If RESPA is used ==> umon and udim correspond to RHF 
     istep = 0
     !if (sys_master) then
     !   call md_kinen
     !   call md_report (istep+istart0, umon, udim, ulr, temp)
     !end if

  end if
  
  !--- begin MD procedure

  do istep = 1, nprod
     
     if ( .not. use_respa ) then 
         call md_integrate (umon, udim, ulr)
     else
         call md_integrate_respa (umon, udim, ulr)
     end if
     
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
  integer :: iat
  real*8  :: dt2
  real*8, dimension(3) :: arr, vrr

  dt2 = 0.5d0*dt

  if (sys_master) then

     call md_nose_chain (dt) ! integrates from t = 0 to t = dt / 2. !
     
     do iat = 1, natom
        arr(1:3)    = -d_rr(1:3,iat) / at_mass(iat)
        vrr(1:3)    = vel(1:3,iat) + arr(1:3)*dt2
        gm_pos(1:3,iat) = gm_pos(1:3,iat)  + dt*vrr(1:3)*bohr2ang
        vel(1:3,iat)= vrr(1:3)
     end do
     
     call md_pbc
  end if

  !--- Velocity Verlet algorithm

  call qc_mpi_barrier ()
  call qc_mpi_real_bcast (gm_pos, 3*natom)

  call qc_bim_pot  (Umon, Udim, Ulr)
  
  if (sys_master) then

     do iat = 1, natom
        arr(1:3) = -d_rr(1:3,iat)/at_mass(iat)
        vel(1:3,iat) = vel(1:3,iat) + dt2*arr(1:3)
     end do
     
     call md_nose_chain (dt)

  end if

end subroutine md_integrate

subroutine md_integrate_respa (Umon, Udim, Ulr)
  use qc_system
  use consts
  use qc_mpi
  use qc_bim
  implicit none
  real*8, intent(out)  :: Umon, Udim, Ulr
  integer :: iat, istep
  real*8  :: dt2, dti, dti2
  real*8, dimension(3) :: arr, vrr

  dt2 = 0.5d0*dt
  dti = dt / nresp
  dti2 = 0.5d0 * dti 

  if (sys_master .and. xorespa) call md_nose_chain (dt)
  
  do istep = 1, nresp
    
    if (sys_master) then

        if (xirespa) call md_nose_chain (dti)
        if (istep .eq. 1) then
            use_respa = .false.
            theory='scf' ! should save the old theory string somewhere
            do iat = 1, natom
                arr(1:3)    = -corr_d_rr(1:3,iat) / at_mass(iat)
                vrr(1:3)    = vel(1:3,iat) + arr(1:3)*dt2
                vel(1:3,iat)= vrr(1:3)
            end do
        end if
        do iat = 1, natom
            arr(1:3)    = -d_rr(1:3,iat) / at_mass(iat)
            vrr(1:3)    = vel(1:3,iat) + arr(1:3)*dti2
            gm_pos(1:3,iat) = gm_pos(1:3,iat)  + dti*vrr(1:3)*bohr2ang
            vel(1:3,iat)= vrr(1:3)
        end do

        call md_pbc

        if (istep .eq. nresp) then
            use_respa = .true.
            theory='mp2'
        end if

    end if 

    call qc_mpi_barrier ()
    call qc_mpi_real_bcast (gm_pos, 3*natom)
    call qc_mpi_logical_bcast1 (use_respa)
    call qc_mpi_char_bcast (theory)

    call qc_bim_pot (Umon, Udim, Ulr)

    if (sys_master) then
        do iat = 1, natom
            arr(1:3)    = -d_rr(1:3,iat) / at_mass(iat)
            vrr(1:3)    = vel(1:3,iat) + arr(1:3)*dti2
            vel(1:3,iat)= vrr(1:3)
        end do
        if (istep .eq. nresp) then
            do iat = 1, natom
                arr(1:3)    = -corr_d_rr(1:3,iat) / at_mass(iat)
                vrr(1:3)    = vel(1:3,iat) + arr(1:3)*dt2
                vel(1:3,iat)= vrr(1:3)
            end do
        end if 
        if (xirespa) call md_nose_chain(dti)
    end if
    
  end do
  
  if (sys_master .and. xorespa) call md_nose_chain (dt)

end subroutine md_integrate_respa

subroutine md_nose_init
  use qc_system
  use consts  ! 
  implicit none
  real(8) :: omega_system, omega2
  integer :: i

  allocate(qmass(nnos))
  allocate(xnh(nnos))
  allocate(vnh(nnos))
  allocate(gnh(nnos))

  Ekin0 = temp0*boltz   ! K-> A.U.

  gfree0 = 3.0d0 * natom
  
  xnh  = 0.0d0
  vnh = 0.0d0
  !time_system  = 8.881d0/au_time ! fs --> A.U.
  time_system  = 20.000d0/au_time ! fs --> A.U. -- why? AK
  omega_system = 2.0d0*pi/time_system
  omega2       = omega_system**2

  qmass(1) = gfree0 * Ekin0 / omega2
  qmass(2:nnos) = Ekin0 / omega2

  gnh   = -1.0d0 * Ekin0 / qmass 

  
  glogv = 0.0d0
  vlogv = 0.0d0
  xlogv = 0.0d0

  pmass = (gfree0 + 3.0d0) * Ekin0 * 0.1d0 !* 100.0d0

  xirespa = .false.
  xorespa = .true.

  return

end subroutine md_nose_init

subroutine md_nose_chain (tau)
  ! Propagates Nose - Hoover chain from t = 0 to t = ___tau/2___
  use qc_system
  use consts
  implicit none
  real*8, intent(in) :: tau
  integer :: iat, iresn, iyosh, inos
  real*8 :: scale, aa

  wdti2 = tau / 2.0d0 / nresn * wyosh
  wdti4 = wdti2 / 2.0d0
  wdti8 = wdti4 / 2.0d0

  call md_kinen ! resets pmom2 and en_kin ==> why do we need pmom2 update at the end of the subroutine?

  scale = 1.0D0
  gnh(1) = (2.0d0*en_kin - gfree0*Ekin0)/qmass(1)

  if (ensemble_nvt) then
      do iresn = 1, nresn
          do iyosh = 1, nyosh
              vnh(nnos) = vnh(nnos) + gnh(nnos) * wdti4(iyosh)

              ! Update thermostat velocities
              do inos = nnos - 1, 1, -1
                  aa = dexp(-wdti8(iyosh) * vnh(inos + 1))
                  vnh(inos) = vnh(inos) * aa * aa + wdti4(iyosh) * gnh(inos) * aa
              end do 
              ! Update thermostat coordinates
              do inos = 1, nnos
                  xnh(inos) = xnh(inos) + vnh(inos) * wdti2(iyosh)
              end do

              scale = scale * dexp(-wdti2(iyosh) * vnh(1))
              gnh(1) = (2.0d0*en_kin*scale*scale - gfree0*Ekin0)/qmass(1)

              do inos = 1, nnos - 1
                  aa = dexp(-wdti8(iyosh) * vnh(inos + 1))
                  vnh(inos) = vnh(inos) * aa * aa + wdti4(iyosh) * gnh(inos) * aa
                  gnh(inos + 1) = (qmass(inos) * vnh(inos) * vnh(inos) - Ekin0)/qmass(inos + 1)
              end do

              vnh(nnos) = vnh(nnos) + gnh(nnos) * wdti4(iyosh)
          end do
      end do
  end if

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
  integer :: iat, i
  
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
     
     do i = 1, nnos
         write (10, *) xnh(i), gnh(i), vnh(i), qmass(i)
     end do

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
     write (45, *) natom
     write (45, '(A,f14.6,A,f10.2,A,f9.3)') &
          '  Upot ', Upot, '  Temp ', temp, &
          '  Time(ps) ', time_ps
     
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
  integer :: iat, i

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

     do i = 1, nnos
         read (10, *) xnh(i), gnh(i), vnh(i), qmass(i)
     end do

     ! Pressure
     read (10, *) xlogv, glogv, vlogv, pmass
     
     close(10)

     call md_kinen

  end if

  call qc_mpi_barrier

  call qc_mpi_int_bcast1 (istart)
  !call qc_mpi_real_bcast (at_chg, natom)
  call qc_mpi_real_bcast (d_rr, 3*natom) ! ?? Propagation is always performed by master thread

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
  real*8 :: vol, qkin, qpot  ! kinetic and potential energy of therm
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
        qkin = 0.5d0 * sum(qmass * vnh * vnh)
        qpot = gfree0* Ekin0 * xnh(1) + Ekin0 * sum(xnh(2:nnos))
        Ten = Ten + qkin + qpot      ! Conserved quantity for N-H dynamics
     end if
     
     !write (6,'(I8,F10.3,5F14.6,F10.2)') &
     write (6,'(I8,F12.5,5F14.6,F10.2,F12.2)') & ! AK : modified format
          istep, time_ps,Upot/Nmol, Umon/Nmol, Udim/Nmol, Ulr/Nmol, &
          Ten/Nmol, temp, press
     !  end if
     
     if (is_temp_tol) then ! Berendsen thermostat -- shouldn't be used except for equilibration runs
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
    !real*8 :: tr


    gm_pos = matmul(lati, gm_pos) ! scale coordinates
    
    do im = 1, nmol ! fold molecules back into central cell
        iat   = mol_iatom(im)
        nsite = mol_nsite(im)
        
        ! X
        !if (gm_pos(1,iat+1) .ge. 1.0d0) then
        !    tr = 1.0d0 * int(gm_pos(1,iat+1))
        !    do is = 1, nsite
        !        gm_pos(1,iat+is) = gm_pos(1,iat+is) - tr
        !    end do
        !end if 
        !if (gm_pos(1,iat+1) .lt. 0.0d0) then
        !    tr = 1.0d0 * int(gm_pos(1,iat+1)) + 1.0
        !    do is = 1, nsite
        !        gm_pos(1,iat+is) = gm_pos(1,iat+is) - tr 
        !    end do
        !end if  
        !
        ! Y
        !if (gm_pos(2,iat+1) .lt. 0.0d0 .or. gm_pos(2,iat+1) .ge. 1.0d0) then
        !    tr = 1.0d0 * int(gm_pos(2,iat+1))
        !    do is = 1, nsite
        !        gm_pos(2,iat+is) = gm_pos(2,iat+is) - tr
        !    end do
        !end if 
        !if (gm_pos(2,iat+1) .lt. 0.0d0) then
        !    tr = 1.0d0 * int(gm_pos(2,iat+1)) + 1.0
        !    do is = 1, nsite
        !        gm_pos(2,iat+is) = gm_pos(2,iat+is) - tr 
        !    end do
        !end if  
        !
        ! Z
        !if (gm_pos(3,iat+1) .lt. 0.0d0 .or. gm_pos(3,iat+1) .ge. 1.0d0) then
        !    tr = 1.0d0 * int(gm_pos(3,iat+1))
        !    do is = 1, nsite
        !        gm_pos(3,iat+is) = gm_pos(3,iat+is) - tr
        !    end do
        !end if 
        !if (gm_pos(3,iat+1) .lt. 0.0d0) then
        !    tr = 1.0d0 * int(gm_pos(3,iat+1)) + 1.0
        !    do is = 1, nsite
        !        gm_pos(3,iat+is) = gm_pos(3,iat+is) - tr 
        !    end do
        !end if  

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
