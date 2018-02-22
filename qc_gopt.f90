subroutine qc_gopt_run (en_pot)
use qc_system
use qc_bim
use qc_mpi
use qc_lattice
use consts
implicit none
real*8, intent(out) :: en_pot

integer, parameter :: maxiter = 50
integer, parameter :: nup = 7
real*8, parameter :: d_max_tol = 1.0d-3
real*8, parameter :: d_rms_tol = 5.0d-4
real*8, parameter :: g_max_tol = 4.5d-4
real*8, parameter :: g_rms_tol = 8.0d-5
real*8, parameter :: x_tol     = 1.0d-16

integer, dimension(2) :: iprint = (/-1, -1/)
character(len=30) :: fname, fname_final

integer :: ia, ii, nvar, nwork, iter, iflag, iout, stallCount
real*8  :: umon, udim, ulr
real*8  :: g_max, g_rms, d_max, d_rms
real*8  :: dx, dy, dz
real*8  :: press, vol, en_pot0
real*8, dimension(:), allocatable :: xx, gxx
real*8, dimension(:), allocatable :: dm, ws
real*8, dimension(:,:), allocatable :: gm_pos0

!--- Print Header Information
iout = 40
fname = 'gopt.xyz'
fname_final = fconfig(1:index(fconfig,'.')-1)//'.opt'

if (sys_master) then
    open (unit=iout,file='gopt.dat', access='append', status='unknown')
    write (iout,'(A)') &
    ' Limited Memory BFGS Quasi-Newton Optimization :'
    write (iout,'(A,A)') &
    '====================================================', &
    '============================'
    write (iout,'(A)') &
    '    step     energy [au]     max displ     rms displ     max force     rms force'
    write (iout,'(24x,4f14.6)') &
    d_max_tol, d_rms_tol, g_max_tol, g_rms_tol
    write(iout,'(A,A)') &
    '----------------------------------------------------', &
    '----------------------------'

    call flush (iout)
end if

!---- (initiate memory)
nvar   = 3*natom
nwork  = nvar*(2*nup + 1) + 2*nup

allocate (xx(nvar))
allocate (gxx(nvar))
allocate (gm_pos0(3, natom))
allocate (dm(nvar))
allocate (ws(nwork))
!----

call qc_mpi_barrier
call qc_mpi_real_bcast (gm_pos, 3*natom)

iter = 0
stallCount = 0
iflag = 0
en_pot = 0.0d0

call qc_bim_pot (umon, udim, ulr)
en_pot = umon + udim + ulr
vol = volume(lat)*ang2bohr**3
press = virt(1,1)+virt(2,2)+virt(3,3) ! positive if virial has been calculated with forces, not grad
press = press/(3.0d0*vol)*au_bar

g_rms = 0.0d0
g_max = 0.0d0
do ia = 1, natom
    g_rms = g_rms + d_rr(1,ia)*d_rr(1,ia)
    g_rms = g_rms + d_rr(2,ia)*d_rr(2,ia)
    g_rms = g_rms + d_rr(3,ia)*d_rr(3,ia)
    !---
    g_max = max (abs (d_rr(1,ia)), g_max)
    g_max = max (abs (d_rr(2,ia)), g_max)
    g_max = max (abs (d_rr(3,ia)), g_max)
end do
g_rms = sqrt(g_rms/nvar)

if (sys_master) then
    write (iout, '(i8,f16.8,28x,2f14.6)') iter, en_pot, g_max, g_rms
    call flush (iout)
    call gopt_save(iter, vol, en_pot, press, fname)
end if

do iter = 1, maxiter
    
    ! back up gm_pos
    gm_pos0 = gm_pos
    do ia = 1, natom
        ii = 3*(ia-1)
        xx(ii+1) = gm_pos(1,ia)*ang2bohr
        xx(ii+2) = gm_pos(2,ia)*ang2bohr
        xx(ii+3) = gm_pos(3,ia)*ang2bohr
        gxx(ii+1)= d_rr(1,ia)
        gxx(ii+2)= d_rr(2,ia)
        gxx(ii+3)= d_rr(3,ia)
    end do

    ! L-BFGS: take a step
    call lbfgs (nvar, nup, xx, en_pot, gxx, .false., dm, & 
    iprint, g_rms_tol, x_tol, ws, iflag)
    if (iflag.lt.0 .and. sys_master) write (*,'(A,i0,A)') "LBFGS ERROR(IFLAG=", iflag, ")"
    if (iflag .lt. 0) exit

    ! update cartesian coordinates
    do ia = 1, natom
        ii = 3*(ia-1)
        gm_pos(1,ia) = xx(ii+1)*bohr2ang
        gm_pos(2,ia) = xx(ii+2)*bohr2ang
        gm_pos(3,ia) = xx(ii+3)*bohr2ang
    end do

    call qc_mpi_barrier
    call qc_mpi_real_bcast (gm_pos, 3*natom)
    
    ! check RMS/MAX shift in coordinates (angstrom)
    d_rms = 0.0d0
    d_max = 0.0d0
    do ia = 1, natom
        dx = gm_pos(1,ia) - gm_pos0(1,ia)
        dy = gm_pos(2,ia) - gm_pos0(2,ia)
        dz = gm_pos(3,ia) - gm_pos0(3,ia)
        d_rms = d_rms + dx*dx + dy*dy + dz*dz
        d_max = max (abs(dx), d_max)
        d_max = max (abs(dy), d_max)
        d_max = max (abs(dz), d_max)
    end do
    d_rms = sqrt(d_rms/nvar)

    if (d_max .lt. d_max_tol) then
        stallCount = stallCount + 1
    else
        stallCount = 0
    end if

    ! Evaluate enthalpy at new geometry
    if ((iflag.gt.0) .or. (d_max.gt.1.0d-5)) then
        call qc_bim_pot (umon, udim, ulr)
        en_pot = umon + udim + ulr
        vol = volume(lat)*ang2bohr**3
        press = virt(1,1)+virt(2,2)+virt(3,3) ! positive if virial has been calculated with forces, not grad
        press = press/(3.0d0*vol)*au_bar
    end if

    ! check RMS/MAX gradient (atomic unit)
    g_rms = 0.0d0
    g_max = 0.0d0
    do ia = 1, natom
        g_rms = g_rms + d_rr(1,ia)*d_rr(1,ia)
        g_rms = g_rms + d_rr(2,ia)*d_rr(2,ia)
        g_rms = g_rms + d_rr(3,ia)*d_rr(3,ia)
        !---
        g_max = max (abs (d_rr(1,ia)), g_max)
        g_max = max (abs (d_rr(2,ia)), g_max)
        g_max = max (abs (d_rr(3,ia)), g_max)
    end do
    g_rms = sqrt(g_rms/nvar)

    ! Report current iteration
    if (sys_master) then
        write (iout, '(i8,f16.8,4f14.6)') iter, en_pot, d_max, d_rms, g_max, g_rms
        call flush (iout)
        call gopt_save(iter, vol, en_pot, press, fname)
    end if

    ! Check convergence
    if ((d_max.lt.d_max_tol) .and. (g_max.lt.g_max_tol) .and. &
        (d_rms.lt.d_rms_tol) .and. (g_rms.lt.g_rms_tol)) then
        iflag = 0
    end if

    ! Exit if tightly converged
    if (iflag .eq. 0) then
        if (sys_master) then
            write(*,*) "    IFLAG=0. L-BFGS is converged."
            if (l_freezecell) call gopt_save(iter, vol, en_pot, press, fname_final)
        end if
        exit
    else if (stallCount .eq. 3) then 
        ! Exit if stalling and normally converged
        if (g_max .lt. atom_gmax) then
            iflag = 0
            if (sys_master) then
                write(*,*) "    Unchanging geometry; acceptable g_max. L-BFGS is converged."
                if (l_freezecell) call gopt_save(iter, vol, en_pot, press, fname_final)
            end if
            exit
        ! Restart if stalling and not converged
        else
            if (sys_master) write(*,*) "    Unchanging geometry; restarting L-BFGS from current geometry"
            iflag = 0
            stallCount = 0
            dm = 0.0d0
            ws = 0.0d0
        end if
    end if

end do
     
if (g_max .lt. atom_gmax) iflag = 0
if (iflag.ne.0 .and. sys_master) write(*,*) "    Warning: GOPT did not converge!"

if (sys_master) then
     write (6,'(A)') "         STRESS TENSOR (BAR)   "
     write (6,'(A)') "        ---------------------  "
     write (6,'(3f12.1)') virt(1,1:3)/vol*au_bar
     write (6,'(3f12.1)') virt(2,1:3)/vol*au_bar
     write (6,'(3f12.1)') virt(3,1:3)/vol*au_bar
     write (6,*) ""
end if

deallocate (xx, gxx, gm_pos0)
deallocate (dm, ws)
if (sys_master) close (iout)

end subroutine qc_gopt_run

subroutine gopt_save (iter, vol, en_pot, press, fname)
    use qc_system
    implicit none
    integer, intent(in) :: iter
    real*8,  intent(in) :: vol, en_pot, press
    character(len=30), intent(in) :: fname
    integer :: ia

    open (45, file=trim(fname), access='append', status='unknown')
    write (45, '(i0)') natom
    write (45, '(f0.8, 5f13.8, i3, A3, f16.8, A3, f16.8, A6, f16.8, A6, i3, A)') lat_a, lat_b, lat_c, lat_alpha, &
    lat_beta, lat_gamma, lat_axis, 'V=', vol, 'E=', en_pot, 'P=', press, '(Iter', iter, ')'
    
    do ia = 1, natom
        write (45, '(A,3f16.8)') at_atnm(ia), gm_pos(1:3,ia)
    end do
    close (45)
end subroutine gopt_save
