! Gradient Descent Cell Optimization
subroutine qc_cell_gopt_run
use qc_system
use qc_mpi
use qc_lattice
use consts
implicit none
integer :: iter, ia, ii, is, im
real*8  :: press, vol, en_pot, enthalpy, pext0au, pext0_tensor(3,3)
real*8  :: alpha, p(6), p0(6), lat0(6), pnorm, pnorm0, en_pot0, enthalpy0
logical :: lowp(3), highp(3), reset
real*8, dimension(:,:), allocatable :: gm_pos0
character(len=30) :: fname, fname_final

allocate(gm_pos0(3,natom))
fname = 'cell_gopt.xyz'
fname_final = fconfig(1:index(fconfig,'.')-1)//'.cellopt'
pext0_tensor = 0.0d0
pext0_tensor(1,1) = pext0
pext0_tensor(2,2) = pext0
pext0_tensor(3,3) = pext0
lowp = .false.
highp = .false.

alpha = 1.0d0
iter = 0

!--Translate center of main atom masses to origin
call center_cell

! Iteration 0
if (sys_master) print *, 'iter', iter

! Relax atom pos, get energy/gradient, enthalpy
en_pot = 0.0d0
call qc_gopt_run (en_pot)
pext0au = pext0/au_bar
vol = volume(lat)*ang2bohr**3
enthalpy = en_pot + pext0au*vol
press = virt(1,1)+virt(2,2)+virt(3,3) ! positive if virial has been calculated with forces, not grad
press = press/(3.0d0*vol)*au_bar

! search direction
call lat_grad(p)
p = -1.0d0 * p
p(1:3) = p(1:3) / bohr2ang
pnorm = norm(p)

! Pressure: E_h/Bohr**3 --> Bar
call report_iter(iter, en_pot, vol, press, enthalpy, fname)

iter = 1
do while ( (maxval(abs(p(1:3))).gt.lat_gmax .or. maxval(abs(p(4:6))).gt.1.0d-2) .and. iter.lt.opt_maxiter )
    if (sys_master) write (*,'(A,i0,A,f8.4,A)') 'Iter ', iter, ' (alpha ', alpha, ')'

    ! backup lattice constants
    lat0(1) = lat_a
    lat0(2) = lat_b
    lat0(3) = lat_c
    lat0(4) = lat_alpha
    lat0(5) = lat_beta
    lat0(6) = lat_gamma

    ! backup lattice gradient
    p0 = p
    pnorm0 = pnorm
    enthalpy0 = enthalpy
    en_pot0 = en_pot
    
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
    if (a_dim) lat_a      = lat_a + 1.0d0*alpha*p(1)
    if (b_dim) lat_b      = lat_b + 1.0d0*alpha*p(2)
    if (c_dim) lat_c      = lat_c + 1.0d0*alpha*p(3)
    if (maxval(abs(p(4:6))).gt. 1.0d-2) then
        if (a_dim .and. b_dim) lat_alpha  = lat_alpha + 10.0d0*alpha*p(4)
        if (a_dim .and. b_dim .and. c_dim) lat_beta   = lat_beta  + 10.0d0*alpha*p(5)
        if (a_dim .and. b_dim .and. c_dim) lat_gamma  = lat_gamma + 10.0d0*alpha*p(6)
    end if
    call update_geometry(gm_pos0)

    ! Evaluate enthalpy/gradients at new geometry
    en_pot = 0.0d0
    call qc_gopt_run (en_pot)
    pext0au = pext0/au_bar
    vol = volume(lat)*ang2bohr**3
    enthalpy = en_pot + pext0au*vol
    press = virt(1,1)+virt(2,2)+virt(3,3) ! positive if virial has been calculated with forces, not grad
    press = press/(3.0d0*vol)*au_bar
    
    call lat_grad(p)
    p = -1.0d0 * p
    p(1:3) = p(1:3) / bohr2ang
    pnorm = norm(p)
    reset = .false.
    do ii = 1, 3
        if (p(ii) < 0.0d0) lowp(ii) = .true.
        if (p(ii) > 0.0d0) highp(ii) = .true.
        reset = reset .or. (lowp(ii) .and. highp(ii))
    end do

    if (reset) then
        alpha = max(alpha*0.6d0, 0.5d0)
        reset = .false.
        lowp = .false.
        highp = .false.
    else if ( (pnorm.lt.pnorm0) .or. (((enthalpy-enthalpy0)/nmol).lt.2.0d-5) .or. (alpha .le.0.51d0)) then
        alpha = min(alpha*1.20d0,5.0d0)
    else
        lat_a = lat0(1)
        lat_b = lat0(2)
        lat_c = lat0(3)
        lat_alpha = lat0(4)
        lat_beta = lat0(5)
        lat_gamma = lat0(6)
        p = p0
        pnorm = pnorm0
        en_pot = en_pot0
        enthalpy = enthalpy0
        alpha = max(alpha*0.5d0, 0.5d0)
        call update_geometry(gm_pos0)
        vol = volume(lat)*ang2bohr**3
    end if
    
    ! Calc new virial press & report
    call report_iter(iter, en_pot, vol, press, enthalpy, fname)

    iter = iter + 1
end do

! Report success/failure
if (maxval(abs(p(1:3))) .lt. lat_gmax) then
    if (sys_master) print *, "CELL OPTIMIZATION CONVERGED"
    if (sys_master) call gopt_save (iter, vol, en_pot, press, fname_final)
else
    if (sys_master) print *, "CELL OPTIMIZATION DID NOT CONVERGE"
end if

deallocate (gm_pos0)

contains
function norm(mat)
    real*8, intent(in) :: mat(6)
    real*8 :: norm
    integer :: i
    norm = 0.0d0
    do i = 1,6
        norm = norm + mat(i)**2
    end do
    norm = sqrt(norm)
end function

end subroutine qc_cell_gopt_run

subroutine report_iter (iter, en_pot, vol, press, enthalpy, fname)
    use qc_system
    use qc_lattice
    use qc_mpi
    use consts
    implicit none
    integer, intent(in) :: iter
    real*8, intent(in) :: en_pot, vol, press, enthalpy
    character(len=30), intent(in) :: fname
    real*8 :: latgrad(6)

    call lat_grad(latgrad)

    if (sys_master) then
        ! Write header
        if (iter .eq. 0) then
            write(*,'(A)')  '@'
            write(*,'(A)') '@ Cell Optimize: Iterative gradient descent, L-BFGS'
            write(*,'(A)') '@ H in a.u., lattice params in Angstroms, degrees'
            write(*,'(A)') '@ ================================================='
            
            write(*,'(A, A5, A14, 6A9, 6A10)') '@', 'step', 'enthalpy', & 
            'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'da', 'db', 'dc', & 
            'dalpha', 'dbeta', 'dgamma'
            
            write(*,'(A,A,A)') '@----------------------------------------------',&
            '----------------------------------------------------------------',&
            '-----------------------'
        end if

        ! Write optimization data
        write(*,'(/,A)') '-----------------------'
        write (*,'(A)') 'CELL PARAMS/GRADIENTS:'
        write (*, '(A, i5, f14.6, 3f9.3, 3f9.2, 6f10.4)') '@', iter, enthalpy, &
        lat_a, lat_b, lat_c, lat_alpha, lat_beta, lat_gamma, latgrad(1)/bohr2ang, &
        latgrad(2)/bohr2ang, latgrad(3)/bohr2ang, latgrad(4:6)
        write(*,'(A,/)') '---------------------'

        ! Write full geometry/data to file
        call gopt_save (iter, vol, en_pot, press, fname)
    end if
end subroutine report_iter

subroutine lat_grad (lat_constants_grad)
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
end subroutine lat_grad

subroutine update_geometry (gm_pos0)
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
end subroutine update_geometry

subroutine center_cell
use qc_system
use qc_mpi
implicit none
real*8 :: center(3)
integer :: im, ia
center = 0.0d0
do im=1,nmol
    ia = mol_iatom(im)
    center(:) = center(:) + gm_pos(:, ia+1)
end do
center = center / nmol
do ia=1,natom
    gm_pos(:,ia) = gm_pos(:,ia) - center(:)
end do
call qc_mpi_real_bcast (gm_pos, 3*natom)
end subroutine center_cell
