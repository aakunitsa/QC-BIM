! LBFGS Cell Optimization
subroutine qc_cell_gopt_run_lbfgs
    use qc_system
    use qc_mpi
    use qc_bim
    use qc_lattice
    use consts
    implicit none
    external LB2
    
    real*8 :: umon, udim, ulr, en_pot, enthalpy, vol 

    integer :: im, ia, is, ii
    
    integer, parameter :: nup = 30
    integer :: nvar, nwork ! nwork = nvar*(2*nup + 1) + 2*nup
    real*8, dimension(:), allocatable :: xx, gxx, diagm, workspace
    real*8, dimension(:,:), allocatable :: gm_pos0
    real*8 :: d_rr0(3), g_lat(6)
    integer, dimension(2) :: iprint = (/-1, -1/)
    integer :: iflag, iter
    character(len=20) :: fname
    
    real*8, parameter :: g_eps     = 1.0d-5
    real*8, parameter :: x_tol     = 1.0d-16
    
    ! atomic gradient convergence criteria
    real*8 ::  g_rms, g_max

    iflag = 0
    nvar = 3*natom + 6
    nwork = nvar*(2*nup + 1) + 2*nup
    fname = 'cell_gopt.xyz'

    allocate(xx(nvar))
    allocate(gxx(nvar))
    allocate(diagm(nvar))
    allocate(workspace(nwork))
    allocate(gm_pos0(3,natom))

    !--Translate center of main atom masses to origin
    call center_cell
    
    ! Iteration 0
    if (sys_master) write(*,'(A)') "CELL OPT: iter 0"
    iter = 0
    
    ! Evaluate enthalpy, gradients
    call qc_bim_pot (umon, udim, ulr)
    en_pot = umon + udim + ulr
    vol = volume(lat)*ang2bohr**3
    enthalpy = en_pot + (pext0/au_bar)*vol
    call lat_grad(g_lat)

    ! check RMS/MAX FORCE
    call grad_criteria(g_max, g_rms)
    
    !--- Report to stdout and write geom to cell_gopt.xyz file
    call report_iter_lbfgs(iter, en_pot, vol, enthalpy, g_max, g_rms, fname)

    do iter = 1, opt_maxiter
        gm_pos0 = gm_pos 
        
        ! Bohr
        xx(1) = lat_a*ang2bohr
        xx(2) = lat_b*ang2bohr
        xx(3) = lat_c*ang2bohr
        gxx(1:3) = g_lat(1:3) ! Bohr^-1

        ! Degrees
        xx(4) = lat_alpha
        xx(5) = lat_beta
        xx(6) = lat_gamma
        gxx(4:6) = g_lat(4:6) ! Degree^-1

        ! scaled/internal coordinates
        ii = 6
        do im = 1,nmol
            is = mol_iatom(im)
            xx(ii+1:ii+3) = matmul(lati,gm_pos(1:3, is+1)) ! dimensionless
            d_rr0 = 0.0d0
            do ia = is+1,is+mol_nsite(im)
                d_rr0 = d_rr0(:) + d_rr(1:3, ia)
            end do
            gxx(ii+1:ii+3) = matmul(lat*ang2bohr,d_rr0) ! dimensionless

            ii = ii + 3
            do ia = is+2,is+mol_nsite(im)
                xx(ii+1:ii+3) = (gm_pos(1:3, ia) - gm_pos(1:3, is+1))*ang2bohr !  bohr
                gxx(ii+1:ii+3) = d_rr(1:3, ia) ! bohr^-1
                ii = ii + 3
            end do
        end do

        call lbfgs(nvar, nup, xx, enthalpy, gxx, .false., diagm, &
                   iprint, g_eps, x_tol, workspace, iflag)

        if (iflag.lt.0 .and. sys_master) print *, "L-BFGS IFLAG<0"
        
        ! Update cell
        lat_a = xx(1)*bohr2ang
        lat_b = xx(2)*bohr2ang
        lat_c = xx(3)*bohr2ang
        lat_alpha = xx(4)
        lat_beta = xx(5)
        lat_gamma = xx(6)
        
        call create_lattice_vectors
        call update_scaling_matrix
    
        ! update cartesian coordinates
        ii = 6
        do im = 1,nmol
            is = mol_iatom(im)
            gm_pos(1:3,is+1) = matmul(lat,xx(ii+1:ii+3))

            ii = ii + 3
            do ia = is+2, is+mol_nsite(im)
                gm_pos(1:3,ia) = gm_pos(1:3,is+1) + xx(ii+1:ii+3)*bohr2ang
                ii = ii + 3
            end do
        end do

        call qc_mpi_real_bcast1 (lat_a)
        call qc_mpi_real_bcast1 (lat_b)
        call qc_mpi_real_bcast1 (lat_c)
        call qc_mpi_real_bcast1 (lat_alpha)
        call qc_mpi_real_bcast1 (lat_beta)
        call qc_mpi_real_bcast1 (lat_gamma)
        call qc_mpi_real_bcast1 (vol)
        call qc_mpi_real_bcast (lat, 9)
        call qc_mpi_real_bcast (lati, 9)
        call qc_mpi_real_bcast (gm_pos, 3*natom)
        
        ! reevaluate enthalpy, gradient at new geometry
        call qc_bim_pot (umon, udim, ulr)
        en_pot = umon + udim + ulr
        vol = volume(lat)*ang2bohr**3
        enthalpy = en_pot + (pext0/au_bar)*vol
        call lat_grad(g_lat)

        ! check convergence: RMS/MAX FORCE
        call grad_criteria(g_max, g_rms)

        !--- Report to stdout and write geom to cell_gopt.xyz file
        call report_iter_lbfgs(iter, en_pot, vol, enthalpy, g_max, g_rms, fname)
        
        ! Quit if converged OR L-BFGS failed
        if (g_max.lt.atom_gmax .and. maxval(abs(g_lat)).lt.lat_gmax) iflag = 0
        if (iflag .le. 0) exit
    end do
    
    if (iflag .ne. 0) then
        if (sys_master) print *, 'DID NOT CONVERGE CELL GEOMETRY!'
    else
        if (sys_master) print *, 'CONVERGED CELL GEOMETRY'
    end if

    deallocate(xx, gxx, diagm, workspace, gm_pos0)
end subroutine qc_cell_gopt_run_lbfgs

subroutine report_iter_lbfgs (iter, en_pot, vol, enthalpy, gmax, grms, fname)
    use qc_system
    use qc_lattice
    use qc_mpi
    use consts
    implicit none
    integer, intent(in) :: iter
    real*8, intent(in) :: en_pot, vol, enthalpy, gmax, grms
    character(len=20), intent(in) :: fname
    real*8 :: latgrad(6)

    call lat_grad(latgrad)

    if (sys_master) then
        ! Write header
        if (iter .eq. 0) then
            write(*,'(A)')  '@'
            write(*,'(A)') '@ Cell Optimize: Full L-BFGS'
            write(*,'(A)') '@ H in a.u., lattice params in Angstroms, degrees'
            write(*,'(A)') '@ ================================================='
            
            write(*,'(A, A5, A14, 6A9, 8A10)') '@', 'step', 'enthalpy', & 
            'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'da', 'db', 'dc', & 
            'dalpha', 'dbeta', 'dgamma', 'max grad', 'rms grad'
            
            write(*,'(A,A,A)') '@----------------------------------------------',&
            '----------------------------------------------------------------',&
            '-------------------------------------------'
        end if

        ! Write optimization data
        write(*,'(/,A)') '-----------------------'
        write (*,'(A)') 'CELL PARAMS/GRADIENTS:'
        write (*, '(A, i5, f14.6, 3f9.3, 3f9.2, 6f10.4,2f10.6)') '@', iter, enthalpy, &
        lat_a, lat_b, lat_c, lat_alpha, lat_beta, lat_gamma, latgrad(1)/bohr2ang, &
        latgrad(2)/bohr2ang, latgrad(3)/bohr2ang, latgrad(4:6), gmax, grms
        write(*,'(A,/)') '---------------------'

        ! Write full geometry/data to file
        call gopt_save (iter, en_pot, vol, enthalpy, fname)
    end if
end subroutine report_iter_lbfgs

subroutine grad_criteria(g_max, g_rms)
    use qc_system
    real*8, intent(out) :: g_max, g_rms
    integer :: ia
    g_max = 0.0d0
    g_rms = 0.0d0
    do ia = 1, natom
        g_rms = g_rms + d_rr(1,ia)*d_rr(1,ia)
        g_rms = g_rms + d_rr(2,ia)*d_rr(2,ia)
        g_rms = g_rms + d_rr(3,ia)*d_rr(3,ia)
        !---
        g_max = max (abs (d_rr(1,ia)), g_max)
        g_max = max (abs (d_rr(2,ia)), g_max)
        g_max = max (abs (d_rr(3,ia)), g_max)
    end do
    g_rms = sqrt(g_rms/(3*natom))
end subroutine grad_criteria
