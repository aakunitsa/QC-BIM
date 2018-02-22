! Mass-weighted Cartesian Force constants matrix 
! Dipole derivatives
module hessian
use qc_system
use consts
use qc_mpi
use qc_neigh
use qc_bim
use chem_interface
implicit none

real*8, dimension(:,:,:,:,:), allocatable :: hess ! Cartesian force constants (a.u.)
real*8, dimension(:,:), allocatable :: d_dipole   ! Dipole derivative
real*8 :: hess0(171) ! only large enough for 3 atom fragments (must increase this)
real*8 :: ddipole0(54)

!-----------------------------------------------
! Finite diff delta parameters in qc_system.f90
! Using defaults from nwchem and center diff
! df/dx = [f(x+delta) - f(x-delta)] / (2*delta)
!-----------------------------------------------

contains

! ------------------------------------------
! qc_hess: main routine to allocate memory,
! calculate lists, ESP, monomer, dimer, LR
! contributions to Hessian.  Then report results
! and clean up
! --------------------------------------------
subroutine qc_hess
    logical :: list_fail
    integer*8 :: cr, cm, c1, c2, elapsed!, im
    real*8 :: rate

    call system_clock(count_rate=cr)
    call system_clock(count_max=cm)
    rate = REAL(cr)

    ! Start timing
    if (sys_master) write(*,'(A)',advance="no") '   BIM_HESS: lists...' 
    call system_clock(c1)

    ! Get neighbor lists (uses gradient subroutine)
    call qc_get_lists(list_fail)
    if (list_fail) then
        call qc_allocate_lists
        call qc_get_lists(list_fail)
        if (list_fail) STOP "UNEXPECTED LIST OVERFLOW"
    end if


    ! Initialize force constants matrix and dipole derivatives
    allocate( hess(3*natom, 3*natom, -n_hess(1):n_hess(1), -n_hess(2):n_hess(2), -n_hess(3):n_hess(3)) )
    allocate( d_dipole(3*natom, 3) )
    hess = 0.0d0
    d_dipole = 0.0d0

    ! Calculate monomer-SCF charges (uses gradient subroutine)
    if (sys_master) write(*,'(A)',advance="no") 'EMBEDDING...' 
    call qc_bim_field_update
    

    ! LR Contribution to Hessian
    if (l_bq) then
        if (sys_master)  write(*,'(A)',advance="no") 'COULOMB...'
        call qc_hess_coulomb
    end if
    
    ! Monomer contributions
    if (sys_master) write(*,'(A)',advance="no") 'mono...' 
    call qc_hess_monomers
    
    ! Dimer contributions
    if (sys_master) write(*,'(i0,A)',advance="no") num_pairs, ' dimers...'
    call qc_hess_dimers
    
    ! End time
    call system_clock(c2)
    elapsed = int(anint((c2-c1)/rate))
    if (sys_master) write(*,'(A,i0,A,i0,A)') "wall: ",elapsed/60,"m ",&
                                        mod(elapsed,60),"s"

    ! MPI Reduce Hessian, Dipole Derivatives and Report
    call qc_mpi_allreduce (hess, (9*natom**2)*(2*n_hess(1)+1)*(2*n_hess(2)+1)*(2*n_hess(3)+1))
    call qc_mpi_allreduce (d_dipole, 3*3*natom)
    if (sys_master) call report_hess
    
    deallocate ( hess )
    deallocate ( d_dipole )
end subroutine qc_hess


! -----------------------------------------------------
! qc_hess_lr: MPI ranks split the monomers of cell 0
! and accumulate the monomer--LR contributions into the 
! global Hessian
! -----------------------------------------------------
subroutine qc_hess_coulomb
    integer :: im

    do im = 1 + sys_me, nmol, sys_nproc
        call get_coulomb_hess(im) ! adds lr contributions to hessian
    end do
    call qc_mpi_barrier
end subroutine qc_hess_coulomb


! --------------------------------------------------------
! qc_hess_monomers: MPI ranks split the monomers of cell 0
! and accumulate the monomer contributions into the global
! Hessian/dipole derivative.  
! Each monomer energy is differentiated wrt ALL
! coordinates in cell 0 -- including both QM and BQ atoms
! --------------------------------------------------------
subroutine qc_hess_monomers

    integer :: imon, bq_idx, ibq, is, is_bq, at_qm, at_bq, n(3)
    real*8 :: k_esp(3,3)
    
    ! MPI run over monomers
    do imon = 1 + sys_me, nmol, sys_nproc

        call get_mono_hess(imon) 

        ! Run over QM atoms
        do is = 1, mol_nsite(imon)
            
            at_qm = mol_iatom(imon)+is
            
            ! Run over BQ molecules
            do bq_idx = 1, nbq(imon)
                ibq  = bq_list(bq_idx, imon)%jm
                n(1) = bq_list(bq_idx, imon)%ja
                n(2) = bq_list(bq_idx, imon)%jb
                n(3) = bq_list(bq_idx, imon)%jc

                ! Run over BQ atoms
                do is_bq = 1, mol_nsite(ibq)
                    at_bq = mol_iatom(ibq)+is_bq
                    call esp_force_consts(at_qm, at_bq, n, k_esp)
                    call get_bqhess_ii(at_qm, -1.0d0, k_esp)
                end do 
            end do 
        end do ! end run over atoms
    
    end do ! end run over monomers
    call qc_mpi_barrier
end subroutine qc_hess_monomers

! -----------------------------------------------------------------------
! qc_hess_dimers: MPI ranks split all the dimers contributing
! to unit cell energy, and accumulate dimer contributions into
! the global Hessian/dipole derivative.  (See comment for get_dimer_hess)
! -----------------------------------------------------------------------
subroutine qc_hess_dimers
    integer :: idx_dimer, ierr, tmp, done
    integer, dimension(MPI_STATUS_SIZE) :: stat

    tmp = 0
    done = -1
    ! PARALLEL
    if (sys_nproc .gt. 1) then
        ! MASTER
        if (sys_master) then
            do idx_dimer = sys_nproc, 3*num_pairs
                call mpi_recv(tmp, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                              0, MPI_COMM_WORLD, stat, ierr)
                call mpi_send(idx_dimer, 1, MPI_INTEGER, stat(MPI_SOURCE), &
                              0, MPI_COMM_WORLD, ierr)
            end do
            
            do idx_dimer = 1, sys_nproc-1
                call mpi_recv(tmp, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                              0, MPI_COMM_WORLD, stat, ierr)
                call mpi_send(done, 1, MPI_INTEGER, stat(MPI_SOURCE), & 
                              0, MPI_COMM_WORLD, ierr)
            end do
        ! WORKER
        else
            idx_dimer = sys_me
            if (idx_dimer .le. 3*num_pairs) call get_dimer_hess(idx_dimer)
            call mpi_send (tmp, 1, MPI_INTEGER, 0, 0, &
                           MPI_COMM_WORLD, ierr)
            call mpi_recv (idx_dimer, 1, MPI_INTEGER, 0, 0, &
                           MPI_COMM_WORLD, stat, ierr)
            do while (idx_dimer .ne. done)
                call get_dimer_hess(idx_dimer)
                call mpi_send (tmp, 1, MPI_INTEGER, 0, 0, &
                               MPI_COMM_WORLD, ierr)
                call mpi_recv (idx_dimer, 1, MPI_INTEGER, 0, 0, &
                               MPI_COMM_WORLD, stat, ierr)
            end do
        end if
    
    ! SERIAL
    else
        do idx_dimer = 1, 3*num_pairs
            call get_dimer_hess(idx_dimer)
        end do
    end if

    call qc_mpi_barrier
end subroutine qc_hess_dimers


!------------------------------------------
! get_coulomb_hess: use BFS to run over all unit 
! cells and  accumulate monomer--LR Hessian
! Note that even if Hess only contains
! interaction force constants with nearest
! neighbor cells, the LR interaction is
! summed out to the full cutoff range (much
! larger).  These terms still contribute to
! the cell0 -- cell0 block of the Hessian.
!------------------------------------------
subroutine get_coulomb_hess (im)
    use lat_queue
    integer, intent(in) :: im
    type(QUEUE_STRUCT), pointer :: lat_vecs
    type(QUEUE_DATA) :: vec, nvec
    logical :: success
    integer :: nhits

    call queue_create( lat_vecs, 200000 )
    allocate ( vis(-100:100,-100:100,-100:100) )
    vis = .false.
    
    vec%a = 0; vec%b = 0; vec%c = 0
    call queue_append_data(lat_vecs, vec, success)
    vis(0,0,0) = .true.
    
    do while (.not. queue_empty(lat_vecs))
        vec = queue_retrieve_data(lat_vecs)
        call accumulate_hess_coulomb(im, vec, nhits)
        if (nhits .gt. 0) then
            if (a_dim) then
                nvec%a = vec%a+1; nvec%b = vec%b; nvec%c = vec%c
                if (.not. visited(nvec)) call queue_append_data(lat_vecs, nvec, success)
                if (.not. success) STOP "***QUEUE OVERFLOW***"
                
                nvec%a = vec%a-1; nvec%b = vec%b; nvec%c = vec%c
                if (.not. visited(nvec)) call queue_append_data(lat_vecs, nvec, success)
                if (.not. success) STOP "***QUEUE OVERFLOW***"
            end if
            
            if (b_dim) then
                nvec%a = vec%a; nvec%b = vec%b+1; nvec%c = vec%c
                if (.not. visited(nvec)) call queue_append_data(lat_vecs, nvec, success)
                if (.not. success) STOP "***QUEUE OVERFLOW***"
                
                nvec%a = vec%a; nvec%b = vec%b-1; nvec%c = vec%c
                if (.not. visited(nvec)) call queue_append_data(lat_vecs, nvec, success)
                if (.not. success) STOP "***QUEUE OVERFLOW***"
            end if
            
            if (c_dim) then
                nvec%a = vec%a; nvec%b = vec%b; nvec%c = vec%c+1
                if (.not. visited(nvec)) call queue_append_data(lat_vecs, nvec, success)
                if (.not. success) STOP "***QUEUE OVERFLOW***"
                
                nvec%a = vec%a; nvec%b = vec%b; nvec%c = vec%c-1
                if (.not. visited(nvec)) call queue_append_data(lat_vecs, nvec, success)
                if (.not. success) STOP "***QUEUE OVERFLOW***"
            end if
        end if
    end do

    deallocate (vis)
    call queue_destroy (lat_vecs)
end subroutine get_coulomb_hess


! ------------------------------------------------------
! get_mono_atom_hess: for a given monomer, perform
! the Hessian calculation and accumulate contributions 
! into the global Hessian and dipole derivative
! -----------------------------------------------------
subroutine get_mono_hess (idx_monomer)
    integer, intent(in) :: idx_monomer
    integer :: row0, row, col, mu, nsite, i

    nsite = mol_nsite(idx_monomer)
    row0 = idx( mol_iatom(idx_monomer)+1, 1 ) - 1 ! row0+1 is x-coordinate of first atom
    
    call system("mkdir -p "//trim(my_scratch))
    call monomer_wrt(idx_monomer)
    call system (trim(jcmd))
    call read_hess(nsite, hess0)
    call read_ddipole(nsite, ddipole0)
    
    i = 0
    do row = row0+1, row0+3*nsite
    do col = row0+1, row
        i = i + 1
        hess(row,col,0,0,0) = hess(row,col,0,0,0) + hess0(i)
        if (col.lt.row) hess(col,row,0,0,0) = hess(col,row,0,0,0) + hess0(i)
    end do
    end do

    i = 0
    do row = row0+1, row0+3*nsite
    do mu = 1,3
        i = i + 1
        d_dipole(row, mu) = d_dipole(row, mu) + ddipole0(i)
    end do
    end do

end subroutine get_mono_hess


! -----------------------------------------------------------------------
! get_dimer_hess: The dimer contribution  [Ei(0)j(n) - Ei(0) - Ej(n)] 
! is differentiated wrt ALL coordinates in cell 0.
! Due to translational symmetry in the energy
! and gradient calculations, only pairs for which i <= j
! are counted and then given a scale factor of 1 if i<j or i==j and n>(0,0,0).
! In contrast, the Hessian calculation scales each i(0)j(n) (i<j, n!=0) by 1.0
! and differentiates wrt coordinates in cell 0 to get the 
! 1.0*[Ei(0)j(n) - Ei(0) - Ej(n)] contribution.  The derivatives wrt coordinates
! in both cells or coordinates in cell n are copied to obtain the 
! 1.0*[Ej(0)i(-n) - Ej(0) - Ei(-n)] contributions.  This effectively "translates"
! each dimer by -n to obtain the contributions which are not explicitly
! included in the pair lists. 
! ------------------------------------------------------------------------
subroutine get_dimer_hess (idx_dimer)
    integer, intent(in) :: idx_dimer
    integer :: im, ip, calc_id, jm, n(3), n0(3), nsite_i, nsite_j
    integer :: is, js, at_i, at_j, mu, i
    integer :: row, col, row0, col0, row00, col00
    real*8 :: k_esp(3,3)

    im = 1
    calc_id = mod(idx_dimer, 3)
    ip = (idx_dimer+2)/3
    do while (ip .gt. npair(im))
        ip = ip - npair(im)
        im = im + 1
    end do

    call qc_bq2_list_get(im, ip, .false.)

    jm = pair_list(ip,im)%jm
    n(1) = pair_list(ip,im)%ja
    n(2) = pair_list(ip,im)%jb
    n(3) = pair_list(ip,im)%jc
    nsite_i = mol_nsite(im)
    nsite_j = mol_nsite(jm)

    ! *** IMPORTANT ***
    ! The scale factor (scalef) should be one for ALL dimers.
    ! These are second derivatives of total energy, NOT the energy
    ! per unit cell. E_cell isn't even defined in terms of coordinates
    ! that allow breaking the periodic symmetry
    
    if (calc_id .eq. 1) then
        ! Eij
        call system("mkdir -p "//trim(my_scratch))
        call dimer_wrt(im, ip)
        call system (trim(jcmd))
        call read_hess(nsite_i + nsite_j, hess0)
        call read_ddipole(nsite_i + nsite_j, ddipole0)

        ! Read dimer Hessian in lower triangle form
        i = 0
        do row00 = 1, 3*(nsite_i+nsite_j)
        do col00 = 1, row00
            i = i + 1
            if (row00 .le. 3*nsite_i) then
                row0 = idx(mol_iatom(im)+1, 1)-1 ! i(0)i(0)
                col0 = row0
                n0(:) = 0
            else
                row0 = idx(mol_iatom(jm)+1,1)-1 - (3*nsite_i) 
                if (col00 .le. 3*nsite_i) then
                    col0 = idx(mol_iatom(im)+1,1)-1 ! j(0)i(-n)
                    n0(:) = -n
                else
                    col0 = row0 ! j(0)j(0)
                    n0(:) = 0
                end if
            end if

            row = row0 + row00
            col = col0 + col00

            if ((abs(n0(1)) .gt. n_hess(1)) .or. (abs(n0(2)) .gt. n_hess(2)) &
                .or. (abs(n0(3)) .gt. n_hess(3))) cycle

            hess(row,col,n0(1),n0(2),n0(3)) = hess(row,col,n0(1),n0(2),n0(3)) + hess0(i)

            if (row.ne.col .or. (.not.cellzero(n0)) ) then
                hess(col,row,-n0(1),-n0(2),-n0(3)) = hess(col,row,-n0(1),-n0(2),-n0(3)) + hess0(i)
            end if
        end do
        end do
        
        i = 0
        do row = idx(mol_iatom(im)+1,1), idx(mol_iatom(im)+nsite_i,3)
        do mu = 1,3
            i = i + 1
            d_dipole(row, mu) = d_dipole(row, mu) + ddipole0(i)
        end do
        end do

        do row = idx(mol_iatom(jm)+1,1), idx(mol_iatom(jm)+nsite_j,3)
        do mu = 1,3
            i = i + 1
            d_dipole(row, mu) = d_dipole(row, mu) + ddipole0(i)
        end do
        end do
        ! ----------------
    else if (calc_id .eq. 2) then
        ! Ei
        call system("mkdir -p "//trim(my_scratch))
        call monomerbq2_wrt(im, ip, 0, .false.)
        call system (trim(jcmd))
        call read_hess(nsite_i, hess0)
        call read_ddipole(nsite_i, ddipole0)

        i = 0 
        do row = idx(mol_iatom(im)+1, 1), idx(mol_iatom(im)+nsite_i, 3)
        do col = idx(mol_iatom(im)+1, 1), row
            i = i + 1
            hess(row,col,0,0,0) = hess(row,col,0,0,0) - hess0(i)
            if (row.ne.col) then 
                hess(col,row,0,0,0) = hess(col,row,0,0,0) - hess0(i)
            end if
        end do
        end do
        
        i = 0
        do row = idx(mol_iatom(im)+1,1), idx(mol_iatom(im)+nsite_i,3)
        do mu = 1,3
            i = i + 1
            d_dipole(row, mu) = d_dipole(row, mu) - ddipole0(i)
        end do
        end do
        !----------
    else if (calc_id .eq. 0) then
        ! Ej
        call system("mkdir -p "//trim(my_scratch))
        call monomerbq2_wrt(im, ip, 1, .false.)
        call system (trim(jcmd))
        call read_hess(nsite_j, hess0)
        call read_ddipole(nsite_j, ddipole0)

        i = 0
        do row = idx(mol_iatom(jm)+1, 1), idx(mol_iatom(jm)+nsite_j, 3)
        do col = idx(mol_iatom(jm)+1, 1), row
            i = i + 1
            hess(row,col,0,0,0) = hess(row,col,0,0,0) - hess0(i)
            if (row.ne.col) then 
                hess(col,row,0,0,0) = hess(col,row,0,0,0) - hess0(i)
            end if
        end do
        end do
        
        i = 0
        do row = idx(mol_iatom(jm)+1,1), idx(mol_iatom(jm)+nsite_j,3)
        do mu = 1,3
            i = i + 1
            d_dipole(row, mu) = d_dipole(row, mu) - ddipole0(i)
        end do
        end do
        !-------
    else
        STOP "bug in modulo"
    end if

    if (calc_id .eq. 1) then
        ! BQ Interaction: all the point-charge interactions cancel in 
        ! considering Eij - Ei - Ej, EXCEPT for those between i and j
        do is = 1, mol_nsite(im)
            at_i = mol_iatom(im)+is
            do js = 1, mol_nsite(jm)
                at_j = mol_iatom(jm)+js
                
                call esp_force_consts(at_i, at_j, n, k_esp)
                
                if (cellzero(n)) then
                    call get_bqhess_ij(at_i, at_j, n, -1.0d0, k_esp)
                else
                    ! translate by -n...BQ interaction of j(0)i(-n)
                    call get_bqhess_ij(at_i, at_j,  n,  -1.0d0, k_esp)
                    call get_bqhess_ij(at_j, at_i, -n, -1.0d0, k_esp)
                end if
            end do
        end do
    end if

end subroutine get_dimer_hess


! -----------------------------------------------------
! accumulate_hess_coulomb: for a given monomer in cell 0
! and cell n, accumulate the second derivatives of the
! Coulomb interaction into the global Hessian
! -----------------------------------------------------
subroutine accumulate_hess_coulomb(im, vec, nhits)
    integer, intent(in) :: im
    TYPE(QUEUE_DATA), intent(in) :: vec
    integer, intent(out) :: nhits
    integer :: ia, is
    integer :: jm, ja, js, n(3)
    real*8 :: ri(3), rj(3), dcel(3), rij(3), rij2
    real*8 :: k_esp(3,3)

    n(1) = vec%a; n(2) = vec%b; n(3) = vec%c
    dcel(:) = n(1)*lat(:,1)+n(2)*lat(:,2)+n(3)*lat(:,3)

    nhits = 0
    ia = mol_iatom(im)
    
    ! Scan over all molecules in unit cell n=(na, nb, nc)
    do jm = 1, nmol
        if (im .eq. jm .and. cellzero(n)) cycle
        
        ja = mol_iatom(jm)
        
        ri = gm_pos(1:3,ia+1)
        rj = gm_pos(1:3,ja+1) + dcel
        
        rij = rj - ri
        rij2 = dot_product(rij, rij)
        
        if (rij2 < rcut_lr2) then
            
            ! Molecule jm(n) is in QM, BQ, or LR region of im(0)
            nhits = nhits + 1

            ! Run over atoms of im(0)
            do is=1,mol_nsite(im)

                ! Run over atoms of jm(n)
                do js = 1,mol_nsite(jm)
                    
                    call esp_force_consts(ia+is, ja+js, n, k_esp)

                    ! -- i(0)i(0) (Upper left K block) --
                    call get_bqhess_ii(ia+is, 1.0d0, k_esp)
                    if (cellzero(n)) then
                        call get_bqhess_ij(ia+is, ja+js, n, 0.5d0, k_esp)
                    else
                        call get_bqhess_ij(ia+is, ja+js, n, 1.0d0, k_esp)
                    end if

                end do !--end run over atoms of j
            end do !--end run over atoms of i
        end if !--end if j-is-in-LR-region
    end do !--end run over molecules j in cell n
end subroutine accumulate_hess_coulomb


! --------------------------------------------------------
! esp_force_consts: Coulomb of point charges i(0) and j(n)
! The second derivatives of qi*qj/r_ij with respect to
! the six Cartesian coordinates form a 6x6 Hessian, which
! can be divided into four 3x3 quadrants with the form
!   [ [K -K], [-K, K] ]
! This subroutine just calculates K in atomic units
! --------------------------------------------------------
subroutine esp_force_consts (at_qm, at_bq, n, k_esp)
    integer, intent(in) :: at_qm, at_bq, n(3)
    real*8, intent(out) :: k_esp(3,3)
    integer :: mu, nu
    real*8 :: ri(3), rj(3), dcel(3), rij(3), qi, qj, rij_norm
    real*8 :: rij_normi, rij_normi3, rij_normi5
    
    if (embed_id .ne. embed_id_esp) then
        k_esp = 0.0d0
        return
    end if

    ri = gm_pos(1:3, at_qm)*ang2bohr
    qi = chg_pos(0, at_qm)
    
    dcel(:) = n(1)*lat(:,1)+n(2)*lat(:,2)+n(3)*lat(:,3)
    rj = (gm_pos(1:3,at_bq) + dcel)*ang2bohr
    qj = chg_pos(0, at_bq)
    
    rij = (ri - rj) ! points from j to i
    rij_norm = sqrt(dot_product(rij,rij))

    ! DEBUG
    if (rij_norm .lt. 0.1d0) stop 'tiny Rij fed to esp_force_consts'

    rij_normi = 1.0d0 / rij_norm ! avoid repeated divison/exponentiation
    rij_normi3 = rij_normi*rij_normi*rij_normi
    rij_normi5 = rij_normi*rij_normi*rij_normi*rij_normi*rij_normi

    do nu=1,3
    do mu=nu,3
        if (mu .eq. nu) then
            !k_esp(mu,nu) = qi*qj/(rij_norm**3) * (3*(rij(mu)/rij_norm)**2 - 1.0d0)
             k_esp(mu,nu) = qi*qj*rij_normi3 * ( 3*(rij(mu)*rij_normi)**2 - 1.0d0)
        else
            !k_esp(mu,nu) = 3.0d0*(qi*qj)*rij(mu)*rij(nu)/(rij_norm**5)
             k_esp(mu,nu) = 3.0d0*qi*qj*rij(mu)*rij(nu)*rij_normi5
             k_esp(nu,mu) = k_esp(mu,nu)
        end if
    end do
    end do
end subroutine esp_force_consts

subroutine get_bqhess_ii(i, scal, k_esp)
    integer, intent(in) :: i
    real*8, intent(in) :: scal, k_esp(3,3)
    integer :: mu, nu, row, col
    do nu=1,3
    do mu=1,3
        row = idx(i,mu)
        col = idx(i,nu)
        hess(row, col, 0, 0, 0) = hess(row, col, 0, 0, 0) + scal*k_esp(mu,nu)
    end do
    end do
end subroutine get_bqhess_ii

subroutine get_bqhess_ij(i, j, n, scal, k_esp)
    integer, intent(in) :: i, j, n(3)
    real*8, intent(in) :: scal, k_esp(3,3)
    integer :: mu, nu, row, col
    if ((abs(n(1)) .le. n_hess(1)) .and. (abs(n(2)) .le. n_hess(2)) &
        .and. (abs(n(3)) .le. n_hess(3))) then
        do nu=1,3
        do mu=1,3
            row = idx(i,mu)
            col = idx(j,nu)
            hess(row, col, n(1), n(2), n(3)) = hess(row, col, n(1), n(2), n(3)) - scal*k_esp(mu,nu)
            if (cellzero(n)) hess(col, row, 0, 0, 0) =  hess(col, row, 0, 0, 0) - scal*k_esp(mu,nu)
        end do
        end do
    end if
end subroutine get_bqhess_ij

subroutine get_bqhess_jj(j, n, scal, k_esp)
    integer, intent(in) :: j, n(3)
    real*8, intent(in) :: scal, k_esp(3,3)
    integer :: mu, nu, row, col
    if (cellzero(n)) then
        do nu=1,3
        do mu=1,3
            row = idx(j,mu)
            col = idx(j,nu)
            hess(row, col, 0, 0, 0) = hess(row, col, 0, 0, 0) + scal*k_esp(mu,nu)
        end do
        end do
    end if
end subroutine get_bqhess_jj

! (atom_i, x) --> index (out of 3N)
function idx(atom, coordinate)
    integer :: idx
    integer, intent(in) :: atom, coordinate
    idx = 3*(atom-1) + coordinate
end function idx


! return true if given (0,0,0)
function cellzero(n)
    logical :: cellzero
    integer, intent(in) :: n(3)
    if (n(1).eq.0 .and. n(2).eq.0 .and. n(3).eq.0) then
        cellzero = .true.
    else
        cellzero = .false.
    end if
end function cellzero


! write the formatted Cartesian force constants and d_dipole
subroutine report_hess
    integer :: na, nb, nc, i, j
    if (sys_master) then 
        open (45, file=fconfig(1:index(fconfig,'.')-1)//'.hess', access='append', status='unknown')
        write (45, '(A)') trim(fconfig)
        write (45, '(A,f0.1,A,f0.2,A,f0.2,A,f0.2,A)') 'BIM Hessian of '//trim(fconfig)//' [',&
           pext0, ' bar / '// trim(theory) // '(' // trim(embedding_field) //')/' // trim(basis_type) // ' / ', &
           sqrt(rcut_qq2), '/', sqrt(rcut_bq2), '/', sqrt(rcut_lr2),']'
        if (cellzero(n_hess)) then
            write(45, '(A)') "lower triangle hessian"
             do i=1,3*natom
             do j=1,i
                write(45,'(f21.10)') hess(i,j,0,0,0)
            end do
            end do
        else
            do na = -n_hess(1), n_hess(1)
            do nb = -n_hess(2), n_hess(2)
            do nc = -n_hess(3), n_hess(3)
                write(45,'(A,3i4)') "cell", na, nb, nc
                do i = 1,3*natom
                    write(45,*) hess(i,1:3*natom,na,nb,nc)
                end do 
            end do
            end do
            end do
        end if
        close(45)
        write(*, '(A,A)') 'Hessian written to ', fconfig(1:index(fconfig,'.')-1)//'.hess'
        
        open (45, file=fconfig(1:index(fconfig,'.')-1)//'.ddipole', access='append', status='unknown')
        write (45, '(A,f0.1,A,f0.2,A,f0.2,A,f0.2,A)') 'BIM Dipole Derivative of '//trim(fconfig)//' [',&
           pext0, ' bar / '// trim(theory) // '/' // trim(basis_type) // ' / ', &
           sqrt(rcut_qq2), '/', sqrt(rcut_bq2), '/', sqrt(rcut_lr2),']'
        do i = 1,3*natom
            write(45,*) d_dipole(i,1:3)
        end do
        close(45)
        write(*, '(A,A)') 'Dipole derivative written to ', fconfig(1:index(fconfig,'.')-1)//'.ddipole'

        print *, "max difference of hess(0) and transpose(hess(0)):", &
            maxval( abs( hess(:,:,0,0,0) - transpose(hess(:,:,0,0,0)) ) )
        if (n_hess(1) .gt. 0) then
            print *, "max difference of hess(1,0,0) and transpose(hess(-1,0,0)):", &
                maxval( abs( hess(:,:,1,0,0) - transpose(hess(:,:,-1,0,0)) ) )
        end if
    end if
end subroutine report_hess

end module hessian
