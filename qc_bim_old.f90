module qc_bim
! Bulk (3D)
use qc_system
use consts
use qc_mpi
use qc_neigh
use chem_interface
implicit none

contains

subroutine qc_bim_pot (en_mon, en_dim, en_lr)
    real*8, intent(out) :: en_mon, en_dim, en_lr
    logical :: list_fail
    integer :: cr, cm, c1, c2, elapsed!, im
    real*8 :: rate

    call system_clock(count_rate=cr)
    call system_clock(count_max=cm)
    rate = REAL(cr)

    if (sys_master) write(*,'(A)',advance="no") '   BIM_POT: lists...' 
    call system_clock(c1)

    call qc_get_lists(list_fail)
    if (list_fail) then
        call qc_allocate_lists
        call qc_get_lists(list_fail)
        if (list_fail) STOP "UNEXPECTED LIST OVERFLOW"
    end if

    if (sys_master) write(*,'(A)',advance="no") 'EMBEDDING...' 
    call qc_bim_field_update
    
    if (sys_master)  write(*,'(A)',advance="no") 'LR...'
    call qc_bim_lr (en_lr, .false.)
    
    if (sys_master) write(*,'(A)',advance="no") 'mono...' 
    call qc_bim_monomers (en_mon)
    
    if (sys_master) write(*,'(i0,A)',advance="no") num_pairs, ' dimers...'
    call qc_bim_dimers (en_dim, .false.)
    
    call system_clock(c2)
    elapsed = int(anint((c2-c1)/rate))
    if (sys_master) write(*,'(A,i0,A,i0,A)') "wall: ",elapsed/60,"m ",&
                                        mod(elapsed,60),"s"
    !if (sys_master) then
    !    do im = 1, nmol
    !        write(*,'(A,i0,A,i0,A)') "  molecule ", im, ": ", nlr(im), " LR's"
    !    end do
    !    write(*,'(A,f12.8)') "Total LR energy: ", en_lr
    !end if
    
end subroutine qc_bim_pot

subroutine qc_bim_pot_oldlists (en_mon, en_dim, en_lr)
    real*8, intent(out) :: en_mon, en_dim, en_lr
    integer :: cr, cm, c1, c2, elapsed!, im
    real*8 :: rate

    call system_clock(count_rate=cr)
    call system_clock(count_max=cm)
    rate = REAL(cr)

    if (sys_master) write(*,'(A)',advance="no") '   BIM_POT: reusing previous lists...' 
    call system_clock(c1)

    if (sys_master) write(*,'(A)',advance="no") 'EMBEDDING...' 
    call qc_bim_field_update
    
    if (sys_master)  write(*,'(A)',advance="no") 'LR...'
    call qc_bim_lr (en_lr, .true.) ! old lists
    
    if (sys_master) write(*,'(A)',advance="no") 'mono...' 
    call qc_bim_monomers (en_mon)
    
    if (sys_master) write(*,'(i0,A)',advance="no") num_pairs, ' dimers...'
    call qc_bim_dimers (en_dim, .true.) ! old lists
    
    call system_clock(c2)
    elapsed = int(anint((c2-c1)/rate))
    if (sys_master) write(*,'(A,i0,A,i0,A)') "wall: ",elapsed/60,"m ",&
                                        mod(elapsed,60),"s"
    !if (sys_master) then
    !    do im = 1, nmol
    !        write(*,'(A,i0,A,i0,A)') "  molecule ", im, ": ", nlr(im), " LR's"
    !    end do
    !    write(*,'(A,f12.8)') "Total LR energy: ", en_lr
    !end if
end subroutine qc_bim_pot_oldlists

subroutine qc_bim_field_update
    select case(embed_id)

    case (embed_id_esp)
        call qc_bim_esp_update
    case (embed_id_dipole)
        call qc_bim_dipole_update
    end select
end subroutine qc_bim_field_update

subroutine qc_bim_esp_update
    integer :: icycl
    integer :: im, is, ia, nsite
    real*8  :: chg_im(5), diff

    mol_ichg = 0
    mol_nchg = 0
    chg_coef = 0.0d0
    chg_pos = 0.0d0
    chg_im = 0.0d0

    do im = 1, nmol
        mol_ichg(im) = mol_iatom(im)
        mol_nchg(im) = mol_nsite(im)
        do is = 1, mol_nsite(im)
            chg_coef(is, mol_ichg(im)+is) = 1.0d0
        end do
    end do

    do icycl = 1, 10
        chg_pos_new = 0.0d0
        do im = 1 + sys_me, nmol, sys_nproc
            ia    = mol_iatom(im)
            nsite = mol_nsite(im)
            call system("mkdir -p "//trim(my_scratch))
            call esp_wrt (im)
            call system (trim(jcmd))
            call esp_read(nsite, chg_im)
            chg_pos_new(0, ia+1:ia+nsite) = chg_im(1:nsite)
            chg_pos_new(1:3, ia+1:ia+nsite) = gm_pos(1:3, ia+1:ia+nsite)
        end do

        call qc_mpi_barrier
        call qc_mpi_allreduce (chg_pos_new, 6*natom*4)

        diff = 0.0d0
        do ia = 1, natom
            diff = diff + (chg_pos_new(0, ia) - chg_pos(0, ia))**2.0d0
            chg_pos(0:3, ia) = chg_pos_new(0:3, ia)
        end do

        diff = sqrt(diff/dble(natom))

        if (diff .lt. 0.001d0) exit
    end do
    if (diff .gt. 0.001d0 .and. sys_master) then 
        print *, "ERROR: ESP DID NOT CONVERGE!"
        call qc_mpi_end
        STOP
    end if

end subroutine qc_bim_esp_update

subroutine qc_bim_dipole_update
    integer :: icycl
    integer :: im, ia, is, nsite, ichg
    real*8  :: dipole0(3), mag_dipole0, diff, nuclear_chg, dpos(3)
    real*8, parameter :: len_dipole = 0.05 ! Bohr
    character(len=2) :: atnm
    
    mol_ichg = 0
    mol_nchg = 0
    chg_coef = 0.0d0
    chg_pos = 0.0d0
    dipole0 = 0.0d0
    
    ichg = 0

    do im= 1, nmol
        
        nuclear_chg = 0.0d0
        mol_ichg(im) = ichg
        if (mol_netcharge(im) .eq. 0) then
            mol_nchg(im) = 2
        else
            mol_nchg(im) = 3
        end if

        ia = mol_iatom(im)
        
        do is = 1, mol_nsite(im)
            atnm = at_atnm(ia+is)

            if (trim(atnm) .eq. 'O') then
                chg_coef( is, ichg+1 : ichg+2 ) = 8.0d0
                if (mol_netcharge(im) .ne. 0) then
                    chg_coef ( is, ichg+3) = 8.0d0
                end if
                nuclear_chg = nuclear_chg + 8.0d0
            
            else if (trim(atnm) .eq. 'H') then
                chg_coef( is, ichg+1 : ichg+2 ) = 1.0d0
                if (mol_netcharge(im) .ne. 0) then
                    chg_coef ( is, ichg+3) = 1.0d0
                end if
                nuclear_chg = nuclear_chg + 1.0d0
            
            else if (trim(atnm) .eq. 'Cl') then
                chg_coef( is, ichg+1 : ichg+2 ) = 17.0d0
                if (mol_netcharge(im) .ne. 0) then
                    chg_coef ( is, ichg+3) = 17.0d0
                end if
                nuclear_chg = nuclear_chg + 17.0d0
            
            else if (trim(atnm) .eq. 'F') then
                chg_coef( is, ichg+1 : ichg+2 ) = 9.0d0
                if (mol_netcharge(im) .ne. 0) then
                    chg_coef(is, ichg+3) = 9.0d0
                end if
                nuclear_chg = nuclear_chg + 9.0d0
            else
                if (sys_master) print *, "need to implement nuclear charge for atom ", atnm
                call qc_mpi_end
                stop
            end if
        end do
        chg_coef(1:mol_nsite(im), ichg+1:ichg+mol_nchg(im)) = chg_coef(1:mol_nsite(im), ichg+1:ichg+mol_nchg(im)) / nuclear_chg
        ichg = ichg + mol_nchg(im)
    end do
    
    do icycl = 1, 10
        chg_pos_new = 0.0d0
        do im = 1 + sys_me, nmol, sys_nproc
            ia    = mol_iatom(im)
            nsite = mol_nsite(im)
            ichg = mol_ichg(im)
            call system("mkdir -p "//trim(my_scratch))
            call dipole_wrt (im)
            call system (trim(jcmd))
            call dipole_read(dipole0) ! Dipole moment vector (au)
            mag_dipole0 = sqrt(dot_product(dipole0, dipole0)) ! au
            
            ! Origin: center of charge (Angstrom)
            dpos = 0.0d0
            do is = 1, mol_nsite(im)
                dpos(:) = dpos(:) + chg_coef( is, ichg+1 ) * gm_pos(:, ia+is)
            end do
            
            ! Place two charges, distance len_dipole apart, along the dipole moment vector
            ! Positions in Angstroms
            ! Charges chosen to reproduce dipole moment
            chg_pos_new(1:3, ichg+1) = dpos(:) + (0.5*len_dipole*bohr2ang/mag_dipole0)*dipole0(:)
            chg_pos_new(1:3, ichg+2)   = dpos(:) - (0.5*len_dipole*bohr2ang/mag_dipole0)*dipole0(:)
            chg_pos_new(0, ichg+1) = mag_dipole0 / len_dipole
            chg_pos_new(0, ichg+2)   = -1.0d0 * mag_dipole0 / len_dipole
            ! If nonzero net charge, place monopole at center-of-charge
            if (mol_netcharge(im) .ne. 0) then
                chg_pos_new(1:3, ichg+3) = dpos(:)
                chg_pos_new(0, ichg+3) = mol_netcharge(im)
            end if
        end do

        call qc_mpi_barrier
        call qc_mpi_allreduce (chg_pos_new, 6*natom*4)
        
        ! RMSD of dipole moment (au)
        diff = 0.0d0
        do im = 1, nmol
            ichg = mol_ichg(im)
            diff = diff + (chg_pos_new(0, ichg+1)*len_dipole - chg_pos(0, ichg+1)*len_dipole)**2
            chg_pos(0:3, ichg+1) = chg_pos_new(0:3, ichg+1)
            chg_pos(0:3, ichg+2) = chg_pos_new(0:3, ichg+2)
            if (mol_netcharge(im) .ne. 0) then
                chg_pos(0:3, ichg+3) = chg_pos_new(0:3, ichg+3)
            end if
        end do

        !if (sys_master) then 
        !    print *, "DIPOLE EMBED ITER", icycl
        !    do im = 1, nmol
        !        print *, "molecule", im
        !        do ichg = mol_ichg(im)+1, mol_ichg(im)+mol_nchg(im)
        !            write(*,'(4f12.6)') chg_pos(0:3, ichg)
        !        end do
        !    end do
        !end if
        diff = sqrt(diff/dble(nmol)) 

        if (diff .lt. 0.0001d0) exit
    end do
    
    if (diff .gt. 0.0001d0 .and. sys_master) then 
        print *, "ERROR: DIPOLES DID NOT CONVERGE!"
        call qc_mpi_end
        STOP
    end if
end subroutine qc_bim_dipole_update

subroutine qc_bim_lr (en_lr, l_old_neigh)
    integer :: im
    real*8 :: u_lr0
    real*8, intent(out) :: en_lr
    logical, intent(in) :: l_old_neigh

    en_lr = 0.0d0
    d_mr = 0.0d0
    m_virt = 0.0d0
    nlr = 0
    
    do im = 1 + sys_me, nmol, sys_nproc
        call get_lr_interaction(u_lr0, im, l_old_neigh) ! modifies d_mr and m_virt
        en_lr = en_lr + u_lr0
    end do

    call qc_mpi_barrier
    call qc_mpi_allreduce1 (en_lr)
    call qc_mpi_allreducei (nlr, nmol)
end subroutine qc_bim_lr

subroutine qc_bim_monomers (u_mon)
    real*8, intent(out) :: u_mon
    integer :: i, j
    integer :: im, ia, is, nsite
    integer :: kb, km, kk, ka, ks, nsite_k, na, nb, nc
    integer :: nchg_i, nchg_k, ichg_i, ichg_k, at_i, at_k
    real*8  :: d_rr0(3,5), upot
    real*8  :: d_ri(3), ri(3), rk(3), dcel(3), rik(3), rik_norm, qi, qk

    u_mon   = 0.0d0
    ! Do not zero out d_mr or m_virt; it contains LR interaction at this point


    do im = 1 + sys_me, nmol, sys_nproc
        nsite = mol_nsite(im)
        ia    = mol_iatom(im)
        
        call system("mkdir -p "//trim(my_scratch))
        call monomer_wrt(im)
        call system (trim(jcmd))
        
        if (l_grd) then
            call pot_grad_read (nsite, upot, d_rr0)
            if (l_bq) call bq_grad_read(d_bq0, im) ! Embedded monomers
        else
            call pot_read(upot)
        end if
    
        u_mon = u_mon + upot

        if (l_grd) then
            do is = 1, nsite
                d_ri = d_rr0(1:3,is)
                d_mr(1:3,ia+is) = d_mr(1:3,ia+is) + d_ri(1:3)
                ri(1:3) = gm_pos(1:3,ia+1)
                do j = 1, 3
                do i = 1, 3
                    m_virt(i,j) = m_virt(i,j) - ri(i)*ang2bohr*d_ri(j) ! force is -d_ri
                end do
                end do
            end do
            if (l_bq) then
                kk = 0
                do kb = 1, nbq(im)
                    km = bq_list(kb,im)%jm
                    na = bq_list(kb,im)%ja
                    nb = bq_list(kb,im)%jb
                    nc = bq_list(kb,im)%jc
                    dcel(:) = na*lat(:,1) + nb*lat(:,2) + nc*lat(:,3)

                    ka = mol_iatom(km)
                    nsite_k = mol_nsite(km)
                    do ks = 1, nsite_k
                        kk = kk + 1
                        d_ri = d_bq0(1:3,kk)
                        d_mr(1:3,ka+ks) = d_mr(1:3,ka+ks) + d_ri(1:3)
                        ri = (gm_pos(1:3,ka+1) + dcel(1:3))
                        do j = 1, 3
                        do i = 1, 3
                            m_virt(i,j) = m_virt(i,j) - ri(i)*ang2bohr*d_ri(j) ! force is -d_ri
                        end do
                        end do
                    end do
                end do
            end if
        end if
        

        ! Correct overcounted QM -- BQ interactions
        nchg_i = mol_nchg(im)
        ichg_i = mol_ichg(im)
        do kb = 1, nbq(im)
            if (bq_list(kb,im)%qm) cycle ! already canceled in QM region
            
            km = bq_list(kb,im)%jm
            na = bq_list(kb,im)%ja
            nb = bq_list(kb,im)%jb
            nc = bq_list(kb,im)%jc
            dcel(:) = na*lat(:,1) + nb*lat(:,2) + nc*lat(:,3)

            ka = mol_iatom(km)
            nsite_k = mol_nsite(km)
            nchg_k = mol_nchg(km)
            ichg_k = mol_ichg(km)
            
            do is = 1, nchg_i
                ri(1:3) = chg_pos(1:3,ichg_i+is)
                qi =      chg_pos(0,  ichg_i+is)
                
                do ks = 1, nchg_k
                    rk(1:3) = chg_pos(1:3,ichg_k+ks) + dcel(1:3)
                    qk = chg_pos(0, ichg_k+ks)
                    rik = (rk - ri)*ang2bohr
                    rik_norm = sqrt(dot_product(rik,rik))

                    ! correct qi--qk energy
                    u_mon = u_mon - 0.5*qi*qk/rik_norm 
                    
                    if (l_grd) then
                        ! correct force on qi
                        d_ri = 0.5d0 * rik * (qi*qk/(rik_norm**3)) ! points from i to k
                        do at_i = 1, nsite
                            d_mr(1:3,ia+at_i) = d_mr(1:3,ia+at_i) - d_ri(1:3)*chg_coef(at_i, ichg_i+is) ! grad of i is +d_ri
                        end do
                        do j = 1,3
                        do i = 1,3
                            m_virt(i,j) = m_virt(i,j) + gm_pos(i,ia+1)*ang2bohr*d_ri(j)
                        end do
                        end do

                        ! correct force on qk, if using BQ_GRAD
                        if (l_bq) then
                            do at_k = 1, nsite_k
                                d_mr(1:3,ka+at_k) = d_mr(1:3,ka+at_k) + d_ri(1:3)*chg_coef(at_k, ichg_k+ks) ! grad of k is -d_ri
                            end do
                            do j = 1,3
                            do i = 1,3
                                m_virt(i,j) = m_virt(i,j) - (gm_pos(i,ka+1)+dcel(i))*ang2bohr*d_ri(j)
                            end do
                            end do
                        end if
                    end if
                end do
            end do
        end do
    end do
            
    call qc_mpi_barrier
    call qc_mpi_allreduce1(u_mon)
    if (l_grd) call qc_mpi_allreduce (d_mr,   3*natom)
    if (l_grd) call qc_mpi_allreduce (m_virt, 9)
    !if (sys_master) call print_mono_gradients

end subroutine qc_bim_monomers

subroutine qc_bim_dimers (u_dim, l_old)
    use qc_lattice
    real*8, intent(out) :: u_dim
    logical, intent(in) :: l_old
    integer :: idx, ierr, tmp, done
    integer, dimension(MPI_STATUS_SIZE) :: stat

    u_dim = 0.0d0
    d_rr  = 0.0d0
    virt  = 0.0d0
    
    tmp = 0
    done = -1
    ! PARALLEL
    if (sys_nproc .gt. 1) then
        ! MASTER
        if (sys_master) then
            do idx = sys_nproc, num_pairs
                call mpi_recv(tmp, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                              0, MPI_COMM_WORLD, stat, ierr)
                call mpi_send(idx, 1, MPI_INTEGER, stat(MPI_SOURCE), &
                              0, MPI_COMM_WORLD, ierr)
            end do
            
            do idx = 1, sys_nproc-1
                call mpi_recv(tmp, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                              0, MPI_COMM_WORLD, stat, ierr)
                call mpi_send(done, 1, MPI_INTEGER, stat(MPI_SOURCE), & 
                              0, MPI_COMM_WORLD, ierr)
            end do
        ! WORKER
        else
            idx = sys_me
            if (idx .le. num_pairs) call calc_dimer(idx, u_dim, l_old)
            call mpi_send (tmp, 1, MPI_INTEGER, 0, 0, &
                           MPI_COMM_WORLD, ierr)
            call mpi_recv (idx, 1, MPI_INTEGER, 0, 0, &
                           MPI_COMM_WORLD, stat, ierr)
            do while (idx .ne. done)
                call calc_dimer(idx, u_dim, l_old)
                call mpi_send (tmp, 1, MPI_INTEGER, 0, 0, &
                               MPI_COMM_WORLD, ierr)
                call mpi_recv (idx, 1, MPI_INTEGER, 0, 0, &
                               MPI_COMM_WORLD, stat, ierr)
            end do
        end if
    
    ! SERIAL
    else
        do idx = 1, num_pairs
            call calc_dimer(idx, u_dim, l_old)
        end do
    end if

    call qc_mpi_barrier
    call qc_mpi_allreduce1 (u_dim)
    call qc_mpi_allreduce  (d_rr, 3*natom)
    call qc_mpi_allreduce  (virt, 9)

    ! --- Total Gradient: Monomer + dimer contributions
    d_rr = d_rr + d_mr
    ! --- Total Virial: Monomer + dimer contributions
    virt = virt + m_virt

end subroutine qc_bim_dimers

subroutine calc_dimer(idx, u_dim, l_old)
    integer, intent(in) :: idx
    real*8, intent(inout) :: u_dim
    logical, intent(in) :: l_old
    integer :: ip, im, ia, jm, ja, is, js, kb, km, ka, ks, kk
    integer :: na, nb, nc, nba, nbb, nbc
    integer :: i, j
    integer :: nsite_i, nsite_ig, nsite_j, nsite_jg, nsite_k
    real*8  :: en_i, en_j, en_ij
    real*8  :: d_rr0(3,12), d_mri(3,6), d_mrj(3,6)
    real*8  :: d_ri(3), ri(3), dcel(3)
    real*8  :: scale

    im = 1
    ip = idx
    do while (ip .gt. npair(im))
        ip = ip - npair(im)
        im = im + 1
    end do

    en_i = 0.0d0
    en_j = 0.0d0
    
    d_mri = 0.0d0
    d_mrj = 0.0d0
    
    d_bq0 = 0.0d0
    
    d_mbqi = 0.0d0
    d_mbqj = 0.0d0

    !--dimer index ip of monomer im
    jm = pair_list(ip,im)%jm
    na = pair_list(ip,im)%ja
    nb = pair_list(ip,im)%jb
    nc = pair_list(ip,im)%jc
    scale = pair_list(ip,im)%scale

    dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

    nsite_i = mol_nsite(im)
    nsite_j = mol_nsite(jm)
    if (l_bsse) then
        nsite_ig = nsite_j
        nsite_jg = nsite_i
    else
        nsite_ig = 0
        nsite_jg = 0
    end if

    !--get union of bq lists
    call qc_bq2_list_get(im, ip, l_old)
    
    !-- dimer: im(0,0,0) -- jm(na, nb, nc)
    !-- store en_ij, d_rr0, d_bq0
    call system("mkdir -p "//trim(my_scratch))
    call dimer_wrt (im, ip)
    call system (trim(jcmd))
    if (l_grd) then
        call pot_grad_read (nsite_i+nsite_j, en_ij, d_rr0)
        if (l_bq) call bq2_grad_read(d_bq0)
    else
        call pot_read(en_ij)
    end if
    
    !print *, ""
    !print *, "IJ (QM)"
    !write(*,'(3f10.4)') d_rr0(1:3, 1)
    !write(*,'(3f10.4)') d_rr0(1:3, 2)
    !write(*,'(3f10.4)') d_rr0(1:3, 3)
    !print *, "IJ (BQ)"
    !write(*,'(3f10.4)') d_bq0(1:3, 1)
    !write(*,'(3f10.4)') d_bq0(1:3, 2)
    !write(*,'(3f10.4)') d_bq0(1:3, 3)

    ! -- monomer im(0,0,0) in bq2 field.
    call system("mkdir -p "//trim(my_scratch))
    call monomerbq2_wrt(im, ip, 0, l_bsse)
    call system (trim(jcmd))
    if (l_grd) then
        call pot_grad_read (nsite_i + nsite_ig, en_i, d_mri)
        if (l_bq) call bq2_grad_read(d_mbqi)
    else
        call pot_read(en_i)
    end if
    
    !print *, "I (QM)"
    !write(*,'(3f10.4)') d_mri(1:3, 1)
    !write(*,'(3f10.4)') d_mri(1:3, 2)
    !write(*,'(3f10.4)') d_mri(1:3, 3)
    !print *, "I (BQ)"
    !write(*,'(3f10.4)') d_mbqi(1:3, 1)
    !write(*,'(3f10.4)') d_mbqi(1:3, 2)
    !write(*,'(3f10.4)') d_mbqi(1:3, 3)

    ! -- monomer jm(na,nb,nc) in bq2 field.
    call system("mkdir -p "//trim(my_scratch))
    call monomerbq2_wrt(im, ip, 1, l_bsse)
    call system (trim(jcmd))
    if (l_grd) then
        call pot_grad_read (nsite_j + nsite_jg, en_j, d_mrj)
        if (l_bq) call bq2_grad_read(d_mbqj)
    else
        call pot_read(en_j)
    end if
    
    !print *, "J (QM)"
    !write(*,'(3f10.4)') d_mrj(1:3, 1)
    !write(*,'(3f10.4)') d_mrj(1:3, 2)
    !write(*,'(3f10.4)') d_mrj(1:3, 3)
    !print *, "J (BQ)"
    !write(*,'(3f10.4)') d_mbqj(1:3, 1)
    !write(*,'(3f10.4)') d_mbqj(1:3, 2)
    !write(*,'(3f10.4)') d_mbqj(1:3, 3)

    ! QM ENERGY
    u_dim = u_dim + scale*(en_ij - en_i - en_j)

    ia = mol_iatom(im)
    ja = mol_iatom(jm)

    ! QM GRADIENT, VIRIAL
    if (l_grd) then

        ! forces on i
        do is = 1, nsite_i
            d_ri = scale*(d_rr0(1:3,is) - d_mri(1:3,is))
            d_rr(1:3,ia+is) = d_rr(1:3,ia+is) + d_ri(1:3)
            ri(1:3) = gm_pos(1:3,ia+1)
            do j = 1, 3
            do i = 1, 3
                virt(i,j) = virt(i,j) - ri(i)*ang2bohr*d_ri(j) ! -d_ri is force
            end do
            end do
        end do
        
        ! forces on ghost j (BSSE)
        do is = 1, nsite_ig
            d_ri = -1.0d0*scale*d_mri(1:3,nsite_i+is)
            d_rr(1:3,ja+is) = d_rr(1:3,ja+is) + d_ri(1:3)
            ri(1:3) = gm_pos(1:3,ja+1)+dcel(1:3)
            do j = 1, 3
            do i = 1, 3
                virt(i,j) = virt(i,j) - ri(i)*ang2bohr*d_ri(j) ! -d_ri is force
            end do
            end do
        end do

        ! gradient on j
        do js = 1, nsite_j
            d_ri = scale*(d_rr0(1:3,nsite_i+js) - d_mrj(1:3,js))
            d_rr(1:3,ja+js) = d_rr(1:3,ja+js) + d_ri(1:3)
            ri(1:3) = (gm_pos(1:3,ja+1)+dcel(1:3))
            do j = 1, 3
            do i = 1, 3
                virt(i,j) = virt(i,j) - ri(i)*ang2bohr*d_ri(j)
            end do
            end do
        end do
        
        ! gradient on ghost i (BSSE)
        do js = 1, nsite_jg
            d_ri = -1.0d0*scale*d_mrj(1:3,nsite_j+js)
            d_rr(1:3,ia+js) = d_rr(1:3,ia+js) + d_ri(1:3)
            ri(1:3) = gm_pos(1:3,ia+1)
            do j = 1, 3
            do i = 1, 3
                virt(i,j) = virt(i,j) - ri(i)*ang2bohr*d_ri(j) ! -d_ri is force
            end do
            end do
        end do

        ! BQ GRAD, VIRIAL
        if (l_bq) then
            kk = 0
            do kb = 1, nbq2
                km = bq2_list(kb)%jm
                nba = bq2_list(kb)%ja
                nbb = bq2_list(kb)%jb
                nbc = bq2_list(kb)%jc
                dcel(:) = nba*lat(:,1)+nbb*lat(:,2)+nbc*lat(:,3)

                ka      = mol_iatom(km)
                nsite_k = mol_nsite(km)

                do ks = 1, nsite_k
                    kk = kk + 1

                    if (km.eq.im .and. nba.eq.0 .and. nbb.eq.0 .and. nbc.eq.0) then
                        d_bq0(1:3,kk) = 0.0d0
                        d_mbqi(1:3,kk) = 0.0d0
                    end if
                    if (km.eq.jm .and. nba.eq.na .and. nbb.eq.nb .and. nbc.eq.nc) then
                        d_bq0(1:3,kk) = 0.0d0
                        d_mbqj(1:3,kk) = 0.0d0
                    end if

                    d_ri(1:3) = scale*(d_bq0(1:3,kk) - d_mbqi(1:3,kk) - d_mbqj(1:3,kk)) !neg
                    d_rr(1:3,ka+ks) = d_rr(1:3,ka+ks) + d_ri
                    ri(1:3) = (gm_pos(1:3,ka+1) + dcel(1:3))
                    do j = 1, 3
                    do i = 1, 3
                        virt(i,j) = virt(i,j) - ri(i)*ang2bohr*d_ri(j) ! neg
                    end do
                    end do
                end do
            end do
        end if
    end if
    !print *, "DIMER GRADIENTS AFTER PAIR", idx 
    !call print_gradients
end subroutine calc_dimer

subroutine print_gradients
    integer :: im, ia, is
    write (6, *) '                   BIM ENERGY GRADIENTS' 
    write (6, *) ' mol(at)           coordinates                    gradient            '
    write (6, *) '----------------------------------------------------------------------'
    do im = 1,nmol
        ia = mol_iatom(im)
        do is=1,mol_nsite(im)
        if (ia+is.lt.10) then 
            write (6, '(1x,i0,A,i0,A,A,3f10.3,3f11.4)') im,'(',ia+is,') ', at_atnm(ia+is),&
            gm_pos(1:3, ia+is), d_rr(1:3,ia+is)
        else
            write (6, '(i0,A,i0,A,A,3f10.3,3f11.4)') im,'(',ia+is,') ', at_atnm(ia+is),&
            gm_pos(1:3, ia+is), d_rr(1:3,ia+is)
        end if
        end do
    end do
end subroutine print_gradients

subroutine print_mono_gradients
    integer :: im, ia, is
    write (6, *) '                   BIM 1-BODY GRADIENTS' 
    write (6, *) ' mol(at)           coordinates                    gradient            '
    write (6, *) '----------------------------------------------------------------------'
    do im = 1,nmol
        ia = mol_iatom(im)
        do is=1,mol_nsite(im)
        if (ia+is.lt.10) then 
            write (6, '(1x,i0,A,i0,A,A,3f10.3,3f11.4)') im,'(',ia+is,') ', at_atnm(ia+is),&
            gm_pos(1:3, ia+is), d_mr(1:3,ia+is)
        else
            write (6, '(i0,A,i0,A,A,3f10.3,3f11.4)') im,'(',ia+is,') ', at_atnm(ia+is),&
            gm_pos(1:3, ia+is), d_mr(1:3,ia+is)
        end if
        end do
    end do
end subroutine print_mono_gradients

end module qc_bim
