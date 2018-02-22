subroutine psi_esp_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im

    open  (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    write (15, '(A)') 'CALC_TYPE ESP'
    call psi_print_geom_one (im, 0, 0, 0)
    call psi_print_esp_bq (im)
    write (15, '(A)') 'theory hf'
    write (15, '(A)') 'basis '//trim(basis_type)

    write (15, *)
    close (15)
end subroutine psi_esp_wrt

subroutine psi_dipole_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im

    stop "not implemented"
end subroutine psi_dipole_wrt

subroutine psi_atomdipole_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im

    stop "not implemented"
end subroutine psi_atomdipole_wrt

subroutine psi_esp_read (nsite,chg0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: chg0 (nsite)
    integer :: iend, is

    chg0 = 0.0d0

    open (unit=15,file=trim(my_scratch)//'/'//'esp.dat',status='old')
    do is = 1, nsite
        read (15, *,iostat=iend) chg0(is)
        if (iend /= 0) stop "did not read esp"
    end do
    close (15)

end subroutine psi_esp_read

subroutine psi_dipole_read (dipole0)
    use qc_system
    implicit none
    real*8, intent(out) :: dipole0(3)
    integer :: eof
    character(len=120) :: buffer
    logical :: found
    stop "not implemented"
end subroutine psi_dipole_read

subroutine psi_atomdipole_read (nsite, charges, dipoles)
    use qc_system
    use qc_mpi
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: charges(nsite)
    real*8, intent(out) :: dipoles(3, nsite)
    integer :: eof, is
    integer :: idx
    character(len=120) :: buffer
    character(len=3) :: symbol
    logical :: found

    charges = 0.0d0
    dipoles = 0.0d0
    found = .false.
    
    stop "not implemented"
end subroutine psi_atomdipole_read

subroutine psi_monomer_wrt (id)
    use qc_system
    implicit none
    integer, intent(in) :: id
    !character(len=500) :: maketempfile

    open (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    write (15, '(A)') 'theory '//trim(theory)
    write (15, '(A)') 'basis '//trim(basis_type)

    if (l_grd) then
        if (.not. use_respa) then
            write (15, '(A)') 'CALC_TYPE GRADIENT'
        else
            write (15, '(A)') 'CALC_TYPE REFGRADIENT'
        end if
    else if (l_frq) then
        write (15, '(A)') 'CALC_TYPE HESSIAN'
    else
        write (15, '(A)') 'CALC_TYPE ENERGY'
    end if

    write (15, *)
    call psi_print_geom_one (id, 0, 0, 0)
    call psi_print_bq (id) ! Monomers embedded

    write (15, *)
    close (15)
    !write(maketempfile, '(a,i0,a)') "cp "//trim(my_scratch)//'/'//trim(jname)//".in mono", id, ".in"
    !call system (maketempfile)
end subroutine psi_monomer_wrt

subroutine psi_dimer_wrt (im, ip)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im, ip
    integer :: nsite_i, nsite_j, jm, n(3)
    !character(len=500) :: maketempfile

    jm = pair_list(ip,im)%jm
    n(1) = pair_list(ip,im)%ja
    n(2) = pair_list(ip,im)%jb
    n(3) = pair_list(ip,im)%jc
    nsite_i = mol_nsite(im)
    nsite_j = mol_nsite(pair_list(ip,im)%jm)

    open (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    write (15, '(A)') 'theory '//trim(theory)
    write (15, '(A)') 'basis '//trim(basis_type)

    if (l_grd) then
        if (.not. use_respa) then
            write (15, '(A)') 'CALC_TYPE GRADIENT'
        else
            write (15, '(A)') 'CALC_TYPE REFGRADIENT'
        end if
    else if (l_frq) then
        write (15, '(A)') 'CALC_TYPE HESSIAN'
    else
        write (15, '(A)') 'CALC_TYPE ENERGY'
    end if

    write (15, *) 
    call psi_print_geom_two (im, ip)
    call psi_print_bq2 (im, ip, 0)

    write (15, *)
    close (15)
    !write(maketempfile, '(a,i0,a,i0,a,i0,a)') "cp "//trim(my_scratch)//'/'//& 
    !trim(jname)//".in dimer", im, "-", jm, "_", n(1), ".in"
    !call system (maketempfile)

end subroutine psi_dimer_wrt

subroutine psi_monomerbq2_wrt (im, ip, mono, ghost)
    ! mono=0: calculate i in union of ij fields
    ! mono!=0: calculate j in union of ij fields
    ! ghost=true: place basis functions at site of the other monomer
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im, ip, mono
    logical, intent(in) :: ghost
    integer :: id, na, nb, nc, jm, jn(3), void
    real*8 :: dcel(3)
    character(len=7) :: mono_str
    !character(len=500) :: maketempfile
    
    jm = pair_list(ip,im)%jm
    jn(1) = pair_list(ip,im)%ja
    jn(2) = pair_list(ip,im)%jb
    jn(3) = pair_list(ip,im)%jc

    !--embed in field of ij
    if (mono .eq. 0) then !--monomer i
        mono_str = '[qm i] '
        id = im
        na = 0
        nb = 0
        nc = 0
        dcel = 0.0d0
        void = 1 !--void charges of i
        !write(maketempfile, '(a,i0,a,i0,a,i0,a)') "cp "//trim(my_scratch)//'/'// &
        !trim(jname)//".in mono", im, "_bq", jm, "_", jn(1), ".in"
    else !--monomer j
        mono_str = '[qm j] '
        id = jm
        na = jn(1)
        nb = jn(2)
        nc = jn(3)
        void = 2 !--void charges of j
        !write(maketempfile, '(a,i0,a,i0,a,i0,a)') "cp "//trim(my_scratch)//'/'// &
        !trim(jname)//".in bq", im, "_mono", jm, "_", na, ".in"
    end if

    open (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    write (15, '(A)') 'theory '//trim(theory)
    write (15, '(A)') 'basis '//trim(basis_type)
    if (l_grd) then
        if (.not. use_respa) then
            write (15, '(A)') 'CALC_TYPE GRADIENT'
        else
            write (15, '(A)') 'CALC_TYPE REFGRADIENT'
        end if
    else if (l_frq) then
        write (15, '(A)') 'CALC_TYPE HESSIAN'
    else
        write (15, '(A)') 'CALC_TYPE ENERGY'
    end if
    
    write(15, *)

    if (ghost) then
        call psi_print_geom_one_ghost (im, ip, mono)
    else
        call psi_print_geom_one (id, na, nb, nc)
    end if
    call psi_print_bq2 (im, ip, void)
    write (15, *)
    close (15)
    !call system (maketempfile)
end subroutine psi_monomerbq2_wrt

subroutine psi_pot_read (upot) 
    use qc_system
    implicit none
    real*8,  intent(out) :: upot
    integer :: iend
    upot = 0.0d0
    open (unit=15,file=trim(my_scratch)//'/'//'en_grad.dat', status='old')
    read (15, *,iostat=iend) upot
    if (iend /= 0) stop "did not read energy"
    close (15)
end subroutine psi_pot_read

subroutine psi_pot_grad_read_respa (natm, upot, dRR0, dRR1)
    use qc_system
    implicit none
    integer, intent(in)  :: natm
    real*8,  intent(out) :: upot
    real*8,  intent(out) :: dRR0(3,natm), dRR1(3,natm)
    integer :: eof, iat, mu, iend
    character(len=120) :: buffer

    real*8 :: upot_corr

    dRR0 = 0.0d0
    dRR1 = 0.0d0
    upot = 0.0d0
    
    open (unit=15,file=trim(my_scratch)//'/'//'en_grad_ref.dat', status='old')
    open (unit=16,file=trim(my_scratch)//'/'//'en_grad_correl.dat', status='old')

    read (15, *,iostat=iend) upot
    if (iend /= 0) stop "did not read reference energy"
    read (16, *,iostat=iend) upot_corr
    if (iend /= 0) stop "did not read correlated energy"

    do iat = 1, natm
        read (15, *,iostat=iend) dRR0(1:3, iat)
        if (iend /= 0) stop "did not read (reference) gradient"
        read (16, *,iostat=iend) dRR1(1:3, iat)
        if (iend /= 0) stop "did not read (correlated) gradient"
        dRR1(1:3, iat) = dRR1(1:3, iat) - dRR0(1:3, iat)
    end do

    close (15)
    close (16)

end subroutine psi_pot_grad_read_respa

subroutine psi_pot_grad_read (natm, upot, dRR0)
    use qc_system
    implicit none
    integer, intent(in)  :: natm
    real*8,  intent(out) :: upot
    real*8,  intent(out) :: dRR0(3,natm)
    integer :: eof, iat, mu, iend
    character(len=120) :: buffer

    dRR0 = 0.0d0
    upot = 0.0d0
    
    open (unit=15,file=trim(my_scratch)//'/'//'en_grad.dat', status='old')
    read (15, *,iostat=iend) upot
    if (iend /= 0) stop "did not read energy"

    do iat = 1, natm
        read (15, *,iostat=iend) dRR0(1:3, iat)
        if (iend /= 0) stop "did not read gradient"
    end do

    close (15)

end subroutine psi_pot_grad_read

subroutine psi_bq_grad_read (d_bq, im)
    ! Read electric field at Bq sites 
    ! multiply by -1.0*charge to get gradient
    ! This gives same results as NW "Bq Force"
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    real*8, intent(out) :: d_bq(3, 2*maxnum_bq*3)
    real*8 :: grad(3)
    integer :: kb, km, kk, ks, is, nsite, iend, nchg, ichg
    
    open (unit=15,file=trim(my_scratch)//'/'//'grid_field.dat',status='old')
    d_bq = 0.0d0
    kk = 0

    do kb = 1, nbq(im) !--for every Bq molecule
        km = bq_list(kb,im)%jm
        nsite = mol_nsite(km)
        nchg = mol_nchg(km)
        ichg = mol_ichg(km)
        do ks = 1, nchg !--for every charge in Bq molecule
            read (15, *,iostat=iend) grad(1:3)
            if (iend /= 0) stop "did not read all bq gradients"
            grad(:) = -1.0d0*chg_pos(0, ichg+ks)*grad(:) ! convert field to grad_bq
            do is = 1, nsite
                d_bq(1:3,kk+is) = d_bq(1:3,kk+is) + grad(:)*chg_coef(is, ichg+ks) !--sum gradient_charge into gradient_nuclei
            end do
        end do
        kk = kk + nsite
    end do

    close (15)

end subroutine psi_bq_grad_read

subroutine psi_bq2_grad_read (d_bq)
    ! Read electric field at Bq sites 
    ! multiply by -1.0*charge to get gradient
    ! Extra logic: some charges are voided (!!!)
    ! The voiding is performed in qc_bim.f90 calc_dimer
    use qc_system
    use qc_neigh
    implicit none
    real*8, intent(out) :: d_bq(3,2*maxnum_bq*3)
    real*8 :: grad(3)
    integer :: is, kb, km, kk, ks, nsite, ichg, nchg, iend

    open (unit=15,file=trim(my_scratch)//'/'//'grid_field.dat',status='old')
    d_bq = 0.0d0
    kk = 0

    do kb = 1, nbq2 !--for every Bq molecule
        km = bq2_list(kb)%jm
        nchg = mol_nchg(km)
        ichg = mol_ichg(km)
        nsite = mol_nsite(km)
        do ks = 1, nchg !--for every charge in Bq molecule
            read (15, *,iostat=iend) grad(1:3)
            if (iend /= 0) stop "did not read all bq gradients"
            grad(:) = -1.0d0*chg_pos(0, ichg+ks)*grad(:) ! convert field to grad_bq
            do is = 1, nsite
                d_bq(1:3,kk+is) = d_bq(1:3,kk+is) + grad(:)*chg_coef(is, ichg+ks) !--sum gradient_charge into gradient_nuclei
            end do
        end do
        kk = kk + nsite
    end do

    close (15)

end subroutine psi_bq2_grad_read

subroutine psi_print_geom_one (im, na, nb, nc)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im, na, nb, nc
    integer :: is, ia, nsite, chg0
    real*8  :: dcel(3)

    ia    = mol_iatom(im)
    nsite = mol_nsite(im)
    dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

    chg0  = mol_netcharge(im)
    write (15, '(A, i0)') 'charge ', chg0
    if (chg0 .eq. 0) then
        write (15, '(A)') 'mult 1'
    else
        write (15, '(A)') 'mult 2'
    end if

    write (15, '(A)') 'GEOM'
    do is = 1, nsite
    write (15, '(A,3f20.10)') at_atnm(ia+is), gm_pos(1:3,ia+is)+dcel
    end do
    write (15, '(A)') 'END'
    write (15, *)
end subroutine psi_print_geom_one

subroutine psi_print_geom_one_ghost (im, ip, mono)
    ! mono=0: geom i, ghost charges j
    ! mono!=0: geom j, ghost charges i
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im, ip, mono
    integer :: mono_id, ghost_id, chg0
    real*8 :: mono_dcel(3), ghost_dcel(3)
    integer :: mono_nsite, mono_start
    integer :: ghost_nsite, ghost_start
    integer :: is

    stop "not implemented"
end subroutine psi_print_geom_one_ghost

subroutine psi_print_geom_two (im, ip)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im, ip
    integer :: jm, na, nb, nc
    integer :: is, ia, ja, nsite_i, nsite_j
    integer :: chg0
    real*8  :: jbox(3)

    jm = pair_list(ip,im)%jm
    na = pair_list(ip,im)%ja
    nb = pair_list(ip,im)%jb
    nc = pair_list(ip,im)%jc
    jbox(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

    ia      = mol_iatom(im)
    nsite_i = mol_nsite(im)

    ja      = mol_iatom(jm)
    nsite_j = mol_nsite(jm)

    chg0 = mol_netcharge(im) + mol_netcharge(jm)
    write (15, '(A, i0)') 'charge ', chg0
    if (chg0 .eq. 0) then
        write (15, '(A)') 'mult 1'
    else
        write (15, '(A)') 'mult 2'
    end if

    write (15, '(A)') 'GEOM'
    do is = 1, nsite_i
    write (15, '(A,3f20.10)') at_atnm(ia+is), gm_pos(1:3,ia+is)
    end do
    do is = 1, nsite_j
    write (15, '(3x,A,3f20.10)') at_atnm(ja+is), (gm_pos(1:3,ja+is) + jbox)
    end do

    write (15, '(A)') 'END'
    write (15, *)
end subroutine psi_print_geom_two

subroutine psi_print_esp_bq (im)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    integer :: ia, kb, km, na, nb, nc, ka, ks, nsite
    real*8  :: dcel(3)

    if ( maxval(abs(chg_pos)) .lt. 1.0d-10 ) then
        write (15, *)
        return
    end if

    write (15, '(A)') 'BQ_CHARGES'
    ia = mol_iatom(im)
    dcel = 0.0d0

    do kb = 1, nbq(im)
        km = bq_list(kb,im)%jm
        na = bq_list(kb,im)%ja
        nb = bq_list(kb,im)%jb
        nc = bq_list(kb,im)%jc
        dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

        ka    = mol_ichg(km)
        nsite = mol_nchg(km)

        ! write charge first, then position
        do ks = 1, nsite
            write (15,'(f16.10,3f20.10)') chg_pos(0,ka+ks), chg_pos(1:3,ka+ks)+dcel
        end do
    end do

    write (15, '(A)') 'END'
    write (15, *)
end subroutine psi_print_esp_bq

subroutine psi_print_bq (im)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    integer :: kb, km, na, nb, nc, ka, ks, nsite
    real*8  :: dcel(3) !, zero

    !zero = 0.0d0

    write (15, '(A)') 'BQ_CHARGES'
    do kb = 1, nbq(im)
        km = bq_list(kb,im)%jm
        na = bq_list(kb,im)%ja
        nb = bq_list(kb,im)%jb
        nc = bq_list(kb,im)%jc
        dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

        ka    = mol_ichg(km)
        nsite = mol_nchg(km)

        do ks = 1, nsite
            write (15,'(f16.10,3f20.10)') chg_pos(0,ka+ks), chg_pos(1:3,ka+ks)+dcel
        end do
    end do

    write (15, '(A)') 'END'
    write (15, *)
end subroutine psi_print_bq

subroutine psi_print_bq2 (im, ip, void)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im, ip, void ! VOID: 0)i and j, 1)i only, 2)j only
    integer :: kidx, km, ka, kb, kc, kat, ks, nsite
    integer :: jm, ja, jb, jc
    real*8  :: dcel(3), zero

    jm   = pair_list(ip,im)%jm
    ja   = pair_list(ip,im)%ja
    jb   = pair_list(ip,im)%jb
    jc   = pair_list(ip,im)%jc

    zero = 0.0d0

    write (15, '(A)') 'BQ_CHARGES'
    do kidx = 1, nbq2
        km    = bq2_list(kidx)%jm
        ka    = bq2_list(kidx)%ja
        kb    = bq2_list(kidx)%jb !MAKE SURE kb isn't simultaneously used for lattice vector and Bq index!
        kc    = bq2_list(kidx)%jc

        dcel(:) = ka*lat(:,1)+kb*lat(:,2)+kc*lat(:,3)
        kat    = mol_ichg(km)
        nsite = mol_nchg(km)

        if (km .eq. im .and. ka .eq. 0 .and. kb .eq. 0 .and. kc .eq. 0 .and. (void .eq. 0 .or. void .eq. 1)) then
            do ks = 1, nsite
                !write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel+1000.0d0, zero
                write (15,'(f16.10,3f20.10)') zero, chg_pos(1:3,kat+ks)+dcel+1000.0d0
            end do
        else if (km .eq. jm .and. ka .eq. ja .and.  kb .eq. jb .and. kc .eq. jc .and. (void .eq. 0 .or. void .eq. 2)) then
            do ks = 1, nsite
                !write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel+1000.0d0, zero
                write (15,'(f16.10,3f20.10)') zero, chg_pos(1:3,kat+ks)+dcel+1000.0d0
            end do
        else
            do ks = 1, nsite
                !write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel, chg_pos(0,kat+ks)
                write (15,'(f16.10,3f20.10)') chg_pos(0,kat+ks), chg_pos(1:3,kat+ks)+dcel
            end do
        end if
    end do
    write (15, '(A)') 'END'

end subroutine psi_print_bq2

subroutine psi_read_hess (nsite, hess0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: hess0(171)
    integer :: row, col, iend, i

    hess0 = 0.0d0
    open (unit=15,file=trim(my_scratch)//'/'//'hess.dat', status='old')
    i = 0
    
    do row = 1, 3*nsite
    do col = 1, row
        i = i + 1
        read (15, *,iostat=iend) hess0(i)
        if (iend /= 0) stop "failed to read hess0(i)"
    end do
    end do
    close (15)
    
    if (i /= 3*nsite + (3*nsite)*(3*nsite-1)/2) stop "BUG didn't read enough into hess0"
end subroutine psi_read_hess

subroutine psi_read_ddipole (nsite, ddipole0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: ddipole0(54)
    integer :: row, mu, i, eof, nread
    character(len=120) :: buffer

    ddipole0 = 0.0d0
    i = 9*nsite
    
    if (i /= 9*nsite) stop "BUG didn't read enough into ddipole0"
end subroutine psi_read_ddipole
