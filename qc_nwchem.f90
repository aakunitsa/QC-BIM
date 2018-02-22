subroutine nw_print_toplevel_directives
    use qc_system
    write (15, '(A)') 'scratch_dir '//trim(my_scratch)
    write (15, '(A)') 'permanent_dir '//trim(my_scratch)
    write (15, '(A,A)') 'start ', trim(jname)
    write (15, '(A)') 'memory total 3800 mb'
    write (15, *)
end subroutine nw_print_toplevel_directives

subroutine qc_nwchem_esp_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im

    open  (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.nw', status='unknown')
    write (15, '(A,I0,A)') 'title "BIM ESP (MONO ', im, ') "'
    call nw_print_toplevel_directives
    call nw_print_geom_one (im, 0, 0, 0)
    call nw_print_basis
    call nw_print_scf
    call nw_print_esp_bq (im)
    write (15, '(A)') 'task scf energy'
    write (15, *)
    call nw_print_esp
    close (15)
end subroutine qc_nwchem_esp_wrt

subroutine qc_nwchem_dipole_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im
    open  (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.nw', status='unknown')
    write (15, '(A,I0,A)') 'title "BIM DIPOLE (MONO ', im, ') "'
    call nw_print_toplevel_directives
    call nw_print_geom_one (im, 0, 0, 0)
    call nw_print_basis
    call nw_print_scf
    call nw_print_esp_bq (im)
    write (15, '(A)') 'task scf energy' ! Dipole is calc'd in multipole analysis of density
    write (15, *)
    close (15)
end subroutine qc_nwchem_dipole_wrt

subroutine qc_nwchem_esp_read (nsite,chg0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: chg0 (nsite)
    integer :: eof, is
    character(len=120) :: buffer
    logical :: found

    chg0 = 0.0d0
    found = .false.

    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out',status='old')

    do 
    read (15, '(A120)', iostat=eof) buffer

    if (eof /= 0) then
        exit
    end if

    if (buffer(39:41) .eq. 'ESP') then
        found = .true.
        read (15, *)
        read (15, *)
        do is = 1, nsite
        read (15, '(A120)') buffer
        read (buffer(46:120),*) chg0(is)
        end do
        exit
    end if
    end do

    close (15)
    if (.not. found) stop "could not find esp charges"

end subroutine qc_nwchem_esp_read

subroutine qc_nwchem_dipole_read (dipole0)
    use qc_system
    implicit none
    real*8, intent(out) :: dipole0(3)
    integer :: eof, i
    character(len=120) :: buffer
    logical :: found

    dipole0 = 0.0d0
    found = .false.
    
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out',status='old')
    do 
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if
        if (buffer(8:25) .eq. 'Multipole analysis') then
            found = .true.
            do i = 1,6
                read (15, *)
            end do
            do i = 1, 3
                read (15, '(A120)') buffer
                read (buffer(15:28),*) dipole0(i) ! NW: read directly in au
            end do
            !write(*,'(A,3f10.4)') "Read Dipole: ", dipole0(1:3)
            exit
        end if
    end do
    
    close (15)
    if (.not. found) stop "could not find dipole moment"
end subroutine qc_nwchem_dipole_read

subroutine qc_nwchem_monomer_wrt (id)
    use qc_system
    implicit none
    integer, intent(in) :: id
    !character(len=500) :: maketempfile

    open (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.nw', status='unknown')
    write (15, '(A,i0,A,i0,A)') 'title "BIM (MONO ', id, '/', nmol, ')"'
    call nw_print_toplevel_directives
    call nw_print_geom_one (id, 0, 0, 0)
    call nw_print_basis
    call nw_print_scf
    call nw_print_theory
    call nw_print_bq (id) ! Monomers embedded
    call nw_print_task

    write (15, *)
    close (15)
    !write(maketempfile, '(a,i0,a)') "cp "//trim(my_scratch)//'/'//trim(jname)//".nw mono", id, ".nw"
    !call system (maketempfile)

end subroutine qc_nwchem_monomer_wrt

subroutine qc_nwchem_dimer_wrt (im, ip)
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

    open (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.nw', status='unknown')
    write (15, '(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A)') &
    'title "BIM DIMER ', &
    im, '(0,0,0)', jm, '(', n(1), ',' ,n(2), ',', n(3), ')  [pair ', &
    ip, ' of ', npair(im), ']"'
    call nw_print_toplevel_directives
    call nw_print_geom_two (im, ip)
    call nw_print_basis
    call nw_print_scf
    call nw_print_theory
    call nw_print_bq2 (im, ip, 0)
    call nw_print_task

    write (15, *)
    close (15)
    !write(maketempfile, '(a,i0,a,i0,a,i0,a)') "cp "//trim(my_scratch)//'/'//& 
    !trim(jname)//".nw dimer", im, "-", jm, "_", n(1), ".nw"
    !call system (maketempfile)

end subroutine qc_nwchem_dimer_wrt

subroutine qc_nwchem_monomerbq2_wrt (im, ip, mono, ghost)
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
        !trim(jname)//".nw mono", im, "_bq", jm, "_", jn(1), ".nw"
    else !--monomer j
        mono_str = '[qm j] '
        id = jm
        na = jn(1)
        nb = jn(2)
        nc = jn(3)
        void = 2 !--void charges of j
        !write(maketempfile, '(a,i0,a,i0,a,i0,a)') "cp "//trim(my_scratch)//'/'// &
        !trim(jname)//".nw bq", im, "_mono", jm, "_", na, ".nw"
    end if

    open (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.nw', status='unknown')
    write (15, '(A,A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A)') &
    'title "BIM MONO-BQ2 ', mono_str, &
    im, '(0,0,0)', jm, '(', jn(1), ',' ,jn(2), ',', jn(3), ')  [pair ', &
    ip, ' of ', npair(im), ']"'
    call nw_print_toplevel_directives
    if (ghost) then
        call nw_print_geom_one_ghost (im, ip, mono)
    else
        call nw_print_geom_one (id, na, nb, nc)
    end if
    call nw_print_basis
    call nw_print_scf
    call nw_print_theory
    call nw_print_bq2 (im, ip, void)
    call nw_print_task
    write (15, *)
    close (15)
    !call system (maketempfile)
end subroutine qc_nwchem_monomerbq2_wrt

subroutine qc_nwchem_pot_read (upot) 
    use qc_system
    implicit none
    real*8,  intent(out) :: upot
    integer :: tstart, tstop, ustart, ustop, eof
    character(len=120) :: buffer
    character(len=20) :: title
    logical :: foundU

    select case(theory_id)
        case (theory_id_scf)
            tstart = 10
            tstop = 18
            ustart = 31
            ustop = 120
            title = 'Total SCF'
        case(theory_id_mp2)
            tstart = 11
            tstop = 19
            ustart = 31
            ustop = 120
            title = 'Total MP2'
        case(theory_id_ccsd)
            tstart = 2
            tstop = 11
            ustart = 31
            ustop = 120
            title = 'Total CCSD'
        case default
            tstart = 10
            tstop = 18
            ustart = 31
            ustop = 120
            title = 'Total SCF'
    end select

    foundU = .false.
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out', status='old')
    do 
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(tstart:tstop) .eq. trim(title)) then
            read (buffer(ustart:ustop), *) upot
            foundU = .true.
            exit
        end if
    end do
    close (15)
    if (.not. foundU) STOP 'ERROR: DID NOT FIND ENERGY'
end subroutine qc_nwchem_pot_read

subroutine qc_nwchem_pot_grad_read (natm, upot, dRR0)
    use qc_system
    implicit none
    integer, intent(in)  :: natm
    real*8,  intent(out) :: upot
    real*8,  intent(out) :: dRR0(3,natm)
    integer :: eof, iat, num
    integer :: gtstart, gtstop, utstart, utstop, ustart, ustop
    character(len=30) :: gtitle, utitle
    character(len=120) :: buffer
    character(len=2)   :: ch2
    real*8 :: xi, yi, zi
    logical :: foundU, foundG

    gtstart = 30
    gtstop = 45
    gtitle = 'ENERGY GRADIENTS'
    select case(theory_id)
        case (theory_id_scf)
            utstart = 10
            utstop = 25
            utitle = 'Total SCF energy'
            ustart = 31
            ustop = 120
        case (theory_id_mp2)
            utstart = 11
            utstop = 26
            utitle = 'Total MP2 energy'
            ustart = 37
            ustop = 120
        case default
            utstart = 10
            utstop = 25
            utitle = 'Total SCF energy'
            ustart = 31
            ustop = 120
    end select

    dRR0 = 0.0d0
    foundU = .false.
    foundG = .false.

    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out', status='old')
    do 
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(utstart:utstop) .eq. trim(utitle)) then
            read (buffer(ustart:ustop), *) upot
            foundU = .true.
        end if

        if (buffer(gtstart:gtstop) .eq. trim(gtitle)) then
            read (15, '(A120)', iostat=eof) buffer
            read (15, '(A120)', iostat=eof) buffer ! atom
            read (15, '(A120)', iostat=eof) buffer ! xyz
            do iat = 1, natm
                read (15, *) num, ch2, xi, yi, zi, dRR0(1:3,iat)
            end do
            foundG = .true.
        end if
        if (foundU .and. foundG) then
            exit
        end if
    end do
    close (15)
    if (.not.foundU) STOP 'WARNING: DID NOT FIND ENERGY'
    if (.not.foundG) STOP 'WARNING: DID NOT FIND GRADIENT'

end subroutine qc_nwchem_pot_grad_read


subroutine qc_nwchem_bq_grad_read (d_bq, im)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    real*8, intent(out) :: d_bq(3, 2*maxnum_bq*3)
    real*8 :: grad(3)
    integer :: kb, km, kk, ks, is, nsite, iend, nchg, ichg
    logical :: valid_file
    
    inquire (FILE=trim(my_scratch)//'/'//trim(jname)//'.bqforce.dat', EXIST=valid_file)
    if (.not. valid_file) then
        STOP "Could not open bqforce.dat"
    end if

    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.bqforce.dat', status='old')
    read (15, *)

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
            do is = 1, nsite
                d_bq(1:3,kk+is) = d_bq(1:3,kk+is) + grad(:)*chg_coef(is, ichg+ks) !--sum gradient_charge into gradient_nuclei
            end do
        end do
        kk = kk + nsite
    end do

    close (15)

end subroutine qc_nwchem_bq_grad_read

subroutine qc_nwchem_bq2_grad_read (d_bq)
    use qc_system
    use qc_neigh
    implicit none
    real*8, intent(out) :: d_bq(3,2*maxnum_bq*3)
    real*8 :: grad(3)
    integer :: kk, kb, km, is, ks, ichg, nsite, nchg
    integer :: iend
    logical :: valid_file
    
    inquire (FILE=trim(my_scratch)//'/'//trim(jname)//'.bqforce.dat', EXIST=valid_file)
    if (.not. valid_file) then
        STOP "Could not open bqforce.dat"
    end if

    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.bqforce.dat', status='old')
    read (15, *)

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
            do is = 1, nsite
                d_bq(1:3,kk+is) = d_bq(1:3,kk+is) + grad(:)*chg_coef(is, ichg+ks) !--sum gradient_charge into gradient_nuclei
            end do
        end do
        kk = kk + nsite
    end do

    close (15)

end subroutine qc_nwchem_bq2_grad_read


subroutine nw_print_geom_one (im, na, nb, nc)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im, na, nb, nc
    integer :: is, ia, nsite
    real*8  :: dcel(3)

    ia    = mol_iatom(im)
    nsite = mol_nsite(im)
    dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

    write (15, '(A,i5)') 'charge', mol_netcharge(im)

    write (15, '(A)') 'geometry units angstroms noautosym noautoz nocenter'
    do is = 1, nsite
    write (15, '(3x,A,3f20.10)') at_atnm(ia+is), gm_pos(1:3,ia+is)+dcel
    end do
    write (15, '(A)') 'end'
    write (15, *)
end subroutine nw_print_geom_one


subroutine nw_print_geom_one_ghost (im, ip, mono)
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

    !--embed in field of ij
    if (mono .eq. 0) then !--monomer i
        mono_id = im
        mono_dcel = 0.0d0
        
        ghost_id = pair_list(ip,im)%jm
        ghost_dcel = pair_list(ip,im)%ja*lat(:,1) + pair_list(ip,im)%jb*lat(:,2) + pair_list(ip,im)%jc*lat(:,3)
    else !--monomer j
        mono_id = pair_list(ip,im)%jm
        mono_dcel = pair_list(ip,im)%ja*lat(:,1) + pair_list(ip,im)%jb*lat(:,2) + pair_list(ip,im)%jc*lat(:,3)
        
        ghost_id = im
        ghost_dcel = 0.0d0
    end if
    
    mono_nsite = mol_nsite(mono_id)
    mono_start = mol_iatom(mono_id)
    ghost_nsite = mol_nsite(ghost_id)
    ghost_start = mol_iatom(ghost_id)
    
    chg0 = mol_netcharge(mono_id)

    write (15, '(A,i5)') 'charge', chg0
    write (15, '(A)') 'geometry units angstroms noautosym noautoz nocenter'

    do is = 1, mono_nsite
        write (15, '(3x,A,3f20.10)') at_atnm(mono_start+is), gm_pos(1:3,mono_start+is)+mono_dcel
    end do

    do is = 1, ghost_nsite
        write (15, '(3x,A,3f20.10)') 'bq'//trim(at_atnm(ghost_start+is)), gm_pos(1:3,ghost_start+is) + ghost_dcel
    end do
    
    write (15, '(A)') 'end'
    write (15, *)

end subroutine nw_print_geom_one_ghost

subroutine nw_print_geom_two (im, ip)
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

    write (15, '(A,i5)') 'charge', chg0

    write (15, '(A)') 'geometry units angstroms noautosym noautoz nocenter'
    do is = 1, nsite_i
    write (15, '(3x,A,3f20.10)') at_atnm(ia+is), gm_pos(1:3,ia+is)
    end do
    do is = 1, nsite_j
    write (15, '(3x,A,3f20.10)') at_atnm(ja+is), (gm_pos(1:3,ja+is) + jbox)
    end do

    write (15, '(A)') 'end'
    write (15, *)
end subroutine nw_print_geom_two

subroutine nw_print_basis
    use qc_system
    implicit none

    write (15, '(A)') 'basis spherical'
!-- 
!    write (15, '(A,A,A)') ' * library ', trim(basis_type), ' file /tmp/kunitsa/aug-cc-pvdz'
!--
    if (basis_type(1:14) .eq. 'aug-cc-pv[dt]z') then
        write (15, '(A)') ' * library aug-cc-pvdz '
    else if (basis_type(1:14) .eq. 'aug-cc-pv[tq]z') then
        write (15, '(A)') ' * library aug-cc-pvtz '
    else
        write (15, '(A,A)') ' * library ', trim(basis_type)
    end if

    ! ADD OTHER FXNS HERE FOR ADDITIONAL GHOST ATOMS
    if (l_bsse) then
        write (15, '(A,A)') ' bqO library O ', trim(basis_type)
        write (15, '(A,A)') ' bqH library H ', trim(basis_type)
    end if
    !--

    write (15, '(A)') 'end'
    write (15, *)
end subroutine nw_print_basis

subroutine nw_print_scf
    write (15, '(A)') "scf "
    write (15, '(A)') "thresh 1e-8"
    !write (15, '(A)') "print low"
    write (15, '(A)') "noprint mulliken"  ! BUG <-- ENSURE THIS DOESNT BREAK ANYTHING
    write (15, '(A)') "end"
    write (15, *)
end subroutine nw_print_scf

subroutine nw_print_theory
    use qc_system 
    if (theory(1:3) .ne. 'scf') then
        write (15, '(A)') theory
        write (15, '(A)') '  freeze atomic'
        if (theory(1:3) .eq. 'mp2') write (15, '(A)') '  tight'
        write (15, '(A)') 'end'
        write (15, *)
    end if
end subroutine nw_print_theory

subroutine nw_print_esp_bq (im)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    integer :: ia, kb, km, na, nb, nc, ka, ks, nsite
    real*8  :: dcel(3)

    if ( maxval(abs(chg_pos)) .lt. 1.0d-10 ) return
    ia = mol_iatom(im)

    write (15, '(A)') 'bq units ang'
    dcel = 0.0d0

    do kb = 1, nbq(im)
        km = bq_list(kb,im)%jm
        na = bq_list(kb,im)%ja
        nb = bq_list(kb,im)%jb
        nc = bq_list(kb,im)%jc
        dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

        ka    = mol_ichg(km)
        nsite = mol_nchg(km)

        do ks = 1, nsite
            write (15,'(3f20.10,f14.8)') chg_pos(1:3,ka+ks)+dcel, chg_pos(0,ka+ks)
        end do
    end do

    write (15, '(A)') 'end'
    write (15, *)
end subroutine nw_print_esp_bq

subroutine nw_print_bq (im)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    integer :: kb, km, na, nb, nc, ka, ks, nsite
    real*8  :: dcel(3) !, zero

    write (15, '(A)') 'bq units ang'
    if (l_bq .and. (.not.l_frq)) write (15, '(A)') ' force '

    !zero = 0.0d0

    do kb = 1, nbq(im)
        km = bq_list(kb,im)%jm
        na = bq_list(kb,im)%ja
        nb = bq_list(kb,im)%jb
        nc = bq_list(kb,im)%jc
        dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

        ka    = mol_ichg(km)
        nsite = mol_nchg(km)

        do ks = 1, nsite
            write (15,'(3f20.10,f14.8)') chg_pos(1:3,ka+ks)+dcel, chg_pos(0,ka+ks)
        end do
    end do

    write (15, '(A)') 'end'
    write (15, *)
end subroutine nw_print_bq


subroutine nw_print_bq2 (im, ip, void)
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

    write (15, '(A)') 'bq units ang'
    if (l_bq) write (15, '(A)') ' force '

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
                write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel, zero
            end do
        else if (km .eq. jm .and. ka .eq. ja .and.  kb .eq. jb .and. kc .eq. jc .and. (void .eq. 0 .or. void .eq. 2)) then
            do ks = 1, nsite
                write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel, zero
            end do
        else
            do ks = 1, nsite
                write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel, chg_pos(0,kat+ks)
            end do
        end if
    end do

write (15, '(A)') 'end'
write (15, *)
end subroutine nw_print_bq2

subroutine nw_print_esp
    write (15, '(A)') 'esp'
    write (15, '(A)') '  recalculate'
    write (15, '(A)') '  restrain hfree'
    write (15, '(A)') 'end'
    write (15, *)
    write (15, '(A)') 'task esp'
    write (15, *)
end subroutine nw_print_esp
    
subroutine nw_read_hess (nsite, hess0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: hess0(171)
    integer :: row, col, iend, i

    hess0 = 0.0d0
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.hess', status='old')
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
end subroutine nw_read_hess

subroutine nw_read_ddipole (nsite, ddipole0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: ddipole0(54)
    integer :: row, mu, iend, i

    ddipole0 = 0.0d0
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.fd_ddipole', status='old')
    i = 0
    
    do row = 1, 3*nsite
    do mu = 1, 3
        i = i + 1
        read (15, *,iostat=iend) ddipole0(i)
        if (iend /= 0) stop "failed to read ddipole0(i)"
    end do
    end do
    close (15)
    if (i /= 9*nsite) stop "BUG didn't read enough into ddipole0"
end subroutine nw_read_ddipole

subroutine nw_print_task
    use qc_system
    if (l_grd) then
        write (15, '(A)') 'task '//theory//' gradient'
    else if (l_frq) then
        write (15, '(A)') 'task '//theory//' hessian numerical'
    else
        write (15, '(A)') 'task '//theory//' energy'
    end if
end subroutine nw_print_task
