subroutine print_route
    use qc_system
    !write (15, '(A)') '%Mem=4000MB'
    write (15, '(A)') '%Chk='//trim(my_scratch)//'/'//trim(jname)//'.chk'
    write (15, *)
end subroutine print_route

subroutine g09_esp_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im

    open  (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    call print_route
    if (maxval(abs(chg_pos)) .lt. 1.0d-10) then
        !write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' NoSymm Density Pop=CHelpG'
        write (15, '(A)') '#n HF/'//trim(basis_type)//' NoSymm Pop=ESP' ! NOTE THAT POP=ESPDipole can improve long range behavior of esp!
    else
        !write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Charge NoSymm Density Pop=CHelpG'
        write (15, '(A)') '#n HF/'//trim(basis_type)//' Charge NoSymm Pop=ESP'
    end if
    write (15, *) 
    write (15, '(A,I0,A)') 'BIM ESP (MONO ', im, ')'
    write (15, *)
    call g09_print_geom_one (im, 0, 0, 0)
    call g09_print_esp_bq (im)
    write (15, *)
    close (15)
end subroutine g09_esp_wrt

subroutine g09_dipole_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im

    open  (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    call print_route
    if (maxval(abs(chg_pos)) .lt. 1.0d-10) then
        !write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' NoSymm Density'
        write (15, '(A)') '#n HF/'//trim(basis_type)//' NoSymm'
    else
        !write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Charge NoSymm Density'
        write (15, '(A)') '#n HF/'//trim(basis_type)//' Charge NoSymm'
    end if
    write (15, *) 
    write (15, '(A,I0,A)') 'BIM DIPOLE (MONO ', im, ')'
    write (15, *)
    call g09_print_geom_one (im, 0, 0, 0)
    call g09_print_esp_bq (im)
    write (15, *)
    close (15)
end subroutine g09_dipole_wrt

subroutine g09_atomdipole_wrt (im)
    use qc_system
    implicit none
    integer, intent(in) :: im

    open  (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    call print_route
    if (maxval(abs(chg_pos)) .lt. 1.0d-10) then
        !write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' NoSymm Density'
        write (15, '(A)') '#n HF/'//trim(basis_type)//' NoSymm Pop=(AtomDipole)'
    else
        !write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Charge NoSymm Density'
        write (15, '(A)') '#n HF/'//trim(basis_type)//' Charge NoSymm Pop=(AtomDipole)'
    end if
    write (15, *) 
    write (15, '(A,I0,A)') 'BIM ATOM-DIPOLE (MONO ', im, ')'
    write (15, *)
    call g09_print_geom_one (im, 0, 0, 0)
    call g09_print_esp_bq (im)
    write (15, *)
    close (15)
end subroutine g09_atomdipole_wrt

subroutine g09_esp_read (nsite,chg0)
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

    if (buffer(2:21) .eq. 'Charges from ESP fit') then
        found = .true.
        read (15, *)
        read (15, *)
        do is = 1, nsite
        read (15, '(A120)') buffer
        read (buffer(12:120),*) chg0(is)
        end do
        exit
    end if
    end do

    close (15)
    if (.not. found) stop "could not read esp charges"

end subroutine g09_esp_read

subroutine g09_dipole_read (dipole0)
    use qc_system
    implicit none
    real*8, intent(out) :: dipole0(3)
    integer :: eof
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
        if (buffer(2:14) .eq. 'Dipole moment') then
            read (15, '(A120)', iostat=eof) buffer
            read (buffer( 7:26),*) dipole0(1) ! G09: read in Debye
            read (buffer(33:52),*) dipole0(2)
            read (buffer(59:78),*) dipole0(3)
            dipole0(:) = dipole0(:) / 2.541766 ! Debye --> a.u.
            !write(*,'(A,3f10.4)') "Read Dipole: ", dipole0(1:3)
            found = .true.
            exit
        end if
    end do
    close (15)
    if (.not. found) stop "could not find dipole moment"
end subroutine g09_dipole_read

subroutine g09_atomdipole_read (nsite, charges, dipoles)
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
    
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out',status='old')
    do
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if
        !if (buffer(2:8).eq.'Charge=') print *, buffer
        if (buffer(12:19).eq.'Charges' .and. buffer(29:42).eq.'Point Dipoles') then
            read (15, '(A120)', iostat=eof) buffer ! 1 2 3 4
            !if (sys_master) print *, ""
            do is = 1, nsite
                read (15, '(A120)', iostat=eof) buffer
                read(buffer, '(i6,A3,f12.6,3f11.6)') idx, symbol, charges(is), dipoles(1:3, is)
                !if (sys_master) then
                !    
                !    write(*,'(A, i0, 2X, A, 2X, A, f6.3, A, 3f6.3)') "Atom ", idx, symbol, "charge=", charges(is), &
                !    " dipole=", dipoles(:,is)
                !end if
            end do
            found = .true.
            exit
        end if
    end do
    close (15)
    if (.not. found) stop "could not find atomic charges/dipoles"
end subroutine g09_atomdipole_read

subroutine g09_monomer_wrt (id)
    use qc_system
    implicit none
    integer, intent(in) :: id
    !character(len=500) :: maketempfile

    open (unit=15, file=trim(my_scratch)//'/'//trim(jname)//'.in', status='unknown')
    call print_route
    if (l_grd) then
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Force Charge NoSymm Density Prop=(Field, Read)'
    else if (l_frq) then
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Freq=Numer Charge NoSymm'
    else
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Charge NoSymm'
    end if
    write (15, *) 
    write (15, '(A,i0,A,i0,A)') 'BIM (MONO ', id, '/', nmol, ')'
    write (15, *) 
    call g09_print_geom_one (id, 0, 0, 0)
    call g09_print_bq (id) ! Monomers embedded

    write (15, *)
    close (15)
    !write(maketempfile, '(a,i0,a)') "cp "//trim(my_scratch)//'/'//trim(jname)//".in mono", id, ".in"
    !call system (maketempfile)
end subroutine g09_monomer_wrt

subroutine g09_dimer_wrt (im, ip)
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
    call print_route
    if (l_grd) then
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Force Charge NoSymm Density Prop=(Field, Read)'
    else if (l_frq) then
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Freq=Numer Charge NoSymm'
    else
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Charge NoSymm'
    end if
    write (15, *) 
    write (15, '(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A)') &
    'BIM DIMER ', &
    im, '(0,0,0)', jm, '(', n(1), ',' ,n(2), ',', n(3), ')  [pair ', &
    ip, ' of ', npair(im), ']'
    write (15, *) 
    call g09_print_geom_two (im, ip)
    call g09_print_bq2 (im, ip, 0)

    write (15, *)
    close (15)
    !write(maketempfile, '(a,i0,a,i0,a,i0,a)') "cp "//trim(my_scratch)//'/'//& 
    !trim(jname)//".in dimer", im, "-", jm, "_", n(1), ".in"
    !call system (maketempfile)

end subroutine g09_dimer_wrt

subroutine g09_monomerbq2_wrt (im, ip, mono, ghost)
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
    call print_route
    if (l_grd) then
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Force Charge NoSymm Density Prop=(Field, Read)'
    else if (l_frq) then
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Freq=Numer Charge NoSymm'
    else
        write (15, '(A)') '#n '//trim(theory)//'/'//trim(basis_type)//' Charge NoSymm'
    end if
    
    write(15, *)
    write (15, '(A,A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A)') &
    'BIM MONO-BQ2 ', mono_str, &
    im, '(0,0,0)', jm, '(', jn(1), ',' ,jn(2), ',', jn(3), ')  [pair ', &
    ip, ' of ', npair(im), ']'
    write (15, *)

    if (ghost) then
        call g09_print_geom_one_ghost (im, ip, mono)
    else
        call g09_print_geom_one (id, na, nb, nc)
    end if
    call g09_print_bq2 (im, ip, void)
    write (15, *)
    close (15)
    !call system (maketempfile)
end subroutine g09_monomerbq2_wrt

subroutine g09_pot_read (upot) 
    use qc_system
    implicit none
    real*8,  intent(out) :: upot
    real*8 :: selfE
    integer :: eof
    character(len=120) :: buffer
    logical :: foundU

    upot = 0.0d0
    selfE = 0.0d0
    foundU = .false.
    ! Self energy of charges
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out', status='old')
    do
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(2:5) .eq. 'Self') then
            read (buffer(30:50), *) selfE
            exit
        end if
    end do
    close (15)

    ! Read Total Energy
    call system('formchk '//trim(my_scratch)//'/'//trim(jname)//'.chk > /dev/null')
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.fchk', status='old')
    do
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(1:12) .eq. 'Total Energy') then
            read (buffer(index(buffer,'R')+1:120), *) upot
            foundU = .true.
            exit
        end if
    end do
    close (15)
    upot = upot - selfE
    if (.not. foundU) STOP 'ERROR: DID NOT FIND ENERGY'
end subroutine g09_pot_read

subroutine g09_pot_grad_read (natm, upot, dRR0)
    use qc_system
    implicit none
    integer, intent(in)  :: natm
    real*8,  intent(out) :: upot
    real*8,  intent(out) :: dRR0(3,natm)
    real*8 :: selfE
    integer :: eof, iat, mu, nread
    character(len=120) :: buffer
    logical :: foundU, foundG


    dRR0 = 0.0d0
    upot = 0.0d0
    selfE = 0.0d0
    foundU = .false.
    foundG = .false.
    
    ! Self energy of charges
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out', status='old')
    do
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(2:5) .eq. 'Self') then
            read (buffer(30:50), *) selfE
            exit
        end if
    end do
    close (15)

    ! Read Total E and Gradients
    call system('formchk '//trim(my_scratch)//'/'//trim(jname)//'.chk > /dev/null')
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.fchk', status='old')
    do 
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(1:12) .eq. 'Total Energy') then
            read (buffer(index(buffer,'R')+1:120), *) upot
            foundU = .true.
        end if

        if (buffer(1:18) .eq. 'Cartesian Gradient') then
            foundG = .true.
            read (15, '(A120)', iostat=eof) buffer
            nread = 0
            do iat = 1, natm
            do mu = 1, 3
                read(buffer(nread*16+1:nread*16+16), '(F16.8)') dRR0(mu, iat)
                nread = nread + 1
                if (nread .eq. 5) then
                    read (15, '(A120)', iostat=eof) buffer
                    nread = 0
                end if
            end do
            end do
        end if
        
        if (foundU .and. foundG) then
            exit
        end if
    end do
    close (15)
    upot = upot - selfE
    if (.not.foundU) STOP 'WARNING: DID NOT FIND ENERGY'
    if (.not.foundG) STOP 'WARNING: DID NOT FIND GRADIENT'

end subroutine g09_pot_grad_read

subroutine g09_bq_grad_read (d_bq, im)
    ! Read electric field at Bq sites 
    ! multiply by -1.0*charge to get gradient
    ! This gives same results as NW "Bq Force"
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    real*8, intent(out) :: d_bq(3, 2*maxnum_bq*3)
    real*8 :: grad(3)
    character(len=120) :: buffer
    integer :: is, kb, km, kk, ks, nsite, ichg, nchg
    integer :: eof 
    logical :: found

    d_bq = 0.0d0
    found = .false.
    
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out',status='old')
    do 
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(42:55) .eq. 'Electric Field') then
            found = .true.
            read (15, '(A120)', iostat=eof) buffer ! XYZ
            read (15, '(A120)', iostat=eof) buffer ! ---
            do is = 1, mol_nsite(im)
                read (15, '(A120)', iostat=eof) buffer ! field at qm site i
            end do
            
            kk = 0
            do kb = 1, nbq(im) !--for every Bq molecule
                km = bq_list(kb,im)%jm
                nsite = mol_nsite(km)
                nchg = mol_nchg(km)
                ichg = mol_ichg(km)
                do ks = 1, nchg !--for every charge in Bq molecule
                    read (15, '(A120)', iostat=eof) buffer 
                    read (buffer(25:66), '(3F14.6)') grad(1:3) ! read E field at bq site
                    grad(:) = -1.0d0*chg_pos(0, ichg+ks)*grad(:) ! convert field to grad_bq
                    do is = 1, nsite
                        d_bq(1:3,kk+is) = d_bq(1:3,kk+is) + grad(:)*chg_coef(is, ichg+ks) ! sum grad_bq into grad_nuclei
                    end do
                end do
                kk = kk + nsite
            end do
        
        end if
    end do

    close (15)
    if (.not. found) STOP "DID NOT FIND ELECTRIC FIELD AT BQ SITE"

end subroutine g09_bq_grad_read

subroutine g09_bq2_grad_read (d_bq)
    ! Read electric field at Bq sites 
    ! multiply by -1.0*charge to get gradient
    ! Extra logic: some charges are voided (!!!)
    ! The voiding is performed in qc_bim.f90 calc_dimer
    use qc_system
    use qc_neigh
    implicit none
    real*8, intent(out) :: d_bq(3,2*maxnum_bq*3)
    real*8 :: grad(3)
    character(len=120) :: buffer
    integer :: is, kb, km, kk, ks, nsite, ichg, nchg, eof
    logical :: found

    d_bq = 0.0d0
    found = .false.

    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.out',status='old')
    do 
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(42:55) .eq. 'Electric Field') then
            found = .true.
            read (15, '(A120)', iostat=eof) buffer ! XYZ
            read (15, '(A120)', iostat=eof) buffer ! ---
            read (15, '(A120)', iostat=eof) buffer ! 1 Atom

            do while (buffer(7:10).eq.'Atom')
                read (15, '(A120)', iostat=eof) buffer ! field at qm site i
            end do
            
            kk = 0
            do kb = 1, nbq2 !--for every Bq molecule
                km = bq2_list(kb)%jm
                nsite = mol_nsite(km)
                nchg = mol_nchg(km)
                ichg = mol_ichg(km)
                do ks = 1, nchg !--for every charge in Bq molecule
                    read (buffer(25:66), '(3F14.6)') grad(1:3)
                    grad(:) = -1.0d0*chg_pos(0, ichg+ks)*grad(:) ! convert field to grad_bq
                    do is = 1, nsite
                        d_bq(1:3,kk+is) = d_bq(1:3,kk+is) + grad(:)*chg_coef(is, ichg+ks) ! sum grad_bq into grad_nuclei
                    end do
                    read (15, '(A120)', iostat=eof) buffer ! field at bq site
                end do
                kk = kk + nsite
            end do
        end if
    end do

    close (15)
    if (.not. found) STOP "DID NOT FIND ELECTRIC FIELD AT BQ SITE"

end subroutine g09_bq2_grad_read

subroutine g09_print_geom_one (im, na, nb, nc)
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

    if (chg0 .eq. 0) then
        write (15, '(A)') '0 1' ! Charge 0 Singlet
    else
        write (15, '(i0, A)'), chg0, ' 1' ! Charge Doublet
    end if

    do is = 1, nsite
    write (15, '(A,3f20.10)') at_atnm(ia+is), gm_pos(1:3,ia+is)+dcel
    end do
    write (15, *)
end subroutine g09_print_geom_one


subroutine g09_print_geom_one_ghost (im, ip, mono)
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
    
    if (chg0 .eq. 0) then
        write (15, '(A)') '0 1' ! Charge 0 Singlet
    else
        write (15, '(i0, A)'), chg0, ' 2' ! Charge Doublet
    end if

    do is = 1, mono_nsite
        write (15, '(A,3f20.10)') at_atnm(mono_start+is), gm_pos(1:3,mono_start+is)+mono_dcel
    end do

    do is = 1, ghost_nsite
        write (15, '(3x,A,3f20.10)') 'bq'//trim(at_atnm(ghost_start+is)), gm_pos(1:3,ghost_start+is) + ghost_dcel
    end do
    
    write (15, *)

end subroutine g09_print_geom_one_ghost

subroutine g09_print_geom_two (im, ip)
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

    if (chg0 .eq. 0) then
        write (15, '(A)') '0 1' ! Charge 0 Singlet
    else
        write (15, '(i0, A)'), chg0, ' 2' ! Charge Doublet
    end if

    do is = 1, nsite_i
    write (15, '(A,3f20.10)') at_atnm(ia+is), gm_pos(1:3,ia+is)
    end do
    do is = 1, nsite_j
    write (15, '(3x,A,3f20.10)') at_atnm(ja+is), (gm_pos(1:3,ja+is) + jbox)
    end do

    write (15, *)
end subroutine g09_print_geom_two

subroutine g09_print_esp_bq (im)
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

        do ks = 1, nsite
            write (15,'(3f20.10,f14.8)') chg_pos(1:3,ka+ks)+dcel, chg_pos(0,ka+ks)
        end do
    end do

    write (15, *)
end subroutine g09_print_esp_bq

subroutine g09_print_bq (im)
    use qc_system
    use qc_neigh
    implicit none
    integer, intent(in) :: im
    integer :: kb, km, na, nb, nc, ka, ks, nsite
    real*8  :: dcel(3) !, zero

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

    if (l_grd) then
        write (15, *)
        do kb = 1, nbq(im)
            km = bq_list(kb,im)%jm
            na = bq_list(kb,im)%ja
            nb = bq_list(kb,im)%jb
            nc = bq_list(kb,im)%jc
            dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

            ka    = mol_ichg(km)
            nsite = mol_nchg(km)

            do ks = 1, nsite
                write (15,'(3f20.10)') chg_pos(1:3,ka+ks)+dcel
            end do
        end do
    end if

    write (15, *)
end subroutine g09_print_bq


subroutine g09_print_bq2 (im, ip, void)
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
                write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel+1000.0d0, zero
            end do
        else if (km .eq. jm .and. ka .eq. ja .and.  kb .eq. jb .and. kc .eq. jc .and. (void .eq. 0 .or. void .eq. 2)) then
            do ks = 1, nsite
                write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel+1000.0d0, zero
            end do
        else
            do ks = 1, nsite
                write (15,'(3f20.10,f14.8)') chg_pos(1:3,kat+ks)+dcel, chg_pos(0,kat+ks)
            end do
        end if
    end do

    if (l_grd) then
        write (15, *)
        do kidx = 1, nbq2
        km    = bq2_list(kidx)%jm
        ka    = bq2_list(kidx)%ja
        kb    = bq2_list(kidx)%jb !MAKE SURE kb isn't simultaneously used for lattice vector and Bq index!
        kc    = bq2_list(kidx)%jc

        dcel(:) = ka*lat(:,1)+kb*lat(:,2)+kc*lat(:,3)
        kat    = mol_ichg(km)
        nsite = mol_nchg(km)

        ! NEW ---
        if (km .eq. im .and. ka .eq. 0 .and. kb .eq. 0 .and. kc .eq. 0 .and. (void .eq. 0 .or. void .eq. 1)) then
            do ks = 1, nsite
                write (15,'(3f20.10)') chg_pos(1:3,kat+ks)+dcel+1000.0d0
            end do
        else if (km .eq. jm .and. ka .eq. ja .and.  kb .eq. jb .and. kc .eq. jc .and. (void .eq. 0 .or. void .eq. 2)) then
            do ks = 1, nsite
                write (15,'(3f20.10)') chg_pos(1:3,kat+ks)+dcel+1000.0d0
            end do
        else
            do ks = 1, nsite
                write (15,'(3f20.10)') chg_pos(1:3,kat+ks)+dcel
            end do
        end if
        end do
    end if

write (15, *)
end subroutine g09_print_bq2

subroutine g09_read_hess (nsite, hess0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: hess0(171)
    integer :: row, col, i, eof, nread
    character(len=120) :: buffer
    logical :: found

    hess0 = 0.0d0
    found = .false.
    i = 0
    call system('formchk '//trim(my_scratch)//'/'//trim(jname)//'.chk > /dev/null')
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.fchk', status='old')
    
    do
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(1:25) .eq. 'Cartesian Force Constants') then
            found = .true.
            
            read (15, '(A120)', iostat=eof) buffer
            nread = 0
            
            i = 0
            do row = 1, 3*nsite
            do col = 1, row
                i = i + 1
                read(buffer(nread*16+1:nread*16+16), '(F16.8)') hess0(i)
                nread = nread + 1
                if (nread .eq. 5) then
                    read (15, '(A120)', iostat=eof) buffer
                    nread = 0
                end if
            end do
            end do
        end if
    end do
    close (15)
    
    if (.not. found) stop "did not find hessian"
    if (i /= 3*nsite + (3*nsite)*(3*nsite-1)/2) stop "BUG didn't read enough into hess0"
end subroutine g09_read_hess

subroutine g09_read_ddipole (nsite, ddipole0)
    use qc_system
    implicit none
    integer, intent(in) :: nsite
    real*8, intent(out) :: ddipole0(54)
    integer :: row, mu, i, eof, nread
    character(len=120) :: buffer
    logical :: found

    ddipole0 = 0.0d0
    found = .false.
    i = 0
    call system('formchk '//trim(my_scratch)//'/'//trim(jname)//'.chk > /dev/null')
    open (unit=15,file=trim(my_scratch)//'/'//trim(jname)//'.fchk', status='old')
    
    do
        read (15, '(A120)', iostat=eof) buffer
        if (eof /= 0) then
            exit
        end if

        if (buffer(1:18) .eq. 'Dipole Derivatives') then
            found = .true.
            
            read (15, '(A120)', iostat=eof) buffer
            nread = 0
            
            i = 0
            do row = 1, 3*nsite
            do mu = 1, 3
                i = i + 1
                read(buffer(nread*16+1:nread*16+16), '(F16.8)') ddipole0(i)
                nread = nread + 1
                if (nread .eq. 5) then
                    read (15, '(A120)', iostat=eof) buffer
                    nread = 0
                end if
            end do
            end do
        end if
    end do
    close (15)
    
    if (.not. found) stop "did not find ddipole"
    if (i /= 9*nsite) stop "BUG didn't read enough into ddipole0"
end subroutine g09_read_ddipole
