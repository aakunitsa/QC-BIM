subroutine read_input (irst, job_id)
    use qc_system
    use consts
    use qc_mpi
    implicit none
    integer, intent(out) :: irst, job_id
    integer :: is, ifg, commentidx
    integer :: iat, im, im0, is0, lerr, natom0, nsite
    integer :: frg_nmol(10)
    integer :: frg_nsite(10)
    integer :: frg_charge(10)
    real*8  :: frg_mass(10), rcut, rcut_bq, rcut_lr
    character(len=256) :: inputFileName
    character(len=10) :: keyword, theory_buf
    character(len=2) :: frg_atnm(10)
    logical :: valid_dir

    irst  = 0 
    job_id = job_id_ener
    nprod = 1000
    nsave = 10
    dt    = 0.5d0
    temp0 = 300.0d0
    pext0 = 0.0d0
    is_temp_tol = .true.

    frg_atnm  = ' '
    frg_mass = 1.0d0
    frg_charge = 0

    l_grd = .false.
    l_frq = .false.
    l_opt = .false.
    l_md  = .false.

    !-- Thermostat parameters
    nnos = 3
    nresn = 1

    opt_maxiter = 50
    atom_gmax = 0.0006d0
    lat_gmax = 0.0006d0
    n_hess = 1

    frg_nmol  = 0
    frg_nsite = 0
    lerr = 0
    l_bq  = .true.
    l_bsse = .false.
    theory = "scf"

    
    if (sys_master) then
        call get_command_argument(1, inputFileName)
        if (len_trim(inputFileName) .gt. 0) then
            open (11, file=trim(inputFileName), status='old')
            if (sys_master) print *, "Reading input file: ", trim(inputFileName)
        else
            open (11, file='qc_bim.inp', status='old')
            if (sys_master) print *, "Reading input file: qc_bim.inp"
        end if
        !---
        read (11, *) keyword, scratch_top
        read (11, *) keyword, code
        read (11, *) keyword, l_bq
        read (11, *) keyword, l_bsse
        read (11, *) keyword, theory_buf
        read (11, *) keyword, embedding_field
        read (11, *) !BEGIN FRAG
        read (11, *) keyword, nfrg
        read (11, *) keyword, (frg_nmol(ifg),  ifg = 1, nfrg)
        read (11, *) keyword, (frg_nsite(ifg), ifg = 1, nfrg)
        read (11, *) keyword, (frg_charge(ifg), ifg = 1, nfrg)
        nsite = sum (frg_nsite)
        do is = 1, nsite
        read (11, *) frg_atnm(is), frg_mass(is)
        frg_mass(is) = frg_mass(is)*amu2au
        end do
        read (11, *) !END FRAG
    end if

    call qc_mpi_barrier
    call qc_mpi_char_bcast (scratch_top)
    call qc_mpi_char_bcast (code)
    
    
    if (scratch_top(1:3) .eq. 'cwd') then
        scratch_top = '.'
    else
        commentidx = index(scratch_top, '!')
        if (commentidx.gt.0) scratch_top = scratch_top(1:commentidx-1)
    end if 

    inquire (FILE=trim(scratch_top)//'/.', EXIST=valid_dir)
    if (.not. valid_dir) then
        call system("mkdir -p "//trim(scratch_top))
    end if

    inquire (FILE=trim(scratch_top)//'/.', EXIST=valid_dir)
    if (.not. valid_dir) then
        if (sys_master) print *, "Invalid scratch directory"
        call qc_mpi_end
        STOP
    end if
    if (valid_dir.and.sys_master) print *, "Validated scratch dir: ", trim(scratch_top)
    if (code(1:3) .eq. 'g09') then 
        if (sys_master) print *, "Using Gaussian 09"
    else if (code(1:2) .eq. 'nw') then
        if (sys_master) print *, "Using NWChem"
    else if (code(1:3) .eq. 'psi') then
        if (sys_master) print *, "Using PSI4"
    else
        if (sys_master) print *, "Error: code must either be 'g09' (Gaussian09) or 'nw' (NWChem 6.*) or 'psi' (Psi4)"
        call qc_mpi_end
        STOP
    end if

    call qc_mpi_logical_bcast1 (l_bq)
    call qc_mpi_logical_bcast1 (l_bsse)
    call qc_mpi_char_bcast (theory_buf)
    call qc_mpi_char_bcast (embedding_field)
    call qc_mpi_int_bcast1 (nfrg)
    call qc_mpi_int_bcast1 (nsite)
    call qc_mpi_int_bcast  (frg_nmol,nfrg)
    call qc_mpi_int_bcast   (frg_nsite,nfrg)
    call qc_mpi_int_bcast   (frg_charge,nfrg)
    call qc_mpi_real_bcast  (frg_mass,nsite)
    do is = 1, nsite
    call qc_mpi_char_bcast (frg_atnm(is))
    end do
    
    theory = trim(theory_buf)
    if (theory(1:3) .eq. 'scf') then
        theory_id = theory_id_scf
    else if (theory(1:3) .eq. 'mp2') then
        theory_id = theory_id_mp2
    else if (theory(1:4) .eq. 'ccsd') then
        theory_id = theory_id_ccsd
    else if (theory(1:7) .eq. 'ccsd(t)') then
        theory_id = theory_id_ccsd_pt
    else
        if (sys_master) write(*,'(A,A)') "This version does not support theory ", theory
        call qc_mpi_end
        stop
    end if

    if (embedding_field(1:3) .eq. 'esp') then
        embed_id = embed_id_esp
    else if (embedding_field(1:3) .eq. 'dip') then
        embed_id = embed_id_dipole
    else if (embedding_field(1:10) .eq. 'atomdipole') then
        embed_id = embed_id_atomdipole
        if (code(1:3).ne.'g09') STOP "atomdipole only implemented with g09"
    else
        if (sys_master) write(*,'(A,A)') "This version does not support embedding option: ", trim(embedding_field)
    end if

    !-- (NATOM) ---
    natom = 0
    nmol  = 0
    do ifg = 1, nfrg
    nmol  = nmol  + frg_nmol(ifg)
    natom = natom + frg_nmol(ifg)*frg_nsite(ifg)
    end do

    allocate (at_atnm(natom))
    allocate (at_mass(natom))
    allocate (mol_iatom(nmol))
    allocate (mol_nsite(nmol))
    allocate (mol_netcharge(nmol))
    allocate (chg_pos(0:3, 6*natom))
    allocate (chg_pos_new(0:3, 6*natom))
    allocate (chg_coef(1:maxval(frg_nsite), 6*natom))
    allocate (mol_ichg(nmol))
    allocate (mol_nchg(nmol))

    iat = 0
    is0 = 0
    im0 = 0

    do ifg = 1, nfrg
       do im = 1, frg_nmol(ifg)
           im0 = im0 + 1
           mol_iatom(im0) = iat
           mol_nsite(im0) = frg_nsite(ifg)
           mol_netcharge(im0) = frg_charge(ifg)

           do is = 1, frg_nsite(ifg)
              iat = iat + 1
              at_atnm(iat) = frg_atnm(is0+is)
              at_mass(iat) = frg_mass(is0+is)
           end do
        end do
        is0 = is0 + frg_nsite(ifg)
    end do

    rcut   = 9.0d0
    rcut_bq= 15.0d0
    rcut_lr= 50.0d0

    if (sys_master) then
        read (11, *) keyword, basis_type
        read (11, *) keyword, rcut, rcut_bq, rcut_lr !M!
        read (11, *) keyword, fconfig ! Name of geom file
        !M! read (11, *) keyword, box(1:3)  ! in A
        !M! read (11, *) keyword, ncelx_qq, ncely_qq, ncelz_qq  
        !M! read (11, *) keyword, ncelx_bq, ncely_bq, ncelz_bq  
        read (11, *) keyword, job_id
        read (11, *) keyword, pstart, pend, pstep     ! in Bar
        pext0 = real(pstart)

        read (11, *) ! blank
        
        read (11, *) ! BEGIN OPTIMIZE INFO
        read (11, *) keyword, l_freezecell
        read (11, *) keyword, optimizer_id
        read (11, *) keyword, atom_gmax
        read (11, *) keyword, lat_gmax
        read (11, *) keyword, opt_maxiter
        read (11, *) ! END OPTIMIZE INFO

        read (11, *) ! blank

        read (11, *) !BEGIN MD
        read (11, *) keyword, irst   ! irst = 0 : initialize, 1 : restart
        read (11, *) keyword, nprod
        read (11, *) keyword, nsave
        read (11, *) keyword, dt
        read (11, *) keyword, temp0, is_temp_tol
        read (11, *) keyword, nnos
        read (11, *) keyword, nresn
        read (11, *) keyword, use_respa
        read (11, *) keyword, nresp ! RESPA splitting parameter
        read (11, *) !END MD

        read (11, *) ! blank

        read (11, *) ! BEGIN FREQ INFO
        read (11, *) keyword, n_hess(1:3)
        read (11, *) ! END FREQ INFO
        
        read (11, *) ! blank

        read (11, *) ! BEGIN RESCALE INFO
        read (11, *) keyword, rescale_param(1:6)
        read (11, *) ! END RESCALE INFO

        if (irst .eq. 0) then
            !--- Read Natom ---
            open (15, file=fconfig,status='old')
            read (15, *) Natom0
            close(15)
        else 
            open (15, file='mdrr.sav',status='old')
            read (15, *) Natom0
            close(15)
        end if
        
        if (natom .ne. natom0) lerr = 1
        close(11)

    end if

    call qc_mpi_barrier
    call qc_mpi_int_bcast1 (lerr)

    if (lerr .eq. 1) then
        if (sys_master) print *, 'Error in Natom '
        call qc_mpi_end
        stop
    end if

    call qc_mpi_char_bcast  (basis_type)
    call qc_mpi_real_bcast1 (rcut)
    call qc_mpi_real_bcast1 (rcut_bq)
    call qc_mpi_real_bcast1 (rcut_lr)
    
    call qc_mpi_char_bcast  (fconfig)
    call qc_mpi_int_bcast1 (job_id)
    call qc_mpi_int_bcast1 (pstart)
    call qc_mpi_int_bcast1 (pend)
    call qc_mpi_int_bcast1 (pstep)
    call qc_mpi_real_bcast1 (pext0)
    
    call qc_mpi_logical_bcast1 (l_freezecell)
    call qc_mpi_int_bcast1 (optimizer_id)
    call qc_mpi_real_bcast1 (atom_gmax)
    call qc_mpi_real_bcast1 (lat_gmax)
    call qc_mpi_int_bcast1 (opt_maxiter)

    call qc_mpi_int_bcast1 (irst)
    call qc_mpi_int_bcast1 (nprod)
    call qc_mpi_int_bcast1 (nsave)
    call qc_mpi_real_bcast1 (dt)
    call qc_mpi_real_bcast1 (temp0)
    call qc_mpi_logical_bcast1 (is_temp_tol)
    call qc_mpi_int_bcast1 (nnos)          ! AK
    call qc_mpi_int_bcast1 (nresn)         ! AK
    call qc_mpi_logical_bcast1 (use_respa) ! AK
    call qc_mpi_int_bcast1 (nresp)         ! AK
    
    call qc_mpi_int_bcast (n_hess, 3)
    call qc_mpi_real_bcast (rescale_param, 6)

    rcut_qq2 = rcut*rcut
    rcut_bq2 = rcut_bq*rcut_bq
    rcut_lr2 = rcut_lr*rcut_lr !M!

    select case (job_id)
    case (job_id_ener)
        l_grd = .false.
        l_frq = .false.
        l_opt = .false.
        l_md = .false.
    case (job_id_engrad)
        l_grd = .true.
        l_frq = .false.
        l_opt = .false.
        l_md = .false.
    case (job_id_opt)
        l_grd = .true.
        l_frq = .false.
        l_opt = .true.
        l_md = .false.
    case (job_id_md)
        l_grd = .true.
        l_frq = .false.
        l_opt = .false.
        l_md = .true.
    case (job_id_freq)
        l_grd = .false.
        l_frq = .true.
        l_opt = .false.
        l_md = .false.
    end select

    if (job_id.eq.job_id_freq .and. l_bq .and. embed_id.ne.embed_id_esp) then
        if (sys_master) print *, 'Classical correction to Hessian is only implemented with atom-centered charges'
        call qc_mpi_end
        stop
    end if

    ensemble_nvt = .true.

    !--
    !temp0 = temp0*boltz  ! K --> A.U.
    !beta  = 1.0d0/temp0 !        
    dt    = dt / au_time
    !print *, 'dt ', dt

    !--
    !--
    if (sys_nproc .gt. 1) then
        if (sys_me .le. 9) then
            write (jname,'(a,i1)') 'job_cpu', sys_me
            write (my_scratch,'(a,i1)') trim(scratch_top)//'/job_cpu', sys_me
            call system("mkdir -p "//trim(my_scratch))
        else if (sys_me .le. 99) then
            write (jname,'(a,i2)') 'job_cpu', sys_me
            write (my_scratch,'(a,i2)') trim(scratch_top)//'/job_cpu', sys_me
            call system("mkdir -p "//trim(my_scratch))
        else if (sys_me .le. 999) then
            write (jname,'(a,i3)') 'job_cpu', sys_me
            write (my_scratch,'(a,i3)') trim(scratch_top)//'/job_cpu', sys_me
            call system("mkdir -p "//trim(my_scratch))
        else if (sys_me .le. 9999) then
            write (jname,'(a,i4)') 'job_cpu', sys_me
            write (my_scratch,'(a,i4)') trim(scratch_top)//'/job_cpu', sys_me
            call system("mkdir -p "//trim(my_scratch))
        else
            print *, "Update code to allow > 9999 procs!"
            stop
        end if
    else
        jname = 'job_cpu0'
        my_scratch = trim(scratch_top)//'/job_cpu0'
        call system("mkdir -p "//trim(my_scratch))
    end if
    
    if (code(1:3) .eq. 'g09') then
        jcmd = '( g09 < '//trim(my_scratch)//'/'//trim(jname)//'.in -scrdir='//trim(my_scratch)//' ) > ' &
        //trim(my_scratch)//'/'//trim(jname)//'.out'
    else if (code(1:2) .eq. 'nw') then
        jcmd = 'nwchem.x '//trim(my_scratch)//'/'//trim(jname)//'.nw  > '//trim(my_scratch)// &
        '/'//trim(jname)//'.out'
    else if (code(1:3) .eq. 'psi') then
        jcmd = 'psi4wrapper.py '//trim(my_scratch)//'/'//trim(jname)//'.in > '//trim(my_scratch)// &
        '/'//trim(jname)//'.out'
    end if

    if (sys_master) then
        write(*,*) '   BIM INPUT   '
        write(*,*) '---------------'
        write (*,'(3x,i0,A)') nmol, " molecules"
        write (*,'(3x,i0,A)') natom, " atoms"
        write(*,*) ""
        write(*,'(3x,A,L2)') 'Use BQ Force?', l_bq
        write(*,'(3x,A,L2)') 'BSSE Correct?',l_bsse
        write(*,'(3x,A,A)')  'Code: ', trim(code)
        write(*,'(3x,A,A)')  'Theory: ', theory
        write(*,'(3x,A,A)') "Basis: ", basis_type
        write(*,'(3x,A,A)') "Embedding: ", trim(embedding_field)
        write(*,'(3x,A,3f10.2)') "QM/BQ/LR cutoffs (Angstrom): ", rcut, rcut_bq, rcut_lr
        write(*,'(3x,A,A)') "Input geometry file: ", fconfig
        write(*,'(3x,A,i3)') "Job ID: ", job_id
        write(*,*) ""
        write(*,'(3x,A,f10.3)') "Ext pressure (bar):", pext0
        write(*,'(3x,A,L2)') "Gradient calc?", l_grd
        write(*,'(3x,A,L2)') "Freq calc?", l_frq
        write(*,'(3x,A,L2)') "Optimization?", l_opt
        write(*,'(3x,A,L2)') "MD?", l_md
        write(*,*) ""
    end if
    return

end subroutine read_input

subroutine read_system 
    use qc_system
    use qc_lattice
    use consts
    use qc_mpi
    implicit none
    character(len=2)   :: ch2
    integer :: im, ia, is, i

    !----------------------------------------------------------------------------!
    !---  read the number of atoms, dimension of box and atomic coordinates -----!
    !----------------------------------------------------------------------------!
    if (sys_master) then
        open(10,file=fconfig,status='old')  ! 
        read(10,*) !Natom
        read(10,*) lat_a, lat_b, lat_c, lat_alpha, lat_beta, lat_gamma, lat_axis

        do i=1, Natom
        read(10,*)ch2, gm_pos(1:3, i)
        enddo
        close (10)
    end if


    call qc_mpi_barrier
    call qc_mpi_real_bcast (gm_pos, 3*natom)
    call qc_mpi_real_bcast1 (lat_a)
    call qc_mpi_real_bcast1 (lat_b)
    call qc_mpi_real_bcast1 (lat_c)
    call qc_mpi_real_bcast1 (lat_alpha)
    call qc_mpi_real_bcast1 (lat_beta)
    call qc_mpi_real_bcast1 (lat_gamma)
    call qc_mpi_int_bcast1 (lat_axis)

    call create_lattice_vectors
    call update_scaling_matrix
    
    call qc_mpi_real_bcast (lat, 9)
    call qc_mpi_real_bcast (lati, 9)

    if (sys_master) then
        write(*,*)   '           INPUT GEOMETRY            '
        write (6, *) ' mol(at)          coordinates        '
        write (6, *) '-------------------------------------'
        do im = 1,nmol
            ia = mol_iatom(im)
            do is=1,mol_nsite(im)
                if (ia+is .lt. 10) then
                    write (6, '(1x,i0,A,i0,A,A,3f10.3)') im,'(',ia+is,') ', at_atnm(ia+is), gm_pos(1:3, ia+is)
                else
                    write (6, '(i0,A,i0,A,A,3f10.3)') im,'(',ia+is,') ', at_atnm(ia+is), gm_pos(1:3, ia+is)
                end if
            end do
        end do
        
        write (*,'(A)') "Lattice vectors:"
        if (a_dim) write(*,'(3x,A,f6.3,A,f6.3,A,f6.3,A)') "(", lat(1,1), ",", lat(2,1), ",", lat(3,1), " )"
        if (b_dim) write(*,'(3x,A,f6.3,A,f6.3,A,f6.3,A)') "(", lat(1,2), ",", lat(2,2), ",", lat(3,2), " )"
        if (c_dim) write(*,'(3x,A,f6.3,A,f6.3,A,f6.3,A)') "(", lat(1,3), ",", lat(2,3), ",", lat(3,3), " )"
        write (*,'(A,f8.3)') 'Volume (Ang**3) ', volume(lat)
        if (.not.(a_dim .or. b_dim .or. c_dim)) write(*,'(3x,A)') "no pbc"
        write (*,*) ""
    end if

end subroutine read_system
