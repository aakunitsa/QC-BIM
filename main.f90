program md_nwchem
    use qc_system
    use consts
    use qc_mpi
    use qc_neigh
    use qc_bim
    use qc_lattice
    use hessian
    implicit none
    real*8 :: Umon, Udim, Ulr, en_pot
    real*8 :: lat0(3)
    integer:: irst, job_id, p, iter
    logical :: valid_dir

    call qc_mpi_start

    call read_input (irst, job_id)

    allocate (gm_pos(3,natom))
    allocate (gm_pos_old(3,natom))

    allocate (d_rr (3,natom))
    !allocate (d_mbq (3,natom))

    if (l_md) then
        allocate (vel (3, natom))
        if (use_respa) allocate (corr_d_rr (3,natom))
    end if

    call read_system 

    !--- assign memory ---
    allocate (npair(nmol))
    allocate (nbq(nmol))
    allocate (nlr(nmol))

    ! pair_list, bq_list, bq2_list, d_bq0, d_mbqi, d_mbqj, d_mr
    call qc_allocate_lists 
    
    select case (job_id)
    
    ! ENERGY
    case (job_id_ener)
        call qc_bim_pot (Umon, Udim, Ulr)
        if (sys_master) then
            write(6,*) "       FINAL BIM RESULTS "
            write(6,*) "      -------------------"
            en_pot = Umon + Udim + Ulr
            write (6, '(A,f20.10)') 'E(mono) ', Umon
            write (6, '(A,f20.10)') 'E(dimer)', Udim
            write (6, '(A,f20.10)') 'E(LR)   ', Ulr
            write (6, '(A,f20.10)') 'Etot    ', en_pot
            write (6, '(A,f20.10)') 'Etot/N  ', en_pot/nmol
        end if

    ! ENERGY AND GRADIENT
    case (job_id_engrad)
        call qc_bim_pot (Umon, Udim, Ulr)
        if (sys_master) then
            write(6,*) "       FINAL BIM RESULTS "
            write(6,*) "      -------------------"
            en_pot = Umon + Udim + Ulr
            write (6, '(A,f20.10)') 'E(mono) ', Umon
            write (6, '(A,f20.10)') 'E(dimer)', Udim
            write (6, '(A,f20.10)') 'E(LR)   ', Ulr
            write (6, '(A,f20.10)') 'Etot    ', en_pot
            write (6, '(A,f20.10)') 'Etot/N  ', en_pot/nmol
            write(6,*) ""
            call print_gradients
            write (6, *) ""
            write (6,'(A)') "          STRESS TENSOR (BAR)  "
            write (6,'(A)') "         ---------------------  "
            write (6,'(3f11.1)') virt(1,1:3)/(volume(lat)*ang2bohr**3)*au_bar
            write (6,'(3f11.1)') virt(2,1:3)/(volume(lat)*ang2bohr**3)*au_bar
            write (6,'(3f11.1)') virt(3,1:3)/(volume(lat)*ang2bohr**3)*au_bar
        end if
    
    ! GEOM/CELL OPT
    case (job_id_opt)
        en_pot = 0.0d0
        if (l_freezecell) then
            call qc_gopt_run (en_pot)
        else if (optimizer_id .eq. opt_id_grad) then

            iter = 0
            do p=pstart, pend, pstep
                iter = iter + 1
                if (sys_master) write(*,'(/,A,i0,A)') "CELL OPT: ", p, " bar"
                
                ! set pressure
                pext0 = real(p)

                ! lat0 <-- previous opt lat
                lat0(1) = lat_a
                lat0(2) = lat_b
                lat0(3) = lat_c

                ! update to extrapolated geometry
                if (iter .gt. 2) then 
                    call rescale
                end if
                
                ! optimize
                call qc_cell_gopt_run
                call qc_mpi_barrier
                
                ! 2-point linear extrapolation of (a,b,c)
                if ( iter .gt. 1) then
                    rescale_param(1) = lat_a + (lat_a - lat0(1))
                    rescale_param(2) = lat_b + (lat_b - lat0(2))
                    rescale_param(3) = lat_c + (lat_c - lat0(3))
                    rescale_param(4) = lat_alpha
                    rescale_param(5) = lat_beta
                    rescale_param(6) = lat_gamma
                end if
            end do
        else
            call qc_cell_gopt_run_lbfgs
        end if
    
    ! MD
    case (job_id_md)
        call qc_md_run (irst)
    
    ! HESSIANS
    case (job_id_freq)
        call qc_hess
    
    case (job_id_rescale)
        call rescale

    end select
    !-----

    deallocate (gm_pos)
    deallocate (gm_pos_old)
    deallocate (d_rr)
    deallocate (d_bq0)
    deallocate (d_mbqi)
    deallocate (d_mbqj)
    deallocate (d_mr)
    deallocate (npair)
    deallocate (nbq)
    deallocate (nlr)
    deallocate (pair_list)
    deallocate (bq_list)
    deallocate (bq2_list)
    deallocate (chg_pos, chg_pos_new)
    deallocate (chg_coef)
    deallocate (mol_ichg)
    deallocate (mol_nchg)
    deallocate (at_mass)
    deallocate (at_atnm)
    if (l_md) then
        deallocate (vel)
    end if
    deallocate (mol_nsite, mol_iatom)

    inquire (FILE=trim(my_scratch)//'/.', EXIST=valid_dir)
    if (valid_dir) call system("rm -r "//trim(my_scratch))
    call qc_mpi_end

    stop

end program md_nwchem
!
!..............................................................................
!
