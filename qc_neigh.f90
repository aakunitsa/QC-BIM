module triple_type
    type triple
        integer :: a, b, c
    end type triple
end module triple_type

module lat_queue
    use triple_type, QUEUE_DATA => triple
    logical, dimension(:,:,:), allocatable :: vis
    include "queues.f90"

    function visited(newvec)
        type(QUEUE_DATA), intent(in) :: newvec
        logical :: visited
        visited = vis(newvec%a, newvec%b, newvec%c)
        vis(newvec%a, newvec%b, newvec%c) = .true.
    end function visited
end module lat_queue

module qc_neigh
    use qc_system
    use consts
    use lat_queue
    implicit none

    type, public :: neigh_typ
        integer :: jm, ja, jb, jc
        logical :: qm
    end type neigh_typ

    type, public :: pair_typ
        integer :: jm, ja, jb, jc
        real*8 :: scale
    end type pair_typ

    integer, dimension(:), allocatable :: npair(:)
    integer, dimension(:), allocatable :: nbq(:), nlr(:)
    integer :: num_pairs
    integer :: nbq2

    type(pair_typ), dimension(:,:), allocatable, target :: pair_list
    type(neigh_typ), dimension(:,:), allocatable, target :: bq_list
    type(neigh_typ), dimension(:), allocatable :: bq2_list

contains

subroutine qc_get_lists (overflow)
    use qc_mpi
    type(QUEUE_STRUCT), pointer :: lat_vecs
    type(QUEUE_DATA) :: vec, nvec
    logical, intent(out) :: overflow
    logical :: success
    integer :: nhits

    num_pairs = 0
    npair = 0
    nbq = 0
    
    !---BREADTH-FIRST-SEARCH: NEIGHBOR CELLS
    call queue_create( lat_vecs, 10000 )
    allocate ( vis(-100:100,-100:100,-100:100) )
    vis = .false.
    
    vec%a = 0; vec%b = 0; vec%c = 0
    call queue_append_data(lat_vecs, vec, success)
    vis(0,0,0) = .true.
    
    do while (.not. queue_empty(lat_vecs))
        vec = queue_retrieve_data(lat_vecs)
        call build_pairs(vec, nhits, overflow)
        if (overflow) exit
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
    
    if (overflow .and. sys_master) then
        write(*,*) "Detected neighbor list overflow"
    end if
end subroutine qc_get_lists

subroutine qc_allocate_lists
    use qc_mpi
    type(QUEUE_STRUCT), pointer :: lat_vecs
    type(QUEUE_DATA) :: vec, nvec
    logical :: success
    integer :: nhits

    if (sys_master) write(*,'(A,/)') "(Re)allocating neighbor lists..."
    num_pairs = 0
    npair = 0
    nbq = 0

    if (allocated(pair_list)) deallocate(pair_list)
    if (allocated(bq_list)) deallocate(bq_list)
    if (allocated(bq2_list)) deallocate(bq2_list)
    if (allocated(d_bq0)) deallocate (d_bq0)
    if (allocated(d_mbqi)) deallocate (d_mbqi)
    if (allocated(d_mbqj)) deallocate (d_mbqj)
    if (allocated(d_mr)) deallocate (d_mr)
    
    !---BREADTH-FIRST-SEARCH: NEIGHBOR CELLS
    call queue_create( lat_vecs, 10000 )
    allocate ( vis(-100:100,-100:100,-100:100) )
    vis = .false.
    
    vec%a = 0; vec%b = 0; vec%c = 0
    call queue_append_data(lat_vecs, vec, success)
    vis(0,0,0) = .true.
    
    do while (.not. queue_empty(lat_vecs))
        vec = queue_retrieve_data(lat_vecs)
        call count_pairs(vec, nhits)
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

    maxnum_pair = maxval(npair) + 10
    maxnum_bq = maxval(nbq) + 10

    allocate (pair_list(maxnum_pair,nmol))
    allocate (bq_list(maxnum_bq,nmol))
    allocate (bq2_list(2*maxnum_bq))
    allocate (d_bq0(3,2*maxnum_bq*3))
    allocate (d_mbqi(3,2*maxnum_bq*3))
    allocate (d_mbqj(3,2*maxnum_bq*3))
    allocate (d_mr (3,natom))
end subroutine qc_allocate_lists

subroutine get_lr_interaction (u_lr, im, l_old_neigh)
    integer, intent(in) :: im
    real*8, intent(out) :: u_lr
    logical, intent(in) :: l_old_neigh
    type(QUEUE_STRUCT), pointer :: lat_vecs
    type(QUEUE_DATA) :: vec, nvec
    logical :: success
    integer :: nhits

    nlr(im) = 0
    u_lr = 0.0d0
    
    call queue_create( lat_vecs, 200000 )
    allocate ( vis(-100:100,-100:100,-100:100) )
    vis = .false.
    
    vec%a = 0; vec%b = 0; vec%c = 0
    call queue_append_data(lat_vecs, vec, success)
    vis(0,0,0) = .true.
    
    do while (.not. queue_empty(lat_vecs))
        vec = queue_retrieve_data(lat_vecs)
        call accumulate_lr(im, u_lr, vec, nhits, l_old_neigh)
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
end subroutine get_lr_interaction

subroutine count_pairs(vec, nhits)
    TYPE(QUEUE_DATA), intent(in) :: vec
    integer, intent(out) :: nhits
    integer :: im, ia
    integer :: jm, na, nb, nc, ja
    real*8 :: ri(3), rj(3), dcel(3), dr(3), r2

    na = vec%a; nb = vec%b; nc = vec%c
    dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)
    
    nhits = 0
    
    do im = 1, nmol
        ia = mol_iatom(im)
        ri = gm_pos(1:3,ia+1)
        do jm = 1, nmol
            if (im .eq. jm .and. na .eq. 0 .and. nb .eq. 0 .and. nc .eq. 0) cycle
            ja = mol_iatom(jm)
            rj = gm_pos(1:3,ja+1) + dcel
            dr = ri - rj
            r2 = dot_product(dr, dr)
            
            !-- within QM and BQ
            if (r2 < rcut_qq2) then
                nhits = nhits + 1
                nbq(im) = nbq(im) + 1
                
                if (jm .gt. im) then
                    npair(im) = npair(im) + 1
                    num_pairs = num_pairs + 1
                else if ( (jm.eq.im) .and. cell_gt_zero(na,nb,nc) ) then
                    npair(im) = npair(im) + 1
                    num_pairs = num_pairs + 1
                end if

            !-- within BQ, outside QM
            else if (r2 < rcut_bq2) then
                nhits = nhits + 1
                nbq(im) = nbq(im) + 1
            end if
        end do 
    end do
end subroutine count_pairs

subroutine build_pairs(vec, nhits, overflow)
    !use qc_mpi
    TYPE(QUEUE_DATA), intent(in) :: vec
    integer, intent(out) :: nhits
    logical, intent(out) :: overflow
    integer :: im, ia
    integer :: jm, na, nb, nc, ja
    real*8 :: ri(3), rj(3), dcel(3), dr(3), r2

    na = vec%a; nb = vec%b; nc = vec%c
    dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)
    
    nhits = 0
    overflow = .false.
    
    do im = 1, nmol
        ia = mol_iatom(im)
        ri = gm_pos(1:3,ia+1)
        do jm = 1, nmol
            if (im .eq. jm .and. na .eq. 0 .and. nb .eq. 0 .and. nc .eq. 0) cycle
            ja = mol_iatom(jm)
            rj = gm_pos(1:3,ja+1) + dcel
            dr = ri - rj
            r2 = dot_product(dr, dr)
            
            !-- within QM and BQ
            if (r2 < rcut_qq2) then
                nhits = nhits + 1
                nbq(im) = nbq(im) + 1
                if (nbq(im) > maxnum_bq) then
                    overflow = .true.
                    exit
                end if
                bq_list(nbq(im),im)%jm = jm
                bq_list(nbq(im),im)%ja = na
                bq_list(nbq(im),im)%jb = nb
                bq_list(nbq(im),im)%jc = nc
                bq_list(nbq(im),im)%qm = .true.
                
                if (jm .gt. im) then
                    npair(im) = npair(im) + 1
                    num_pairs = num_pairs + 1
                    !if (sys_master) write(*,*) im, 'has', npair(im)
                    if (npair(im) > maxnum_pair) then
                        overflow = .true.
                        exit
                    end if
                    pair_list(npair(im),im)%jm = jm
                    pair_list(npair(im),im)%ja = na
                    pair_list(npair(im),im)%jb = nb
                    pair_list(npair(im),im)%jc = nc
                    pair_list(npair(im),im)%scale = 1.0d0
                else if ( (jm.eq.im) .and. cell_gt_zero(na,nb,nc) ) then
                    npair(im) = npair(im) + 1
                    num_pairs = num_pairs + 1
                    !if (sys_master) write(*,*) im, 'has', npair(im)
                    if (npair(im) > maxnum_pair) then
                        overflow = .true.
                        exit
                    end if
                    pair_list(npair(im),im)%jm = jm
                    pair_list(npair(im),im)%ja = na
                    pair_list(npair(im),im)%jb = nb
                    pair_list(npair(im),im)%jc = nc
                    pair_list(npair(im),im)%scale = 1.0d0
                end if

            !-- within BQ, outside QM
            else if (r2 < rcut_bq2) then
                nhits = nhits + 1
                nbq(im) = nbq(im) + 1
                if (nbq(im) > maxnum_bq) then
                    overflow = .true.
                    exit
                end if
                bq_list(nbq(im),im)%jm = jm
                bq_list(nbq(im),im)%ja = na
                bq_list(nbq(im),im)%jb = nb
                bq_list(nbq(im),im)%jc = nc
                bq_list(nbq(im),im)%qm = .false.
            end if
        end do 
        if (overflow) exit
    end do
end subroutine build_pairs

subroutine qc_bq2_list_get(im, ip, l_old_neigh)
    integer, intent(in) :: im, ip
    logical, intent(in) :: l_old_neigh
    integer :: jm, ka, kb, na, nb, nc
    real*8 :: ri(3), rk(3), dr(3), r2, dcel(3)

    ! monomer im and its bq's (this will include monomers i and j)
    if (l_old_neigh) then
        ri = gm_pos_old(1:3,mol_iatom(im)+1)
    else
        ri = gm_pos(1:3,mol_iatom(im)+1)
    end if

    bq2_list(1)%jm = im
    bq2_list(1)%ja = 0
    bq2_list(1)%jb = 0
    bq2_list(1)%jc = 0
    nbq2 = 1

    do kb=1,nbq(im)
        nbq2 = nbq2 + 1
        bq2_list(nbq2)%jm = bq_list(kb,im)%jm
        bq2_list(nbq2)%ja = bq_list(kb,im)%ja
        bq2_list(nbq2)%jb = bq_list(kb,im)%jb
        bq2_list(nbq2)%jc = bq_list(kb,im)%jc
    end do

    ! bq's of jm not within bq cutoff of im (will NOT include monomers i or j)
    jm = pair_list(ip,im)%jm
    na = pair_list(ip,im)%ja
    nb = pair_list(ip,im)%jb
    nc = pair_list(ip,im)%jc
    if (l_old_neigh) then
        dcel(:) = na*lat_old(:,1)+nb*lat_old(:,2)+nc*lat_old(:,3)
    else
        dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)
    end if

    do kb=1,nbq(jm)
        ka = mol_iatom(bq_list(kb,jm)%jm)
        if (l_old_neigh) then
            rk = gm_pos_old(1:3,ka+1) + bq_list(kb,jm)%ja*lat_old(:,1) + bq_list(kb,jm)%jb*lat_old(:,2) &
                                      + bq_list(kb,jm)%jc*lat_old(:,3) 
        else
            rk = gm_pos(1:3,ka+1) + bq_list(kb,jm)%ja*lat(:,1) + bq_list(kb,jm)%jb*lat(:,2) & 
                                  + bq_list(kb,jm)%jc*lat(:,3) 
        end if
        rk = rk + dcel
        dr = rk - ri
        r2 = dot_product(dr, dr)
        if (r2 < rcut_bq2) then 
            cycle !-- already got from im
        else
            nbq2 = nbq2 + 1
            bq2_list(nbq2)%jm = bq_list(kb,jm)%jm
            bq2_list(nbq2)%ja = bq_list(kb,jm)%ja + na
            bq2_list(nbq2)%jb = bq_list(kb,jm)%jb + nb
            bq2_list(nbq2)%jc = bq_list(kb,jm)%jc + nc
        end if
    end do
end subroutine qc_bq2_list_get

subroutine accumulate_lr(im, u_lr, vec, nhits, l_old_neigh)
    integer, intent(in) :: im
    real*8, intent(inout) :: u_lr
    TYPE(QUEUE_DATA), intent(in) :: vec
    integer, intent(out) :: nhits
    logical, intent(in) :: l_old_neigh
    logical :: completeLRcell, outOfRange
    integer :: ia, is, i, j, nchg_i, ichg_i, nchg_j, ichg_j
    integer :: jm, na, nb, nc, ja, js, at_i, at_j
    real*8 :: ri(3), rj(3), qi, qj, dcel(3), rij(3), rij2, rij_norm
    real*8 :: d_ri(3)

    completeLRcell = .false.
    outOfRange = .false.
    na = vec%a; nb = vec%b; nc = vec%c
    dcel(:) = na*lat(:,1)+nb*lat(:,2)+nc*lat(:,3)

    nhits = 0
    ia = mol_iatom(im)
    nchg_i = mol_nchg(im)
    ichg_i = mol_ichg(im)

    do jm = 1, nmol
        if (im .eq. jm .and. na .eq. 0 .and. nb .eq. 0 .and. nc .eq. 0) cycle
        
        ja = mol_iatom(jm)
        nchg_j = mol_nchg(jm)
        ichg_j = mol_ichg(jm)
        
        if (l_old_neigh) then
            ri = gm_pos_old(1:3,ia+1)
            rj = gm_pos_old(1:3,ja+1) + na*lat_old(:,1)+nb*lat_old(:,2)+nc*lat_old(:,3)
        else
            ri = gm_pos(1:3,ia+1)
            rj = gm_pos(1:3,ja+1) + dcel
        end if
        
        rij = rj - ri
        rij2 = dot_product(rij, rij)
        
        if (rij2 < rcut_qq2) then
            nhits = nhits + 1
        else if (rij2 < rcut_bq2) then
            nhits = nhits + 1
        else if (rij2 < rcut_lr2) then
            nhits = nhits + 1
            nlr(im) = nlr(im) + 1
            completeLRcell = .true.
            do is=1, nchg_i
                ri = chg_pos(1:3,ichg_i+is)
                qi = chg_pos(0,  ichg_i+is)
                do js = 1, nchg_j
                    rj = chg_pos(1:3,ichg_j+js) + dcel
                    qj = chg_pos(0,  ichg_j+js)
                    rij = (rj - ri)*ang2bohr ! points from i to j
                    rij_norm = sqrt(dot_product(rij,rij))
                    u_lr = u_lr + 0.5d0*qi*qj/rij_norm
                    if (l_grd) then

                        d_ri = 0.5d0 * rij * (qi*qj/(rij_norm**3))
                        
                        do at_i = 1, mol_nsite(im)
                            d_mr(1:3,ia+at_i) = d_mr(1:3,ia+at_i) + d_ri(1:3)*chg_coef(at_i, ichg_i+is) ! why are we doing this???
                        end do
                        do at_j = 1, mol_nsite(jm)
                            d_mr(1:3,ja+at_j) = d_mr(1:3,ja+at_j) - d_ri(1:3)*chg_coef(at_j, ichg_j+js)
                        end do

                        !--LR Contribution to Virial Tensor
                        do j = 1, 3
                        do i = 1, 3
                            m_virt(i,j) = m_virt(i,j) - gm_pos(i,ia+1)*ang2bohr*d_ri(j) ! force on i is -d_ri
                            m_virt(i,j) = m_virt(i,j) + (gm_pos(i,ja+1)+dcel(i))*ang2bohr*d_ri(j) ! force on k is +d_ri
                        end do
                        end do
                    end if
                end do
            end do
        else
            outOfRange = .true.
        end if
    end do 

    if (completeLRcell .and. outOfRange) then
        do jm = 1, nmol
            if (im .eq. jm .and. na .eq. 0 .and. nb .eq. 0 .and. nc .eq. 0) cycle
            
            ja = mol_iatom(jm)
            nchg_j = mol_nchg(jm)
            ichg_j = mol_ichg(jm)
            
            if (l_old_neigh) then
                ri = gm_pos_old(1:3,ia+1)
                rj = gm_pos_old(1:3,ja+1) + na*lat_old(:,1)+nb*lat_old(:,2)+nc*lat_old(:,3)
            else
                ri = gm_pos(1:3,ia+1)
                rj = gm_pos(1:3,ja+1) + dcel
            end if
            
            rij = rj - ri
            rij2 = dot_product(rij, rij)
            
            if (rij2 >= rcut_lr2) then
                nhits = nhits + 1
                nlr(im) = nlr(im) + 1
                do is=1, nchg_i
                    ri = chg_pos(1:3,ichg_i+is)
                    qi = chg_pos(0,  ichg_i+is)
                    do js = 1, nchg_j
                        rj = chg_pos(1:3,ichg_j+js) + dcel
                        qj = chg_pos(0,  ichg_j+js)
                        rij = (rj - ri)*ang2bohr ! points from i to j
                        rij_norm = sqrt(dot_product(rij,rij))
                        u_lr = u_lr + 0.5d0*qi*qj/rij_norm
                        if (l_grd) then

                            d_ri = 0.5d0 * rij * (qi*qj/(rij_norm**3))
                            
                            do at_i = 1, mol_nsite(im)
                                d_mr(1:3,ia+at_i) = d_mr(1:3,ia+at_i) + d_ri(1:3)*chg_coef(at_i, ichg_i+is)
                            end do
                            do at_j = 1, mol_nsite(jm)
                                d_mr(1:3,ja+at_j) = d_mr(1:3,ja+at_j) - d_ri(1:3)*chg_coef(at_j, ichg_j+js)
                            end do

                            !--LR Contribution to Virial Tensor
                            do j = 1, 3
                            do i = 1, 3
                                m_virt(i,j) = m_virt(i,j) - gm_pos(i,ia+1)*ang2bohr*d_ri(j) ! force on i is -d_ri
                                m_virt(i,j) = m_virt(i,j) + (gm_pos(i,ja+1)+dcel(i))*ang2bohr*d_ri(j) ! force on k is +d_ri
                            end do
                            end do
                        end if
                    end do
                end do
            end if
        end do 
    end if
end subroutine accumulate_lr

function cell_gt_zero (na, nb, nc)
    integer, intent(in) :: na, nb, nc
    logical :: cell_gt_zero
    cell_gt_zero = (na.gt.0).or.(na.eq.0 .and. nb.gt.0).or.(na.eq.0 .and. nb.eq.0 .and. nc.gt.0)
end function cell_gt_zero
end module qc_neigh
