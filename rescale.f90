subroutine rescale
use qc_system
use qc_lattice
use qc_mpi
integer :: ii, im, is, ia
real*8, dimension(:,:), allocatable :: gm_pos0

allocate(gm_pos0(3,natom))

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

! change lattice parameters
lat_a = rescale_param(1)
lat_b = rescale_param(2)
lat_c = rescale_param(3)
lat_alpha = rescale_param(4)
lat_beta = rescale_param(5)
lat_gamma = rescale_param(6)

call update_geometry(gm_pos0)

! Write if invoked outside of optimization
if (sys_master .and. (.not. l_opt)) then
    write(*, '(i0)') natom
    write (*, '(f0.8, 5f14.8, i3)') lat_a, lat_b, lat_c, lat_alpha, &
    lat_beta, lat_gamma, lat_axis
    do im = 1,nmol
        ia = mol_iatom(im)
        do is=1,mol_nsite(im)
            write (*, '(A,3f16.8)') at_atnm(ia+is), gm_pos(1:3, ia+is)
        end do
    end do
end if

deallocate(gm_pos0)
end subroutine rescale
