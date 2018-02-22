module qc_lattice
use qc_system
use consts
implicit none
contains

function cross(a, b)
    real*8, dimension(3) :: cross
    real*8, dimension(3), intent(in) :: a, b
    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross

function volume(lat)
    real*8 :: volume
    real*8, dimension(3,3), intent(in) :: lat
    volume = abs(dot_product( lat(:,1), cross(lat(:,2), lat(:,3)) ))
end function volume

subroutine create_lattice_vectors
    real*8 :: ab, ac, bc, cx, cy, cz
    
    a_dim = abs(lat_a) > 1.0d-5
    b_dim = abs(lat_b) > 1.0d-5
    c_dim = abs(lat_c) > 1.0d-5

    ab = lat_gamma / RADIAN
    ac = lat_beta / RADIAN
    bc = lat_alpha / RADIAN
    
    if (lat_axis .eq. 0) then
        lat(:, 1) = lat_a * (/     1.0d0,     0.0d0, 0.0d0 /)
        lat(:, 2) = lat_b * (/ dcos(ab), dsin(ab), 0.0d0 /)
        cx = dcos(ac)
        cy = ( dcos(bc) - dcos(ac)*dcos(ab) ) / dsin(ab)
        cz = dsqrt( 1.0d0 - cx*cx - cy*cy )
    else if (lat_axis .eq. 1) then
        lat(:, 1) = lat_a * (/ 0.0d0,     1.0d0,     0.0d0 /)
        lat(:, 2) = lat_b * (/ 0.0d0, dcos(ab), dsin(ab) /)
        cy = dcos(ac)
        cz = ( dcos(bc) - dcos(ab)*dcos(ac) ) / dsin(ab)
        cx = dsqrt( 1.0d0 - cy*cy - cz*cz )
    else if (lat_axis .eq. 2) then
        lat(:, 1) = lat_a * (/ 0.0d0,     0.0d0,     1.0d0 /)
        lat(:, 2) = lat_b * (/ 0.0d0, dsin(ab), dcos(ab) /)
        cz = dcos(ac)
        cy = ( dcos(bc) - dcos(ab)*dcos(ac) ) / dsin(ab)
        cx = dsqrt( 1.0d0 - cy*cy - cz*cz )
    else
        STOP "Invalid choice for a axis"
    end if
    lat(:, 3) = lat_c * (/ cx, cy, cz /)
end subroutine create_lattice_vectors 

subroutine update_scaling_matrix
    real*8 :: lat0(3,3), lat2d(2,2), lat2d0(2,2), lat2di(2,2)
    
    ! 3D
    if (a_dim .and. b_dim .and. c_dim) then
        lat0 = lat
        call inverse(lat0, lati, 3)
    ! 2D: c = 0; a,b vecs in xy-plane
    else if (a_dim .and. b_dim .and. (.not. c_dim) .and. lat_axis.eq.0) then
        lat2d(1:2,1) = lat(1:2,1)
        lat2d(1:2,2) = lat(1:2,2)
        lat2d0 = lat2d
        call inverse(lat2d0, lat2di, 2)
        lati(1:2,1) = lat2di(1:2,1)
        lati(1:2,2) = lat2di(1:2,2)
        lati(:,3) = 0.0d0
        lati(3,:) = 0.0d0
        lati(3,3) = 1.0d0
        lat(:,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
    ! 1D: b=c=0; a vec along x-axis
    else if (a_dim .and. (.not. b_dim) .and. (.not. c_dim) .and. lat_axis.eq.0) then
        lati(:,1) = (/ 1.0d0/lat_a, 0.0d0, 0.0d0 /)
        lati(:,2) = (/ 0.0d0, 1.0d0, 0.0d0 /)
        lat(:,2) = (/ 0.0d0, 1.0d0, 0.0d0 /)
        lati(:,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
        lat(:,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
    ! 0D: no PBC
    else if ( .not. (a_dim .or. b_dim .or. c_dim) ) then
        lat(:,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
        lati(:,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
        lat(:,2) = (/ 0.0d0, 1.0d0, 0.0d0 /)
        lati(:,2) = (/ 0.0d0, 1.0d0, 0.0d0 /)
        lat(:,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
        lati(:,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
    else
        STOP "INVALID 1D or 2D LATTICE VECTOR FORMAT"
    end if
end subroutine update_scaling_matrix

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse


subroutine calc_lattice_parameters
    lat_a = sqrt(dot_product(lat(:,1),lat(:,1)))
    lat_b = sqrt(dot_product(lat(:,2),lat(:,2)))
    lat_c = sqrt(dot_product(lat(:,3),lat(:,3)))
    lat_gamma = dacos( dot_product(lat(:,1),lat(:,2))/(lat_a*lat_b) ) * RADIAN
    lat_beta = dacos( dot_product(lat(:,1),lat(:,3))/(lat_a*lat_c) ) * RADIAN
    lat_alpha = dacos( dot_product(lat(:,2),lat(:,3))/(lat_b*lat_c) ) * RADIAN
end subroutine calc_lattice_parameters

subroutine lat_angle_grad(d_al, d_be, d_ga)
    real*8, dimension(3,3), intent(out) :: d_al, d_be, d_ga
    real*8, parameter :: delta = 1.0d-4
    real*8, dimension(3,3) :: lat0
    real*8 :: angle0
    
    lat0 = lat
    
    !dlat/dalpha (units Angstrom/degree)
    angle0 = lat_alpha

    lat_alpha = lat_alpha + delta
    call create_lattice_vectors
    d_al = (lat - lat0) / delta

    lat_alpha = angle0
    lat = lat0
    !dlat/dbeta (units Angstrom/degree)
    angle0 = lat_beta

    lat_beta = lat_beta + delta
    call create_lattice_vectors
    d_be = (lat - lat0) / delta

    lat_beta = angle0
    lat = lat0
    
    !dlat/dgamma (units Angstrom/degree)
    angle0 = lat_gamma

    lat_gamma = lat_gamma + delta
    call create_lattice_vectors
    d_ga = (lat - lat0) / delta

    lat_gamma = angle0
    lat = lat0
end subroutine lat_angle_grad

end module qc_lattice
