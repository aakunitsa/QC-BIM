module qc_mpi
  use mpi
  implicit none
  
  integer :: sys_nproc
  integer :: sys_me
  logical :: sys_master

contains

  subroutine qc_mpi_start
    integer :: ierr
    integer :: comm = MPI_COMM_WORLD  ! global communicator

    sys_me = 0
    sys_nproc = 1
    
    call mpi_init (ierr)
    call mpi_comm_rank (comm, sys_me, ierr)
    call mpi_comm_size (comm, sys_nproc, ierr)
    
    sys_master = sys_me .eq. 0

    !if (sys_master .and. sys_nproc .gt. 1) then
    !   write (6, '(2x,a,10x,i5,a)') 'mpi is running with ', sys_nproc, ' cpus'
    !end if

  end subroutine qc_mpi_start


  subroutine qc_mpi_end
    integer :: ierr

    call mpi_finalize (ierr)

  end subroutine qc_mpi_end


  subroutine qc_mpi_barrier
    integer :: ierr
    
    call mpi_barrier (MPI_COMM_WORLD, ierr)

  end subroutine qc_mpi_barrier


  subroutine qc_mpi_allreduce1c (val)
    complex*8, intent(inout) :: val
    complex*8 :: sum_val
    integer  :: ierr
    
    call mpi_allreduce (val, sum_val, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
         MPI_COMM_WORLD, ierr)

    val = sum_val

  end subroutine qc_mpi_allreduce1c


  subroutine qc_mpi_allreduce1 (val)
    real*8, intent(inout) :: val
    real*8 :: sum_val
    integer  :: ierr
    
    call mpi_allreduce (val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)

    val = sum_val

  end subroutine qc_mpi_allreduce1




  subroutine qc_mpi_allreduce1i (val)
    integer, intent(inout) :: val
    integer :: sum_val
    integer  :: ierr

    call mpi_allreduce (val, sum_val, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD, ierr)

    val = sum_val

  end subroutine qc_mpi_allreduce1i



  subroutine qc_mpi_allreduce (val, ndim)
    integer, intent(in)       :: ndim
    real*8, dimension(ndim), intent(inout) :: val
    real*8, dimension(:), allocatable  :: sum_val
    integer :: ierr
    
    allocate (sum_val(ndim))

    sum_val = val

    call mpi_allreduce (val, sum_val, ndim, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)

    val = sum_val

    deallocate (sum_val)

  end subroutine qc_mpi_allreduce
  
  subroutine qc_mpi_allreducei (val, ndim)
    integer, intent(in)       :: ndim
    integer, dimension(ndim), intent(inout) :: val
    integer, dimension(:), allocatable  :: sum_val
    integer :: ierr
    
    allocate (sum_val(ndim))

    sum_val = val

    call mpi_allreduce (val, sum_val, ndim, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD, ierr)

    val = sum_val

    deallocate (sum_val)

  end subroutine qc_mpi_allreducei

  
  subroutine qc_mpi_reduce1c (val)
    complex*8, intent(inout) :: val
    complex*8 :: sum_val
    integer  :: ierr
    
    call mpi_reduce (val, sum_val, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
         0, MPI_COMM_WORLD, ierr)

    val = sum_val

  end subroutine qc_mpi_reduce1c


  subroutine qc_mpi_reduce1 (val)
    real*8, intent(inout) :: val
    real*8 :: sum_val
    integer  :: ierr
    
    call mpi_reduce (val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         0, MPI_COMM_WORLD, ierr)

    val = sum_val

  end subroutine qc_mpi_reduce1



  subroutine qc_mpi_reduce (val, ndim)
    integer, intent(in)       :: ndim
    real*8, dimension(ndim), intent(inout) :: val
    real*8, dimension(:), allocatable  :: sum_val
    integer :: ierr
    
    allocate (sum_val(ndim))

    sum_val = val

    call mpi_reduce (val, sum_val, ndim, MPI_DOUBLE_PRECISION, MPI_SUM, &
         0, MPI_COMM_WORLD, ierr)

    val = sum_val

    deallocate (sum_val)

  end subroutine qc_mpi_reduce




  subroutine qc_mpi_real_bcast (buffer, ndim)
    integer, intent(in)    :: ndim
    real*8, dimension(ndim), intent(inout) :: buffer
    integer :: ierr, root

    root = 0 ! Master process; see mpi_init for details
    call mpi_bcast (buffer, ndim, MPI_DOUBLE_PRECISION, &
         root, MPI_COMM_WORLD, ierr)

  end subroutine qc_mpi_real_bcast



  subroutine qc_mpi_int_bcast (buffer, ndim)
    integer, intent(in)   :: ndim
    integer, dimension(ndim), intent(inout) :: buffer
    integer :: ierr, root

    root = 0
    call mpi_bcast (buffer, ndim, MPI_INTEGER, &
         root, MPI_COMM_WORLD, ierr)

  end subroutine qc_mpi_int_bcast



  subroutine qc_mpi_real_bcast1 (buffer) 
    real*8, intent(inout) :: buffer
    integer  :: ierr, root, ndim

    root = 0
    ndim = 1
    call mpi_bcast (buffer, ndim, MPI_DOUBLE_PRECISION, &
         root, MPI_COMM_WORLD, ierr)

  end subroutine qc_mpi_real_bcast1


  
  subroutine qc_mpi_int_bcast1 (buffer) 
    integer, intent(inout) :: buffer
    integer :: ierr, root, ndim

    root = 0
    ndim = 1
    call mpi_bcast (buffer, ndim, MPI_INTEGER, &
         root, MPI_COMM_WORLD, ierr)

  end subroutine qc_mpi_int_bcast1

    
  subroutine qc_mpi_logical_bcast1 (buffer) 
    logical, intent(inout) :: buffer
    integer :: ierr, root, ndim

    root = 0
    ndim = 1
    call mpi_bcast (buffer, ndim, MPI_LOGICAL, &
         root, MPI_COMM_WORLD, ierr)

  end subroutine qc_mpi_logical_bcast1

  subroutine qc_mpi_char_bcast (buffer) 
    character(*), intent(inout) :: buffer
    integer :: ierr, root, ndim

    root = 0
    ndim = 1
    call mpi_bcast (buffer, len(buffer), MPI_CHARACTER, &
         root, MPI_COMM_WORLD, ierr)

  end subroutine qc_mpi_char_bcast

end module qc_mpi
