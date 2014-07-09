!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : mmpi
!!! source  : m_mpi.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 08/09/2006 by li huang
!!!           02/27/2010 by li huang
!!!           07/09/2014 by li huang
!!! purpose : define my own mpi calls, inspired by famous quantum espresso
!!!           code. we note that the original mpi interfaces/subroutines 
!!!           are rather complicated for newbies, thus we try to wrap the
!!!           most important and useful (not all) mpi subroutines (in our
!!!           opinion) in this module to facilite the usage of mpi. these
!!!           subroutines are
!!!               MPI_INIT(),
!!!               MPI_FINALIZE(),
!!!               MPI_WTIME(),
!!!               MPI_WTICK(),
!!!               MPI_BARRIER(),
!!!               MPI_DIMS_CREATE(),
!!!               MPI_CART_CREATE(),
!!!               MPI_CART_COORDS(),
!!!               MPI_COMM_SPLIT(),
!!!               MPI_COMM_RANK(),
!!!               MPI_COMM_SIZE(),
!!!               MPI_GET_PROCESSOR_NAME(),
!!!               MPI_BCAST(),
!!!               MPI_GATHER(),
!!!               MPI_GATHERV(),
!!!               MPI_ALLGATHER(),
!!!               MPI_ALLGATHERV(),
!!!               MPI_REDUCE(),
!!!               MPI_ALLREDUCE(),
!!!           etc. in the module, we also try to implement a light-weight
!!!           error handler. enjoy it!
!!! status  : unstable
!!! comment : this module has been tested under the following environment:
!!!               mpich1    1.2.7p1
!!!               mpich2    1.2.1p1
!!!               mvapich2  1.2.0p1
!!!               intel mpi 3.2.0
!!!-----------------------------------------------------------------------

  module list
  end module list
