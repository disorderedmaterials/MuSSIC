!******************************************************************
!     parallel.h
!*****************************************************************
!  These variable that are needed in MPI.
!******************************************************************
      integer*8 ierr,myid,numprocs
      integer*8 natoms_proc(20),mol_proc(3,20),nmol_proc(20)
      integer*8 histogram_proc(6,6,5000),bound_hist_proc(6,6,5000),
     &unbound_hist_proc(6,6,5000),io400

      common/mpi/ierr,myid,numprocs,natoms_proc,mol_proc
      common/mpi/nmol_proc,histogram_proc,bound_hist_proc,
     &unbound_hist_proc
      common/mpiio/io400

