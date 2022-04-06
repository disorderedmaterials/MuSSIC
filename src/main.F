      Program scattering_cg

!    -computes the neutron scattering from an atomistic trajectory.
!    -reads the input files 'input.dat', 'model.dat' and the trajectory file 'trajectory' 
!    -writes the output to the files named as provided in  the input file
!
!     input files:  - input.dat  (main input file containing data like atom types,
!     box dimensions, isotopologues, scattering lengths and binwidth and
!     binrange)
!                   - model.dat (model file with system description like
!                    number of molecules, sequence of atoms in each
!                    molecule ) 
!                   - trajectory (trajectory in xyz format)
!     Output files: 
!                    ***-unweighted.gr   (unweighted partial pair distribuion functions)
!                    ***-unweighted.sq   (unweighted partial pair structure factor)
!                    unweighted-total.sq (total unweighted structure factor) 
!                    ***-weighted.sq     (weighted partial pair strucutre factor)
!                    ###-weighted-total.sq   (total weighted structure factor)
!                    ***-weighted.gr     (weighted partial pair distribution function) 
!                    weighted-total.gr   (total weighted radial distribution function)
!
!                 (*** indicate all pairs of atomtypes in the system description)
!                 (### indicate the isotopologue given in the input.dat)
!                               

                     
      implicit none
      include 'system.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'
      include 'omp_lib.h'


      double precision elapsetime1,elapsetime2
      integer*8 istep,id

cc$omp parallel
cc$omp& private ( id )
cc
cc  Have each thread say hello.
cc
c      id = omp_get_thread_num ( )
c      write ( *, '(a,i8)' ) '  Hello from process ', id

cc$omp end parallel

c      stop

      istop=0

      t1 = 0.00
      t2 = 0.00
      t3 = 0.00

      call CPU_TIME(t1)

c  opens input and output files
      call ofiles

c     stops if one of the files needed does not exist

      if (istop.gt.0) then
          stop
      endif
c  reads data from the input file
      call getdata

      if (istop.gt.0) then
          stop
      endif
        
               
c  reads data from the model file      

      call rdmodel

        
c     stops if one of the input lines is incomprehensible
      if (istop.gt.0) then
          stop
      endif
c initialize values and vectors needed for the calculation
      call initialize

             
      call CPU_TIME(t2)

c  Main loop over number of frames start here

      write(io7,*)'Scattering calculation for atomistic system
     & (06-04-2022) '

     
      do istep = 1,nframes

c read configuration from the trajectory file      
      call rdtraj

c     stops if the trajectory file is unavailable/incorrect
      if (istop.gt.0) then
          stop
      endif

c   populate the histogram bin for radial distribution function

      call rdf

      enddo

c   calculates the partial pair radial distribution 

      call partial_rdf

c   calculates the partial pair structure factors 
      
      call sq

c   calculates the total structure factor and G(r) for isotpologue mixtures     

      call neutronsq

c      call finalize

      call CPU_TIME(t3)

        elapsetime1 = t2-t1
        elapsetime2 = t3-t2

      write(io7,50)  elapsetime1, elapsetime1/3600
      write(io7,100) elapsetime2, elapsetime2/3600

 50   format(5x,'initial set up time',g15.6,'s',g15.6,'h')
 100  format(5x,'running time',g15.6,'s',g15.6,'h')

      stop
      end



