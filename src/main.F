c      SPDX-License-Identifier: GPL-3.0-or-later

      Program scattering_cg

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!    DESCRIPTION:
!    computes the neutron scattering either from a coarse-grain
!    trajectory or from an atomistic trajectory which is then compared
!    against CG models
!              
!    INPUT FILES: 
!        'input.dat'  (main input file containing data like atom types,
!                         box dimensions, isotopologues, scattering lengths, 
!                         binwidth and binrange)
! 
!        'model.dat'  (model file with system description like
!                         number of molecules, sequence of atoms in each molecule )

!        'trajectory' or  (atomistic trajectory in xyz format)
!        'CG_trajectory'  (CG-trajectory in xyz format)
              

!
!    OUTPUT FILES:              
!    -writes the output to the files named as provided in  the input file
!
!                   
!                    ***-CG-unweighted.gr   (unweighted partial pair distribuion functions)
!                    ***-CG-unweighted.sq   (unweighted partial pair structure factor)
!                    CG-unweighted-total.sq (total unweighted structure factor) 
!                    ***-CG-weighted.sq     (weighted partial pair strucutre factor)
!                    ###-CG-weighted-total.sq   (total weighted structure factor)
!                    ***-CG-weighted.gr     (weighted partial pair distribution function) 
!                    CG-weighted-total.gr   (total weighted radial distribution function)
!
!                 (*** indicate all pairs of atomtypes in the system description)
!                 (### indicate the isotopologue given in the input.dat)
!                 filename without CG indicate the same files from atomistic scattering
                               
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                               
                     
      implicit none
      include 'system.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'
      include 'omp_lib.h'


      double precision elapsetime1,elapsetime2
      integer*8 istep,id,i,fu_1

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*)'Scattering calculation: version 1.0 (26-01-2023) '
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
cccccccccccccccccccccccccccccccccccccccccccccc
c     check on  available parallel threads
cccccccccccccccccccccccccccccccccccccccccccccc

c$omp parallel
c$omp& private ( id )
c
c  Have each thread say hello.
c
c        id = omp_get_thread_num ( )
c       write ( *, '(a,i8)' ) '  thread num ', id

c$omp end parallel

ccccccccccccccccccccccccccccccccccccccccccccccccc


      istop=0

      t1 = 0.00
      t2 = 0.00
      t3 = 0.00

      call CPU_TIME(t1)
ccccccccccccccccccccccccccccccccc
c  opens input and output files
ccccccccccccccccccccccccccccccccc

      call ofiles

ccccccccccccccccccccccccccccccccccccccccccccccccc
c stops if one of the files needed does not exist
ccccccccccccccccccccccccccccccccccccccccccccccccc      

      if (istop.gt.0) then
          stop
      endif

ccccccccccccccccccccccccccccccccccccccccccccc      
c  reads data from the input file 'input.dat'
cccccccccccccccccccccccccccccccccccccccccccc

      call getdata

      write(*,*)'read input.dat file successfully'

cccccccccccccccccccccccccccccccccccccccccccccc
c stops if there is a missing/wrong input data 
cccccccccccccccccccccccccccccccccccccccccccccc     

      if (istop.gt.0) then
          stop
      endif
        
ccccccccccccccccccccccccccccccccccccccccccccc               
c  reads data from the model file 'model.dat'
ccccccccccccccccccccccccccccccccccccccccccccc      

      call rdmodel



      write(*,*)'read model.dat file successfully'

ccccccccccccccccccccccccccccccccccccccccccccccccc      
c stops if there is any missing/wrong model data
ccccccccccccccccccccccccccccccccccccccccccccccccc


      if (istop.gt.0) then
          stop
      endif

cccccccccccccccccccccccc
c note the cpu time       
      call CPU_TIME(t2)
cccccccccccccccccccccccc      


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  check whether the run type is 'default' or 'CG_compare'
c
c   if  'default' computes  F_CG(Q) for the given CG_trajectory  (refer documentation for details)
c
c   if  'compare_CG' compute F(Q) for refernce atomisic simulation and
c   compare with F_CG(Q) by mapping atoms in the atomistic trajectory
c   onto a CG bead   (refer documentation for deatils)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      if(run_type.eq.'compare_CG')then


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Scattering calculation from atomistic simulations starts here
      write(io7,*)'Computing F(Q) using the atomistic trajectory '
      write(*,*)'******started F(Q) calculation using atomistic trajecto
     &ry*******'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialize number of  histogram bins and vectors needed for F(Q) calculation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       call initialize
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   if isotope susbstitution is manually done, write the trajectory
c   and reinitialize the arrays
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if(isubstitute.eq.1) then
         call wrtraj
          
       
         call reinitialize
        endif

       
       
      write(*,*)'Initialization of vectors for F(Q) completed'

      write(*,*)'Loop begins over number of frames to populate 
     &histogram for g(r) calculation'


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Loop over number of frames in atomisitc trajectory start here
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do istep = 1,nframes


     
      write(*,*)'frame:',istep


ccccccccccccccccccccccccccccccccccccccccccccc
c read configuration from the trajectory file
ccccccccccccccccccccccccccccccccccccccccccccc

      call rdtraj

      

     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     stops if the trajectory file is unavailable/incorrect
      if (istop.gt.0) then
          stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccc
c   populate the histogram bin for rdf calculation
cccccccccccccccccccccccccccccccccccccccccccccccccc

      call rdf

      enddo

             
      write(*,*)'Loop ends over number of frames; histogram computed'

cccccccccccccccccccccccccccccccccccccccc
c  Loop  over number of frames ends here 
cccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c   calculates the partial pair correlation for all possible pairs of atom types
c   using the histogram bins             

c             g_{ij}(r)   where i and j are the indices of atom types
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      call partial_rdf

       
      write(*,*)'partial g(r) computed'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   calculates the partial pair structure factors by performing the
c   Fourier transfrom on g_ij(r)

c    S_{ij} (Q)  where i and j are the indices of atom types       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sq


      write(*,*)'partial S(Q) computed'
      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   calculates the total structure factor F(Q) for different isotopologues given in input.dat

c   For a given isotpologue,  F(Q) = sum of  F_{ij}(Q)  where i and j are the indices of atom types 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call neutronsq


      write(*,*)'Weighted S(Q) computed'
 
        fu_1 = 120+num_iso       
 
       close(io2)

      if(isubstitute.eq.1) close(fu_1)
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    perform the CG mapping following the mapping scheme given in model.dat

c   - generates CG trajectory by grouping atoms in the atomistic trajectory into a CG bead 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        call cg_mapping

       open (unit = io3, file = 'CG_trajectory',  status = 'UNKNOWN')


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Scattering calculation from atomistic simulations ends  here
      write(io7,*)'F(Q) calculation from atomistic trajectory completed'
      write(*,*)'F(Q) calculation from atomistic trajectory completed'


      write(*,*)'Generated CG trajectory following the mapping scheme
     &given in model file'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
           
            endif     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Scattering calulcation on the coarse grain trajectory obtained from CG_mapping
      write(io7,*)'*****Started F_CG(Q) calculation using CG trajectory'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

        

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialize values and vectors needed for F_CG(Q)  calculation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      call initialize_cg


        if(isubstitute.eq.1.and.run_type.eq.'default') then
         call wrtraj_cg                 
         call reinitialize_cg
        endif

      
      write(*,*)'Initialization of vectors for F(Q) completed'
      write(*,*)'Loop begins over number of frames to populate
     & histogram for g(r) calculation'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Loop over number of frames in the CG trajectory starts here
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do istep = 1,nframes

      write(*,*)'frame:',istep
ccccccccccccccccccccccccccccccccccccccccccccc
c read configuration from the trajectory file
ccccccccccccccccccccccccccccccccccccccccccccc      

      call rdtraj_cg

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     stops if the trajectory file is unavailable/incorrect
      if (istop.gt.0) then
          stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   populate the histogram bin for radial distribution function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call rdf_cg

  
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  loop over number of frames in CG trajectory ends here
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
     

      write(*,*)'Loop ends over number of frames; histogram computed'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c   calculates the partial pair correaltion function for all possible
c   pairs of bead types


c             g_{ij}(r)   where i and j are the indices of bead types
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call partial_rdf_cg

      
      write(*,*)'partial g(r) computed'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   calculates the bead radius from the CG trajectory

c     - computes the distance between then consective beads in a chain and
c        average over the number of frames to get the radius of the bead 
c       (refer to documention for more details)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       call bead_size


      write(*,*)'bead radius computed'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c   calculates the partial pair structure factors by performing the
c   Fourier transfrom on g_ij(r)

c    S_{CG,ij} (Q)  where i and j are the indices of bead types       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       

      call sq_cg


      write(*,*)'partial S(Q) computed'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   calculates the total structure factor F(Q) for different isotopologues given in input.dat

c   For a given isotpologue,  F_CG(Q) = sum of  F_{CG,ij}(Q)  where i and j are the indices of bead types 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call neutronsq_cg


      write(*,*)'Weighted S(Q) computed'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Scattering calculation on CG trajectory  ends  here
      write(io7,*)'calculation of F_CG(Q) on CG trajectory completed'
      write(io7,*)'*****calculation of F_CG(Q) on CG trajectory 
     &completed*******'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  note the cpu time and compute the time taken to run the whole calculation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call CPU_TIME(t3)

        elapsetime1 = t2-t1
        elapsetime2 = t3-t2

      write(io7,50)  elapsetime1/nproc, elapsetime1/(3600*nproc)
      write(io7,100) elapsetime2/nproc, elapsetime2/(3600*nproc)

      write(*,50)  elapsetime1/nproc, elapsetime1/(3600*nproc)
      write(*,100) elapsetime2/nproc, elapsetime2/(3600*nproc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

 50   format(5x,'initial set up time',g15.6,'seconds',g15.6,'hrs')
 100  format(5x,'running time',g15.6,'seconds',g15.6,'hrs')

      stop
      end



