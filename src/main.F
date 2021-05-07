      Program scattering_calc
      implicit none
      include 'system.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'



      double precision elapsetime1,elapsetime2
      integer istep

      istop=0

      t1 = 0.00
      t2 = 0.00
      t3 = 0.00

      call CPU_TIME(t1)

c  open input and out files
      call ofiles

c     stops if one of the files needed does not exist

      if (istop.gt.0) then
          stop
      endif
c  read data from the input file
      call getdata

      if (istop.gt.0) then
          stop
      endif
        
        
c  read data from the model file      

      call rdmodel

        
c     stops if one of the input lines is incomprehensible
      if (istop.gt.0) then
          stop
      endif
c initialize values and vectors needed for the calculation
      call initialize

               

      call CPU_TIME(t2)

      write(io7,*)'Scattering_calculation_CG version: 11-03-2021 '

     

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

c   calculate the partial pair radial distribution 

      call partial_rdf

c   calculate the partial pair structure factors 
      
      call sq

c   calculate the weighted sum of the structure factors and G(r)       

      call neutronsq

c   calculate the total structure factor for  isotpologue mixtures     

      call neutronsq_deuterated


c      call finalize

      call CPU_TIME(t3)

        elapsetime1 = t2-t1
        elapsetime2 = t3-t2

      write(io2,50)  elapsetime1, elapsetime1/3600
      write(io2,100) elapsetime2, elapsetime2/3600

 50   format(5x,'initial set up time',g15.6,'s',g15.6,'h')
 100  format(5x,'running time',g15.6,'s',g15.6,'h')

      stop
      end



