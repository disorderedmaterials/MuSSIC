      subroutine ofiles
      implicit none

c     start_doc
c     Name
c       ofiles
c
c     Usage
c       call ofiles
c
c     Description
c       Open all input and output files 
c       are there
c       includes 
c       file `iochan.h'.
c       file `control.inc'.
c
c     
c     Parameters
c     istop  passed to the main to stop the program
c     
c
c     end_doc

      include 'iochan.h'
      include 'control.inc'

      integer ios

c     ... 0 stderr (preconnected)
c     checks on the files supposed to be already there

c     check if  input file is there      



      open (unit = io7, file = 'fort.7',  status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c   ... 0  neutron scaterring lengths 

c      open(unit = io0, file = 'weights.dat', status= 'OLD',
c     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')
c      if (ios.ne.0) then 
c         write (io7,*) '* FATAL (ofile) weigths file does not exist'
c         istop=1
c      endif



c   .... 1 input data       

      open (unit = io1, file = 'input.dat',  status = 'OLD',
     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')
      if (ios.ne.0) then 
         write (io7,*) '* FATAL (ofile) input data file does not exist'
         istop=1
      endif

       open (unit = io00, file = 'model.dat',  status = 'OLD',
     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')
      if (ios.ne.0) then
         write (io7,*) '* FATAL (ofile) input data file does not exist'
         istop=1
      endif


c     ... 2 trajectory file 

      open (unit = io2, file = 'trajectory',  status = 'OLD',
     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')
      if (ios.ne.0) then 
         write (io7,*) '* FATAL (ofile) trajectory file does not exist'
         istop=1
      endif
      if(ierror.eq.1) then
      write(io7,*)'error input file does not exist'
      endif



c     ... 3 partial unweighted g(r)

c      open (unit = io3, file = 'unweighted.gr',  status = 'UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 4 total unweighted g(r) 
      open (unit = io4,file = 'unweighted-total.gr',status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 5 partial unweighted S(Q)
c      open (unit = io5, file = 'unweighted.sq',  status = 'UNKNOWN',
c     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 6 total unweighted S(Q)  
      open (unit=io66, file = 'unweighted-total.sq',status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 8 partial weighted S(Q)  
c      open (unit = io8, file = 'weighted.sq',status = 'UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 9 total weighted S(Q)  
c      open (unit = io9, file = 'weighted_total.sq',status = 'UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 10 total weighted G(r)
c      open (unit=io10,file='weighted_total.gr',status='UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 11 partial weighted S(Q) - deuterated
c      open (unit=io11,file = 'deut_weighted.sq',status='UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 12 total weighted S(Q)  - deuterated
c      open (unit=io12,file='deut_weighted_total.sq',status='UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 13 total weighted G(r)  - deuterated
c      open (unit=io13,file='deut_weighted_total.gr',status='UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 14 partial weighted S(Q)  - mixed
c      open (unit=io14,file='mix_weighted.sq',status='UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 15 total weighted S(Q) - mixed
c      open (unit=io15,file = 'mix_weighted_total.sq',status='UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

cc     ... 16 total weighted G(r)  - mixed
c      open (unit=io16,file='mix_weighted_total.gr',status='UNKNOWN',
c     $      access = 'SEQUENTIAL', form = 'FORMATTED')

       end
