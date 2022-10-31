      subroutine ofiles
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       ofiles
c
c     Usage
c       call ofiles
c
c     Description
c       Open all input and output files 

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
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'iochan.h'
      include 'control.inc'

      integer ios

c     ... 0 stderr (preconnected)
c     checks on the files supposed to be already there

c     check if  input file is there      



      open (unit = io7, file = 'fort.7',  status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')


c   .... 1 input.dat       

      open (unit = io1, file = 'input.dat',  status = 'OLD',
     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')
      if (ios.ne.0) then 
         write (io7,*) '* FATAL (ofile) input data file does not exist'
         istop=1
      endif

c   .... 2  model.dat

       open (unit = io00, file = 'model.dat',  status = 'OLD',
     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')
      if (ios.ne.0) then
         write (io7,*) '* FATAL (ofile) input data file does not exist'
         istop=1
      endif


c     ... 3 trajectory file 

      open (unit = io2, file = 'trajectory',  status = 'OLD',
     $      iostat=ios, access = 'SEQUENTIAL', form = 'FORMATTED')
      if (ios.ne.0) then 
         write (io7,*) '* atomisitc trajectory does not exist'
c         istop=1
      endif

      if(ierror.eq.1) then
      write(io7,*)'error input file does not exist'
      endif

c     ... 4 CG trajectory file

      open (unit = io3, file = 'CG_trajectory',  status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')



       end
