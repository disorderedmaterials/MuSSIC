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


c     ... 0 stderr (preconnected)
c     checks on the files supposed to be already there

c     check if  input file is there      

          itraj = 1
          icgtraj = 1

      open (unit = io7, file = 'out.log',  status = 'UNKNOWN',
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
         write (io7,*) '* FATAL (ofile) model data file does not exist'
         istop=1
      endif




       end
