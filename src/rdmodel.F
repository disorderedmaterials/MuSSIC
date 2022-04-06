      subroutine rdmodel

c     start_doc
c     Name
c       rdmodel
c     Usage
c       call rdmodel
c     Description
c       reads model parameters 
c       skips comment lines starting with #
c*******************************************
c       elseif (line(1:1) .eq. '#') then
c       goto 1 
c******************************************* 

c     then writes the header of the output file
c     Imported
c     
c     end_doc 

      implicit none
      include 'system.inc'
      include 'iochan.h'
      include 'control.inc'
      include 'model.inc' 
     
      integer*8 i,j,inum,inat,inmol,ilabel

        inum = 0
        inat = 0
        inmol = 0
        ilabel=0
    
2     read (io00, '(a80)', end = 200, 
     $err = 800) line

c     types of molecules in the  mixture
      if (line(1:11) .eq. 'num_species') then
         read(io00,*) num
         write(io7,*)'number of species',num
          inum = 1
         goto 2
c  number of molecules of each species
      elseif (line(1:12) .eq. 'nmol_species') then
            write(io7,*)'number of molecules of each species'
            do i = 1,num
            read(io00,*)nmol(i)
            write(io7,*)nmol(i)
            enddo
            inmol = 1
            goto 2


c     size of molecules
      elseif (line(1:16) .eq. 'size_per_species') then
          write(io7,*)'size of the molecule of each species'
          do i = 1,num
          read (io00, *) nat(i)
          write(io7,*)nat(i)
          enddo
         inat=1
         goto 2
      
c      labels of the atoms in molecules in sequence
      elseif (line(1:14) .eq. 'sequence_atoms') then
              write(io7,*)'label of atoms in sequence'
               do i = 1,num
                 do j = 1,nat(i)
                     read(io00,*)label(i,j)
                     write(io7,*)label(i,j)
                  enddo
               enddo
         ilabel=1
         goto 2


c     comment line     
      elseif (line(1:1) .eq. '#') then
      goto 2 
c     stop
      elseif (line(1:9) .eq. 'endinput') then

c     unrecognised line
      else
         write(io7,*) 'FATAL (getdata)'
         write(io7,*) 'bad input line'
         write(io7,*) line                                      
         istop=1
         stop
      endif
200    continue

      return

 800  write (io7, *) '*** can''t read input file'
      istop=1
      end 

