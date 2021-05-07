      subroutine getdata

c     start_doc
c     Name
c       getdata
c     Usage
c       call getdata
c     Description
c       read in from input data file
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
      include 'model.inc'
      include 'iochan.h'
      include 'control.inc'
     
      integer intypes,inframes,icutoff,ibox,ibin_width,isizemol,iqmin,
     &iqdelta,iqmax,ibroadening,iFWHM,iwindowfn,ilabel,iex,ibc,i,j,
     &inum_iso
     
  
      intypes = 0
      inframes = 0 
      icutoff = 0
      ibox = 0
      ibin_width = 0
      isizemol = 0
      iqmin = 0
      iqdelta = 0
      iqmax = 0
      ibroadening = 0
      iFWHM = 0
      iwindowfn = 0
      ilabel = 0
      iex = 0
      ibc = 0
      inum_iso = 0
      

1     read (io1, '(a80)', end = 100, 
     $err = 900) line

c     title
      if (line(1:5) .eq. 'title') then
         read(io1,*) title
         goto 1

c     number of atom types
      elseif (line(1:6) .eq. 'ntypes') then
         read (io1, *) ntypes
         intypes=1
         goto 1
      
c     number of frames 
      elseif (line(1:7) .eq. 'nframes') then
         read (io1, *) nframes
         inframes=1
         goto 1
c     length of the simulation box
      elseif (line(1:10) .eq. 'box_length') then
         read (io1, *) box
         ibox=1
         goto 1

c     histogram bin width
      elseif (line(1:9) .eq. 'bin_width') then
         read (io1, *) bin_width
         ibin_width=1
         goto 1
 

c     cutoff for the distance calculation
      elseif (line(1:6) .eq. 'cutoff') then
         read (io1, *) cutoff
         icutoff=1
         goto 1

c     label of the atoms
              
      elseif (line(1:10) .eq. 'label_atom') then
         do i = 1,ntypes
         read (io1, *) label_type(i)
         enddo
         ilabel = 1
         goto 1

c     neutron weights         
      elseif (line(1:18) .eq. 'scattering_lengths') then
         do i = 1,ntypes
         read (io1, *) bc(i)
         bc(i) = bc(i)*0.1d0
         enddo
         read(io1,*)bc(ntypes+1)
         bc(ntypes+1) = bc(ntypes+1)*0.1d0
         ibc = 1
         goto 1

c     exchangeableatom type

      elseif (line(1:18) .eq. 'exchangeable_atoms') then
         read(io1,*)noex
         do i = 1,noex
         read (io1, *)label_ex(i)
         enddo
         iex = 1
         goto 1

c     number of isotopologues         

       elseif(line(1:20).eq.'number_isotopologues') then
         read(io1,*)num_iso
         inum_iso = 1
          goto 1
c    isotpe details in each isotopologue

       elseif(line(1:13).eq.'isotopologues') then
    
          do i = 1,num_iso
          read(io1,*)        
          read(io1,*)
          
          do j = 1,ntypes
          read(io1,*)label_iso(i,j), natural(i,j),isotope(i,j)
          enddo
          enddo
          goto 1

c     minimum value of Q range
      elseif (line(1:4) .eq. 'qmin') then
         read (io1, *) qmin
         iqmin=1
         goto 1

c     step size of Q axis            
      elseif (line(1:6) .eq. 'qdelta') then
         read (io1, *) qdelta
         iqdelta=1
         goto 1

c     maximum value of Q range
      elseif (line(1:4) .eq. 'qmax') then
         read (io1, *) qmax
         iqmax=1
         goto 1

c     broadening function
      elseif (line(1:10) .eq. 'broadening') then
         read (io1, *) broadfn
         ibroadening = 1
         goto 1

C    FWHM value
      elseif (line(1:4) .eq.'FWHM') then 
              read(io1, *) FWHM
              iFWHM = 1
         goto 1

C     window function 
      elseif (line(1:8) .eq.'windowfn') then
         read(io1, *) windowfn
         iwindowfn=0
         goto 1



c     comment line     
      elseif (line(1:1) .eq. '#') then
      goto 1 
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
100    continue

      return

 900  write (io7, *) '*** can''t read input file'
      istop=1
      end 

