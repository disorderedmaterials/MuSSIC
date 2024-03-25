      subroutine rdtraj_cg 
ccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       rdtraj  
c     Usage
c       call rdtraj 
c     Description
c       read in trajectory  file 
c       number of atoms
c       comment line: frame number
c       atom_label       x  y  z  in Angstroms
c        label_atom(l) x(l) y(l) z(l) 
c     end_doc
ccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'model.inc'

      integer*8 temp,i,j,is,fu_1,l,ii,jj

      double precision minimum


C********************************
C**** Read- trajectory file *****
C********************************
C**** Read configuration data ***



c read number of beads
        read(io3,*,iostat=ierr)natoms

        if(ierr.ne.0) goto 676

        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        read(io3,*)boxx,boxy,boxz
        else
        read(io3,*)
        endif

      
c loop over number of beads to read x, y and z cordinates
        do i = 1,natoms
        
        read(io3,*)label_atom(i),x(i),y(i),z(i)
        
        enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialize the nr to zeros ( matrix giving the number of beads per  bead type)  

        do i = 1,maxtyp
        nr(i) = 0
        enddo
          
           if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt') then
ccccccccccccccccccccccc
c  computing the cutoff

          minimum = boxy

          if(Boxx.lt.minimum) minimum =boxx
          if(Boxz.lt.minimum) minimum = boxz

          minimum = minimum/2.0d0


        cutoff =  minimum -  mod(minimum,bin_width)

cccccccccccccccccccccccccccccccccccccc
c   number of bins for rdf calculation

          nbins = int(cutoff/bin_width)


          if(nbins.ge.max_binr) max_binr = nbins

          do i = 1, max_binr
          frame(i)= frame(i) + 1
          enddo

         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prepare coordinate matrix according to the bead type

         do i = 1,natoms

              do j = 1,bead_types

              if(label_atom(i).eq.label_bead(j)) then
                   itype(i)=j
                                
               nr(j) = nr(j) + 1

               temp = nr(j)
c matrix rr(bead type 'i',nunber of beads of type 'i', position coordinate)
              rr(j,temp,1) = x(i)
              rr(j,temp,2) = y(i)
              rr(j,temp,3) = z(i)

              nmol_label(j,temp) = mol_num(i)

              endif

            enddo

        enddo


c        do i = 1,bead_types
c        write(io7,*)'natoms:',natoms,'bead type',i,'no. of beads',nr(i)
c        enddo
        
676      continue 
      end 
