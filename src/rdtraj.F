      subroutine rdtraj 
ccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'model.inc'
      include 'manual.inc'

      integer*8 temp,i,j,l,ii,jj,fu_1,is

      double precision minimum
       
C********************************
C**** Read- trajectory file *****
C********************************
C**** Read configuration data ***

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c read number of atoms



        if(isubstitute.eq.1) then

        do is = 1,num_iso

        fu_1 = 120 + is

        read(fu_1,*,iostat=ierr)natoms

                if (ierr.ne.0) goto 678

        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        read(fu_1,*)boxx,boxy,boxz
        else
        read(fu_1,*)
        endif

        l = 0
        do ii = 1,num

             do jj = 1,nmol(ii)

c loop over number of atoms to read x, y and z cordinates
        do i = 1,nat(ii)
        l = l + 1
        read(fu_1,*)label_atom(l),x(l),y(l),z(l)

        enddo

        enddo

        enddo

        enddo

         else


        read(io2,*,iostat=ierr) natoms

        if (ierr.ne.0) goto 678

        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        read(io2,*)boxx,boxy,boxz
        else
        read(io2,*)
        endif

        l = 0
        do ii = 1,num

             do jj = 1,nmol(ii)

c loop over number of atoms to read x, y and z cordinates
        do i = 1,nat(ii)
        l = l + 1
        read(io2,*)label_atom(l),x(l),y(l),z(l)

        enddo

        enddo

        enddo


        endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
c initialize the nr to zeros ( matrix giving the number of atoms per atom type)  

        do i = 1,maxtyp
        nr(i) = 0
        enddo

        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt') then

c  computing the cutoff

          minimum = boxy

          if(Boxx.lt.minimum) minimum =boxx
          if(Boxz.lt.minimum) minimum = boxz

          minimum = minimum/2.0d0


        cutoff =  minimum -  mod(minimum,bin_width)

ccccccccccccccccccccccccccccccccccccccc
c   number of bins for rdf calculation

          nbins = int(cutoff/bin_width)

          if(nbins.gt.max_binr) max_binr = nbins

          do i = 1, max_binr
          frame(i)= frame(i) + 1
          enddo

         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prepare coordinate matrix according to the atom type

         do i = 1,natoms

              do j = 1,ntypes

              if(label_atom(i).eq.label_type(j)) then

                   itype(i)=j
                                
               nr(j) = nr(j) + 1

               temp = nr(j)
c matrix rr(atom type 'i',nunber of atoms of type 'i', position coordinate)
              rr(j,temp,1) = x(i)
              rr(j,temp,2) = y(i)
              rr(j,temp,3) = z(i)

              endif

            enddo

        enddo

c          do i = 1,ntypes
c         do j = 1,nr(i)
c          write(io7,*)i,j,rr(i,j,1),rr(i,j,2),rr(i,j,3)
c          enddo
c          enddo

         
              



c        write(io7,*)'reading atomistic trajectory'
678       continue
    
          end 
