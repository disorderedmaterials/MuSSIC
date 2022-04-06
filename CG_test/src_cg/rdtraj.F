      subroutine rdtraj 
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

      implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'model.inc'

      integer temp,i,j,l,ii,jj


C********************************
C**** Read- trajectory file *****
C********************************
C**** Read configuration data ***


c read number of atoms
        read(io3,*)natoms
        read(io3,*)

      
c loop over number of atoms to read x, y and z cordinates
        do i = 1,natoms
        
        read(io3,*)label_atom(i),x(i),y(i),z(i)
c        label_atom(l) = label_bead(ii,i)
c        write(io7,*)label_atom(l),x(l),y(l),z(l)
        enddo

            

    


c initialize the nr to zeros ( matrix giving the number of atoms per atom type)  

        do i = 1,bead_types
        nr(i) = 0
        enddo
          
   
c prepare coordinate matrix according to the atom type

         do i = 1,natoms

              do j = 1,bead_types

              if(label_atom(i).eq.label_bead(j)) then
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

        do i = 1,bead_types
        write(io7,*)natoms,i,nr(i)
        enddo
        
 
      end 
