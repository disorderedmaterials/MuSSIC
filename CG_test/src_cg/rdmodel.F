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
     
      integer i,j,k,kk,inum,inat,inmol,ilabel
      character*10 temp_char

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
      



       elseif(line(1:20).eq.'avg_num_atoms_CGbead') then
        read(io00,*)avg_noa_bead
        write(io7,*)'average number of atoms per bead:', avg_noa_bead
       goto 2

       elseif(line(1:16).eq.'num_CGbead_types') then
       read(io00,*)bead_types
       write(io7,*)'number of CG bead types:',bead_types
       goto 2 
          

       elseif(line(1:21).eq.'num_atom_types_CGbead')then
       write(io7,*)'number of atoms types per CG bead'
       do i = 1,bead_types
       noa_bead(i) = 0
       read(io00,*)num_atmtyp(i)
       write(io7,*)num_atmtyp(i)
       write(io7,*)'label of the atom type and number of those atoms in 
     & the bead'
       do j = 1,num_atmtyp(i)
       read(io00,*)atmtyp(i,j),num_atm_bead(i,j)
       write(io7,*)atmtyp(i,j),num_atm_bead(i,j)
       noa_bead(i) = noa_bead(i) + num_atm_bead(i,j)
       enddo
       enddo
       goto 2

       elseif(line(1:9).eq.'bead_size')then
        do i = 1,bead_types
        read(io00,*)lbl_bdtyp(i),label_bead(i),size_bead(i)
        enddo
        goto 2 

       elseif(line(1:10).eq.'cg_mapping')then
       do i = 1,num
       read(io00,*)
       write(io7,*)'molecule:',i
       read(io00,*)bead_sec(i)
       write(io7,*)'number of sections:',bead_sec(i)
       write(io7,*)'no. of beads, no. of repetitions, atom step'
       do j = 1,bead_sec(i)
       read(io00,*)bead_num(i,j),rep_num(i,j),atom_step(i,j)
       write(io7,*)bead_num(i,j),rep_num(i,j),atom_step(i,j)
c       write(io7,*)'bead_label,atom list size, atoms for com'
       do k = 1,bead_num(i,j)
       read(io00,*)lbl_bead(i,j,k),list_size(i,j,k)
       write(io7,*)lbl_bead(i,j,k),list_size(i,j,k)
       do kk = 1,list_size(i,j,k)
       read(io00,*)atm_list(i,j,k,kk)
       write(io7,*)atm_list(i,j,k,kk)
       enddo
        enddo
       enddo
       enddo



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

