      subroutine rdmodel
cccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccc


      implicit none
      include 'system.inc'
      include 'iochan.h'
      include 'control.inc'
      include 'model.inc' 
     
      integer*8 i,j,k,kk,inum,inat,inmol,ilabel

        inum = 0
        inat = 0
        inmol = 0
        ilabel=0
        iradius = 0
        ibonded = 0
    
2     read (io00, '(a80)', end = 200, 
     $err = 800) line

ccccccccccccccccccccccccccccccccccccccccccccccc
c     types of molecules in the  mixture
      if (line(1:11) .eq. 'num_species') then
         read(io00,*) num
         write(io7,*)'number of species',num
          inum = 1
         goto 2

cccccccccccccccccccccccccccccccccccccccccccccccc
c  number of molecules of each species
      elseif (line(1:12) .eq. 'nmol_species') then
            write(io7,*)'number of molecules of each species'
            do i = 1,num
            read(io00,*)nmol(i)
            write(io7,*)nmol(i)
            enddo
            inmol = 1
            goto 2

cccccccccccccccccccccccccccccccccccccccccccccccc
c     size of molecules
      elseif (line(1:16) .eq. 'size_per_species') then
          write(io7,*)'size of the molecule of each species'
          do i = 1,num
          read (io00, *) nat(i)
          write(io7,*)nat(i)
          enddo
         inat=1
         goto 2

ccccccccccccccccccccccccccccccccccccccccccccccccc      
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

cccccccccccccccccccccccccccccccccccccccccccccccccc         
c     number of CG bead types  
       elseif(line(1:16).eq.'num_CGbead_types') then
       read(io00,*)bead_types
       write(io7,*)'number of CG bead types:',bead_types
       goto 2


ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     label of the bead type
       elseif(line(1:10).eq.'label_bead') then
       write(io7,*)'label of the bead types'
       do i = 1,bead_types
       read(io00,*)lbl_bdtyp(i),label_bead(i)
       write(io7,*)lbl_bdtyp(i),label_bead(i)
       enddo
       goto 2


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     raidus of the bead
       elseif(line(1:11).eq.'bead_radius') then
       write(io7,*)'raidus of the bead types'
       iradius = 1
       do i = 1,bead_types
       read(io00,*)sigma(1,i)
       write(io7,*)sigma(1,i)
       enddo
       goto 2

ccccccccccccccccccccccccccccccccccccccccccccccccccc       
c      number of beads per molecule
       elseif(line(1:18).eq.'beads_per_molecule') then
       write(io7,*)'number of beads per molecule'
       do i = 1,num
       read(io00,*)nbd(i)
       write(io7,*)i,nbd(i)
       enddo
       goto 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      bead pair connectivity
       elseif(line(1:17).eq.'bead_connectivity') then

       read(io00,*) ((bead_pair(i,j), j=1,bead_types), i=1,bead_types)
       write(io7,*) ((bead_pair(i,j), j=1,bead_types), i=1,bead_types)

         ibonded = 1 
         do i = 1,bead_types
 
          bonded(i) =  0

           do j = 1,bead_types

           if(bead_pair(i,j).eq.1) bonded(i) = 1

           

           enddo

         enddo

       goto 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c      number of atom types in CG bead  and the number of atoms of each atom type    
       elseif(line(1:21).eq.'num_atom_types_CGbead')then
       write(io7,*)'number of atoms types per CG bead'
       do i = 1,bead_types
       noa_bead(i) = 0.0d0
       read(io00,*)num_atmtyp(i)
       write(io7,*)num_atmtyp(i)
       write(io7,*)' atom type and the number of those atoms in
     & the bead'
       do j = 1,num_atmtyp(i)
       read(io00,*)atmtyp(i,j),num_atm_bead(i,j)
       write(io7,*)atmtyp(i,j),num_atm_bead(i,j)
       noa_bead(i) = noa_bead(i) + num_atm_bead(i,j)
       enddo
       enddo
       goto 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c      CG mapping model to convert atomistic to CG trajectory

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

        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     comment line     
      elseif (line(1:1) .eq. '#') then
      goto 2 
c     stop
      elseif (line(1:9) .eq. 'endinput') then

c     unrecognised line
      else
         write(io7,*) 'Error: wrong input in input.dat'
         write(io7,*) line                                      
         istop=1
         stop
      endif
200    continue

        if(iradius.eq.0.and.ibonded.eq.0) then

         write(io7,*)'Error: no input for bead radius'
      write(io7,*)'either radii of the beads or bead connectivity
     &matrix has to be provided'

         stop
         endif


      return

 800  write (io7, *) '*** can''t read input file'
      istop=1
      end 

