       subroutine wrtraj
ccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       wrtraj  
c     Usage
c       call wrtraj 
c     Description
c       write in trajectory with isotopologue substitution 
c       number of atoms
c       comment line: box dimensions in case of NPT
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

      integer*8 i,j,k,l,ii,jj,jjj,temp_int

      integer*8 is,flag_cnt,mm,nn,cnt,iflag,
     & natoms1,natoms2,imol,iframe,nmol1(num),fu_1

      double precision iso_ratio,ran_num(num,natoms)
      double precision x1(natoms),y1(natoms),z1(natoms),
     &x2(natoms),y2(natoms),z2(natoms)
       character(len=10) label1(natoms),label2(natoms),file_id

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

            
          do is = 1,num_iso

         file_id = text2(is)


         filename1= 'trajectory' // trim(adjustl(file_id))

          fu_1 = 120+is

          open(unit=fu_1,file=filename1,status='unknown')

                  
C********************************
C**** Read- trajectory file *****
C********************************
C**** Read configuration data ***
         ierr = 0

         do while(ierr.eq.0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read number of atoms
         
        
        read(io2,*,iostat=ierr)natoms

        if(ierr.ne.0) exit

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

       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ntypes1 = ntypes
            bead_types1 = bead_types
            noex1 = noex
            num1 = num
           
            do i = 1,num
            nat2(i) = nat(i)
            nmol2(i) = nmol(i)
            nmol1(i) = nmol1(i)
            enddo

            do imol = 1,num

          if(deut_flag(is,imol).eq.0) then

          iso_ratio = 0.0d0

         elseif(deut_flag(is,imol).eq.1) then

          iso_ratio = 1.0d0

            do i = 1,ntypes

            if(label_type(i)(1:1).eq.'H') label_type(i)(1:1) = 'D'

            enddo


         elseif(deut_flag(is,imol).eq.2) then

                

          iso_ratio = iso_fraction(is,imol)

          do i = 1,ntypes

         
          if(label_type(i)(1:1).eq.'H'.and.mol(i).eq.imol) then

                 ntypes1 = ntypes1 + 1
                  num1 = num1 + 1
          label_type(ntypes1) = label_type(i)
          label_type(ntypes1)(1:1) = 'D'

           nmol2(num1) = int(nmol(imol)*iso_fraction(is,imol))
         nmol2(imol)=(nmol(imol)-int(nmol(imol)*iso_fraction(is,imol)))
           deut_flag(is,num1) = 2

           
                  nat2(num1) = nat(imol)

                    do j = 1,noex

                   if(label_ex(j).eq.label_type(i)) then
                     noex1 = noex1+1
                    ex(ntypes1) = 1
                   endif

                   enddo

                   label_iso(is,ntypes1) = label_type(ntypes1)
                   natural(is,ntypes1) = 1.0
                   isotope(is,ntypes1) = 0.0
                   natural(is,i) = 1.0
                   isotope(is,i) = 0.0

           endif

           enddo

               bead_sec(num1)= bead_sec(imol)
               nbd(num1) = nbd(imol)

                 do j = 1,bead_sec(num1)
                 bead_num(num1,j) = bead_num(imol,j)
                 rep_num(num1,j) = rep_num(imol,j)
                 atom_step(num1,j) = atom_step(imol,j)
                 enddo

                 do j = 1,bead_sec(num1)                

                 do jj = 1,bead_num(num1,j)
                 lbl_bead(num1,j,jj) = lbl_bead(imol,j,jj)
                 lbl_bead(num1,j,jj)(1:1) = 'D'
                 list_size(num1,j,jj) = list_size(imol,j,jj)

                 enddo

                 enddo

                 do j = 1,bead_sec(num1)

                 do jj = 1,bead_num(num1,j)

                 do jjj = 1,list_size(num1,j,jj)
                  atm_list(num1,j,jj,jjj) = atm_list(imol,j,jj,jjj)

                 enddo

                 enddo

                 enddo

                 bead_types1 = bead_types1+bead_types
                 do j = 1,bead_types
                 label_bead(j+bead_types) = label_bead(j)
                 label_bead(j+bead_types)(1:1) = 'D'
                  
                  
                 sigma_rmsd(1,j+bead_types) = sigma_rmsd(1,j)

                 if(iradius.eq.1) then
                         sigma(1,j+bead_types) = sigma(1,j)

                   endif

                   if (ibonded.eq.1) then

                    do jj = 1,bead_types

                  bead_pair(j+bead_types,jj)  = 0
                  bead_pair(j,jj+bead_types) = 0

                  bead_pair(j+bead_types,jj+bead_types)=bead_pair(j,jj)

          if(bead_pair(j+bead_types,jj+bead_types).eq.1)
     & bonded(j+bead_types)=1

                    enddo
                    endif

                    num_atmtyp(j+bead_types) = num_atmtyp(j)
                    

                    do jj = 1,num_atmtyp(j)
                    atmtyp(j+bead_types,jj) = atmtyp(j,jj)
                    temp_int = atmtyp(j,jj)
           if(label_type(temp_int)(1:1).eq.'H')atmtyp(j+bead_types,jj)
     & =ntypes1
                    

                    num_atm_bead(j+bead_types,jj) = num_atm_bead(j,jj)

                    enddo


                   
                 enddo


                 
                  
          endif

         
       
          
         nmol1(imol) = int(nmol(imol)*iso_ratio)

         do j = 1,nmol1(imol)

88       ran_num(imol,j) = int(rand()*nmol(imol))+1


          do k = 1,j-1

          if(ran_num(imol,k).eq.ran_num(imol,j)) goto 88

          enddo

          enddo

          enddo

                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     

            flag_cnt = 0 
            cnt = 0
            mm = 0
            nn = 0


        do ii = 1, num

        do i = 1,nmol(ii)

        iflag = 0 

         do k = 1,nmol1(ii)

         if(ran_num(ii,k).eq.i)iflag = 1

         enddo


         flag_cnt = flag_cnt + iflag

         if(iflag.eq.0) then

           do j = 1,nat(ii)

           cnt = cnt + 1

           mm = mm + 1

           label1(mm) = label_atom(cnt)
           x1(mm) = x(cnt)
           y1(mm) = y(cnt)
           z1(mm) = z(cnt)

           enddo


          elseif(iflag.eq.1) then


           do j = 1,nat(ii)

           cnt = cnt + 1

            nn = nn + 1

            if(label_atom(cnt)(1:1).eq.'H') then

            label2(nn) = label_atom(cnt)

            label2(nn)(1:1) = 'D'

            else
            label2(nn) = label_atom(cnt)
             endif


             x2(nn) = x(cnt)
             y2(nn) = y(cnt)
             z2(nn) = z(cnt)

             enddo

             endif

             enddo

             enddo

            

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           
        write(fu_1,*)natoms
       

        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        write(fu_1,*)boxx,boxy,boxz
        else
        write(fu_1,*)
        endif

         natoms1 = mm
         natoms2 = nn

        do i = 1,natoms1

        write(fu_1,*)label1(i),x1(i),y1(i),z1(i)
        enddo

        do i = 1,natoms2

        write(fu_1,*)label2(i),x2(i),y2(i),z2(i)
        enddo

        
          enddo

         
          

        close(io2)
       open (unit = io2, file = 'trajectory',  status = 'unknown')
         
          enddo

        end
