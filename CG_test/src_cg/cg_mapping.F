        subroutine cg_mapping
        implicit none
        
      include 'system.inc'
      include 'iochan.h'
      include 'control.inc'
      include 'config.inc'
      include 'model.inc'


      integer i,ii,j,jj,k,kk,iii,l,ll,total_beads,ndxaux,sel,temp1,
     & temp,atmcnt
         double precision cen_x(10), cen_y(10),cen_z(10),Dx,Dy,Dz,
     &cox,coy,coz

C**** Read- trajectory file *****
C********************************
C**** Read configuration data ***

         total_beads = 0

          do iii = 1,nframes

c read number of atoms
        read(io2,*)natoms
        read(io2,*)

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

        

           l = 0
           atmcnt = 0
           

           do ii = 1,num

              
           if(ii.gt.1) atmcnt = atmcnt + nmol(ii-1)*nat(ii-1)

                  
                 nbd(ii) = 0

                if(ii.le.num-1) then

                do jj = 1,nmol(ii)

                 do i = 1,bead_sec(ii)

                 do j = 1,rep_num(ii,i)

                 ndxaux = (j-1)*atom_step(ii,i)+((jj-1)*nat(ii))+atmcnt


                 do k = 1,bead_num(ii,i)

                 if(jj.eq.1) nbd(ii) = nbd(ii) + 1

                  l = l + 1

                      bead_lbl(l) = lbl_bead(ii,i,k)

                      if(jj.eq.1) then

                              temp = nbd(ii)

                              label(ii,temp) = bead_lbl(l)

                      endif

                      temp1 = ndxaux + atm_list(ii,i,k,1)


                                           
 
                      
                     cen_x(1) = x(temp1)
                     cen_y(1) = y(temp1)
                     cen_z(1) = z(temp1)

                     cox = cen_x(1)
                     coy = cen_y(1)
                     coz = cen_z(1)





                  do kk = 2,list_size(ii,i,k)

                  sel =  ndxaux + atm_list(ii,i,k,kk)

                       
             
                 Dx = x(sel) - cen_x(kk-1)
                 Dy = y(sel) - cen_y(kk-1)
                 Dz = z(sel) - cen_z(kk-1)



Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)



                 cen_x(kk) = cen_x(kk-1) + Dx
                 cen_y(kk) = cen_y(kk-1) + Dy
                 cen_z(kk) = cen_z(kk-1) + Dz


                 cox = cox + cen_x(kk)
                 coy = coy + cen_y(kk)
                 coz = coz + cen_z(kk)

                 enddo

                 

                 cox = cox/dble(list_size(ii,i,k))
                 coy = coy/dble(list_size(ii,i,k))
                 coz = coz/dble(list_size(ii,i,k))

        cox = cox - boxx* dble(Idint((cox/boxx) +9999.5d0) - 9999)
        coy = coy - boxy* dble(Idint((coy/boxy) +9999.5d0) - 9999)
        coz = coz - boxz* dble(Idint((coz/boxz) +9999.5d0) - 9999)

                  bead_x(l) = cox
                  bead_y(l) = coy
                  bead_z(l) = coz

                              
                 enddo

                 enddo

                 enddo

                 enddo

         else



                do jj = 1,nmol(ii),1

                 do i = 1,bead_sec(ii)

                 do j = 1,rep_num(ii,i)

                 ndxaux = (j-1)*atom_step(ii,i)+((jj-1)*nat(ii))+atmcnt


                 do k = 1,bead_num(ii,i)

                 if(jj.eq.1) nbd(ii) = nbd(ii) + 1

                  l = l + 1

                      bead_lbl(l) = lbl_bead(ii,i,k)

                      if(jj.eq.1) then

                              temp = nbd(ii)

                              label(ii,temp) = bead_lbl(l)

                      endif

                      temp1 = ndxaux + atm_list(ii,i,k,1)


                                           
 
                      
                     cen_x(1) = x(temp1)
                     cen_y(1) = y(temp1)
                     cen_z(1) = z(temp1)

                     cox = cen_x(1)
                     coy = cen_y(1)
                     coz = cen_z(1)





                  do kk = 2,list_size(ii,i,k)

                  sel =  ndxaux + atm_list(ii,i,k,kk)

                       
             
                 Dx = x(sel) - cen_x(kk-1)
                 Dy = y(sel) - cen_y(kk-1)
                 Dz = z(sel) - cen_z(kk-1)



Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)



                 cen_x(kk) = cen_x(kk-1) + Dx
                 cen_y(kk) = cen_y(kk-1) + Dy
                 cen_z(kk) = cen_z(kk-1) + Dz


                 cox = cox + cen_x(kk)
                 coy = coy + cen_y(kk)
                 coz = coz + cen_z(kk)

                 enddo

                 

                 cox = cox/dble(list_size(ii,i,k))
                 coy = coy/dble(list_size(ii,i,k))
                 coz = coz/dble(list_size(ii,i,k))

        cox = cox - boxx* dble(Idint((cox/boxx) +9999.5d0) - 9999)
        coy = coy - boxy* dble(Idint((coy/boxy) +9999.5d0) - 9999)
        coz = coz - boxz* dble(Idint((coz/boxz) +9999.5d0) - 9999)

                  bead_x(l) = cox
                  bead_y(l) = coy
                  bead_z(l) = coz

                              
                 enddo

                 enddo

                 enddo

                 enddo

               endif                 

                 enddo

                 total_beads = l


                  write(io3,'(I6)')total_beads
                  write(io3,*) 
                  do k = 1,total_beads
        write(io3,*)bead_lbl(k),bead_x(k),
     & bead_y(k),bead_z(k)

                  enddo

                  

           enddo

           write(io7,*)'molecules and bead labels'
           do i = 1,num

           write(io7,*)i,nbd(i)
           do j = 1,nbd(i)
           write(io7,*)label(i,j)
           enddo

           enddo
           

           close(io3)

         end 
