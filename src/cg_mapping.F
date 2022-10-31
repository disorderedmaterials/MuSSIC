        subroutine cg_mapping
        implicit none
        
      include 'system.inc'
      include 'iochan.h'
      include 'control.inc'
      include 'config.inc'
      include 'model.inc'


      integer*8 i,ii,j,jj,k,kk,iii,l,ll,total_beads,ndxaux,sel,temp1,
     & temp,atmcnt
         double precision cen_x(100), cen_y(100),cen_z(100),Dx,Dy,Dz,
     &cox,coy,coz,box,boy,boz,mass,total_mass,bondx(100),bondy(100),
     &bondz(100)

C**** Read- trajectory file *****
C********************************
C**** Read configuration data ***
       open (unit = io2, file = 'trajectory',  status = 'unknown')

         total_beads = 0

          do iii = 1,nframes

c read number of atoms
        read(io2,*)natoms
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

       

            
c        write(io7,*)'cg_mapping: read atomistic trajectory'
        

           l = 0
           atmcnt = 0
           

           do ii = 1,num

              
           if(ii.gt.1) atmcnt = atmcnt + nmol(ii-1)*nat(ii-1)

                  
                 temp = 0

                if(ii.le.num) then

                do jj = 1,nmol(ii)

                 do i = 1,bead_sec(ii)

                 do j = 1,rep_num(ii,i)

                 ndxaux = (j-1)*atom_step(ii,i)+((jj-1)*nat(ii))+atmcnt


                 do k = 1,bead_num(ii,i)

c                if(jj.eq.1) temp = temp + 1

                  l = l + 1

                      bead_lbl(l) = lbl_bead(ii,i,k)

                      if(jj.eq.1) then

                              temp = temp + 1

                              label(ii,temp) = bead_lbl(l)

                      endif

                      temp1 = ndxaux + atm_list(ii,i,k,1)


c                   if(label_atom(temp1) .eq.'C')  mass = 12.011d0
c                   if(label_atom(temp1) .eq.'H')  mass = 1.000784d0
c                   if(label_atom(temp1) .eq.'O')  mass = 15.999d0
c                   if(label_atom(temp1) .eq.'N')  mass = 14.0067d0
                                       
 
                      
                     cen_x(1) = x(temp1)
                     cen_y(1) = y(temp1)
                     cen_z(1) = z(temp1)

c                     cox = mass*cen_x(1)
c                     coy = mass*cen_y(1)
c                     coz = mass*cen_z(1)

                     cox = cen_x(1)
                     coy = cen_y(1)
                     coz = cen_z(1)

c                     total_mass = mass


                  do kk = 2,list_size(ii,i,k)

                  sel =  ndxaux + atm_list(ii,i,k,kk)

                                    
                 Dx = x(sel) - cen_x(kk-1)
                 Dy = y(sel) - cen_y(kk-1)
                 Dz = z(sel) - cen_z(kk-1)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

                Dx = Dx - boxx*nint(Dx/boxx)
                Dy = Dy - boxy*nint(Dy/boxy)
                Dz = Dz - boxz*nint(Dz/boxz)

c                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
c                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
c                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)


                 cen_x(kk) = cen_x(kk-1) + Dx
                 cen_y(kk) = cen_y(kk-1) + Dy
                 cen_z(kk) = cen_z(kk-1) + Dz


c                   if(label_atom(sel) .eq.'C')  mass = 12.011d0
c                   if(label_atom(sel) .eq.'H')  mass = 1.000784d0
c                   if(label_atom(sel) .eq.'O')  mass = 15.999d0
c                   if(label_atom(sel) .eq.'N')  mass = 14.0067d0
                 
c                 cox = cox + mass*cen_x(kk)
c                 coy = coy + mass*cen_y(kk)
c                 coz = coz + mass*cen_z(kk)

c                 total_mass = total_mass + mass

                 
                 cox = cox + cen_x(kk)
                 coy = coy + cen_y(kk)
                 coz = coz + cen_z(kk)

                 enddo

c                 box = 0.0d0
c                 boy = 0.0d0
c                 boz = 0.0d0

c                 do kk = 1,list_size(ii,i,k)-1

c                  bondx(kk) = (cen_x(kk)+cen_x(kk+1))/2.0d0
c                  bondy(kk) = (cen_y(kk)+cen_y(kk+1))/2.0d0
c                  bondz(kk) = (cen_z(kk)+cen_z(kk+1))/2.0d0

c                  box = box + bondx(kk)
c                  boy = boy + bondy(kk)
c                  boz = boz + bondz(kk)

c                  enddo

                   
c                 box = box/dble(list_size(ii,i,k)-1)
c                 boy = boy/dble(list_size(ii,i,k)-1)
c                 boz = boz/dble(list_size(ii,i,k)-1)

                 cox = cox/dble(list_size(ii,i,k))
                 coy = coy/dble(list_size(ii,i,k))
                 coz = coz/dble(list_size(ii,i,k))

c                 cox = cox/total_mass
c                 coy = coy/total_mass
c                 coz = coz/total_mass

                 if(cox.gt.boxx) cox = cox -boxx
                 if(cox.lt.0.0d0) cox = cox +boxx
                 if(coy.gt.boxy) coy = coy -boxy
                 if(coy.lt.0.0d0) coy = coy +boxy
                 if(coz.gt.boxz) coz = coz -boxz
                 if(coz.lt.0.0d0) coz = coz +boxz
                 
c                 if(box.gt.boxx) box = box -boxx
c                 if(box.lt.0.0d0) box = box +boxx
c                 if(boy.gt.boxy) boy = boy -boxy
c                 if(boy.lt.0.0d0) boy = boy +boxy
c                 if(boz.gt.boxz) boz = boz -boxz
c                 if(boz.lt.0.0d0) boz = boz +boxz
c        cox = cox - boxx* dble(Idint((cox/boxx) +9999.5d0) - 9999)
c        coy = coy - boxy* dble(Idint((coy/boxy) +9999.5d0) - 9999)
c        coz = coz - boxz* dble(Idint((coz/boxz) +9999.5d0) - 9999)

                  bead_x(l) = cox
                  bead_y(l) = coy
                  bead_z(l) = coz


                 do kk = 1,list_size(ii,i,k)-1

                  bondx(kk) = (cen_x(kk)+cen_x(kk+1))/2.0d0
                  bondy(kk) = (cen_y(kk)+cen_y(kk+1))/2.0d0
                  bondz(kk) = (cen_z(kk)+cen_z(kk+1))/2.0d0

                  box = box + bondx(kk)
                  boy = boy + bondy(kk)
                  boz = boz + bondz(kk)

                  enddo



c                  bead_x(l) = box
c                  bead_y(l) = boy
c                  bead_z(l) = boz
                              
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


                  l = l + 1

                      bead_lbl(l) = lbl_bead(ii,i,k)


                 if(jj.eq.1) then
                 temp = temp + 1

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
                Dx = Dx - boxx*nint(Dx/boxx)
                Dy = Dy - boxy*nint(Dy/boxy)
                Dz = Dz - boxz*nint(Dz/boxz)

c                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
c                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
c                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)



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


                 if(cox.gt.boxx) cox = cox -boxx
                 if(cox.lt.0.0d0) cox = cox +boxx
                 if(coy.gt.boxy) coy = coy -boxy
                 if(coy.lt.0.0d0) coy = coy +boxy
                 if(coz.gt.boxz) coz = coz -boxz
                 if(coz.lt.0.0d0) coz = coz +boxz
c        cox = cox - boxx* dble(Idint((cox/boxx) +9999.5d0) - 9999)
c        coy = coy - boxy* dble(Idint((coy/boxy) +9999.5d0) - 9999)
c        coz = coz - boxz* dble(Idint((coz/boxz) +9999.5d0) - 9999)

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
        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        write(io3,*)boxx,boxy,boxz
        else
        write(io3,*)
        endif
                  do k = 1,total_beads
        write(io3,*)bead_lbl(k),bead_x(k),
     & bead_y(k),bead_z(k)

                  enddo

                  

           enddo

           write(io7,*)'sequence of beads in the molecule'
           do i = 1,num

           write(io7,*)'molecule',i,'with beads',nbd(i)
           do j = 1,nbd(i)
           write(io7,*)label(i,j)
           enddo

           enddo

           
            write(io7,*)'cg mapping completed'           
           

           close(io2)
         
           close(io3)

         end 