        subroutine cg_mapping
        implicit none
        
      include 'system.inc'
      include 'iochan.h'
      include 'control.inc'
      include 'config.inc'
      include 'model.inc'


      integer*8 i,ii,j,jj,k,kk,iii,l,total_beads,ndxaux,sel,temp1,
     & temp,atmcnt,btyp,bl,el,bl_init,fu_1
         double precision cen_x(100), cen_y(100),cen_z(100),Dx,Dy,Dz,
     &cox,coy,coz,rmsd


         double precision cen_x1(100), cen_y1(100),cen_z1(100),Dx1,Dy1,
     &Dz1,cox1,coy1,coz1

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




ccccccccccccccccccccccccccccccccccccccccccc

             do i = 1,num
                do j = 1,nmol(i)

                bl = ((j-1)*nat(i))+1

                el = bl+nat(i)-1

                xpbc(bl) = x(bl)
                ypbc(bl) = y(bl)
                zpbc(bl) = z(bl)

                do jj = bl+1,el

                     Dx = x(jj) - xpbc(jj-1)
                     Dy = y(jj) - ypbc(jj-1)
                     Dz = z(jj) - zpbc(jj-1)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

c               Dx = Dx - boxx*nint(Dx/boxx)
c               Dy = Dy - boxy*nint(Dy/boxy)
c               Dz = Dz - boxz*nint(Dz/boxz)

                xpbc(jj) = x(jj) - boxx*nint(Dx/boxx)
                ypbc(jj) = y(jj) - boxy*nint(Dy/boxy)
                zpbc(jj) = z(jj) - boxz*nint(Dz/boxz)

                enddo

             enddo

             enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc             
c        write(io7,*)'cg_mapping: read atomistic trajectory'
        

           l = 0
           atmcnt = 0
           

           do ii = 1,num

              
           if(ii.gt.1) atmcnt = atmcnt + nmol(ii-1)*nat(ii-1)

                  
                 temp = 0


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

                                       
 
                      
                    cen_x(1) = xpbc(temp1)
                    cen_y(1) = ypbc(temp1)
                    cen_z(1) = zpbc(temp1)

                    cen_x1(1) = x(temp1)
                    cen_y1(1) = y(temp1)
                    cen_z1(1) = z(temp1)


                     cox = cen_x(1)
                     coy = cen_y(1)
                     coz = cen_z(1)

                     
                     cox1 = cen_x1(1)
                     coy1 = cen_y1(1)
                     coz1 = cen_z1(1)



                  do kk = 2,list_size(ii,i,k)

                  sel =  ndxaux + atm_list(ii,i,k,kk)

                                    
                 Dx = xpbc(sel) - cen_x(kk-1)
                 Dy = ypbc(sel) - cen_y(kk-1)
                 Dz = zpbc(sel) - cen_z(kk-1)

                 Dx1 = x(sel) - cen_x1(kk-1)
                 Dy1 = y(sel) - cen_y1(kk-1)
                 Dz1 = z(sel) - cen_z1(kk-1)
Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

c                Dx = Dx - boxx*nint(Dx/boxx)
c                Dy = Dy - boxy*nint(Dy/boxy)
c                Dz = Dz - boxz*nint(Dz/boxz)

           Dx1 = Dx1 - boxx* dble(Idint((Dx1/boxx) +9999.5d0) - 9999)
           Dy1 = Dy1 - boxy* dble(Idint((Dy1/boxy) +9999.5d0) - 9999)
           Dz1 = Dz1 - boxz* dble(Idint((Dz1/boxz) +9999.5d0) - 9999)


                 cen_x(kk) = cen_x(kk-1) + Dx
                 cen_y(kk) = cen_y(kk-1) + Dy
                 cen_z(kk) = cen_z(kk-1) + Dz


                 cen_x1(kk) = cen_x1(kk-1) + Dx1
                 cen_y1(kk) = cen_y1(kk-1) + Dy1
                 cen_z1(kk) = cen_z1(kk-1) + Dz1


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

                 cox1 = cox1 + cen_x1(kk)
                 coy1 = coy1 + cen_y1(kk)
                 coz1 = coz1 + cen_z1(kk)
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


                 cox1 = cox1/dble(list_size(ii,i,k))
                 coy1 = coy1/dble(list_size(ii,i,k))
                 coz1 = coz1/dble(list_size(ii,i,k))


                 rmsd = 0.0d0

                   if(list_size(ii,i,k).gt.1) then

                 do kk = 1,list_size(ii,i,k)

          rmsd = rmsd + (cen_x1(kk)-cox1)**2+(cen_y1(kk)-coy1)**2+
     & (cen_z1(kk)-coz1)**2 

           

c          write(*,*) rmsd,cen_x(kk),cen_y(kk),cen_z(kk),cox,coy,coz
              
                  enddo

                   else

                   rmsd = rmsd+ 0.5d0

                  endif

                  
                  rmsd = rmsd/dble(list_size(ii,i,k))


                  do kk = 1,bead_types
                  if(bead_lbl(l).eq.label_bead(kk)) btyp = kk
                  enddo 

                  

                  sigma_rmsd(1,btyp) = sigma_rmsd(1,btyp) + dsqrt(rmsd) 

                  sigma_rmsd_cnt(btyp) = sigma_rmsd_cnt(btyp) + 1.0d0


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

        cox1 = cox1 - boxx* dble(Idint((cox1/boxx) +9999.5d0) - 9999)
        coy1 = coy1 - boxy* dble(Idint((coy1/boxy) +9999.5d0) - 9999)
        coz1 = coz1 - boxz* dble(Idint((coz1/boxz) +9999.5d0) - 9999)

                  bead_x(l) = cox1
                  bead_y(l) = coy1
                  bead_z(l) = coz1


c                 do kk = 1,list_size(ii,i,k)-1

c                  bondx(kk) = (cen_x(kk)+cen_x(kk+1))/2.0d0
c                  bondy(kk) = (cen_y(kk)+cen_y(kk+1))/2.0d0
c                  bondz(kk) = (cen_z(kk)+cen_z(kk+1))/2.0d0

c                  box = box + bondx(kk)
c                  boy = boy + bondy(kk)
c                  boz = boz + bondz(kk)

c                  enddo



c                  bead_x(l) = box
c                  bead_y(l) = boy
c                  bead_z(l) = boz
                              
                 enddo

                 enddo

                 enddo

                 enddo




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

           do i = 1,bead_types

           sigma_rmsd(1,i) = sigma_rmsd(1,i)/sigma_rmsd_cnt(i)

           write(io7,*)'radius of the bead from atomic distribution'
     &,i,sigma_rmsd(1,i)
           
           enddo
           
           
            write(io7,*)'cg mapping completed'           

          

           close(io2)
           close(io3)

         end 
