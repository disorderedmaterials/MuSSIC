        subroutine rdf_cg
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       rdf_cg
c     Usage
c       call rdf_cg
c     Description
c       performs histogram binning for rdf calculation
c
c   histogram(bead types,bead types, nbins) = self terms + cross terms
c       bound_hist(bead types, atom types, nbins)
c       unbound_hist(bead types, atom types, nbins)
c
c     end_doc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'model.inc'

      integer*8 temp,temp1,temp2,ii,jj,kk,i,j,k,bin,bl,el,
     &typeJ,typeK,bl_init
      double precision Dx,Dy,Dz,distance
      double precision hist_local(bead_types,bead_types,max_binr),
     &bound_hist_local(bead_types,bead_types,max_binr)


         do i = 1,bead_types
           do j = 1,bead_types
            do k = 1,max_binr

            

        hist_local(i,j,k) = 0.0d0
        bound_hist_local(i,j,k) = 0.0d0

             enddo

            enddo

          enddo

ccccccccccccccc
c   self terms
ccccccccccccccc
c       write(io7,*)'rdf calculation'
c       write(io7,*)'self terms'

      do i = 1,bead_types



         temp = nr(i)


c$omp parallel
c$omp& shared ( i,rr, temp, boxx, boxy, boxz, bin_width, nbins )
c$omp& shared (max_binr,cutoff)
c$omp& private ( j,k,Dx,Dy,Dz,distance,bin )

c$omp do reduction ( + : hist_local )

             do j = 1,temp-1


                       do k = j+1, temp


                 Dx = rr(i,k,1) - rr(i,j,1)
                 Dy = rr(i,k,2) - rr(i,j,2)
                 Dz = rr(i,k,3) - rr(i,j,3)



Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc
                Dx = Dx - boxx*nint(Dx/boxx)
                Dy = Dy - boxy*nint(Dy/boxy)
                Dz = Dz - boxz*nint(Dz/boxz)

c                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
c                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
c                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)



                distance = dsqrt(Dx**2+Dy**2+Dz**2)


                bin = ceiling(distance/bin_width)


               if(bin.le.max_binr.and.distance.le.cutoff) then
               hist_local(i,i,bin) = hist_local(i,i,bin)+1.0d0
               endif


                      enddo

             enddo
c$omp end do

c$omp end parallel
      
        enddo

        
ccccccccccccccc
c  cross terms
ccccccccccccccc

c          write(io7,*)'calculating cross terms'

            do ii = 1,bead_types-1

                   

                do jj = ii+1,bead_types


                  if(ii.eq.jj) cycle



                   temp1 = nr(ii)
                   temp2 = nr(jj)


c$omp parallel
c$omp& shared ( ii,jj,rr, temp, boxx, boxy, boxz, bin_width, nbins)
c$omp& shared (max_binr,cutoff)
c$omp& private ( j,k,Dx,Dy,Dz,distance,bin )

c$omp do collapse(2)  reduction ( + : hist_local )

         do j = 1,temp1


          do k = 1, temp2

                     Dx = rr(ii,j,1) - rr(jj,k,1)
                     Dy = rr(ii,j,2) - rr(jj,k,2)
                     Dz = rr(ii,j,3) - rr(jj,k,3)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc
                Dx = Dx - boxx*nint(Dx/boxx)
                Dy = Dy - boxy*nint(Dy/boxy)
                Dz = Dz - boxz*nint(Dz/boxz)



c                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
c                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
c                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)

               distance = dsqrt(Dx**2+Dy**2+Dz**2)

             bin = ceiling(distance/bin_width)


       if(bin.le.max_binr.and.distance.le.cutoff)
     &hist_local(ii,jj,bin)=hist_local(ii,jj,bin)+1.0d0

             enddo

             enddo

c$omp end do

c$omp end parallel


             enddo

             enddo




cccccccccccccccccccccccccccccccccccccccccccccccc
c    histogram from intramolecular pairs
cccccccccccccccccccccccccccccccccccccccccccccccc

             bl_init = 0

             do i = 1,num

               if(i.le.1) then

                bl_init = 0

                else

                bl_init = bl_init + nmol(i-1)*nbd(i-1)

                endif


             do j = 1,nmol(i)

                bl = ((j-1)*nbd(i))+1+bl_init

                el = bl+nbd(i)-1
 

                xcgpbc(bl) = bead_x(bl)
                ycgpbc(bl) = bead_y(bl)
                zcgpbc(bl) = bead_z(bl)

                do jj = bl+1,el



                     Dx = bead_x(jj) - xcgpbc(jj-1)
                     Dy = bead_y(jj) - ycgpbc(jj-1)
                     Dz = bead_z(jj) - zcgpbc(jj-1)


Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc


                Dx = Dx - boxx*nint(Dx/boxx)
                Dy = Dy - boxy*nint(Dy/boxy)
                Dz = Dz - boxz*nint(Dz/boxz)

                xcgpbc(jj) = xcgpbc(jj-1) +Dx
                ycgpbc(jj) = ycgpbc(jj-1) +Dy
                zcgpbc(jj) = zcgpbc(jj-1) +Dz

c      write(io7,*)jj,xcgpbc(jj),ycgpbc(jj),zcgpbc(jj),
c     & dsqrt(Dx**2+Dy**2+Dz**2)

                enddo

                enddo

                enddo


c            write(io7,*)'intra molecular pairs'

             bl_init = 0

             do i = 1,num

             if(i.le.1) then

                bl_init = 0

                else

                bl_init = bl_init + nmol(i-1)*nbd(i-1)

                endif

                if(nbd(i).gt.1) then

c$omp parallel
c$omp& shared ( i,nmol,nat,itype,boxx,boxy,boxz,bin_width,nbins)
c$omp& shared (xcgpbc,ycgpbc,zcgpbc,max_binr,cutoff,bl_init)
c$omp& private ( j,bl,el,jj,kk,typeJ,typeK,Dx,Dy,Dz,distance,bin )

c$omp do  reduction ( + : bound_hist_local )
            
                  do j = 1,nmol(i)

                       bl = ((j-1)*nbd(i))+1+bl_init
                       el = bl+nbd(i)-1


c                       write(io7,*)i,nmol(i),bl,el

                       

                      do jj = bl,el-1

                       typeJ = itype(jj)

                       do kk = jj+1,el

                       typeK = itype(kk)

c                     Dx = bead_x(kk) - bead_x(jj)
c                     Dy = bead_y(kk) - bead_y(jj)
c                     Dz = bead_z(kk) - bead_z(jj)

c                     Dx = xcgpbc(kk) - xcgpbc(jj)
c                     Dy = ycgpbc(kk) - ycgpbc(jj)
c                     Dz = zcgpbc(kk) - zcgpbc(jj)

                     Dx = xcgpbc(kk) - xcgpbc(jj)
                     Dy = ycgpbc(kk) - ycgpbc(jj)
                     Dz = zcgpbc(kk) - zcgpbc(jj)
Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc
                  

c                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
c                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
c                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)

c                Dx = Dx - boxx*nint(Dx/boxx)
c                Dy = Dy - boxy*nint(Dy/boxy)
c                Dz = Dz - boxz*nint(Dz/boxz)

             distance = dsqrt(Dx**2+Dy**2+Dz**2)


             bin = ceiling(distance/bin_width)


          if(distance.le.cutoff) then
          bound_hist_local(typeJ,typeK,bin)=
     & bound_hist_local(typeJ,typeK,bin)+1.0d0
        
          
           if(iradius.eq.0) then
          if(bead_pair(typeJ,typeK).eq.1.and.(kk.eq.jj+1)) then


c           write(io7,*)typeJ,typeK,distance,'bead size'

           sigma(1,typeJ) = sigma(1,typeJ) + (distance)
           sigma(2,typeJ) = sigma(2,typeJ) + (distance)**2


           sigma(1,typeK) = sigma(1,typeK) + (distance)
           sigma(2,typeK) = sigma(2,typeK) + (distance)**2

           sigma_cnt(typeJ) = sigma_cnt(typeJ) + 1.0d0
           sigma_cnt(typeK) = sigma_cnt(typeK) + 1.0d0


           endif

           endif

           endif

          enddo

          enddo

          enddo
c$omp end do

c$omp end parallel

          endif

          enddo
ccccccccccccccccccccccccccccccccc
c   normalize the histogram bins
ccccccccccccccccccccccccccccccccc          



            do i = 1,bead_types

                do j = 1,bead_types

                 
             if(i.gt.j) then

               do k =1,max_binr
               hist_local(i,j,k) = hist_local(j,i,k)
               bound_hist_local(i,j,k) = bound_hist_local(j,i,k)

               enddo

                    endif
               

                    do k = 1,max_binr

        histogram(i,j,k)=histogram(i,j,k)+ hist_local(i,j,k)/
     &((nr(i)*nr(j))/(boxx*boxy*boxz))

       bound_hist(i,j,k)=bound_hist(i,j,k)+ bound_hist_local(i,j,k)/
     &((nr(i)*nr(j))/(boxx*boxy*boxz))


       unbound_hist(i,j,k) = histogram(i,j,k)-bound_hist(i,j,k)

       
                 
                     enddo

                 enddo
            
             enddo

c             write(io7,*)'histogram for rdf calculation completed'


             
              end

