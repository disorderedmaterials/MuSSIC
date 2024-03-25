        subroutine rdf
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                
c     start_doc
c     Name
c       rdf
c     Usage
c       call rdf
c     Description
c       performs histogram binning for rdf calculation
c
c   histogram(atom types,atom types, nbins) = self terms + cross terms 
c       bound_hist(atom types, atom types, nbins)
c       unbound_hist(atom types, atom types, nbins)
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

      double precision hist_local(ntypes,ntypes,max_binr),
     &bound_hist_local(ntypes,ntypes,max_binr),box_vol
c      double precision,dimension(:,:,:),allocatable :: hist_local
c      double precision,dimension(:,:,:),allocatable :: bound_hist_local

c      allocate(hist_local(ntypes,ntypes,max_binr))
c      allocate(bound_hist_local(ntypes,ntypes,max_binr))



         do i = 1,ntypes
           do j = 1,ntypes
            do k = 1,max_binr

        hist_local(i,j,k) = 0.0d0
        bound_hist_local(i,j,k) = 0.0d0

             enddo

            enddo

          enddo


cccccccccccccccccccccccccccccccccccc      
c   self terms in histogram binning


      do i = 1,ntypes

         temp = nr(i)

          
c$omp parallel
c$omp& shared (i,rr, temp, boxx, boxy, boxz, bin_width, nbins)
c$omp& shared (max_binr, cutoff)
c$omp& private (j,k,Dx,Dy,Dz,distance,bin )

c$omp do reduction ( + : hist_local)


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

c        write(io7,*)'histogram for self terms taken'

         
ccccccccccccccccccccccc      
c  cross terms
ccccccccccccccccccccccc

c   write(*,*)'calculating cross terms'

            do ii = 1,ntypes-1

                do jj = ii+1,ntypes


                   temp1 = nr(ii)
                   temp2 = nr(jj)

c$omp parallel 
c$omp& shared ( ii,jj,rr, temp, boxx, boxy, boxz, bin_width, nbins ) 
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

       if(bin.le.max_binr.and.distance.le.cutoff) then
          hist_local(ii,jj,bin)=hist_local(ii,jj,bin)+1.0d0

       endif

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

                bl_init = bl_init + nmol(i-1)*nat(i-1)

                endif


                do j = 1,nmol(i)

                bl = ((j-1)*nat(i))+1+bl_init

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


                Dx = Dx - boxx*nint(Dx/boxx)
                Dy = Dy - boxy*nint(Dy/boxy)
                Dz = Dz - boxz*nint(Dz/boxz)
 

                xpbc(jj) = xpbc(jj-1) +Dx
                ypbc(jj) = ypbc(jj-1) +Dy
                zpbc(jj) = zpbc(jj-1) +Dz


                enddo

             enddo

             enddo

                bl_init = 0

             do i = 1,num

               if(i.le.1) then

                bl_init = 0

                else

                bl_init = bl_init + nmol(i-1)*nat(i-1)

                endif

c$omp parallel
c$omp& shared ( i,nmol,nat,itype,boxx,boxy,boxz,bin_width,nbins)
c$omp& shared (xpbc,ypbc,zpbc,max_binr,cutoff,bl_init)
c$omp& private ( j,bl,el,jj,kk,typeJ,typeK,Dx,Dy,Dz,distance,bin )

c$omp do  reduction ( + : bound_hist_local )

                  do j = 1,nmol(i)

                       bl = ((j-1)*nat(i))+1+bl_init
                       el = bl+nat(i)-1

                      do jj = bl,el

                       typeJ = itype(jj)

                       do kk = bl,el

                            if(jj.ne.kk) then

                       typeK = itype(kk)

                     Dx = xpbc(kk) - xpbc(jj)
                     Dy = ypbc(kk) - ypbc(jj)
                     Dz = zpbc(kk) - zpbc(jj)


             distance = dsqrt(Dx**2+Dy**2+Dz**2)


             bin = ceiling(distance/bin_width)

          if(distance.le.cutoff) then
          bound_hist_local(typeJ,typeK,bin)=
     &  bound_hist_local(typeJ,typeK,bin)+1.0d0

          endif

          endif

          enddo

          enddo

           enddo


c$omp end do

c$omp end parallel

          enddo

cccccccccccccccccccccccccccccccccccccc
c normalization of the histogram bins
cccccccccccccccccccccccccccccccccccccc          

c        write(io7,*)'histogram for inter molecular terms taken'
                
            do i = 1,ntypes

                do j = 1,ntypes


         if(i.eq.j) then

         do k =1,max_binr
         bound_hist_local(i,j,k) = 0.5d0*bound_hist_local(i,j,k)
         enddo

         endif
                  if(i.gt.j) then

                          do k =1,max_binr
                      hist_local(i,j,k) = hist_local(j,i,k)
                      bound_hist_local(i,j,k) = bound_hist_local(j,i,k)
                          enddo

                    endif




                       

                    do k = 1,max_binr


       histogram(i,j,k)=histogram(i,j,k)+ (hist_local(i,j,k))
     &/((nr(i)*nr(j))/(boxx*boxy*boxz))

       bound_hist(i,j,k)=bound_hist(i,j,k)+(bound_hist_local(i,j,k))
     &/((nr(i)*nr(j))/(boxx*boxy*boxz))


        unbound_hist(i,j,k) = histogram(i,j,k)-bound_hist(i,j,k)
                 
                     enddo

                 enddo
            
             enddo



              end
