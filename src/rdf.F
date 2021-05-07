        subroutine rdf
        implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'model.inc'

      integer temp,temp1,temp2,ii,jj,kk,i,j,k,l,bin,beginloop,endloop,
     &typeJ,typeK
      double precision Dx,Dy,Dz,distance


c   self terms


      do i = 1,ntypes

         temp = nr(i)

             do j = 1,temp

                       do k = j+1, temp

                 Dx = rr(i,k,1) - rr(i,j,1)
                 Dy = rr(i,k,2) - rr(i,j,2)
                 Dz = rr(i,k,3) - rr(i,j,3)


Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

                Dx = Dx - box* dble(Idint((Dx/box) +9999.5d0) - 9999)
                Dy = Dy - box* dble(Idint((Dy/box) +9999.5d0) - 9999)
                Dz = Dz - box* dble(Idint((Dz/box) +9999.5d0) - 9999)



                distance = dsqrt(Dx**2+Dy**2+Dz**2)


                bin = floor(distance/bin_width)


               if(bin.lt.nbins) then
               histogram(i,i,bin) = histogram(i,i,bin)+1.0d0
               endif

                      enddo

             enddo
      
        enddo

        

c  cross terms

            do ii = 1,ntypes

                do jj = 1,ntypes

                  if(ii.eq.jj) cycle

c                 if(nr(ii).gt.nr(jj)) cycle

c                 if((nr(ii).eq.nr(jj)).and.(ii.gt.jj)) cycle


                   temp1 = nr(ii)
                   temp2 = nr(jj)

         do j = 1,temp1

          do k = 1, temp2

                     Dx = rr(ii,j,1) - rr(jj,k,1)
                     Dy = rr(ii,j,2) - rr(jj,k,2)
                     Dz = rr(ii,j,3) - rr(jj,k,3)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc


                Dx = Dx - box* dble(Idint((Dx/box) +9999.5d0) - 9999)
                Dy = Dy - box* dble(Idint((Dy/box) +9999.5d0) - 9999)
                Dz = Dz - box* dble(Idint((Dz/box) +9999.5d0) - 9999)

               distance = dsqrt(Dx**2+Dy**2+Dz**2)


             bin = floor(distance/bin_width)

       if(bin.lt.nbins)histogram(ii,jj,bin)=histogram(ii,jj,bin)+1.0d0

             enddo

             enddo

             enddo

             enddo


cccccccccccccccccccccccccccccccccccccccccccccccc
c    histogram from intramolecular pairs
cccccccccccccccccccccccccccccccccccccccccccccccc

              l = 0
             do i = 1,num

                  do j = 1,nmol(i)

                      beginloop = l+1
                      endloop = l + nat(i) 

                      l = l + nat(i)

                      do jj = beginloop,endloop-1

                       typeJ = itype(jj)

                       do kk = jj+1,endloop

                       typeK = itype(kk)

                     Dx = x(kk) - x(jj)
                     Dy = y(kk) - y(jj)
                     Dz = z(kk) - z(jj)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc
                  

                Dx = Dx - box* dble(Idint((Dx/box) +9999.5d0) - 9999)
                Dy = Dy - box* dble(Idint((Dy/box) +9999.5d0) - 9999)
                Dz = Dz - box* dble(Idint((Dz/box) +9999.5d0) - 9999)


             distance = dsqrt(Dx**2+Dy**2+Dz**2)


             bin = floor(distance/bin_width)


          if(bin.lt.nbins) then
          bound_hist(typeJ,typeK,bin)=bound_hist(typeJ,typeK,bin)+1.0d0
        
          endif

          enddo

          enddo


          enddo

          enddo


            do i = 1,ntypes

                do j = 1,ntypes


                                

                    do k = 1,nbins

               unbound_hist(i,j,k) = histogram(i,j,k)-bound_hist(i,j,k)
                 
                     enddo

                 enddo
            
             enddo

              end

