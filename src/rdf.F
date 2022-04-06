        subroutine rdf
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

        implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'model.inc'

      integer*8 temp,temp1,temp2,ii,jj,kk,i,j,k,l,bin,bl,el,
     &typeJ,typeK
      double precision Dx,Dy,Dz,distance


c   self terms

c       write(*,*)'rdf calculation started'
c       write(*,*)'self terms'

      do i = 1,ntypes



         temp = nr(i)

c$omp parallel 
c$omp& shared ( i,rr, temp, boxx, boxy, boxz, bin_width, nbins ) 
c$omp& private ( j,k,Dx,Dy,Dz,distance,bin ) 

c$omp do reduction ( + : histogram )
             do j = 1,temp-1


                       do k = j+1, temp

                 Dx = rr(i,k,1) - rr(i,j,1)
                 Dy = rr(i,k,2) - rr(i,j,2)
                 Dz = rr(i,k,3) - rr(i,j,3)


            
Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)



                distance = dsqrt(Dx**2+Dy**2+Dz**2)


                bin = ceiling(distance/bin_width)


               if(bin.le.nbins) then
               histogram(i,i,bin) = histogram(i,i,bin)+1.0d0
               endif


                      enddo

             enddo
      
      

        
c$omp end do

c$omp end parallel

       enddo

      
c  cross terms

c              write(*,*)'calculating cross terms'

            do ii = 1,ntypes-1

                   

                do jj = ii+1,ntypes




                   temp1 = nr(ii)
                   temp2 = nr(jj)


c$omp parallel 
c$omp& shared ( ii,jj,rr, temp, boxx, boxy, boxz, bin_width, nbins ) 
c$omp& private ( j,k,Dx,Dy,Dz,distance,bin ) 

c$omp do collapse(2)  reduction ( + : histogram )
         do j = 1,temp1


          do k = 1, temp2

                     Dx = rr(ii,j,1) - rr(jj,k,1)
                     Dy = rr(ii,j,2) - rr(jj,k,2)
                     Dz = rr(ii,j,3) - rr(jj,k,3)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc


                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)

               distance = dsqrt(Dx**2+Dy**2+Dz**2)


             bin = ceiling(distance/bin_width)

       if(bin.le.nbins)histogram(ii,jj,bin)=histogram(ii,jj,bin)+1.0d0

             enddo

             enddo


c$omp end do

c$omp end parallel

             enddo

             enddo


cccccccccccccccccccccccccccccccccccccccccccccccc
c    histogram from intramolecular pairs
cccccccccccccccccccccccccccccccccccccccccccccccc

c            write(*,*)'intra molecular pairs'


             do i = 1,num


c$omp parallel 
c$omp& shared ( i,nmol,nat,itype,x,y,z,boxx,boxy,boxz,bin_width,nbins) 
c$omp& private ( j,bl,el,jj,kk,typeJ,typeK,Dx,Dy,Dz,distance,bin ) 

c$omp do  reduction ( + : bound_hist )
                  do j = 1,nmol(i)


                       bl = ((j-1)*nat(i))+1
                       el = bl+nat(i)

                      do jj = bl,el-1

                       typeJ = itype(jj)

                       do kk = jj+1,el

                       typeK = itype(kk)

                     Dx = x(kk) - x(jj)
                     Dy = y(kk) - y(jj)
                     Dz = z(kk) - z(jj)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc
                  

                Dx = Dx - boxx* dble(Idint((Dx/boxx) +9999.5d0) - 9999)
                Dy = Dy - boxy* dble(Idint((Dy/boxy) +9999.5d0) - 9999)
                Dz = Dz - boxz* dble(Idint((Dz/boxz) +9999.5d0) - 9999)


             distance = dsqrt(Dx**2+Dy**2+Dz**2)


             bin = ceiling(distance/bin_width)


          if(bin.le.nbins) then
          bound_hist(typeJ,typeK,bin)=bound_hist(typeJ,typeK,bin)+1.0d0
        
          endif

          enddo

          enddo


          enddo


c$omp end do

c$omp end parallel

          enddo


            do i = 1,ntypes

                do j = 1,ntypes


                                

                    do k = 1,nbins

               unbound_hist(i,j,k) = histogram(i,j,k)-bound_hist(i,j,k)
                 
                     enddo

                 enddo
            
             enddo

              end

