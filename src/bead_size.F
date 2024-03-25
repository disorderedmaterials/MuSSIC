      subroutine bead_size
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       bead_size
c     Usage
c       call bead_size
c     Description
c       Computes bead size  and  the form factor for all bead types
c
c     end_doc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              
      implicit none

      include 'system.inc'
      include 'iochan.h'
      include 'control.inc'
      include 'config.inc'
      include 'model.inc'
      include 'rdf.inc'

        integer*8 i,j
        double precision tempr(max_binr)

cccccccccccccccccccccccccccccccccccccccccccccc
c    normalize the bead size        


         do i = 1,bead_types

          if(iradius.eq.0) then

         if(bonded(i).ne.0) then

           sigma(1,i) = sigma(1,i)/sigma_cnt(i)

           sigma(1,i) = sigma(1,i)/2.0d0

c           write(*,*)i,sigma(1,i),sigma_cnt(i)

c         sigma(1,i) = sigma(1,i)/(2.0d0*sigma_cnt(i))
c         sigma(2,i) = sigma(2,i)/(2.0d0*sigma_cnt(i))
            sigma(2,i) = sigma(1,i)**2
         elseif(bonded(i).eq.0) then

c                 do j = 1,max_binr
c                 tempr(j) = partialgr(i,i,j)
c                 enddo

c          sigma(1,i) = binr(maxloc(tempr,1))/2.0d0
c          sigma(2,i) = sigma(1,i)**2

            sigma(1,i) = noa_bead(i)*0.2d0 

          endif

           endif

         
         if(run_type.eq.'compare_CG_atomistic') then

          sigma(1,i) = sigma_rmsd(1,i)

         endif

         std_dev(i) = dsqrt(sigma(2,i) - sigma(1,i)**2)

         write(io7,*)'bead radius for  ',label_bead(i),'is',sigma(1,i)

         enddo


         
ccccccccccccccccccccccccccccccccccccc
c   initialisation of form factor

c             write(io7,*)'Form factor'

            if(FF_label.eq.'Gaussian')then

             do i = 1, bead_types

c                write(io7,*)'bead radius',i,size_bead(i)

                 do j = 1, wbins

          FF(i,j) = dexp(-0.5d0*(binw(j)*(0.51d0)*sigma(1,i))**2.0d0)

                 enddo

             enddo

             elseif(FF_label.eq.'uniform') then

             do i = 1, bead_types

c                write(io7,*)'bead radius',i,size_bead(i)

                 do j = 1, wbins

          FF(i,j) = (3.0d0/(binw(j)*sigma(1,i))**3)*
     &(dsin(binw(j)*sigma(1,i))-(dcos(binw(j)*
     & sigma(1,i))*binw(j)*sigma(1,i)))

                 enddo

             enddo

            elseif(FF_label.eq.'no_FF')then

                   do i = 1,bead_types

                   do j = 1,wbins

                   FF(i,j) = 1.0d0

                   enddo

                   enddo

                endif

c                do i = 1,bead_types

c                do j = 1,wbins

c                write(999,*)i,j,FF(i,j)

c                enddo
c                enddo

      end


