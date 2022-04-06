        subroutine form_total(data1,data2,n_bins)
        implicit none
         
        include 'system.inc'
        include 'config.inc'
        include 'control.inc'
        include 'iochan.h'
        include 'model.inc'


        double precision multiplier,total_factor
        double precision data1(maxtyp,maxtyp,maxbins),
     & data2(maxbins)

        integer i,j,k,n_bins

         


            do i = 1,bead_types

                do j = 1,bead_types


        if(i .eq.j) then
        multiplier = 1.0d0
        else
        multiplier = 2.0d0
        endif


        if(i.lt.j)cycle


         total_factor = (real(nr(i)*nr(j))/real(natoms)**2)*multiplier

              do k = 1,n_bins

               data2(k) = data2(k) + data1(i,j,k)*total_factor

               enddo


               enddo

                enddo


           return
          end
