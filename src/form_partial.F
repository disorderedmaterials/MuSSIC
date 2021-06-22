        subroutine  form_partial(data1,data2)
           implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'

        double precision lower_shell_limit,delta,number_density,
     & multiplier,shell_volume,factor
        double precision data1(maxtyp,maxtyp,maxbins),
     & data2(maxtyp,maxtyp,maxbins)

        integer i,j,k


              do i = 1,ntypes

                do j = 1,ntypes

                     lower_shell_limit = 0.0d0

                     number_density = nr(j)/(boxx*boxy*boxz)

        if(i .eq.j) then
        multiplier = 2.0d0
        else
        multiplier = 1.0d0
        endif


c        if(i.lt.j)cycle

                    do k = 1,nbins


        shell_volume =(4.0d0/3.0d0)*PI*((lower_shell_limit+bin_width)**3
     & -lower_shell_limit**3)

        factor = nr(i)*shell_volume*number_density

        data2(i,j,k)=(data1(i,j,k)*(multiplier/factor))/nframes


        lower_shell_limit = lower_shell_limit + bin_width


                   enddo
                 enddo
            enddo

           return
          end
