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

        integer*8 i,j,k


              do i = 1,ntypes

                do j = 1,ntypes

                     lower_shell_limit = 0.0d0


      if(i .eq.j) then
        multiplier = 2.0d0
        else
        multiplier = 1.0d0
        endif


        if(i.gt.j)cycle

                    do k = 1,max_binr


        shell_volume =(4.0d0/3.0d0)*PI*((lower_shell_limit+bin_width)**3
     & -lower_shell_limit**3)


         if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then

     
        data2(i,j,k)=(data1(i,j,k)*(multiplier/shell_volume))/frame(k)

        else

        data2(i,j,k)=(data1(i,j,k)*(multiplier/shell_volume))/nframes

        endif


        lower_shell_limit = lower_shell_limit + bin_width



                   enddo
                 enddo
            enddo

           return
          end
