       subroutine sinefFT(data1,data2,normfactor,broadening,window)
        implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
    

        double precision data1(max_binr),data2(wbins+1),normfactor
        double precision ft,deltaX,omega,oldX(max_binr),oldY(max_binr),
     &param,broaden,wind,wstep,omegamax

        integer*8 m,cnt,n_bins
        character*40 broadening
        character*40 window

              

                n_bins = max_binr-1
                wstep = qdelta
                omega = qmin
                omegamax = qmax

               do m = 1,n_bins
               oldX(m) = binr(m)
               oldY(m) = data1(m)

               enddo

               

    
        param= FWHM/(2.0d0 * dsqrt(2.0d0 * dlog(2.0d0)))

            cnt = 0

            do while( omega.le. omegamax)
            cnt = cnt+1
            ft = 0.0d0
            if (omega .gt. 0.0d0)then

            do m = 1,n_bins-1

                deltaX = oldx(m + 1) - oldx(m)


c Calculate broadening
        if(broadening.eq.'OmegaDependentGaussian')then


      broaden=dexp(-(0.5d0*oldx(m)*oldx(m)*(param*omega)*(param*omega)))

       

         elseif(broadening.eq.'none') then

               broaden = 1.0d0


       endif


c apply window function

        if(window.eq.'Lorch0')then
        wind = dsin(oldx(m) * param) / (oldx(m) * param)

        elseif(window.eq.'none') then
        wind = 1.0d0
        endif


       ft=ft+dsin(oldX(m)*omega)*oldX(m)*oldY(m)
     &*broaden*wind*deltaX

           
             enddo


c Normalise w.r.t. omega
            ft = ft/omega

             else


         do m = 1,n_bins-1

         deltaX = oldx(m + 1) - oldx(m)


c Calculate broadening
        if(broadening.eq.'OmegaDependentGaussian')then


      broaden=dexp(-(0.5d0*oldx(m)*oldx(m)*(param*omega)*(param*omega)))


         elseif(broadening.eq.'none') then

               broaden = 1.0d0


       endif


c apply window function

        if(window.eq.'Lorch0')then
        wind = dsin(oldx(m) * param) / (oldx(m) * param)

        elseif(window.eq.'none') then
        wind = 1.0d0
        endif


                ft = ft + (oldx(m) * oldy(m)*broaden *wind* deltaX)

                enddo

        endif
c Add point
       
        data2(cnt) = ft*normfactor

        if (bin_type.eq.'log_linear') then
        omega = omega + binq(cnt)
        elseif(bin_type.eq.'linear') then
         omega = omega + qdelta
        endif


        enddo

         

                return
                end

