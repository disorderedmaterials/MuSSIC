       subroutine sineFT(data1,data2,normfactor,broadening,window,mode)
        implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      character*20 mode

      double precision ft, deltaX,data1(maxbins),data2(maxbins),
     &normfactor,omega,oldX(maxbins),oldY(maxbins),param,broaden,
     &newX(maxbins),newY(maxbins),wind,wstep,omegamax

        integer m,nX,cnt,n_bins
        character*40 broadening
        character*40 window

              
        if(mode.eq.'forward')then

                n_bins = nbins-1
                wstep = qdelta
                omega = qmin
                omegamax = qmax

               do m = 1,n_bins
               oldX(m) = binr(m)
               oldY(m) = data1(m)

               enddo

        param= FWHM/(2.0d0 * dsqrt(2.0d0 * dlog(2.0d0)))

        elseif(mode.eq.'backward') then

                n_bins = wbins
                wstep = bin_width
                omega = binr(1)
                omegamax = binr(n_bins)


                do m = 1,n_bins
                oldX(m) = binw(m)
                oldY(m) = data1(m)
                enddo
      
        param = PI/binr(n_bins)

        endif

            cnt = 0

            do while( omega.le. omegamax)
            cnt = cnt+1
            ft = 0.0d0
            if (omega .gt. 0.0d0)then

            do m = 1,n_bins-1

                deltaX = oldx(m + 1) - oldx(m)


c Calculate broadening
        if(broadening.eq.'OmegaDependentGaussian')then

                if(mode.eq.'forward') then

      broaden=dexp(-(0.5d0*oldx(m)*oldx(m)*(param*omega)*(param*omega)))

                 elseif(mode.eq.'backward')then
                         
      broaden=dexp(-(oldx(m)*oldx(m))/(2.0d0*(param*omega)*
     &(param*omega)))

                endif
         
         elseif(broadening.eq.'none') then

               broaden = 1.0d0


       endif


c apply window function

        if(window.eq.'Lorch0')then
        wind = dsin(oldx(m) * param) / (oldx(m) * param)

        elseif(window.eq.'none') then
        wind = 1.0d0
        endif


         ft = ft + sin(x(m)*omega) * x(m) * y(m)*broaden*wind*deltaX

             enddo

c Normalise w.r.t. omega
            ft = ft/omega

             else


         do m = 1,n_bins-1

         deltaX = oldx(m + 1) - oldx(m)


c Calculate broadening
        if(broadening.eq.'OmegaDependentGaussian')then

                if(mode.eq.'forward') then

      broaden=dexp(-(0.5d0*oldx(m)*oldx(m)*(param*omega)*(param*omega)))

                 elseif(mode.eq.'backward')then

      broaden=dexp(-(oldx(m)*oldx(m))/(2.0d0*(param*omega)*
     &(param*omega)))

                endif

         elseif(broadening.eq.'none') then

               broaden = 1.0d0


       endif


c apply window function

        if(window.eq.'Lorch0')then
        wind = dsin(oldx(m) * param) / (oldx(m) * param)

        elseif(window.eq.'none') then
        wind = 1.0d0
        endif


                ft = ft + x(m) * y(m)*broaden *wind* deltaX

                enddo

        endif
c Add point
       
        data2(cnt) = ft*normfactor

        omega = omega + wstep

        enddo

                return
                end

