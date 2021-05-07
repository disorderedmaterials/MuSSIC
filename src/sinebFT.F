       subroutine sinebFT(data1,data2,normfactor,window)
        implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'

      double precision ft, deltaX,data1(wbins),data2(nbins),
     &normfactor,omega,oldX(wbins),oldY(wbins),param,broaden,
     &newX(nbins),newY(nbins),wind,wstep,omegamax

        integer m,nX,cnt,n_bins
        character*40 broadening
        character*40 window

              

                n_bins = wbins
                wstep = bin_width
                omega = binr(1)
                omegamax = binr(nbins-1)




                do m = 1,n_bins
                oldX(m) = binw(m)
                oldY(m) = data1(m)

                enddo

      
            param = PI/binw(wbins)



            cnt = 0

            do while( omega.le. omegamax)
            cnt = cnt+1
            ft = 0.0d0
            if (omega .gt. 0.0d0)then

            do m = 1,n_bins-1

                deltaX = oldx(m + 1) - oldx(m)




c apply window function

        if(window.eq.'Lorch0')then
        wind = dsin(oldx(m)*param)/(oldx(m)*param)
        
        elseif(window.eq.'none') then
        wind = 1.0d0
        endif
        

       ft = ft+dsin(oldx(m)*omega)*oldx(m)*oldy(m)*wind*deltaX

             enddo

             

c Normalise w.r.t. omega

            ft = ft/omega

             else


         do m = 1,n_bins-2

         deltaX = oldx(m + 1) - oldx(m)



c apply window function

        if(window.eq.'Lorch0')then
        wind = dsin(oldx(m) * param) / (oldx(m) * param)

        elseif(window.eq.'none') then
        wind = 1.0d0
        endif

             ft = ft+oldx(m)*oldy(m) *wind* deltaX

                enddo

        endif
c Add point
       
        data2(cnt) = ft*normfactor


        omega = omega + wstep


        enddo

                return
                end

