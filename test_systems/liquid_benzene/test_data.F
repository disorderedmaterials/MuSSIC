        program test_data
        implicit none

        double precision   errmax,err,binw1(1000),binw2(1000),
     & data1(1000),data2(1000),tempr

        integer  wbins,i,exit_1,exit_0


        open(unit=1,file='H-weighted-total.sq', status='unknown')

        open(unit=2,file='C6H6-weighted-totalsq.dat',status = 'unknown')

         wbins = 600
         errmax = 0.0d0
         exit_1 = 1
         exit_0 = 0


         do i = 1,wbins

         read(1,*)binw1(i),data1(i)
         read(2,*)binw2(i),data2(i)

         if(abs(data1(i)).gt.0.1d0) then

         if(NINT(binw1(i)).eq.NINT(binw2(i))) then

         err = ((data1(i)-data2(i))/data1(i))*100
c         write(*,*)i,data1(i),data2(i),abs(err)
         endif
          
         if (err.gt.errmax) errmax = err

         endif

         enddo

           if(abs(errmax).gt.30.0d0) then
             write(*,*)exit_1
           else
             write(*,*)exit_0
           endif


           stop
           end program




