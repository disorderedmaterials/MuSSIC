        subroutine sq
        implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'

      integer i,j,k,file_cnt,fu
      character*10 label1,label2
      character*100 file_name
      double precision data_1(nbins),data_2(nbins),data_3(nbins),
     & data_4(wbins),data_5(wbins),data_6(wbins),norm


        norm  = 4.0d0*PI*(dble(natoms)/(boxx*boxy*boxz))

            do i = 1,ntypes
               do j = 1,ntypes

c              if(nr(i).gt.nr(j)) cycle

c                 if((nr(i).eq.nr(j)).and.(i.gt.j)) cycle


                  do k = 1,nbins-1

         partialgr(i,j,k) = partialgr(i,j,k) - 1.0d0
         unbound_partialgr(i,j,k) = unbound_partialgr(i,j,k) - 1.0d0

             data_1(k) = partialgr(i,j,k)
             data_2(k) = bound_partialgr(i,j,k)
             data_3(k) = unbound_partialgr(i,j,k)


                  enddo


         call sinefFT(data_1,data_4,norm,broadfn,windowfn)
         call sinefFT(data_2,data_5,norm,broadfn,windowfn)
         call sinefFT(data_3,data_6,norm,broadfn,windowfn)


                 do k = 1,wbins

                 partialsq(i,j,k) = data_4(k)
                 bound_partialsq(i,j,k) = data_5(k)
                 unbound_partialsq(i,j,k) = data_6(k)

                 enddo

                  enddo
                  enddo




       call form_total(partialsq,totalsq,wbins)


        file_cnt = 0

       do i = 1,ntypes

          do j = 1,ntypes

          if(i.gt.j) cycle

          label1 = label_type(i)
          label2 = label_type(j)

       file_name = 'unweighted-' // trim(adjustl(label1)) //
     &  '-'  //  trim(adjustl(label2)) // '.sq'

                file_cnt = file_cnt + 1

                fu = 200 + file_cnt


                open(unit = fu,file=trim(file_name),status='unknown')



              do k = 1,wbins

       write(fu,*)binw(k),partialsq(i,j,k),
     & bound_partialsq(i,j,k),
     & unbound_partialsq(i,j,k)

               enddo

        enddo
        enddo


      do i = 1,wbins

      write(io66,*)binw(i),totalsq(i)

      enddo

                  end

