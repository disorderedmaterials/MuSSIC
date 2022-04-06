        subroutine partial_rdf
        implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'model.inc'

      integer i,j,k,file_cnt,fu
      character*10 label1,label2
      character*100 file_name


        call form_partial(histogram,partialgr)

       call form_partial(bound_hist,bound_partialgr)
       call form_partial(unbound_hist,unbound_partialgr)

       call form_total(partialgr,totalgr,nbins)

        file_cnt = 0

       do i = 1,bead_types

          do j = 1,bead_types 

          if(i.gt.j) cycle

          label1 = label_bead(i)
          label2 = label_bead(j)

       file_name = 'unweighted-' // trim(adjustl(label1)) // 
     &  '-'  //  trim(adjustl(label2)) // '.gr'

                file_cnt = file_cnt + 1

                fu = 100+ file_cnt


                open(unit = fu,file=trim(file_name),status='unknown')



              do k = 1,nbins

       write(fu,*)binr(k),partialgr(i,j,k),bound_partialgr(i,j,k),
     & unbound_partialgr(i,j,k)

               enddo

        enddo
        enddo

      do i = 1,nbins

      write(io4,*)binr(i),totalgr(i)

      enddo

         

       end

