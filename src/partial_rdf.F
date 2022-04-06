        subroutine partial_rdf
c     start_doc
c     Name
c       partial_rdf
c     Usage
c       call partial_rdf
c     Description
c       Computes partial and total  g(r) from histogram by normalization
c       and summation
c       
c     end_doc

        implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'

      integer*8 i,j,k,l,ii,file_cnt,fu
      character*10 label1,label2
      character*100 file_name
      double precision xCentre,newY(nbins),param,broaden,xtemp,norm

c    normalize histogram values to from partial g(r)
c    sum all the partial g(r) to form total g(r)      

        call form_partial(histogram,partialgr)

       call form_partial(bound_hist,bound_partialgr)
       call form_partial(unbound_hist,unbound_partialgr)

       call form_total(partialgr,totalgr,nbins)


c  writes the partial g(r) values to the files named after the atom types         

        file_cnt = 0

       do i = 1,ntypes

          do j = 1,ntypes 

          if(i.gt.j) cycle

          label1 = label_type(i)
          label2 = label_type(j)

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

c   writes the total g(r) to the file

      do i = 1,nbins

      write(io4,*)binr(i),totalgr(i)

      enddo

         

       end

