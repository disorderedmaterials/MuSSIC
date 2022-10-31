        subroutine partial_rdf_cg
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             
        implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'model.inc'

      integer*8 i,j,k,file_cnt,fu
      character*10 label1,label2
      character*100 file_name

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    normalize histogram values to from partial g(r)
c    sum all the partial g(r) to form total g(r)


       call form_partial_cg(histogram,partialgr)

       call form_partial_cg(bound_hist,bound_partialgr)
       call form_partial_cg(unbound_hist,unbound_partialgr)

       call form_total_cg(partialgr,totalgr,max_binr)

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c writes the partial g(r) and total g(r) to the files        

        file_cnt = 0

       do i = 1,bead_types

          do j = 1,bead_types 

          if(i.gt.j) cycle

          label1 = label_bead(i)
          label2 = label_bead(j)

       file_name = 'CG-unweighted-' // trim(adjustl(label1)) // 
     &  '-'  //  trim(adjustl(label2)) // '.gr'

                file_cnt = file_cnt + 1

                fu = 100+ file_cnt


                open(unit = fu,file=trim(file_name),status='unknown')



              do k = 1,max_binr

       write(fu,*)binr(k),partialgr(i,j,k),bound_partialgr(i,j,k),
     & unbound_partialgr(i,j,k)

               enddo

        enddo
        enddo

      do i = 1,max_binr

      write(io24,*)binr(i),totalgr(i)

      enddo

         

       end

