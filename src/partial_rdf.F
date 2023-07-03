        subroutine partial_rdf
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

      integer*8 i,j,k,file_cnt,fu
      character*10 label1,label2
      character*100 file_name

cccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c    normalize histogram values to from partial g(r)
c    sum all the partial g(r) to form total g(r) 


       call form_partial(histogram,partialgr)
       call form_partial(bound_hist,bound_partialgr)
       call form_partial(unbound_hist,unbound_partialgr)

       
       call form_total(partialgr,totalgr,nbins)



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
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
           open(unit = 14,file='unweighted-total.gr',status='unknown')

      write(fu,*)'# r(Angstroms)              partial g(r)         part
     &ial g(r) - intra         partial  g(r) - inter '


              do k = 1,nbins

       write(fu,*)binr(k),partialgr(i,j,k),bound_partialgr(i,j,k),
     & unbound_partialgr(i,j,k)

               enddo

        enddo
        enddo
ccccccccccccccccccccccccccccccccccccc
c   writes the total g(r) to the file
ccccccccccccccccccccccccccccccccccccc        

      write(14,*)'# r(Angstroms)            total g(r)'    
      do i = 1,nbins

      write(14,*)binr(i),totalgr(i)

      enddo

         

       end

