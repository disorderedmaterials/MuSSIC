        subroutine sq

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                
c     start_doc
c     Name
c       sq
c     Usage
c       call sq
c     Description
c       performs the Fourier transfom on partial g(r) to obtain partial S(Q)
c     end_doc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'

      integer*8 i,j,k,file_cnt,fu
      character*10 label1,label2
      character*100 file_name
      double precision data_1(nbins),data_2(nbins),data_3(nbins),
     & data_4(wbins+1),data_5(wbins+1),data_6(wbins+1),norm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  normalization factor      


        norm  = 4.0d0*PI*(dble(natoms)/(boxx*boxy*boxz))

            do i = 1,ntypes
               do j = 1,ntypes

                if(i.gt.j) cycle

                  do k = 1,nbins-1

         partialgr(i,j,k) = partialgr(i,j,k) - 1.0d0
         unbound_partialgr(i,j,k) = unbound_partialgr(i,j,k) - 1.0d0

             data_1(k) = partialgr(i,j,k)
             data_2(k) = bound_partialgr(i,j,k)
             data_3(k) = unbound_partialgr(i,j,k)


                  enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  performs Fourier transform of partial g(r) to compute partial S(Q)


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
ccccccccccccccccccccccccccccccccccccccccccccccccc
c summation of partial S(Q) to form total S(Q)


       call form_total(partialsq,totalsq,wbins,ntypes)
       call form_total(bound_partialsq,total_boundsq,wbins,ntypes)
       call form_total(unbound_partialsq,total_unboundsq,wbins,ntypes)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c writes the partial S(Q) to a file named after the atom types       


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

      write(fu,*)'# Q(1/Angstroms)     partial unweighted S(Q)  
     &partial unweighted S(Q)-intra  partial unweighted S(Q)-inter '


              do k = 1,wbins

       write(fu,*)binw(k),partialsq(i,j,k),
     & bound_partialsq(i,j,k),
     & unbound_partialsq(i,j,k)

               enddo

        enddo
        enddo

ccccccccccccccccccccccccccc
c write the values to file 
      write(io16,*)'# Q(1/Angstroms)        unweighted S(Q)-total   
     &unweighted S(Q)-intra         unweighted S(Q)-inter '

      do i = 1,wbins

       write(io16,*)binw(i),totalsq(i),total_boundsq(i),
     & total_unboundsq(i)

      enddo

                  end

