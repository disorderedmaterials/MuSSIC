        subroutine sq_cg
        implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       sq_cg
c     Usage
c       call sq_cg
c     Description
c       performs the Fourier transfom on partial g(r) to obtain partial S(Q)
c     end_doc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'model.inc'
      include 'rdf.inc'
      include 'sq.inc'

      integer*8 i,j,k,file_cnt,fu
      character*10 label1,label2
      character*100 file_name
      double precision data_1(max_binr),data_2(max_binr),
     &data_3(max_binr),data_4(maxbins),data_5(maxbins),data_6(maxbins),
     &norm

cccccccccccccccccccccccc
c   normalization factor 


        norm  = 4.0d0*PI*(dble(natoms)/(boxx*boxy*boxz))

  
            do i = 1,bead_types
               do j = 1,bead_types



                  do k = 1,max_binr-1

         partialgr(i,j,k) = partialgr(i,j,k) - 1.0d0
         unbound_partialgr(i,j,k) = unbound_partialgr(i,j,k) - 1.0d0

             data_1(k) = partialgr(i,j,k)
             data_2(k) = bound_partialgr(i,j,k)
             data_3(k) = unbound_partialgr(i,j,k)


                  enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Fourier transform of partial g(r)  to compute partial S(Q)                  

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


cccccccccccccccccccccccccccccccccccccccccccccccc
c   summing the partial S(Q) to get total F_CG(Q)                  

       call form_total_cg(partialsq,totalsq,wbins)
       call form_total_cg(bound_partialsq,total_boundsq,wbins)
       call form_total_cg(unbound_partialsq,total_unboundsq,wbins)

ccccccccccccccccccccccccccccccccccccccccccccccccc
c   write partial S(Q) and total F_CG(Q) to files       

        file_cnt = 0

       do i = 1,bead_types

          do j = 1,bead_types

          if(i.gt.j) cycle

          label1 = label_bead(i)
          label2 = label_bead(j)

       file_name = 'CG-unweighted-' // trim(adjustl(label1)) //
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

      write(io26,*)binw(i),totalsq(i),total_boundsq(i),
     & total_unboundsq(i)

      enddo

                  end

