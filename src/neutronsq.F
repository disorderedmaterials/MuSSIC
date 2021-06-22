         subroutine neutronsq
           implicit none
      include 'system.inc'
      include 'model.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'

           integer i,j,k

       double precision multiplier,norm,data_1(wbins),data_2(nbins),
     & bcd(ntypes),bj,bk

           character(len=10) file_id,label1,label2
           character(len=100) file_name,file_name1,file_name2
           integer fu,fu1,fu2,deu_step,temp,typej,typek,file_cnt
            double precision wt_factor, total_num


          do i = 1,ntypes
          c(i) = dble(nr(i))/dble(natoms)
          enddo

               do deu_step = 1,num_iso

                 file_id = text2(deu_step)

                  


      file_name1= trim(adjustl(file_id)) // '-weighted-total'// '.sq'

      file_name2= trim(adjustl(file_id)) // '-weighted-total'// '.gr'

                 fu1 = 400+deu_step
                 fu2 = 500+deu_step


                open(unit = fu1,file=trim(file_name1),status='unknown')
                open(unit = fu2,file=trim(file_name2),status='unknown')


                 do i = 1,wbins
                 d_totalweightedsq(i) = 0.0d0
                 enddo

                 do i = 1,nbins
                 d_totalweightedgr(i) = 0.0d0
                 enddo

                 do i = 1,ntypes

                 do j = 1,ntypes

                 intra_weight(i,j)= 0.0d0
                 intra_norm(i,j) = 0.0d0

                 enddo

                 enddo


                 do i = 1,ntypes

            if(natural(deu_step,i)+isotope(deu_step,i).eq.1.0d0) then
      bcd(i)=natural(deu_step,i)*bc(i)+isotope(deu_step,i)*bc(ntypes+1)
            else
        write(io7,*)'natural+isotope conc. not equal to 1.0 for sample',
     & deu_step,'and for the type', label_iso(deu_step,i)
            stop
            endif
                 enddo


                 wt_factor = 0.0d0
                  total_num = 0.0d0

                   do i = 1,ntypes
                   if(ex(i).eq.1) then

                          wt_factor = wt_factor + dble(nr(i))*bcd(i)
                          total_num = total_num + dble(nr(i))

                    endif
                    enddo

                    do i = 1,ntypes
                     
                    if(ex(i).eq.1) bcd(i) = wt_factor/total_num

                    enddo
                  


           do i = 1,num


           if(deut_flag(deu_step,i).ne.2) then

             do j = 1,nat(i)

             typej = model_type(i,j)

             bj = bcd(typej)

             do k = j,nat(i)

             typek = model_type(i,k)

             bk = bcd(typek)


             if(intra_flag(typej,typek).eq.0) cycle


             intra_weight(typej,typek) = intra_weight(typej,typek)+ 
     & (nmol(i)*iso_fraction(deu_step,i)*bj*bk)
             intra_norm(typej,typek) = intra_norm(typej,typek) + 
     & (iso_fraction(deu_step,i)*nmol(i))


            enddo

            enddo

            else


             do j = 1,nat(i)

             typej = model_type(i,j)

             if(ex(typej).eq.1) then
             bj = bcd(typej)
               else
              bj = bc(typej)
              endif
            

             do k = j,nat(i)

             typek = model_type(i,k)

             if(ex(typek).eq.1) then

             bk = bcd(typek)
             else
             bk = bc(typek)
             endif



             if(intra_flag(typej,typek).eq.0) cycle


             intra_weight(typej,typek) = intra_weight(typej,typek)+
     & (nmol(i)*(1.0d0-iso_fraction(deu_step,i))*bj*bk)
             intra_norm(typej,typek) = intra_norm(typej,typek) +
     & ((1.0d0-iso_fraction(deu_step,i))*nmol(i))


            enddo

            enddo

             do j = 1,nat(i)

             typej = model_type(i,j)
             
             if(ex(typej).eq.1) then
             
              bj = bcd(typej)
               else

              if(isotope(deu_step,typej).gt.0.0d0) then
              bj = bc(ntypes+1)
              else
              bj = bc(typej)
              endif
              endif


             do k = j,nat(i)

             typek = model_type(i,k)

             if(ex(typek).eq.1) then

             bk = bcd(typek)
             else
             if(isotope(deu_step,typek).gt.0.0d0)then
             bk = bc(ntypes+1)
             else
             bk = bc(typek)
             endif
             endif



             if(intra_flag(typej,typek).eq.0) cycle


             intra_weight(typej,typek) = intra_weight(typej,typek)+
     & (nmol(i)*iso_fraction(deu_step,i)*bj*bk)
             intra_norm(typej,typek) = intra_norm(typej,typek) +
     & (iso_fraction(deu_step,i)*nmol(i))


            enddo

            enddo


            endif

            enddo


            do i = 1,ntypes

            do j = 1,ntypes

            if(intra_flag(i,j).eq.1) then


            if(intra_norm(i,j).eq.0.0d0) cycle

         intra_weight(i,j) = intra_weight(i,j)/intra_norm(i,j)
         intra_weight(i,j) = intra_weight(i,j) *c(i)*c(j)

             endif

             enddo

             enddo




         do i = 1,ntypes

           do j = 1,ntypes

                 
                 if(i.eq.j)then
                         multiplier = 1.0d0
                 else
                         multiplier = 2.0d0
                  endif
               
                  do k = 1,wbins

        d_weightedsq(i,j,k)=(unbound_partialsq(i,j,k)*c(i)*c(j)*bcd(i)
     &*bcd(j)+bound_partialsq(i,j,k)*intra_weight(i,j))*multiplier
        d_bound_wsq(i,j,k)=bound_partialsq(i,j,k)*intra_weight(i,j)
     & *multiplier

        d_unbound_wsq(i,j,k)=d_weightedsq(i,j,k)-d_bound_wsq(i,j,k)

        if(i.ge.j) then
       d_totalweightedsq(k) = d_totalweightedsq(k) + d_weightedsq(i,j,k)
         endif
                  enddo



                  do k = 1,nbins


      d_weightedgr(i,j,k) = (partialgr(i,j,k)-1.0d0)*c(i)*c(j)
     &*bcd(i)*bcd(j)
    

      d_totalweightedgr(k) = d_totalweightedgr(k) + d_weightedgr(i,j,k)
                   enddo
       

                  enddo

                  enddo

                 file_cnt = 0

         do i = 1,ntypes
         do j = 1,ntypes


          if(i.gt.j) cycle

                 label1 = label_type(i)
                 label2 = label_type(j)

                 file_id = text2(deu_step)

                  
      file_name = trim(adjustl(file_id)) // '-weighted-' //
     &trim(adjustl(label1)) // '-' // trim(adjustl(label2)) // '.sq'

                file_cnt = file_cnt + 1

                fu = 300+file_cnt


                open(unit = fu,file=trim(file_name),status='unknown')


              do k = 1,wbins

       write(fu,*)binw(k),d_weightedsq(i,j,k),d_bound_wsq(i,j,k),
     & d_unbound_wsq(i,j,k)

               enddo
        enddo
        enddo

        do i = 1,wbins

        write(fu1,*)log(binw(i)),d_totalweightedsq(i)

        enddo
              
           norm = 1.0d0/(2.0d0*PI*PI*(dble(natoms)/(boxx*boxy*boxz)))


           do k = 1,wbins

          data_1(k) = d_totalweightedsq(k)
         

          enddo


         call sinebFT(data_1,data_2,norm,windowfn)

       

         do k = 1,nbins
           
            d_Gr(k) = data_2(k)

         enddo


         do i = 1,nbins-1

         write(fu2,*)binr(i),d_Gr(i)

         enddo

          enddo

             end


