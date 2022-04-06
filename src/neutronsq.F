      subroutine neutronsq
c     start_doc
c     Name
c       neutronsq
c     Usage
c       call neutronsq
c     Description
c       
c      multiplies partial S(Q) with appropriate weighting factor 
c      Weighting factors are computed according to the isotopologue
c      The bound partial S(Q) are multiplied by those weights seperately
c      and then combined with unbound partial S(Q) to form total
c      weighted  S(Q)
c
c     end_doc

      implicit none
      include 'system.inc'
      include 'model.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'

           integer*8 i,j,k

       double precision multiplier,norm,data_1(wbins),data_2(nbins),
     & bcd(ntypes),bj,bk,print_weight(ntypes,ntypes),
     & print_intra(ntypes,ntypes)

           character(len=10) file_id,label1,label2
           character(len=100) file_name,file_name1,file_name2
           integer*8 fu,fu1,fu2,deu_step,temp,typej,typek,file_cnt
            double precision wt_factor, total_num

c   concentration of the atom types            

          do i = 1,ntypes
          c(i) = dble(nr(i))/dble(natoms)
          enddo

c    opens the files according to the isotopologue name

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

                 print_weight(i,j) = 0.0d0
                 print_intra(i,j) = 0.0d0
                 intra_weight(i,j)= 0.0d0
                 intra_norm(i,j) = 0.0d0

                 enddo

                 enddo

c    compute the scattering length density of the atom type from isotope ratio              


                 do i = 1,ntypes

            if(natural(deu_step,i)+isotope(deu_step,i).eq.1.0d0) then
      bcd(i)=natural(deu_step,i)*bc(i)+isotope(deu_step,i)*bc(ntypes+1)
            else
        write(io7,*)'natural+isotope conc. not equal to 1.0 for sample',
     & deu_step,'and for the type', label_iso(deu_step,i)
            stop
            endif
                 enddo


c change the scattering length density if there are exchangeables                 

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

c  compute the weight for intra molecular contribution                     


           do i = 1,num


           if(deut_flag(deu_step,i).ne.2) then

             do j = 1,nat(i)

             typej = model_type(i,j)

             bj = bcd(typej)

             do k = j,nat(i)

             typek = model_type(i,k)

             bk = bcd(typek)

c              write(io7,*)'intra',text2(deu_step),j,k,bj,bk

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

c             write(io7,*)'intra',text2(deu_step),j,k,bj,bk
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


c             write(io7,*)'intra',text2(deu_step),j,k,bj,bk
             intra_weight(typej,typek) = intra_weight(typej,typek)+
     & (nmol(i)*iso_fraction(deu_step,i)*bj*bk)
             intra_norm(typej,typek) = intra_norm(typej,typek) +
     & (iso_fraction(deu_step,i)*nmol(i))


            enddo

            enddo


            endif

            enddo


c  normalize the intra_molecular weighting factor


            do i = 1,ntypes

            do j = 1,ntypes

            if(intra_flag(i,j).eq.1) then


            if(intra_norm(i,j).eq.0.0d0) cycle

         intra_weight(i,j) = intra_weight(i,j)/intra_norm(i,j)
         intra_weight(i,j) = intra_weight(i,j) *c(i)*c(j)

             endif

             enddo

             enddo


c   multipliy the partial S(Q)s with the weighting factors seperately

        write(io7,*)'ci,bi'

         do i = 1,ntypes

           do j = 1,ntypes

                 
                 if(i.eq.j)then
                         multiplier = 1.0d0
                 else
                         multiplier = 2.0d0
                  endif
               
                  print_weight(i,j) = print_weight(i,j) +
     & c(i)*c(j)*bcd(i)*bcd(j)*multiplier

                  print_intra(i,j) = print_intra(i,j) + 
     & intra_weight(i,j)*multiplier


           write(io7,*)text2(deu_step),i,j,c(i),c(j),bcd(i),bcd(j)

                  do k = 1,wbins

 

        d_weightedsq(i,j,k)=(unbound_partialsq(i,j,k)*c(i)*c(j)*bcd(i)
     &*bcd(j)+bound_partialsq(i,j,k)*intra_weight(i,j))*multiplier
        d_bound_wsq(i,j,k)=bound_partialsq(i,j,k)*intra_weight(i,j)
     & *multiplier

        d_unbound_wsq(i,j,k)=d_weightedsq(i,j,k)-d_bound_wsq(i,j,k)

        if(i.le.j) then
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

c  write the weighting matrix to log file     
c write the weighted partial S(Q) to a file              

                 file_cnt = 0

                 write(io7,*)'weight matrix(i,j)'

         do i = 1,ntypes
         do j = 1,ntypes


         write(io7,*)i,j,print_weight(i,j),print_intra(i,j)

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

        write(fu1,*)binw(i),d_totalweightedsq(i)

        enddo

c   back Fourier transform to obtain G(r) from total S(Q)        
              
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


