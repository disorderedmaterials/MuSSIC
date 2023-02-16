         subroutine neutronsq_cg
           implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       neutronsq_cg
c     Usage
c       call neutronsq_cg
c     Description
c
c      multiplies partial S(Q) with appropriate weighting factor
c      Weighting factors are computed according to the isotopologue
c      The bound partial S(Q) are multiplied by those weights seperately
c      and then combined with unbound partial S(Q) to form total
c      weighted  S(Q)
c
c     end_doc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           
      include 'system.inc'
      include 'model.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'

           integer*8 i,j,k,int1,int2,ii

      double precision multiplier,norm,data_1(maxbins),data_2(max_binr),
     &bcd(ntypes),bj,bk,bc1,bc2,bc3,bc4,ratio1,ratio2,ratio,
     &mix_ratio(bead_types),wt_factor,cnt
     

           character(len=10) file_id,label1,label2
           character(len=100) file_name,file_name1,file_name2
           integer*8 fu,fu1,fu2,deu_step,typej,typek,file_cnt,
     & mix_int,mix_flag(bead_types)
            double precision total_num,temp(maxtyp),temp1(maxtyp),
     & temp2(maxtyp),self_term(wbins)

            avg_noa_bead = 0.0d0

cccccccccccccccccccccccccccccc
c   concentration of the beads
cccccccccccccccccccccccccccccc

          do i = 1,bead_types
          c(i) = dble(nr(i))/dble(natoms)
          write(io7,*)c(i),'c(i)',i
          avg_noa_bead = avg_noa_bead + noa_bead(i)*c(i)
          enddo

cccccccccccccccccccccccccccccccccccc
c  loop over number of isotopologues
cccccccccccccccccccccccccccccccccccc

               do deu_step = 1,num_iso

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   open the files according to the isotopologue name
ccccccccccccccccccccccccccccccccccccccccccccccccccccc               

                 file_id = text2(deu_step)


      file_name1= trim(adjustl(file_id)) // '-CG-weighted-total'// '.sq'

      file_name2= trim(adjustl(file_id)) // '-CG-weighted-total'// '.gr'

                 fu1 = 400+deu_step
                 fu2 = 500+deu_step


                open(unit = fu1,file=trim(file_name1),status='unknown')
                open(unit = fu2,file=trim(file_name2),status='unknown')

ccccccccccccccccccccccccccccccccccc
c    initialize the arrays to zero
ccccccccccccccccccccccccccccccccccc

                 do i = 1,wbins
                 d_totalweightedsq(i) = 0.0d0
                 d_totalboundsq(i) = 0.0d0
                 d_totalunboundsq(i) = 0.0d0
                 self_term(i) = 0.0d0
                 enddo

                 do i = 1,max_binr
                 d_totalweightedgr(i) = 0.0d0
                 enddo

                 do i = 1,bead_types

                 weight(i) = 0.0d0

                 do j = 1,bead_types

                 weight1(i,j) = 0.0d0
                 intra_weight(i,j)= 0.0d0
                 intra_norm(i,j) = 0.0d0

                 enddo

                 enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute the scattering length densities considering the isotpe ratio

            do i = 1,ntypes

         if(natural(deu_step,i)+isotope(deu_step,i).eq.1.0d0) then
       bcd(i)=natural(deu_step,i)*bc(i)+isotope(deu_step,i)*bc(ntypes+1)
            else
        write(io7,*)'natural+isotope conc. not equal to 1.0 for sample',
     & deu_step,'and for the type', label_iso(deu_step,i)
          stop
         endif
          enddo                 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c change the scattering length density for exchangeables


                 wt_factor = 0.0d0
                  total_num = 0.0d0

                   do i = 1,ntypes
                   if(ex(i).eq.1) then

                           cnt = 0.0d0 

                           do j = 1,bead_types

                           do k = 1,num_atmtyp(j)

            if(i.eq.atmtyp(j,k)) cnt =cnt+dble(nr(j))*num_atm_bead(j,k)

                           enddo
                           enddo 

                          wt_factor = wt_factor + cnt*bcd(i)
                          total_num = total_num + cnt

                    endif
                    enddo

                    write(io7,*)'scattering lengths'

                    do i = 1,ntypes

                    if(ex(i).eq.1) bcd(i) = wt_factor/total_num

                    write(io7,*)i,bcd(i),bc(i)

                    enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute the effective scattering length density for each bead
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                    
                    
              do i = 1,bead_types

                           temp(i) = 0.0d0
                           temp1(i) = 0.0d0
                           temp2(i) = 0.0d0

                     do k = 1,num_atmtyp(i)
                     int1 = atmtyp(i,k)
           
           temp(i) = temp(i) + num_atm_bead(i,k) *bcd(int1)
           temp1(i) = temp1(i) + num_atm_bead(i,k)
     & *(bc(int1))

            if(isotope(deu_step,int1).gt.0.0d0) then
        temp2(i) = temp2(i) + (num_atm_bead(i,k)*bc(ntypes+1))
       

     
            else
         temp2(i) = temp2(i) + (num_atm_bead(i,k)*bc(int1))
   
            endif


          enddo


          enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    add the flag to the exchangables and beads with isotope ratios          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                 do i = 1,bead_types

                   mix_flag(i) = 0
                   mix_ratio(i) = 1.0d0

                 do ii = 1,num_atmtyp(i)

                    int1 = atmtyp(i,ii)

       if((isotope(deu_step,int1).gt.0.0).and.(isotope(deu_step,int1)
     & .lt.1.0).and.(ex(int1).eq.0)) then
          mix_int = int1
          mix_flag(i) = 1
          mix_ratio(i) = natural(deu_step,int1)
      elseif((isotope(deu_step,int1).gt.0.0).and.(isotope(deu_step,int1)
     & .lt.1.0).and.(ex(int1).eq.1)) then
          mix_int = int1
          mix_flag(i) = 0
          mix_ratio(i) = natural(deu_step,int1)
        endif

         enddo

         enddo



ccccccccccccccccccccccccccccccccccccccccccc
c compute weights for self scattering term         
ccccccccccccccccccccccccccccccccccccccccccc

                do i = 1,bead_types


                if(mix_flag(i).eq.0) then
       


                 do j = 1,num_atmtyp(i)

                 int1 = atmtyp(i,j)
                
                   bc1 = bcd(int1)                    

                  do k = j,num_atmtyp(i)

                   int2 = atmtyp(i,k)
                   bc2 = bcd(int2)


                   if(k.eq.j) then
                           multiplier = 1.0d0
                   else
                           multiplier = 2.0d0
                   endif


                 
         weight(i) = weight(i) + (c(i)*
     & num_atm_bead(i,j)* num_atm_bead(i,k)*bc1*bc2*multiplier)



                 enddo

                 enddo


                  do j = 1, num_atmtyp(i)

                  int1 = atmtyp(i,j)

           bc1 = bcd(int1)

          weight(i) =weight(i) - (c(i)*
     & num_atm_bead(i,j)*bc1**2)



                  enddo


          else


                 do j = 1,num_atmtyp(i)

                 int1 = atmtyp(i,j)
                 bc1 = bcd(int1)
                 bc3 = bcd(int1)
                 ratio1 = 1.0

                 if(int1.eq.mix_int) then
                   bc1 = bc(int1)
                   bc3 = bc(ntypes+1)
                   ratio1 = mix_ratio(i)
                  endif
                   
                    
                  do k = j,num_atmtyp(i)

                   int2 = atmtyp(i,k)

                   bc2 = bcd(int2)
                   bc4 = bcd(int2)
                   ratio2 = 1.0

                   if(int2.eq.mix_int) then
                     bc2 = bc(int2)
                     bc4 = bc(ntypes+1)
                     ratio2 = mix_ratio(i)
                   endif



                   if(k.eq.j) then
                           multiplier = 1.0d0
                   else
                           multiplier = 2.0d0
                   endif


       if(ratio1.eq.mix_ratio(i).or.ratio2.eq.mix_ratio(i)) 
     & ratio = mix_ratio(i)
                 
         weight(i) = weight(i) + (ratio*c(i)*
     & num_atm_bead(i,j)* num_atm_bead(i,k)*bc1*bc2*multiplier)


         weight(i) = weight(i) + ((1.0-ratio)*c(i)*
     & num_atm_bead(i,j)* num_atm_bead(i,k)*bc3*bc4*multiplier)



                 enddo

                 enddo


                  do j = 1, num_atmtyp(i)

                  int1 = atmtyp(i,j)

           bc1 = bcd(int1)
           bc2 = bcd(int1)
           ratio = 1.0

           if(int1.eq.mix_int) then
                   bc1 = bc(int1)
                   bc2 = bc(ntypes+1)
                  ratio = mix_ratio(i)
          endif 

          weight(i) =weight(i) - (ratio*c(i)*
     & num_atm_bead(i,j)*bc1**2)

          

          weight(i) =weight(i) - ((1.0-ratio)*c(i)*
     & num_atm_bead(i,j)*bc2**2)


                  enddo

           endif

                             
             enddo

                      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute weights for the inter molecular scattering term                      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                
                      do i = 1,bead_types

                      do j = i,bead_types

                      if(i.eq.j) then
                      multiplier = 1.0d0
                      else
                      multiplier = 2.0d0
                       endif


c        if(mix_ratio(i).lt.1.0)then

c        term1 = (mix_ratio(i)*temp1(i)+(1.0-mix_ratio(i))*temp2(i))
c        endif


c        if(mix_ratio(j).lt.1.0)then
c
c        term2 = (mix_ratio(j)*temp1(j)+(1.0-mix_ratio(j))*temp2(j))
c        endif
        

         weight1(i,j)=weight1(i,j)+  temp(i)*temp(j)



   

                    enddo
                    enddo

ccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute weights for intra molecular scattering

            do i = 1,num

            if(deut_flag(deu_step,i).ne.2) then

            do j = 1,nbd(i)


            typej = model_type(i,j)

            
            do k = j,nbd(i)


            typek = model_type(i,k)

            
            if(intra_flag(typej,typek).eq.0) cycle

            if(typej.eq.typek) then
                    multiplier = 1.0d0
            else
                    multiplier = 2.0d0
             endif

            

      intra_weight(typej,typek) = intra_weight(typej,typek)+
     &(temp(typej)*temp(typek)*nmol(i))
       intra_norm(typej,typek) = intra_norm(typej,typek) +
     & (nmol(i))


        enddo

        enddo

        else


            do j = 1,nbd(i)


            typej = model_type(i,j)


            if(ex_cg(deu_step,typej).eq.1) then
             bj = temp(typej)
               else
             bj = temp1(typej)

               endif


               do k = j,nbd(i)


            typek = model_type(i,k)



            if(ex_cg(deu_step,typek).eq.1) then

                    bk = temp(typek)
             else
                     bk = temp1(typek)
             endif

             if(intra_flag(typej,typek).eq.0) cycle

             intra_weight(typej,typek) = intra_weight(typej,typek)+
     & (nmol(i)*(1.0d0-iso_fraction(deu_step,i))*bj*bk)
             intra_norm(typej,typek) = intra_norm(typej,typek) +
     & ((1.0d0-iso_fraction(deu_step,i))*nmol(i))

        
              enddo

              enddo

            do j = 1,nbd(i)

             typej = model_type(i,j)

             if(ex_cg(deu_step,typej).eq.1) then

                     bj = temp(typej)
               else

              if(isotope_cg(deu_step,typej).gt.0.0d0) then
                      bj = temp2(typej)
              else
                      bj = temp1(typej)
              endif
              endif


             do k = j,nbd(i)

             typek = model_type(i,k)

             if(ex_cg(deu_step,typek).eq.1) then

                     bk = temp(typek)
             else
             if(isotope_cg(deu_step,typek).gt.0.0d0)then
                     bk = temp2(typek)
             else
                     bk = temp1(typek)
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
         
                   
            
           write(io7,*)'INTRA_WEIGHT matrix'

           do i = 1,bead_types
           do j = 1,bead_types

           if (i.eq.j) then
              multiplier = 1.0d0
           else
              multiplier = 2.0d0
           endif

           if(intra_flag(i,j).eq.1) then


c      write(io7,*)'intra',deu_step,i,j,intra_weight(i,j),intra_norm(i,j)
c           if(intra_norm(i,j).eq.0.0d0) cycle

      intra_weight(i,j) = intra_weight(i,j)/intra_norm(i,j)
c       write(io7,*)deu_step,i,j,intra_weight(i,j)
c        intra_weight(i,j) = intra_weight(i,j)*multiplier *c(i)*c(j)

         if(i.le.j) then
        write(io7,*)deu_step,i,j,intra_weight(i,j)*c(i)*c(j)*multiplier
         endif

            endif

            enddo

            enddo


            
           write(io7,*)'INTER_WEIGHT matrix'

           do i = 1,bead_types
           do j = 1,bead_types

           

            if(i.eq.j) then
                    multiplier = 1.0d0
            else
                    multiplier = 2.0d0
              endif

              if(i.le.j) then

        write(io7,*)deu_step,i,j,weight1(i,j)*c(i)*c(j)*multiplier
               endif


            enddo

            enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c weighted sum of partial S(Q) using the weights calculated above            
         

         do i = 1,bead_types

        do k = 1,wbins

       d_totalweightedsq(k)=d_totalweightedsq(k)+((weight(i)*FF(i,k)**2)
     &*(1.0d0/avg_noa_bead))

       self_term(k) = self_term(k) + ((weight(i)*FF(i,k)**2)
     &*(1.0d0/avg_noa_bead))

 
c       d_totalunboundsq(k)=d_totalunboundsq(k)+((weight(i)*FF(i,k)**2)
c     &*(1.0d0/avg_noa_bead))
        enddo
        

   
           do j = i,bead_types

                 
            if(i.eq.j) then
                    multiplier = 1.0d0
            else
                    multiplier = 2.0d0

              endif
               
                  do k = 1,wbins



c      d_weightedsq(i,j,k)= (partialsq(i,j,k)*weight1(i,j))
c    &*c(i)*c(j)*multiplier*FF(i,k)*FF(j,k)*(1.0d0/avg_noa_bead)

       d_unbound_wsq(i,j,k)= (unbound_partialsq(i,j,k)*weight1(i,j))
     &*c(i)*c(j)*multiplier*FF(i,k)*FF(j,k)*(1.0d0/avg_noa_bead)

       d_bound_wsq(i,j,k)= (bound_partialsq(i,j,k)*intra_weight(i,j)
     &*c(i)*c(j)*multiplier)*FF(i,k)*FF(j,k)*(1.0d0/avg_noa_bead)

       d_weightedsq(i,j,k) = d_unbound_wsq(i,j,k) + d_bound_wsq(i,j,k) 
        
       d_totalweightedsq(k) = d_totalweightedsq(k)+d_bound_wsq(i,j,k)+
     & d_unbound_wsq(i,j,k)
       d_totalboundsq(k) = d_totalboundsq(k) + d_bound_wsq(i,j,k)
       d_totalunboundsq(k)=d_totalunboundsq(k)+ d_unbound_wsq(i,j,k)
         
                  enddo

                  enddo

                  enddo


cccccccccccccccccccccccccccccccccccccccccccccccc
c  write partial S(Q) and F(Q) to the files                  

                 file_cnt = 0

         do i = 1,bead_types
         do j = i,bead_types


                 label1 = label_bead(i)
                 label2 = label_bead(j)

                 file_id = text2(deu_step)

                  
      file_name = trim(adjustl(file_id)) // '-CG-weighted-' //
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

        write(fu1,*)binw(i),d_totalweightedsq(i),d_totalboundsq(i),
     & d_totalunboundsq(i),self_term(i)
c        write(153,*)deu_step,binw(i),test1(i)

        enddo
              
           norm = 1.0d0/(2.0d0*PI*PI*(dble(natoms)/(boxx*boxy*boxz)))


           do k = 1,wbins

          data_1(k) = d_totalweightedsq(k)
         

          enddo


c         call sinebFT(data_1,data_2,norm,windowfn)

       

         do k = 1,max_binr
           
            d_Gr(k) = data_2(k)

         enddo


         do i = 1,max_binr

         write(fu2,*)binr(i),d_Gr(i)

         enddo

          enddo

             end


