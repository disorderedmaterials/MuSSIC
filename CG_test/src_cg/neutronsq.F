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

           integer i,j,k,int1,int2,ii,normal_factor

       double precision multiplier,norm,data_1(maxbins),data_2(nbins),
     &bcd(ntypes),bj,bk,bc1,bc2,bc3,bc4,ratio1,ratio2,ratio,
     &mix_ratio(bead_types),wt_factor,test1(wbins),test2(wbins)
     

           character(len=10) file_id,label1,label2
           character(len=100) file_name,file_name1,file_name2
           integer fu,fu1,fu2,deu_step,typej,typek,file_cnt,
     & mix_int,mix_flag(bead_types)
            double precision total_num,temp(maxtyp),temp1(maxtyp),
     & temp2(maxtyp)


          do i = 1,bead_types
          c(i) = dble(nr(i))/dble(natoms)
          write(io7,*)c(i),'c(i)',i
          enddo


                

               do deu_step = 1,num_iso


          do i = 1,wbins
          test1(i) = 0.0d0
          test2(i) =0.0d0
          enddo

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

                 do i = 1,bead_types

                 weight(i) = 0.0d0

                 do j = 1,bead_types

                 weight1(i,j) = 0.0d0
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


c                 write(io7,*)bcd(1),bcd(2)
                 wt_factor = 0.0d0
                  total_num = 0.0d0

                   do i = 1,bead_types
                   if(ex(i).eq.1) then

                          wt_factor = wt_factor + dble(nr(i))*bcd(i)
                          total_num = total_num + dble(nr(i))

                    endif
                    enddo

                    do i = 1,ntypes

                    if(ex(i).eq.1) bcd(i) = wt_factor/total_num

                    enddo




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




                 do i = 1,bead_types

                   mix_flag(i) = 0

                 do ii = 1,num_atmtyp(i)

                    int1 = atmtyp(i,ii)

       if((isotope(deu_step,int1).gt.0.0).and.(isotope(deu_step,int1)
     & .lt.1.0).and.(ex(int1).eq.0)) then
          mix_int = int1
          mix_flag(i) = 1
          mix_ratio(i) = natural(deu_step,int1)
           endif

         enddo


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
                  ratio = mix_int
          endif 

          weight(i) =weight(i) - (ratio*c(i)*
     & num_atm_bead(i,j)*bc1**2)

          

          weight(i) =weight(i) - ((1.0-ratio)*c(i)*
     & num_atm_bead(i,j)*bc2**2)


                  enddo

           endif

                             
             enddo

                      



                
                      do i = 1,bead_types

                      do j = i,bead_types

                      if(i.eq.j) then
                      multiplier = 1.0d0
                      else
                      multiplier = 2.0d0
                       endif

         weight1(i,j)=weight1(i,j)+(multiplier*c(i)*c(j)*
     & temp(i)*temp(j))

                    enddo
                    enddo

c  weights for bound partialsq




            do j = 1,bead_types

            ratio1 = 1.0

            if(mix_flag(j).eq.1) then
             ratio1 = mix_ratio(j)
            endif

            
            do k = j,bead_types

            ratio2 = 1.0

              if(mix_flag(k).eq.1) then
                      ratio2 = mix_ratio(k)
              endif
            
            if(intra_flag(j,k).eq.0) cycle

            if(j.eq.k) then
                    multiplier = 1.0d0
            else
                    multiplier = 2.0d0
             endif

            
             if((ratio1.eq.1.0).and.(ratio2.eq.1.0)) then

      intra_weight(j,k) = intra_weight(j,k)+ (
     &temp(j)*temp(k)*c(j)*c(k)*multiplier)
c           intra_norm(typej,typek) = intra_norm(typej,typek) +
c     & (iso_fraction(deu_step,i)*nmol(i))

        elseif((ratio1.ne.1.0).and.(ratio2.eq.1.0)) then


      intra_weight(j,k) = intra_weight(j,k)+ ((ratio1*temp1(j)+
     & (1.0-ratio1)*temp2(j))*temp(k)*c(j)*c(k)*multiplier)
c             intra_norm(typej,typek) = intra_norm(typej,typek) +
c     & (iso_fraction(deu_step,i)*nmol(i))

           elseif((ratio1.eq.1.0).and.(ratio2.ne.1.0)) then

      intra_weight(j,k) = intra_weight(j,k)+ ((ratio2*temp1(k)+
     & (1.0-ratio2)*temp2(k))*temp(j)*c(j)*c(k)*multiplier)
c             intra_norm(typej,typek) = intra_norm(typej,typek) +
c     & (iso_fraction(deu_step,i)*nmol(i))

          elseif((ratio1.ne.1.0).and.(ratio2.ne.1.0)) then

      intra_weight(j,k) = intra_weight(j,k)+((ratio1*temp1(j)*
     & temp1(k)+(1.0-ratio1)*temp2(j)*temp2(k))
     &*c(j)*c(k)*multiplier)
c             intra_norm(typej,typek) = intra_norm(typej,typek) +
c     & (iso_fraction(deu_step,i)*nmol(i))

           endif

           enddo

           enddo
                   
            


            
c  weights for bound partialsq

c           do i = 1,num


c           if(deut_flag(deu_step,i).ne.2) then  

c             do j = 1,nbd(i)

c             typej = model_type(i,j)

             

c             do k = j,nbd(i)

c             typek = model_type(i,k)

             


c             if(intra_flag(typej,typek).eq.0) cycle


c       intra_weight(typej,typek) = intra_weight(typej,typek)+ (nmol(i)*
c     &iso_fraction(deu_step,i)*temp(typej)*temp(typek))
c             intra_norm(typej,typek) = intra_norm(typej,typek) +
c     & (iso_fraction(deu_step,i)*nmol(i))


c       intra_weight(typej,typek) = intra_weight(typej,typek)+ 
c     & (temp(typej)*temp(typek))
c             intra_norm(typej,typek) = intra_norm(typej,typek) + 1.0d0
c     & (iso_fraction(deu_step,i)*nmol(i))
c            enddo

c            enddo
                   
c            else


c             do j = 1,nbd(i)

c             typej = model_type(i,j)

c             if(ex_cg(deu_step,typej).eq.1) then
c             bj = temp(typej)
c               else
c                 bj = temp1(typej)
c              endif


          
c             do k = j,nbd(i)

c             typek = model_type(i,k)

c             if(ex_cg(deu_Step,typek).eq.1) then

c                     bk = temp(typek)
c             else
c                     bk = temp1(typek)
c             endif



c             if(intra_flag(typej,typek).eq.0) cycle


c             intra_weight(typej,typek) = intra_weight(typej,typek)+
c     & (nmol(i)*(1.0d0-iso_fraction(deu_step,i))*bj*bk)
c             intra_norm(typej,typek) = intra_norm(typej,typek) +
c     & ((1.0d0 - iso_fraction(deu_step,i))*nmol(i))


c             intra_weight(typej,typek) = intra_weight(typej,typek)+
c     & (bj*bk)
c             intra_norm(typej,typek) = intra_norm(typej,typek) + 1.0d0
c     & ((1.0d0 - iso_fraction(deu_step,i))*nmol(i))
c            enddo

c            enddo

c             do j = 1,nbd(i)

c             typej = model_type(i,j)

c             if(ex_cg(deu_step,typej).eq.1) then

c              bj = temp(typej)
c               else

c              bj = temp2(typej)
c              endif




c             do k = j,nbd(i)

c             typek = model_type(i,k)

c             if(ex_cg(deu_step,typek).eq.1) then

c             bk = temp(typek)
c             else
c             bk = temp2(typek)
c             endif




c             if(intra_flag(typej,typek).eq.0) cycle


c             intra_weight(typej,typek) = intra_weight(typej,typek)+
c     & (nmol(i)*iso_fraction(deu_step,i)*bj*bk)
c             intra_norm(typej,typek) = intra_norm(typej,typek) +
c     & (nmol(i)*iso_fraction(deu_step,i))


c             intra_weight(typej,typek) = intra_weight(typej,typek)+
c     & (bj*bk)
c            intra_norm(typej,typek) = intra_norm(typej,typek) + 1.0d0
c    & (nmol(i)*iso_fraction(deu_step,i))
c            enddo
c            enddo
c            endif

c            enddo



           

c           write(io7,*)'INTRA_WEIGHT matrix'

           do i = 1,bead_types
           do j = 1,bead_types

c           if (i.eq.j) then
c              multiplier = 1.0d0
c           else
c              multiplier = 2.0d0
c           endif

           if(intra_flag(i,j).eq.1) then


c      write(io7,*)'intra',deu_step,i,j,intra_weight(i,j),intra_norm(i,j)
c           if(intra_norm(i,j).eq.0.0d0) cycle

c       intra_weight(i,j) = intra_weight(i,j)/intra_norm(i,j)
c       write(io7,*)deu_step,i,j,intra_weight(i,j)
c        intra_weight(i,j) = intra_weight(i,j)*multiplier *c(i)*c(j)

        write(io7,*)deu_step,i,j,intra_weight(i,j)

            endif

            enddo

            enddo


            do i = 1,bead_types

        do k = 1,wbins

       test2(k)=test2(k)+((weight(i)*FF(i,k)**2)
     &*(1.0d0/avg_noa_bead))
 
        enddo
        enddo

        do i =1,wbins
        write(143,*)deu_step,binw(i),test2(i)
        enddo

         do i = 1,bead_types

        do k = 1,wbins

       d_totalweightedsq(k)=d_totalweightedsq(k)+((weight(i)*FF(i,k)**2)
     &*(1.0d0/avg_noa_bead))
c     &*(1.0d0/noa_bead(i)))
 
        enddo


   
           do j = i,bead_types

                 
c                 if(i.eq.j)then
c                         multiplier = 1.0d0
c                 else
c                         multiplier = 2.0d0
c                  endif
               
                  do k = 1,wbins

c         if(i.eq.j) then


       d_weightedsq(i,j,k)= ((unbound_partialsq(i,j,k)*weight1(i,j))+
     &(bound_partialsq(i,j,k)*intra_weight(i,j)))*
     &FF(i,k)*FF(j,k)*(1.0d0/avg_noa_bead)
c     &FF(i,k)*FF(j,k)*((dble(nr(i))+dble(nr(j)))/
c     &dble((nr(i)*noa_bead(i))+(nr(j)*noa_bead(j))))

c     &FF(i,k)*FF(j,k)*(1.0d0/avg_noa_bead)
c     &FF(i,k)*FF(j,k)*(((dble(noa_bead(i))+dble(noa_bead(j))))/
c     &(2.0d0*dble(noa_bead(i))*dble(noa_bead(j))))


       write(253,*)deu_step,i,j,k,unbound_partialsq(i,j,k),weight1(i,j),
     & bound_partialsq(i,j,k),intra_weight(i,j)
c         else
c        d_weightedsq(i,j,k)= ((unbound_partialsq(i,j,k)*weight1(i,j))+
c     &(bound_partialsq(i,j,k)*intra_weight(i,j)))*
c     &FF(i,k)*FF(j,k)*((2.0d0*(dble(noa_bead(i))+dble(noa_bead(j))))/
c     &(2.0d0*dble(noa_bead(i))*dble(noa_bead(j))))
c        endif

c    &FF(i,k)*FF(j,k)*((noa_bead(i)*nr(i)*noa_bead(j)*nr(j))/
c    &dble(175413))


c     &FF(i,k)*FF(j,k)*(1.0d0/avg_noa_bead)
c     &FF(i,k)*FF(j,k)*(1.0d0/((noa_bead(i)+noa_bead(j))*0.5d0))

c     &FF(i,k)*FF(j,k)*(1.0d0/((noa_bead(i)+noa_bead(j))*0.5d0))


c        d_weightedsq(i,j,k)= (unbound_partialsq(i,j,k)*weight1(i,j))*
c     &FF(i,k)*FF(j,k)*(1.0d0/avg_noa_bead)

        
       d_totalweightedsq(k) = d_totalweightedsq(k) + d_weightedsq(i,j,k)
       test1(k) = test1(k) + d_weightedsq(i,j,k)
         
                  enddo

                  enddo

                  enddo

                 file_cnt = 0

         do i = 1,bead_types
         do j = i,bead_types


                 label1 = label_bead(i)
                 label2 = label_bead(j)

                 file_id = text2(deu_step)

                  
      file_name = trim(adjustl(file_id)) // '-weighted-' //
     &trim(adjustl(label1)) // '-' // trim(adjustl(label2)) // '.sq'

                file_cnt = file_cnt + 1

                fu = 300+file_cnt


                open(unit = fu,file=trim(file_name),status='unknown')


              do k = 1,wbins

       write(fu,*)binw(k),d_weightedsq(i,j,k)
       

               enddo
        enddo
        enddo

        do i = 1,wbins

        write(fu1,*)binw(i),d_totalweightedsq(i)
        write(153,*)deu_step,binw(i),test1(i)

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


