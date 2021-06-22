      subroutine initialize
      implicit none
      include 'system.inc'
      include 'model.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'

      integer i,j,k,typej,typek,temp
         double precision delta



          nbins = int(cutoff/bin_width)


          wbins = idnint((qmax-qmin)/qdelta)

            delta = bin_width/2.0d0

            do i = 1,nbins
            binr(i) =  delta
            delta = delta + bin_width
            enddo

           

            delta = qmin
            do i = 1,wbins
            binw(i) = delta
            delta = delta + qdelta
            enddo


c initialisation of arrays

       do i = 1,ntypes


          do j = 1,ntypes

          intra_flag(i,j) = 0

          do k = 1,nbins

          histogram(i,j,k) = 0.0d0
          bound_hist(i,j,k) = 0.0d0
          unbound_hist(i,j,k) = 0.0d0
          totalgr(k) = 0.0d0
          enddo

           do k = 1,wbins
           totalsq(k) = 0.0d0
           enddo

          enddo
         enddo



              do i = 1,ntypes

                   ex(i) = 0

                   do j = 1,noex

                   if(label_ex(j).eq.label_type(i)) then
                    ex(i) = 1
                   endif

                   enddo

             enddo



              do i = 1,wbins
             totalweightedsq(i) = 0.0d0
              enddo

              do i = 1,nbins
              totalweightedgr(i) = 0.0d0
              enddo



                
             


                   do i = 1,num

                     

                     do j = 1,nat(i)


                        do k = 1,ntypes
   
                     if(label(i,j).eq.label_type(k)) model_type(i,j)= k

                     
                        enddo

                        temp = model_type(i,j)

                        mol(temp) = i



                       enddo

                       enddo





                       do i = 1,num

                        do j = 1,nat(i)

                         
                         typej = model_type(i,j)

                         

                            do k = 1,nat(i)

                            typek = model_type(i,k)


                            intra_flag(typej,typek) = 1



                            enddo

                            enddo

                            enddo

                

                do i = 1,num_iso

                     do j = 1,num

                       deut_flag(i,j) = 0
                       iso_fraction(i,j)= 1.0

                  do k = 1, ntypes

                      if(mol(k).eq.j) then

                         if(isotope(i,k).eq.1.0) then

                             deut_flag(i,j) = 1

              elseif (isotope(i,k).gt.0.0.and.isotope(i,k).lt.1.0) then

                      deut_flag(i,j) = 2

                      iso_fraction(i,j) = isotope(i,k)

                       endif

               endif

                    enddo

                enddo

                enddo






              
         end

