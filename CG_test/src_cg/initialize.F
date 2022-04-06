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
         double precision delta,radius(maxtyp)



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


c   initialisation of form factor

c             write(io7,*)'Form factor'

            if(FF_label.eq.'Gaussian')then

             do i = 1, bead_types

c                write(io7,*)'bead radius',i,size_bead(i)

                 do j = 1, wbins

          FF(i,j) = dexp(-0.5d0*(binw(j)*(0.51d0)*size_bead(i))**2.0d0)

                 enddo

             enddo            

             elseif(FF_label.eq.'uniform') then

             do i = 1, bead_types

c                write(io7,*)'bead radius',i,size_bead(i)

                 do j = 1, wbins

          FF(i,j) = (3.0d0/(binw(j)*size_bead(i))**3)*
     &(dsin(binw(j)*size_bead(i))-(dcos(binw(j)*
     & size_bead(i))*binw(j)*size_bead(i)))

                 enddo

             enddo

           elseif(FF_label.eq.'no_FF')then

                   do i = 1,bead_types

                   do j = 1,wbins

                   FF(i,j) = 1.0d0

                   enddo

                   enddo

                endif

c             do i = 1,wbins
c             write(io7,*)binw(i),FF(1,i),FF(2,i),FF(3,i)
c             enddo

            

             


c initialisation of arrays

       do i = 1,bead_types


          do j = 1,bead_types

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

                     

                     do j = 1,nbd(i)


                        do k = 1,bead_types
   
                if(label(i,j).eq.label_bead(k)) model_type(i,j)= k

                     
                        enddo

                       temp = model_type(i,j)

                        mol(temp) = i



                       enddo

                       enddo





                       do i = 1,num

                        do j = 1,nbd(i)

                         
                         typej = model_type(i,j)

                         

                            do k = 1,nbd(i)

                            typek = model_type(i,k)


                            intra_flag(typej,typek) = 1


                            enddo

                            enddo

                            enddo


                            do i = 1,num_iso

                            do j = 1, bead_types

                            ex_cg(i,j) = 0
                            isotope_cg(i,j) = 0.0
                            natural_cg(i,j) = 1.0

                            do k = 1,num_atmtyp(j)


                              temp = atmtyp(j,k)

                              if(ex(temp).eq.1) ex_cg(i,j) = 1
         if(isotope(i,temp).gt.0.0) then
        isotope_cg(i,j) = isotope(i,temp)
        natural_cg(i,j) = 1.0 - isotope(i,temp)
         endif

        

                           enddo

                           enddo

                           enddo



                

                do i = 1,num_iso

                     do j = 1,num

                       deut_flag(i,j) = 0
                       iso_fraction(i,j)= 1.0

                  do k = 1, bead_types

                      if(mol(k).eq.j) then

                         if(isotope_cg(i,k).eq.1.0) then

                             deut_flag(i,j) = 1

         elseif (isotope_cg(i,k).gt.0.0.and.isotope_cg(i,k).lt.1.0) then

                      deut_flag(i,j) = 2

                      iso_fraction(i,j) = isotope_cg(i,k)

                       endif

               endif

                    enddo

                enddo

                enddo



         end

