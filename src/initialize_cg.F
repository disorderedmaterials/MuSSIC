      subroutine initialize_cg
cccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       initialize_cg
c     Usage
c       call initialize_cg
c     Description
c       initializes  parameters and arrays for RDF caluclation
c       defines the bins for histogram binning
c
c     end_doc
ccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'system.inc'
      include 'model.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'

      integer*8 ii,jj,l,p,i,j,k,typej,typek,temp
         double precision delta,minimum,boxx1,boxy1,boxz1

c     ... 5 total unweighted g(r)
      open (unit=io24,file='CG-unweighted-total.gr',status ='UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 6 total unweighted S(Q)
      open (unit=io26,file='CG-unweighted-total.sq',status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')
ccccccccccccccccccccccccc
c  computing the cutoff

          minimum = boxy

          if(Boxx.lt.Boxy) minimum =boxx
          if(Boxz.lt.minimum) minimum = boxz

          minimum = minimum/2.0d0


        cutoff =  minimum -  mod(minimum,bin_width)
ccccccccccccccccccccccccccccccccccccccccccc
c   initialize the bead size arrays to zero        

         do i = 1,bead_types

          if(iradius.eq.0..and.ibonded.eq.1) then

        
         sigma(1,i) = 0.0d0
         sigma(2,i) = 0.0d0
       
         sigma_cnt(i) = 0.0d0
       

         endif

         enddo 

ccccccccccccccccccccccccccccccccccccccc
c   number of bins for rdf calculation

          nbins = int(cutoff/bin_width)
          max_binr = nbins
ccccccccccccccccccccccccccccccccccccccccc
c   initialize the frame cnt for the bins          

          do i = 1,max_binr
          frame(i) = 0
          enddo
ccccccccccccccccccccccccccccc          
c   number of bins in Q-space

          qmin = (2.0d0*pi)/cutoff

          wbins = idnint((qmax-qmin)/qdelta)

ccccccccccccccccccccccccccccccccccccccccccccc
c   define the bin values for rdf calculation

          delta = bin_width/2.0d0

            do i = 1,maxbins
            binr(i) =  delta
            delta = delta + bin_width
            enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   define the bin values for Q-space depending on the binning regime

          if(bin_type.eq.'linear') then

           
            delta = qmin
            do i = 1,wbins
            binw(i) = delta
            delta = delta + qdelta
            enddo

           elseif(bin_type.eq.'log_linear') then
            delta = qmin
             do i = 1,maxbins

             if(delta.lt.qmax) then

             binq(i) = qdelta * dtanh((log_width/qdelta)*i)


             binw(i) = delta + binq(i)/2.0d0


            delta = delta + binq(i)

             wbins = i

             endif

             enddo

             endif

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c  finding the sequence of the beads in the molecule             


            boxx1 = boxx
            boxy1 = boxy
            boxz1 = boxz

        write(io7,*)' sequence of the beads in the molecule'
c read number of atoms
        read(io3,*)natoms
        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        read(io3,*)boxx,boxy,boxz
        else
        read(io3,*)
        endif
        l = 0
        p = 0

        do ii = 1,num

             do jj = 1,nmol(ii)

        
c        write(io7,*)'molecule',ii,'with beads:', nbd(ii)
c loop over number of atoms to read x, y and z cordinates
        do i = 1,nbd(ii)
        l = l + 1
        mol_num(l) = jj + p 
        read(io3,*)label_atom(l),x(l),y(l),z(l)

        if (jj.eq.1)then
         label(ii,i) = label_atom(l)


c        write(io7,*)label(ii,i)

         endif
        enddo

             enddo

             p = p + nmol(ii)

        enddo

        close(io3)
      open (unit = io3, file = 'CG_trajectory',  status = 'unknown')

        boxx = boxx1
        boxy = boxy1
        boxz = boxz1

c        endif

             

ccccccccccccccccccccccccccc
c initialisation of arrays

       do i = 1,bead_types


          do j = 1,bead_types

          intra_flag(i,j) = 0

          do k = 1,maxbins

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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  finding the exchangable atom types  and tagging them         


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

              do i = 1,maxbins
              totalweightedgr(i) = 0.0d0
              enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Label the beads that belong to same molecule. This helps computing
c  the rdf correlations due to intramolecular beads

                
                    do i = 1,num

                     do j = 1,nbd(i)


                        do k = 1,bead_types
   
                if(label(i,j).eq.label_bead(k)) model_type(i,j)= k

                        enddo

                       temp = model_type(i,j)

                        mol(temp) = i

                       

                       enddo

                       enddo

c                       write(io7,*)'molecule number for bead type'

c                       do i = 1,bead_types
c                       write(io7,*)i,mol(i)
c                       enddo

c                       stop



                       do i = 1,num

                        do j = 1,nbd(i)

                         
                         typej = model_type(i,j)

                            do k = 1,nbd(i)

                            typek = model_type(i,k)


                            intra_flag(typej,typek) = 1


                            enddo

                            enddo

                            enddo


                            write(io7,*)'intra flag'

                            do i = 1,bead_types

                            do j = 1,bead_types

             write(io7,*)label_bead(i),label_bead(j),intra_flag(i,j)

                            enddo
                            enddo


                            do i = 1,num_iso

                            do j = 1, bead_types

                            ex_cg(i,j) = 0
                            isotope_cg(i,j) = 0.0
                            natural_cg(i,j) = 1.0


                            do k = 1,num_atmtyp(j)


                              temp = atmtyp(j,k)
                                         

           if(ex(temp).eq.1) then
            ex_cg(i,j) = 1
           endif


         if(isotope(i,temp).gt.0.0) then
        isotope_cg(i,j) = isotope(i,temp)
        natural_cg(i,j) = 1.0 - isotope(i,temp)
         endif


        
                           enddo

                           enddo

                           enddo

                

               write(io7,*)'deuteration flag'
               
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

c               write(io7,*)i,j,deut_flag(i,j)

                    enddo

                enddo

                enddo



         end

