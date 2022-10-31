      subroutine initialize
cccccccccccccccccccccccccccccccccccccccccccccccc
c     start_doc
c     Name
c       initialize
c     Usage
c       call initialize
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

      integer*8 ii,jj,l,i,j,k,typej,typek,temp
         double precision delta,minimum

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         

c     ... 5 total unweighted g(r)
      open (unit = io14,file = 'unweighted-total.gr',status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')

c     ... 6 total unweighted S(Q)
      open (unit=io16, file = 'unweighted-total.sq',status = 'UNKNOWN',
     $      access = 'SEQUENTIAL', form = 'FORMATTED')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  computing the cutoff 

          minimum = boxy

          if(Boxx.lt.minimum) minimum =boxx
          if(Boxz.lt.minimum) minimum = boxz

          minimum = minimum/2.0d0


        cutoff =  minimum -  mod(minimum,bin_width) 
                   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   number of bins for rdf calculation

          nbins = int(cutoff/bin_width)
          max_binr = nbins


          do i = 1,max_binr
          frame(i) = 0
          enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   number of bins in Q-space     

         qmin = (2.0d0*pi)/cutoff

          wbins = idnint((qmax-qmin)/qdelta)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   define the bin values for rdf calculation          

            delta = bin_width/2.0d0

            do i = 1,maxbins
            binr(i) =  delta
            delta = delta + bin_width
            enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc           
c   define the bin values for Q-space

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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialisation of arrays for histogram binning

       do i = 1,ntypes


          do j = 1,ntypes

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


              do i = 1,wbins
             totalweightedsq(i) = 0.0d0
              enddo

              do i = 1,maxbins
              totalweightedgr(i) = 0.0d0
              enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    check for the exchangeables and mark the array 'ex(atom type)' as '1'


              do i = 1,ntypes

                   ex(i) = 0

                   do j = 1,noex

                   if(label_ex(j).eq.label_type(i)) then
                    ex(i) = 1
                   endif

                   enddo

              enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read number of atoms
        read(io2,*)natoms
        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        read(io2,*)boxx,boxy,boxz
        else
        read(io2,*)
        endif

        l = 0
        do ii = 1,num

             do jj = 1,nmol(ii)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop over number of atoms to read x, y and z cordinates
        do i = 1,nat(ii)
        l = l + 1
        read(io2,*)label_atom(l),x(l),y(l),z(l)

        if (jj.eq.1) label_at(ii,i) = label_atom(l)

        enddo

             enddo

        enddo

        close(io2)
      open (unit = io2, file = 'trajectory',  status = 'unknown')


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Label the atoms that belong to same molecule. This helps computing
c  the rdf correlations due to intramolecular atoms    

                   do i = 1,num

                     

                     do j = 1,nat(i)


                        do k = 1,ntypes
   
                 if(label_at(i,j).eq.label_type(k)) model_type(i,j)= k

                     
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Assign  the  values of ratio of isotope to that atom type                 

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

