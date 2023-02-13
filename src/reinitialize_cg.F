      subroutine reinitialize_cg
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
      include 'manual.inc'

      integer*8 is,ii,jj,l,i,j,k,typej,typek,temp,fu_1
      double precision boxx1,boxy1,boxz1



           ntypes = ntypes1

           bead_types = bead_types1

c           deallocate(rr)
c           allocate(rr(ntypes,natoms,3))

           num = num1
           noex = noex1
c           nmol = nmol2


           do i = 1,num_iso

           do j = 1,num

           if(deut_flag(i,j).eq.2) then

c           ntypes = ntypes1
c           num = num1
c           noex = noex1
           nmol(j) = nmol2(j)



           endif

           enddo

           enddo

           do i = 1,num
           nat(i) = nat2(i)
           nbd(i) = nbd2(i)
           enddo

           



           write(io7,*)'after manual substitution of the isotope'

           write(io7,*)'number of atom types', ntypes
           write(io7,*)'number of species',num
           write(io7,*)'number of molecules of each species'
           do i = 1,num
           write(io7,*)nmol(i)
           enddo
           
           write(io7,*)'number of atoms in each molecule'
           do i = 1,num
           write(io7,*)nat(i)
           enddo


           write(io7,*)'number of echangeables',noex
           do i = 1,noex
           write(io7,*)ex(i)
           enddo

           write(io7,*)'number of bead types',bead_types
           do i = 1,bead_types
           write(io7,*)'label of the bead',i,label_bead(i)
           enddo

           write(io7,*)'number of beads in each molecule'
           do i = 1,num
           write(io7,*)nbd(i)
           enddo

           if(iradius.eq.1) then
           write(io7,*)'bead radius provided in model file'
         do i = 1,bead_types
         write(io7,*)sigma(1,i)

         enddo
           endif


           write(io7,*)'bead rradius from atomic distribution'

         do i = 1,bead_types
         write(io7,*)sigma_rmsd(1,i)
         enddo

           if(ibonded.eq.1) then

        write(io7,*) ((bead_pair(i,j), j=1,bead_types), i=1,bead_types)

            endif

        write(io7,*)'number of atoms types per CG bead'
       do i = 1,bead_types
       noa_bead(i) = 0
       write(io7,*)num_atmtyp(i)
       write(io7,*)' atom type and the number of those atoms in
     & the bead'
       do j = 1,num_atmtyp(i)
       write(io7,*)atmtyp(i,j),num_atm_bead(i,j)
       noa_bead(i) = noa_bead(i) + num_atm_bead(i,j)
       enddo
       enddo
       

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialisation of arrays for histogram binning

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


c              do i = 1,wbins
c             totalweightedsq(i) = 0.0d0
c              enddo

c             do i = 1,maxbins
c              totalweightedgr(i) = 0.0d0
c              enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    check for the exchangeables and mark the array 'ex(atom type)' as '1'

c             do i = 1,ntypes

c                  ex(i) = 0

c                  do j = 1,noex

c                  if(label_ex(j).eq.label_type(i)) then
c                   ex(i) = 1
c                  endif

c                  enddo
c             enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read number of atoms

             boxx1 = boxx
             boxy1 = boxy
             boxz1 = boxz

           do is = 1, num_iso
          fu_1 = 120 + is 

         close(fu_1)

          open(unit=fu_1,file=filename1,status='unknown')

        read(fu_1,*)natoms
        if(ensemble_type.eq.'NPT'.or.ensemble_type.eq.'npt')then
        read(fu_1,*)boxx,boxy,boxz
        else
        read(fu_1,*)
        endif

        l = 0
        do ii = 1,num

             do jj = 1,nmol(ii)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop over number of atoms to read x, y and z cordinates
        do i = 1,nbd(ii)
        l = l + 1
        read(fu_1,*)label_atom(l),x(l),y(l),z(l)

        if (jj.eq.1) label(ii,i) = label_atom(l)

        enddo

             enddo

        enddo

        close(fu_1)
      open(unit=fu_1,file =filename1,status = 'unknown')

      boxx = boxx1
      boxy = boxy1
      boxz = boxz1


     
       enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Label the atoms that belong to same molecule. This helps computing
c  the rdf correlations due to intramolecular atoms    

                   do i = 1,num

                      do j = 1,nbd(i)


                        do k = 1,bead_types
   
          if(label(i,j).eq.label_bead(k)) model_type(i,j)= k

                        enddo

                        temp = model_type(i,j)

                        mol(temp) = i


c        write(io7,*)'mol num',i,'atm num',j,'type',temp,'mol',mol(temp)

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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Assign  the  values of ratio of isotope to that atom type                 


                do i = 1,num_iso

                     do j = 1,num

                       deut_flag(i,j) = 0
                       iso_fraction(i,j)= 1.0

                  do k = 1, ntypes

                      if(mol(k).eq.j) then

                         if(isotope(i,k).ge.0.0) then

                             deut_flag(i,j) = 1

                       endif

               endif

                    enddo



                enddo

                enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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


              
               end

