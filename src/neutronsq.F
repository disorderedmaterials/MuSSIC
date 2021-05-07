         subroutine neutronsq
         implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'

           integer i,j,k

           double precision multiplier,norm,data_1(wbins),data_2(nbins)
           double precision intraweight(ntypes,ntypes)

          do i = 1,ntypes
          c(i) = dble(nr(i))/dble(natoms)
          enddo


        do i = 1,ntypes

         do j = 1,ntypes


           intraweight(i,j) = c(i)*c(j)*
     &  bc(i)*bc(j)




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

        weightedsq(i,j,k)=partialsq(i,j,k)*c(i)*c(j)*bc(i)*bc(j)
     & *multiplier
        bound_wsq(i,j,k)=bound_partialsq(i,j,k)*intraweight(i,j)
     & *multiplier
        unbound_wsq(i,j,k)=weightedsq(i,j,k)-bound_wsq(i,j,k)

        if(i.ge.j) then
      totalweightedsq(k) = totalweightedsq(k) + weightedsq(i,j,k)
         endif
                  enddo



                  do k = 1,nbins


      weightedgr(i,j,k) = (partialgr(i,j,k)-1.0d0)*c(i)*c(j)
     &*bc(i)*bc(j)
    

      totalweightedgr(k) = totalweightedgr(k) + weightedgr(i,j,k)
                   enddo
       

                  enddo

                  enddo

         do i = 1,ntypes
         do j = 1,ntypes

              do k = 1,wbins

       write(io8,*)i,j,k,weightedsq(i,j,k),bound_wsq(i,j,k),
     & unbound_wsq(i,j,k),binw(k)

               enddo
        enddo
        enddo

        do i = 1,wbins

        write(io9,*)binw(i),totalweightedsq(i)

        enddo
              
           norm = 1.0d0/(2.0d0*PI*PI*(dble(natoms)/box**3))


           do k = 1,wbins

          data_1(k) = totalweightedsq(k)
         

          enddo


         call sinebFT(data_1,data_2,norm,windowfn)

       

         do k = 1,nbins
           
            Gr(k) = data_2(k)

         enddo


         do i = 1,nbins-1

         write(io10,*)binr(i),Gr(i)

         enddo


             end


