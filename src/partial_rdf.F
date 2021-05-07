        subroutine partial_rdf
        implicit none
      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'

      integer i,j,k


        call form_partial(histogram,partialgr)

       call form_partial(bound_hist,bound_partialgr)
       call form_partial(unbound_hist,unbound_partialgr)

       call form_total(partialgr,totalgr,nbins)


        do i = 1,ntypes
        do j = 1,ntypes

              do k = 1,nbins

       write(io3,*)i,j,k,partialgr(i,j,k),bound_partialgr(i,j,k),
     & unbound_partialgr(i,j,k),binr(k)

               enddo

        enddo
        enddo

      do i = 1,nbins

      write(io4,*)binr(i),totalgr(i)

      enddo

         

       end

