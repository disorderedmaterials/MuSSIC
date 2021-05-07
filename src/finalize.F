      subroutine finalize
      implicit none

      include 'system.inc'
      include 'config.inc'
      include 'control.inc'
      include 'iochan.h'
      include 'rdf.inc'
      include 'sq.inc'
      include 'neutronsq.inc'


      integer i,j,k

c        do i = 1,ntypes
c        do j = 1,ntypes

c              do k = 1,nbins

c       write(io3,*)i,j,k,partialgr(i,j,k),bound_partialgr(i,j,k),
c     & unbound_partialgr(i,j,k),binr(k)

c               enddo

c        enddo
c        enddo

c      do i = 1,nbins

c      write(io4,*)binr(i),totalgr(i)

c      enddo


c        do i = 1,ntypes
c        do j = 1,ntypes

c              do k = 1,wbins

c       write(io5,*)i,j,k,partialsq(i,j,k),bound_partialsq(i,j,k),
c     & unbound_partialsq(i,j,k),binw(k)

c               enddo

c        enddo
c        enddo

c      do i = 1,wbins

c      write(io6,*)binw(i),totalsq(i)

c      enddo


c        do i = 1,ntypes
c        do j = 1,ntypes

c              do k = 1,wbins

c       write(io8,*)i,j,k,weightedsq(i,j,k),binw(k)

c               enddo

c        enddo
c        enddo

c      do i = 1,wbins

c      write(io9,*)binw(i),totalweightedsq(i)

c      enddo

      do i = 1,nbins

      write(io10,*)binr(i),Gr(i)

      enddo



      end
