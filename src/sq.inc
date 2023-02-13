c   sq.inc
c
c
       double precision partialsq(maxtyp,maxtyp,maxbins),
     & bound_partialsq(maxtyp,maxtyp,maxbins),
     & unbound_partialsq(maxtyp,maxtyp,maxbins)

     
       double precision totalsq(maxbins)
       double precision total_boundsq(maxbins)
       double precision total_unboundsq(maxbins)

      common/sq1/ partialsq,bound_partialsq,unbound_partialsq
      common/sq2/ totalsq,total_boundsq,total_unboundsq
