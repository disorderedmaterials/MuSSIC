c  rdf.inc
c
	
      double precision histogram(maxtyp,maxtyp,maxbins),
     & bound_hist(maxtyp,maxtyp,maxbins),
     & unbound_hist(maxtyp,maxtyp,maxbins)

       double precision partialgr(maxtyp,maxtyp,maxbins),
     & bound_partialgr(maxtyp,maxtyp,maxbins),
     & unbound_partialgr(maxtyp,maxtyp,maxbins)

       double precision totalgr(maxbins)
     
       common/rdf1/ histogram,bound_hist,unbound_hist 
       common/rdf2/ partialgr,bound_partialgr,unbound_partialgr
       common/rdf3/ totalgr
