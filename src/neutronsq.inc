C neutronsq.inc

	double precision weightedsq(maxtyp,maxtyp,maxbins),
     & weightedgr(maxtyp,maxtyp,maxbins),totalweightedsq(maxbins),
     & totalweightedgr(maxbins),Gr(maxbins)

       double precision bound_wsq(maxtyp,maxtyp,maxbins),
     & unbound_wsq(maxtyp,maxtyp,maxbins)



       double precision d_weightedsq(maxtyp,maxtyp,maxbins),
     & d_weightedgr(maxtyp,maxtyp,maxbins),d_totalweightedsq(maxbins),
     & d_totalweightedgr(maxbins),d_Gr(maxbins)

       double precision d_bound_wsq(maxtyp,maxtyp,maxbins),
     & d_unbound_wsq(maxtyp,maxtyp,maxbins)


        double precision intra_weight(maxtyp,maxtyp),
     & weight(maxtyp,maxtyp),intra_norm(maxtyp,maxtyp)



       common/neu1/ weightedsq,weightedgr,totalweightedsq,
     & totalweightedgr,Gr
       common/ neu2/ bound_wsq,unbound_wsq

       common/neu3/ d_weightedsq,d_weightedgr,d_totalweightedsq,
     & d_totalweightedgr,d_Gr
       common/ neu4/ d_bound_wsq,d_unbound_wsq

       common/ neu5/ intra_weight,weight,intra_norm



