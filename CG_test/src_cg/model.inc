C     model.inc


      integer num, nmol(maxnum),nat(maxnum),intra_flag(800,800),
     & model_type(800,800),nbd(800),noa_bead(10)
      character*10 label(maxnum,800)

      integer bead_types,num_atmtyp(maxtyp),
     &num_atm_bead(maxtyp,maxtyp),atmtyp(maxtyp,maxtyp),
     &ex_cg(maxtyp,maxtyp)


     
      double precision size_bead(maxtyp),natural_cg(maxtyp,maxtyp),
     &isotope_cg(maxtyp,maxtyp)

      double precision avg_noa_bead,FF(maxtyp,maxbins)

      character*10 lab_atm_bead(maxtyp,maxtyp,70),
     &lbl_bdtyp(maxtyp),label_bead(maxtyp)


      character*10 lbl_bead(maxtyp,maxtyp,maxtyp)

      integer bead_num(10,10),rep_num(10,10),atom_step(10,10),
     &bead_sec(10),list_size(10,10,10),atm_list(10,10,10,10)

      common/ mod1/    num,nmol,nat,intra_flag,model_type,nbd,noa_bead
      common/ mod2/    label
      common/ mod3/    lbl_bead
      common/ mod4/ bead_types,num_atmtyp,num_atm_bead,atmtyp,ex_cg
      common/ mod5/ avg_noa_bead,FF
      common/ mod6/ lab_atm_bead,lbl_bdtyp,label_bead
      common/ mod7/ size_bead,natural_cg,isotope_cg
      common/ mod8/ bead_num,rep_num,atom_step,bead_sec,
     &list_size,atm_list
