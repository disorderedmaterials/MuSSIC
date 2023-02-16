C     model.inc


      integer*8 num, nmol(maxnum),nat(maxnum),intra_flag(maxtyp,maxtyp),
     & model_type(maxnum,3000),nbd(maxtyp)
      character*10 label_at(maxnum,3000),label(maxnum,3000)

      integer*8 bead_types,num_atmtyp(maxtyp),
     &atmtyp(maxtyp,maxtyp),ex_cg(maxtyp,maxtyp),
     &bead_pair(maxtyp,maxtyp),bonded(maxtyp)


     
      double precision sigma(2,maxtyp),sigma_cnt(maxtyp),
     &std_dev(maxtyp),natural_cg(maxtyp,maxtyp),
     & isotope_cg(maxtyp,maxtyp),num_atm_bead(maxtyp,maxtyp)

      double precision avg_noa_bead,FF(maxtyp,maxbins),noa_bead(maxtyp)

      character*10 lab_atm_bead(maxtyp,maxtyp,70),
     &lbl_bdtyp(maxtyp),label_bead(20)


      character*10 lbl_bead(20,20,20)

      integer*8 bead_num(20,20),rep_num(20,20),atom_step(20,20),
     &bead_sec(20),list_size(20,20,20),atm_list(20,20,20,100),iradius,
     &ibonded

      
      double precision sigma_rmsd(2,maxtyp),sigma_rmsd_cnt(maxtyp)

      common/ mod1/    num,nmol,nat,intra_flag,model_type,nbd
      common/ mod2/    label_at,label
      common/ mod3/    lbl_bead
      common/ mod4/ bead_types,num_atmtyp,atmtyp,ex_cg
      common/ mod5/ bead_pair,bonded,avg_noa_bead,FF,noa_bead
      common/ mod6/ lab_atm_bead,lbl_bdtyp,label_bead
      common/ mod7/ sigma,sigma_cnt,std_dev,natural_cg,isotope_cg,
     &num_atm_bead
      common/ mod8/ bead_num,rep_num,atom_step,bead_sec,
     &list_size,atm_list,iradius,ibonded
      common/ mod9/ sigma_rmsd,sigma_rmsd_cnt
