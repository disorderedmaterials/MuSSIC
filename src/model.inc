C     model.inc


      integer num, nmol(maxnum),nat(maxnum),intra_flag(maxtyp,maxtyp),
     & model_type(maxtyp,maxtyp)
      character*10 label(maxnum,maxtyp)

      common/ mod1/    num,nmol,nat,intra_flag,model_type
      common/ mod2/    label

