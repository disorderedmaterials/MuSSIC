C     model.inc


      integer num, nmol(maxnum),nat(maxnum),intra_flag(70,70),
     & model_type(70,70)
      character*10 label(maxnum,70)

      common/ mod1/    num,nmol,nat,intra_flag,model_type
      common/ mod2/    label

