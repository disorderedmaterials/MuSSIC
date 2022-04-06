C     model.inc


      integer*8 num, nmol(maxnum),nat(maxnum),intra_flag(1000,1000),
     & model_type(1000,1000)
      character*10 label(maxnum,1000)

      common/ mod1/    num,nmol,nat,intra_flag,model_type
      common/ mod2/    label

