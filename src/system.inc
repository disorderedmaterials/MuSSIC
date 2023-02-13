C     System.Inc

      character*40 title, broadfn,windowfn,line,FF_label,run_type
      double precision cutoff,boxx,boxy,boxz,bin_width,qmin,qdelta,
     &qmax,FWHM
           
       integer*8 natoms,ntypes,nframes,sizemol,maxnatoms,maxtyp,maxbins,
     & nbins,wbins,maxnum,nproc

      common/ sys1/    cutoff,boxx,boxy,boxz,bin_width,qmin,qdelta,qmax
      common/ sys2/    FWHM
      common/ sys3/    natoms,ntypes,nframes,sizemol,nbins,wbins,nproc
      common/ sys4/    title,broadfn,windowfn,line,FF_label,run_type
     

C vectors and matrixes sizes
      parameter (maxnatoms = 1000000)   ! max number of atoms
      parameter (maxtyp = 10)   ! max number of atom types
      parameter (maxbins = 10000) !max number of histogram bins     
      parameter (maxnum = 10) ! max number of species
C physical constants
      double precision PI
      parameter (PI = 3.14159 26535 89793 23846 26433 8d00)

      double precision bc(maxtyp),c(maxtyp)
      common/ sys4/ bc,c

      double precision binr(maxbins),binw(maxbins),binq(maxbins),
     &log_width

      common/ sys5/ binr,binw,binq,log_width

      integer*8 num_iso,ex(maxtyp),noex,mol(maxtyp),
     &deut_flag(maxtyp,maxtyp),rmsd_cnt(maxtyp),isubstitute

      common/ sys6/ num_iso,ex,noex,mol,deut_flag,rmsd_cnt,isubstitute
      
      character*10 label_iso(maxtyp,maxtyp),label_ex(maxtyp),
     &text1(maxtyp),text2(maxtyp),bin_type,ensemble_type,iso_sub
      double precision natural(maxtyp,maxtyp),isotope(maxtyp,maxtyp),
     &iso_fraction(maxtyp,maxtyp)
      

      common/ sys7/ label_iso,label_ex,text1,text2,bin_type
      common/ sys10/ ensemble_type,iso_sub

      common/ sys8/ natural,isotope,iso_fraction

      character*10 label_type(maxtyp)      
      double precision t1,t2,t3
      common/ sys9/ label_type
      common/tt/ t1,t2,t3

      integer*8 max_binr

       common/sys11/ max_binr
