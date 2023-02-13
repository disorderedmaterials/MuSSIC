C     Config.Inc


C     Information about configuration 
C
C
      integer*8 itype(maxnatoms),nr(maxtyp),frame(maxbins)
      character*10 label_atom(maxnatoms)
      double precision x(maxnatoms),y(maxnatoms),z(maxnatoms)
     &,rr(maxtyp,maxnatoms,3)

c      double precision,dimension(:,:,:),allocatable :: rr
c      double precision,dimension(:,:,:),allocatable :: x,y,z


      double precision xpbc(maxnatoms),ypbc(maxnatoms),zpbc(maxnatoms)
      double precision xcgpbc(maxnatoms),ycgpbc(maxnatoms),
     & zcgpbc(maxnatoms)
       character*10 bead_lbl(maxnatoms)
       double precision bead_x(maxnatoms),bead_y(maxnatoms),
     &bead_z(maxnatoms)

      

      common/ conf1/ x,y,z,xpbc,ypbc,zpbc,xcgpbc,ycgpbc,zcgpbc
      common/ conf2/ rr
      common/ conf3/ itype,nr,frame,label_atom
      common/ conf4/ bead_lbl
      common/ conf5/ bead_x,bead_y,bead_z      
 
