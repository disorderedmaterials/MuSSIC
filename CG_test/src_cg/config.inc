C     Config.Inc


C     Information about configuration 
C
C
      integer itype(maxnatoms),nr(maxtyp)
      character*10 label_atom(maxnatoms)
      double precision x(maxnatoms),y(maxnatoms),z(maxnatoms)
     &,rr(maxtyp,maxnatoms,3)

       character*10 bead_lbl(maxnatoms)
       double precision bead_x(maxnatoms),bead_y(maxnatoms),
     &bead_z(maxnatoms)

      

      common /conf1/   x,y,z
      common/ conf2/ rr
      common/ conf3/ itype,nr,label_atom
      common/ conf4/ bead_lbl
      common/ conf5/ bead_x,bead_y,bead_z      
 
