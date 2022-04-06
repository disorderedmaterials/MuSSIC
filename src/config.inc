C     Config.Inc


C     Information about configuration 
C
C
      integer*8 itype(maxnatoms),nr(maxtyp)
      character*10 label_atom(maxnatoms)
      double precision x(maxnatoms),y(maxnatoms),z(maxnatoms)
     &,rr(maxtyp,maxnatoms,3)
      

      common /conf1/   x,y,z
      common/ conf2/ rr
      common/ conf3/ itype,nr,label_atom      
 
