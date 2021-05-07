      module scfstuff
c     f90 module for all variable allocatable vectors 
c     useful for hybrid scf calculations       
   

c     total number of particles
      integer npart, nvert,cnvert

c     debbugging variables
      integer npoint     

c     number of cells in x  y  z  direction
      integer mx, my, mz , mxy,ncell

c     cell length in x  y  and z directio
      real*8 lx, ly, lz, vcell, vbox
      real*8 ilx, ily, ilz, ivcell
      real*8 lxh, lyh, lzh
      real*8 lxhsq, lyhsq, lzhsq
      real*8 shrt_fctr, lng_fctr
      real*8 self_energy_factor


c     some constants used in the f_mfieldtp.F
      real*8 constant1
      real*8 xkomp
      real*8 xtype

c     some constants used in the dens_xyz.F   
      real*8  lx1,lx2
      real*8  ly1,ly2
      real*8  lz1,lz2

      real*8  ilx1,ilx2
      real*8  ily1,ily2
      real*8  ilz1,ilz2

c     vertex indexes vector v(ncells+1,8), list vector of linked list
      integer, allocatable :: v(:,:), list(:), head(:), vertex(:,:,:)
      integer, allocatable :: cvertex(:,:,:)

      integer, allocatable :: dlist(:), dhead(:)

c     particles  coordinates
      real*8  komp, vscf_old, vconv
      real*8, allocatable :: d(:), d_pos(:,:), dprev(:), mass_cell(:)
      real*8, allocatable :: d_ac(:),d_old(:)
      real*8, allocatable :: g(:,:), gdx(:,:) , gdy(:,:), gdz(:,:)
      real*8, allocatable :: f_mfx(:), f_mfy(:), f_mfz(:), chi(:,:) 

      real*8, allocatable :: c_d(:), c_pos(:)
      real*8, allocatable :: c_d_ac(:),c_d_old(:),d_factor(:)
      real*8, allocatable :: eps_d(:)
      real*8, allocatable :: c_g(:,:)
      real*8, allocatable :: grad_x(:),grad_y(:),grad_z(:)
      real*8, allocatable :: cgdx_new(:), cgdy_new(:), cgdz_new(:)

      real*8, allocatable :: f_mx_new(:), f_my_new(:), f_mz_new(:)
      real*8, allocatable :: feng(:),CD(:),f_eng(:)

      real*8, allocatable :: tot_ene(:),elec_energy(:),self_energy(:)
      real*8 tot_ele,etol,omega

C     Sigbjorn 4/04/2018
      real*8, allocatable :: P_eng(:),P_eng_x(:),P_eng_y(:),P_eng_z(:)
      real*8, allocatable :: eps_x(:),eps_y(:),eps_z(:)
      real*8, allocatable :: P_g(:,:)
      real*8, allocatable :: polx_new(:), poly_new(:), polz_new(:)

C     

      
c     cell parameters nrow=no. of vertexes in a row ncol=no. of vert. in a col
c     nlay= no. of vert. in a layer

      integer nrow, ncol, nlay, nfloor
      save
      end module

