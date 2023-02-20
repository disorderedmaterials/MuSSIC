        subroutine cg_init
        implicit none
        
      include 'system.inc'
      include 'iochan.h'
      include 'control.inc'
      include 'config.inc'
      include 'model.inc'


      integer*8 i,ii,j,k,temp

       

        

           

           do ii = 1,num

              

                  
                 temp = 0



                 do i = 1,bead_sec(ii)

                 do j = 1,rep_num(ii,i)



                 do k = 1,bead_num(ii,i)



                              temp = temp + 1

                              label(ii,temp) = lbl_bead(ii,i,k)



                 enddo

                 enddo

                 enddo


                 enddo


         end 
