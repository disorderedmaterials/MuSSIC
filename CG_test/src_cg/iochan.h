c *** include file
c *** define I/O units 

c     added type defs for implicit none modules
      integer io00,io0,io1,io2,io3,io4,io5,io6,io7,io8, io9
      integer io10, io11, io12,io13,io14,io15,io16,io66 


      parameter (io0 = 100)
      parameter (io00 = 200)
        
      parameter (io1 = 1)
c     input : data file
c     formatted, sequential
        
      parameter (io2 = 2)
c     input : trajectory file
c     formatted, sequential
        
      parameter (io3 = 3)
c     output : unweighted partial  g(r)               
c     formatted, sequential
        
        parameter (io4 = 4)
c     output : unweighted total  g(r)               
c     formatted, sequential
        
        parameter (io5 = 5)
c     output : unweighted partial  S(Q)               
c     formatted, sequential
        
        parameter (io66 = 66)
c     output : unweighted total  S(Q)             
c     formatted, sequential
        
        parameter (io7 = 7)
c     output :     print output
c     formatted, sequential
        
        parameter (io8 = 8)
c     output : weighted   S(Q)             
c     formatted, sequential
        
        parameter (io9 = 9)
c     output : weighted total F(Q)             
c     formatted, sequential

        parameter (io10 = 10)
c     output : weighted total  G(r)             
c     formatted, sequential
        
        parameter (io11 = 11)
c     output : weighted   S(Q) - deuterated            
c     formatted, sequential
        
        parameter (io12 = 12)
c     output : weighted total F(Q) - deuterated    
c     formatted, sequential

        parameter (io13 = 13)
c     output : weighted total  G(r) - deuterated 
c     formatted, sequential
        
        parameter (io14 = 14)
c     output : weighted   S(Q) - mixed
c     formatted, sequential
        
        parameter (io15 = 15)
c     output : weighted total F(Q) - mixed         
c     formatted, sequential

        parameter (io16 = 16)
c     output : weighted total  G(r) - mixed      
c     formatted, sequential
        
