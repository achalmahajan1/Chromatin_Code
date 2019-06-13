MODULE Variables
	IMPLICIT NONE
	        INTEGER*8, PARAMETER :: N = 5, Nc = 2, Ntr = 0 !N = number of beads to be simulated, Nc = Number of chains, Ntr = Number of tracer particles
                                                                !Can be either bubbles of tracer particles		
		REAL*8, PARAMETER :: PI = 4.d0*ATAN(1.d0), Rs = 8.d0, aoc = 0.8d0,boc = 0.8d0
		
		REAL*8, PARAMETER :: Fm = 150.d0, Hc = 120.d0, Bpar = 3.d0! Active force magnitude (negative for contractile)
		
		REAL*8, PARAMETER :: kon = 200.d0, koff = 500.d0 !Rate of switch on/off for the active forces
		
		REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Fdpx,Fdpy,Fdpz,Fdpxa,Fdpya,Fdpza,trel,tswitch,Fljx,Fljy,Fljz,u_1,vna,&
		                                        Qx,Qy,Qz,Qo,Fbx,Fby,Fbz,x_temp1,y_temp1,z_temp1!Spring force variables
		
		INTEGER*8, ALLOCATABLE, DIMENSION(:,:) :: dip! dipole (either 1 or 0) based on the motor dynamics 
		
		REAL*8, ALLOCATABLE, DIMENSION(:) :: gamma1,alpha,beta
		
		INTEGER*8 :: nn(1280,6)
		
		REAL*8 :: p(2562,3)
		
END MODULE Variables
