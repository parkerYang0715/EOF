module LBMset
implicit none
integer, parameter:: opposite(0:8) = (/0,3,4,1,2,7,8,5,6/)
real(kind=8),allocatable :: f(:,:,:),f_potential(:,:,:)  ! f: for flow momentum, f:_potential: electripotential
real(kind=8),allocatable :: f_nega(:,:,:),f_posi(:,:,:) ! f_nega: negative ion concen, f:_potential: positive ion concen
real(kind=8),allocatable :: rho(:,:)
real(kind=8),allocatable :: u(:,:),v(:,:),psi(:,:),psi_old(:,:),con_n(:,:),con_p(:,:)  ! flow velocity psi: electric potential
real(kind=8),allocatable :: weight(:),weight_p(:) ! weight_p: weighting coefficient for electripotential
real(kind=8),allocatable :: cx(:),cy(:)  ! grid velocity
real(kind=8),allocatable :: gradP_x(:,:),gradP_y(:,:)  ! gradient of electric potential
real(kind=8),allocatable :: streamfun_s(:,:),streamfun_p(:,:),streamfun_T(:,:)
real(kind=8),allocatable :: u_analytic_s(:,:),v_analytic_s(:,:),u_analytic_p(:,:),v_analytic_p(:,:)
real(kind=8),allocatable :: u_analytic_T(:,:),v_analytic_T(:,:),psi_analytic(:,:)
integer :: mstep
real*8 :: dt
real*8 :: ck ! Lattice unit per time step
real*8 :: omega
real*8 :: omega_potential
real*8 :: omega_con
real*8 :: tol
end module LBMset
module grid
implicit none
real(kind=8),allocatable::x(:),y(:)
REAL(kind=8) :: xmin = 0.0D0
REAL(kind=8) :: xmax
REAL(kind=8) :: ymin = 0.0D0
REAL(kind=8) :: ymax
integer :: nxp
integer :: nyp
real*8 :: dx
real*8 :: dy
end module grid
module physical_constant ! nondimensional
implicit none
real*8,  parameter :: PI_8 = 4.0D0 * atan (1.0_8)
real*8 :: visco
real*8 :: visco_dim
real*8 :: diffusity ! m^2/s
real*8 :: Dref ! m^2/s
real*8 :: D_lb
real*8 :: epi ! epi= epi_r*epi_0
real*8 :: epi_ref  ! set nondimen(epi)=0.25  ,  nondimen(epi)=epi/epi_ref
real*8 :: Faraday_con
real*8 :: kb ! Boltzmann constant  (J/K)
real*8 :: Temperature
real*8 :: e ! unit charge
real*8 :: KbTovere  ! psi = KbTovere(~0.0258) * nondimensional(psi)
real*8 :: Lref  ! m
real*8 :: Uref     ! m/s
real*8 :: Cref   !mol/m^3
real*8 :: Ex   
real*8 :: Ex_dim
real*8 :: dP  
real*8 :: dP_dim  
real*8 :: Re
real*8 :: Pe
real*8 :: L_debye ! Debye length
real*8 :: kappa
real*8 :: nondimen_coe1 ! ~  0.0013  for source term electric Poisson equation
real*8 :: nondimen_coe2 ! ~2.4965 for source term of NS 1000.0D0:density of fluid
real*8 :: zeta  ! 
real*8 :: zeta_dim  ! 
real*8 :: sigma
real*8 :: sigma_dim
real*8 :: bulk_con = 0.10D0
real(kind=8),allocatable:: bottomBC(:),topBC(:)
end module physical_constant
MODULE FileIndex
  IMPLICIT NONE
  
  ! IO1 : input.dat
  ! IO2 : Convergence DATA
  ! IO3 : Concentration Record

  INTEGER,PARAMETER :: IO1 = 101
  INTEGER,PARAMETER :: IO2 = 102
  INTEGER,PARAMETER :: IO3 = 103  
END MODULE FileIndex
program PNPcoupleNSequation_Benchmark
use LBMset
use grid
use physical_constant
use FileIndex
implicit none
real*8 :: Ma   ! Mach number for incompressible Ma should be small
real*8:: es,er,er_ref ! converge test for poisson equation
real :: start,finish
real :: CPU_Poi_s,CPU_Poi_f,CPU_Poi_T
real*8:: conservation_p,conservation_n ! check if the concentration is conservative 
real*8 :: INITIAL_CP,INITIAL_CN 
real*8 :: dt_NP   !physical time scale of NernstPlanck
real*8 :: dt_NS   !physical time scale of NavierStokes
integer :: NP_NS  !time scale ratio of NP to NS
integer:: i,j,k,kk,kkk,k_NS    ! kkk: time step for poisson equation 
INTEGER :: nb_ticks_initial ! initial value of the clock tick counter
INTEGER ::nb_ticks_final  ! final value of the clock tick counter
INTEGER ::nb_ticks_max  ! maximum value of the clock counter
INTEGER ::nb_ticks_sec  ! number of clock ticks per second
INTEGER ::nb_ticks  ! number of clock ticks of the code
REAL :: elapsed_time ! real time in seconds
! Finish Announcement
OPEN (UNIT=IO2,FILE='Convergence DATA_7_2mum_case1.dat')
OPEN (UNIT=IO3,FILE='7_2mum_case1 Concentration Record.dat')

! input physical constant and computational setting
CALL FileIn()
dt_NP=Lref**2/Dref
dt_NS=Lref/Uref
NP_NS=dt_NP*1.0000001/dt_NS
!--------
write(*,*)'------Check FILE input------'
write(*,*)'TOTAL TIME STEP',mstep
write(*,*)'epi=',epi
write(*,*)'Relaxation_TIME_potential=',1.0D0/omega_potential
write(*,*)'diffusity=',diffusity
write(*,*)'Relaxation_TIME_ion=',1.0D0/omega_con
write(*,*)'Lref=',Lref
write(*,*)'Uref=',Uref
write(*,*)'visco_dim=',visco_dim
write(*,*)'visco=',visco
write(*,*)'relaxation_time for NS eq:',1./omega
write(*,*)'Cref=',Cref
write(*,*)'(dimension) zeta=',zeta*KbTovere
write(*,*)'(dimension) Ex=',Ex*KbTovere/Lref
write(*,*)'DebyeLength',L_debye!=dsqrt(0.50D0*KbTovere*epi/Cref/bulk_con/Faraday_con) ! Debye length
write(*,*)'kappa=',kappa!=1.0D0/L_debye 
write(*,*)'(dimension) sigma_dim',	sigma_dim
write(*,*)'(dimensionless) sigma',	sigma
write(*,*)'(dimension) dP=', dP_dim
write(*,*)'(dimensionless) dP',	dP
write(*,*)'ratio of dp/E=',dP_dim/(Ex*KbTovere/Lref)
write(*,*)'nondimen_coe1',nondimen_coe1!=Lref**2*Cref*Faraday_con/epi_ref/KbTovere  ! ~  0.0013  for source term electric Poisson equation
write(*,*)'nondimen_coe2',nondimen_coe2!=Faraday_con*Cref*KbTovere/Uref**2
write(*,*)'NernstPlanck time scale=',dt_NP, '(s)'
write(*,*)'NavierStokes time scale=',dt_NS, '(s)'
write(*,*)'ratio of NernstPlanck time scale / NavierStokes time scale=',NP_NS
write(*,*)'------FINISH Check FILE input------'
!------------
! Initialization
call gridgen()
call allocateLBMset()
CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
!CALL constructImage()  
allocate(topBC(0:nxp),bottomBC(0:nxp))
topBC=0.
bottomBC=0.

do i=0,nxp
	topBC(i)=zeta
	bottomBC(i)=zeta
end do
topBC(nxp)=topBC(0)
bottomBC(nxp)=bottomBC(0)
k_NS=0
do j=0,nyp
	do i=0,nxp
		rho(i,j)=5.
	end do
end do
u=0.0
v=0.0
psi=0. !initialization of electric potential
psi(:,0)=bottomBC
psi(:,nyp/2)=0. ! the artificial electripotential condition
psi(:,nyp)=topBC
do j=1,nyp-1
	do i=0,nxp
		psi(i,j)=psi(i,0)+j*(psi(i,nyp)-psi(i,0))/nyp ! Initial Guess!!!!
	end do
end do
psi_old=psi
con_n=bulk_con
con_p=bulk_con
psi_old=psi
con_n(:,0)=bulk_con*dexp(psi(:,0))
con_p(:,0)=bulk_con*dexp(-psi(:,0))
con_n(:,nyp)=bulk_con*dexp(psi(:,nyp))
con_p(:,nyp)=bulk_con*dexp(-psi(:,nyp))
INITIAL_CN=0.
INITIAL_CP=0.
do j=0,nyp
	do i=0,nxp 
		INITIAL_CN=INITIAL_CN+con_n(i,j)
		INITIAL_CP=INITIAL_CP+con_p(i,j)
		do k=0,4
			f_nega(k,i,j)=weight_p(k)*con_n(i,j)
			f_posi(k,i,j)=weight_p(k)*con_p(i,j)
		end do
	end do
end do
write(*,*) 'initiall concentration_positive =',INITIAL_CP
write(*,*) 'initiall concentration_negative =',INITIAL_CN
write(*,*) 'psi(nxp/8,0)=',psi(nxp/8,0)
write(*,*) 'psi(nxp/8,nyp)=',psi(nxp/8,nyp)
write(*,*) 'psi(nxp/4,nyp)=',psi(nxp/4,nyp)
write(*,*)' con_n(nxp/2,nyp/2)=',con_n(nxp/2,nyp/2),'con_p(nxp/2,nyp/2)=',con_p(nxp/2,nyp/2)
write(*,*) 'con_n(nxp-10,10)=',con_n(nxp-10,10), 'con_p(nxp-10,nyp-10)=',con_p(nxp-10,nyp-10)
write(*,*)'nondimen_coe1=',nondimen_coe1
write(*,*)'==================================================================='
call cpu_time(start)
CPU_Poi_T=0.
CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
mainloop: do kk=1,mstep
psi_old=0.
call cpu_time(CPU_Poi_s)
ELECTRIC_potential: do kkk=1,6999000
call collision_pot
call streamingD2Q5_p
call f_potbound
call potentialcalcu() 
psi(:,nyp/2)=0.  ! the artificial electripotential condition
con_n(:,nyp/2)=bulk_con
con_p(:,nyp/2)=bulk_con
er=0.
er_ref=0.
do j=0,nyp
	do i=0,nxp
		er=er+dabs(psi(i,j)-psi_old(i,j))
		er_ref=er_ref+dabs(psi_old(i,j))
	end do
end do
er=er/er_ref
if(mod(kkk,80000)==0)then
write(*,*)'Poisson step',kkk,'iterative er=',er
write(*,*)'psi(nxp/2-5,nyp/2-5)=',psi(nxp/2-5,nyp/2-5)
end if
if (er<tol) then
	write(*,*)'er/er_ref=',er/er_ref
	write(*,*)'NernstPlanckLoop=',kk,'converge when kkk step=',kkk
	write(IO2,*)kk,er/er_ref,kkk
end if
if(er<tol) exit
psi_old=psi
end do ELECTRIC_potential
call cpu_time(CPU_Poi_f)
CPU_Poi_T=CPU_Poi_T+CPU_Poi_f-CPU_Poi_s
call gradientP_x()
call gradientP_y()
do k_NS=1,NP_NS  ! NavierStokes solver 
	call collision()  
	call streaming_NS()
	call sfbound()  
	call rhouv()
end do
call collision_con()
call streamingD2Q5_con
call f_conbound
call concen_calcu
!======================
write(*,*)'====step',kk,'===='
write(*,*) 'MAX(u)=',maxval(abs(u))
write(*,*)'con_n(nxp/2,nyp/2)=',con_n(nxp/2,nyp/2),'con_p(nxp/2,nyp/2)=',con_p(nxp/2,nyp/2)
write(*,*)'con_n(nxp-5,8)=',con_n(nxp-5,8), 'con_p(nxp-5,nyp-8)=',con_p(nxp-5,nyp-8)
20 format(F14.9,F14.9,F14.9,F14.9,F14.9)
conservation_n=0.
conservation_p=0.
do j=0,nyp
do i=0,nxp 
conservation_n=conservation_n+con_n(i,j)
conservation_p=conservation_p+con_p(i,j)
end do
end do
!==========================
write(*,*) 'the total concentration_positive =',conservation_p
write(*,*) 'the total concentration_negative =',conservation_n
write(IO3,*)conservation_p,conservation_n

! ====================  OUTPUT VTK FILE  ==============================
!if(kk .LE. 480 .and. mod(kk,40)==0)then
!call outvtk_visualization(kk)
!call outputMATLAB(kk)
!end if
if(kk .gt. 480 .and. kk .lt. 12000)then
if(mod(kk,400)==0) then
!call outvtk_visualization(kk)
call outputMATLAB(kk)
end if
end if
END DO mainloop
call cpu_time(finish)
CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
nb_ticks = nb_ticks_final - nb_ticks_initial
IF (nb_ticks_final < nb_ticks_initial) then
	nb_ticks = nb_ticks + nb_ticks_max
endif
elapsed_time  = REAL(nb_ticks) / nb_ticks_sec
write(*,*)'elapsed_time:',elapsed_time   ! elapsed_time
write(*,*)'cputime:',finish-start
OPEN (UNIT=689,FILE='CPU TIME.dat')
write(689,*)'elapsed_time:',elapsed_time
write(689,*)'CPU_time:',finish-start
OPEN (UNIT=777,FILE='7_2mum_case1 EP UV cp cn.dat')
do j=0,nyp
do i=0,nxp
	write(777,20) KbTovere*psi(i,j),Uref*u(i,j),Uref*v(i,j), Cref*con_p(i,j),Cref*con_n(i,j)
end do
end do
write(*,*)'ratio: Poisson/Total=',CPU_Poi_T/(finish-start)*100.0,'%'
end program

subroutine outputMATLAB(kk)
!======
use LBMset, ONLY:u,v,psi,con_n,con_p
use grid, ONLY:nxp, nyp, x, y
use physical_constant, ONLY: Uref, Cref, KbTovere, Faraday_con
implicit none
integer,intent(in)::kk
integer::i,j
character(len=38)::filename
WRITE(FileName ,"(I5.5)") kk
FileName = "7_2mum_case1_step" // TRIM(FileName) // ".dat"
OPEN( UNIT = 1900 , FILE = TRIM(FileName) , STATUS = 'REPLACE' )
do j=0,nyp
	do i=0,nxp
		write(1900,20) KbTovere*psi(i,j),Uref*u(i,j),Uref*v(i,j), Cref*con_p(i,j),Cref*con_n(i,j)
	end do
end do
20 format(F14.9,F14.9,F14.9,F14.9,F14.9)
end subroutine outputMATLAB

subroutine outvtk_visualization(kk)  !output VTK file for PARAVIEW
!======
use LBMset, ONLY:u,v,psi,con_n,con_p
use grid, ONLY:nxp, nyp, x, y
use physical_constant, ONLY: Uref, Cref, KbTovere, Faraday_con
implicit none
integer,intent(in)::kk
integer::i,j
character(len=38)::filename
character(len=42)::filename2
WRITE(FileName ,"(I5.5)") kk
FileName = "7_2mum_case1_step_" // TRIM(FileName) // ".vtk"
WRITE(FileName2 ,"(I5.5)") kk
filename2='case4_step_'// TRIM(FileName2)
OPEN( UNIT = 1000 , FILE = TRIM(FileName) , STATUS = 'REPLACE' )
write(1000,'(A)')'# vtk DataFile Version 4.0'
write(1000,'(A)') filename2
write (1000, '(A)' ) 'ASCII'
write (1000,'(A)')'DATASET POLYDATA'
write(1000,9) 'POINTS ' ,(nxp+1)*(nyp+1) ,'double'
9 format(A8,I12,A8)
do j=0,nyp
    do i=0,nxp
        write(1000,*) x(i),y(j),0  
    end do
end do
write(1000,*)'POLYGONS',    nxp*nyp  ,    5*nxp*nyp
do j=0,nyp-1
    do i=0,nxp-1
        write(1000,7) 4,i+(nxp+1)*j,i+1+(nxp+1)*j,i+(nxp+2)+(nxp+1)*j,i+(nxp+1)+(nxp+1)*j
    end do
end do
7 format (I4,I12,I12,I12,I12)
write(1000,'(A10,I13)') 'POINT_DATA ',(nxp+1)*(nyp+1)
write(1000,'(A)') 'FIELD FieldData    5'
write(1000,'(A25,I13,A9)') 'Concentration_cation 1 ',  (nxp+1)*(nyp+1) , ' double'  !positive ion
do j=0,nyp
	do i=0,nxp
		write(1000,'(f15.8)') Cref*con_p(i,j)
	end do
end do
write(1000,'(A25,I13,A9)') 'Concentration_anion 1 ',  (nxp+1)*(nyp+1) , ' double'  !negative ion
do j=0,nyp
	do i=0,nxp
		write(1000,'(f15.8)') Cref*con_n(i,j)
	end do
end do
write(1000,'(A16,I13,A9)') 'NetCharge 1 ',  (nxp+1)*(nyp+1) , ' double'
do j=0,nyp
	do i=0,nxp
		write(1000,'(f15.8)') Faraday_con*(Cref*con_p(i,j)-Cref*con_n(i,j))
	end do
end do
write(1000,'(A24,I13,A9)') 'ElectricPotential 1 ',  (nxp+1)*(nyp+1) , ' double'
do j=0,nyp
	do i=0,nxp
		write(1000,'(f15.8)') KbTovere*psi(i,j)
	end do
end do
write(1000,'(A10,I13,A9)') 'Velocity 3 ',  (nxp+1)*(nyp+1) , ' double'
do j=0,nyp
	do i=0,nxp
		write(1000,6) u(i,j)*Uref,v(i,j)*Uref,0
	end do
end do
6 format (F14.9,F14.9,F5.2)
end subroutine outvtk_visualization

SUBROUTINE FileIn()
  USE FileIndex
  USE LBMset
  USE grid
  USE physical_constant
  IMPLICIT NONE
  ! open input.dat =================================================================
  OPEN( UNIT = IO1 , FILE = 'input.dat' , STATUS = 'OLD' )
  ! start the program ==============================================================
    ! read in input.dat ==================================================
	READ(IO1,*) mstep
	READ(IO1,*) nxp
	xmax=nxp
	READ(IO1,*) nyp
	ymax=nyp
	READ(IO1,*) tol
	READ(IO1,*) visco_dim
	READ(IO1,*) diffusity
	READ(IO1,*) Dref
	D_lb=diffusity/Dref
	READ(IO1,*) epi
	READ(IO1,*) epi_ref
	READ(IO1,*) Faraday_con
	READ(IO1,*) kb ! Boltzmann constant  (J/K)
	READ(IO1,*) Temperature
	READ(IO1,*) e ! unit charge
	KbTovere=kb*Temperature/e    ! psi = KbTovere(~0.0258) * nondimensional(psi)
	READ(IO1,*) Lref  ! m
	READ(IO1,*) Uref     ! m/s
	visco=visco_dim/Uref/Lref
	READ(IO1,*) Cref   !mol/m^3
	READ(IO1,*) Ex   ! V/m
	READ(IO1,*) dP_dim  !1.250D0*10.0D0**(-7)*(16.0D0/80.0D0)**3*3.0D0**2 	
	READ(IO1,*) zeta ! (V)
    CLOSE(IO1)
    
	! nondimensionization
	Re=abs(epi*Ex_dim*zeta/1000.0D0/visco_dim)*(ymax*Lref)/visco_dim
	Pe=Uref*Lref/Dref
	L_debye=dsqrt(0.50D0*KbTovere*epi/Cref/bulk_con/Faraday_con) ! Debye length
	kappa=1.0D0/L_debye 
	sigma=zeta*epi*kappa ! build the relationship between zeta and sigma
	nondimen_coe1=Lref**2*Cref*Faraday_con/epi_ref/KbTovere  ! ~  0.0013  for source term electric Poisson equation
	nondimen_coe2=Faraday_con*Cref*KbTovere/Uref**2/1000.0D0 ! ~2.4965 for source term of NS 1000.0D0:density of fluid
	Ex_dim=Ex  ! dimensional
	Ex=Ex/KbTovere*Lref    ! dimensionless
	dP=dP_dim*Lref/Uref**2/1000.0D0
	zeta_dim=zeta ! dimensional
	zeta=zeta/KbTovere ! dimensionless
	sigma_dim=sigma
	sigma=sigma/KbTovere*Lref/epi
    ! compute dx, dy =====================================================
    dx = dble( ( xmax - xmin ) /nxp)
    dy = dble( ( ymax - ymin ) /nyp)
    ! compute dt =========================================================
	dt = dx**2
	ck=dx/dt ! Lattice unit per time step
	omega=1.0D0/(3.0D0*visco*dt/dx**2+0.50D0)
	omega_potential=1.0D0/(3.0D0*epi/epi_ref*dt/dx**2+0.50D0)
	omega_con=1.0D0/(3.0D0*D_lb*dt/dx**2+0.50D0)  
  RETURN
END SUBROUTINE FileIn


subroutine allocateLBMset()
use grid
use LBMset
use physical_constant
implicit none
integer::i
allocate(f(0:8,0:nxp,0:nyp),f_potential(0:4,0:nxp,0:nyp),f_nega(0:4,0:nxp,0:nyp),f_posi(0:4,0:nxp,0:nyp))
allocate(psi(0:nxp,0:nyp),psi_old(0:nxp,0:nyp),rho(0:nxp,0:nyp),weight(0:8),weight_p(0:4),u(0:nxp,0:nyp),v(0:nxp,0:nyp),con_n(0:nxp,0:nyp),con_p(0:nxp,0:nyp))
allocate(cx(0:8),cy(0:8))
allocate(gradP_x(0:nxp,0:nyp),gradP_y(0:nxp,0:nyp))
allocate(streamfun_s(0:nxp,0:nyp),streamfun_p(0:nxp,0:nyp),streamfun_T(0:nxp,0:nyp))
allocate(u_analytic_s(0:nxp,0:nyp),v_analytic_s(0:nxp,0:nyp),u_analytic_p(0:nxp,0:nyp),v_analytic_p(0:nxp,0:nyp))
allocate(u_analytic_T(0:nxp,0:nyp),v_analytic_T(0:nxp,0:nyp),psi_analytic(0:nxp,0:nyp))
f=0.
f_potential=0.
f_nega=0.
f_posi=0.
rho=0
psi=0
streamfun_s=0
u_analytic_s=0
u_analytic_s=0
weight(0)=4.0D0/9.0D0
do i=1,4
	weight(i)=1.0D0/9.0D0
	weight(i+4)=1.0D0/36.0D0
end do
weight_p(0)=1.0D0/3.0D0
do i=1,4
	weight_p(i)=1.0D0/6.0D0
end do
cx(0)=0.0D0
cx(1)=1.0D0
cx(2)=0.0D0
cx(3)=-1.0D0
cx(4)=0.0D0
cx(5)=1.0D0
cx(6)=-1.0D0
cx(7)=-1.0D0
cx(8)=1.0D0
cy(0)=0.0D0
cy(1)=0.0D0
cy(2)=1.0D0
cy(3)=0.0D0
cy(4)=-1.0D0
cy(5)=1.0D0
cy(6)=1.0D0
cy(7)=-1.0D0
cy(8)=-1.0D0

cx=cx*dx/dt
cy=cy*dx/dt
gradP_x=0.
gradP_y=0.

end subroutine allocateLBMset

subroutine gridgen
use grid , ONLY: x, y, nxp, nyp, dx, dy
implicit none
integer ::i,j
allocate (x(0:nxp),y(0:nyp))
do i=0,nxp
	x(i)=i*dx
end do
do j=0,nyp
	y(j)=j*dy
end do
return
end subroutine gridgen

subroutine collision()
use LBMset , ONLY: f, u, v, cx, cy, rho, con_n, con_p, gradP_x, gradP_y, weight, omega, ck, dt
use grid, ONLY: nxp, nyp
use physical_constant , ONLY: dP, Ex, nondimen_coe2
implicit none
real*8::t1,t2,feq
real*8::force_x,force_y,g,Pressureforce
integer::i,j,k
real*8 :: c_1  !constant of equilibrium distribution function
real*8 :: c_2  !constant of equilibrium distribution function
real*8 :: c_3  !constant of equilibrium distribution function
real*8 :: c_4  !constant of equilibrium distribution function
!---------
c_1=3.0D0*dP/ck**2  
c_2=3.0D0/ck**2  
c_3=4.50D0/ck**4 
c_4=1.50D0/ck**2 
!---------
!$omp parallel
!$omp do private(i,j,k,t1,t2,Pressureforce,g,feq)
DO j=0,nyp
DO i=0,nxp
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
force_x=-nondimen_coe2*(con_p(i,j)-con_n(i,j))*(gradP_x(i,j)+Ex)
force_y=-nondimen_coe2*(con_p(i,j)-con_n(i,j))*gradP_y(i,j)
DO k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
Pressureforce=c_1*weight(k)*cx(k)  ! 3/ck^2 = 1/Cs^2, Cs=c/sqrt(3)
g=c_2*rho(i,j)*weight(k)*(cx(k)*force_x+cy(k)*force_y)
if(j==0 .or. j==nyp) g=0.
if(j==0 .or. j==nyp) Pressureforce=0.
feq=rho(i,j)*weight(k)*(1.0D0+c_2*t2+c_3*t2*t2-c_4*t1)
f(k,i,j)=omega*feq+(1.0D0-omega)*f(k,i,j)+dt*g+dt*Pressureforce
END DO
END DO
END DO
!$omp end do
!$omp end parallel
return
end subroutine collision

subroutine collision_pot()
use LBMset , ONLY: f_potential, psi, con_n, con_p, weight_p, omega_potential, dt
use grid, ONLY: nxp, nyp
use physical_constant, ONLY: nondimen_coe1
implicit none
real*8::g,feq
integer::i,j,k
!$omp parallel
!$omp do private(i,j,k,g,feq)
DO j=0,nyp
	DO i=0,nxp
		DO k=0,4
			g=weight_p(k)*nondimen_coe1*(con_p(i,j)-con_n(i,j))  ! source term
			if(j==0 .or. j==nyp) g=0.
			feq=psi(i,j)*weight_p(k)
			f_potential(k,i,j)=omega_potential*feq+(1.0D0-omega_potential)*f_potential(k,i,j)+dt*g
		END DO
	END DO
END DO
!$omp end do
!$omp end parallel
return
end subroutine collision_pot

subroutine collision_con()
use LBMset,ONLY: f_nega, f_posi, gradP_x, gradP_y, con_n, con_p, weight_p, omega_con, u, v, cx, cy, ck
use grid, ONLY: nxp, nyp
use physical_constant, ONLY: Ex, Pe, D_lb
implicit none
real*8::t1,t2,feq
integer::i,j,k
!$omp parallel
!$omp do private(i,j,k,t1,t2,feq)
DO j=0,nyp
DO i=0,nxp
DO k=0,4
	t1=Pe*u(i,j)-(-1.0D0)*(Ex+gradP_x(i,j))*D_lb
	if(j==0 .or. j==nyp)then
		t1=0
	end if
	t2=Pe*v(i,j)-(-1.0D0)*gradP_y(i,j)*D_lb
	feq=con_n(i,j)*weight_p(k)*(1.0+3.0D0*(t1*cx(k)+t2*cy(k))/ck**2)
	f_nega(k,i,j)=omega_con*feq+(1.0D0-omega_con)*f_nega(k,i,j)
END DO
DO k=0,4
	t1=Pe*u(i,j)-(1.0D0)*(Ex+gradP_x(i,j))*D_lb
	if(j==0 .or. j==nyp)then
		t1=0
	end if
	t2=Pe*v(i,j)-(1.0D0)*gradP_y(i,j)*D_lb
	feq=con_p(i,j)*weight_p(k)*(1.0+3.0D0*(t1*cx(k)+t2*cy(k))/ck**2)
	f_posi(k,i,j)=omega_con*feq+(1.0D0-omega_con)*f_posi(k,i,j)
END DO
END DO
END DO
!$omp end do
!$omp end parallel
return
end subroutine collision_con

subroutine streaming_NS()
use grid, ONLY: nxp, nyp
use LBMset, ONLY: f
implicit none
integer::i,j
! streaming
DO j=0,nyp
DO i=nxp,1,-1 !RIGHT TO LEFT
f(1,i,j)=f(1,i-1,j)
END DO
DO i=0,nxp-1 !LEFT TO RIGHT
f(3,i,j)=f(3,i+1,j)
END DO
END DO
DO j=nyp,1,-1 !TOP TO BOTTOM
DO i=0,nxp
f(2,i,j)=f(2,i,j-1)
END DO
DO i=nxp,1,-1
f(5,i,j)=f(5,i-1,j-1)
END DO
DO i=0,nxp-1
f(6,i,j)=f(6,i+1,j-1)
END DO
END DO
DO j=0,nyp-1 !BOTTOM TO TOP
DO i=0,nxp
f(4,i,j)=f(4,i,j+1)
END DO
DO i=0,nxp-1
f(7,i,j)=f(7,i+1,j+1)
END DO
DO i=nxp,1,-1
f(8,i,j)=f(8,i-1,j+1)
END DO
END DO
return
end subroutine streaming_NS

subroutine streamingD2Q5_p
use LBMset, ONLY: f_potential
use grid, ONLY: nxp, nyp
implicit none
!real*8,parameter::Faraday_con=96485.0D0
integer::i,j
! streaming
DO j=0,nyp
DO i=nxp,1,-1 !RIGHT TO LEFT
f_potential(1,i,j)=f_potential(1,i-1,j)
END DO
DO i=0,nxp-1 !LEFT TO RIGHT
f_potential(3,i,j)=f_potential(3,i+1,j)
END DO
END DO
DO j=nyp,1,-1 !TOP TO BOTTOM
DO i=0,nxp
f_potential(2,i,j)=f_potential(2,i,j-1)
END DO
END DO
DO j=0,nyp-1 !BOTTOM TO TOP
DO i=0,nxp
f_potential(4,i,j)=f_potential(4,i,j+1)
END DO
END DO
return
end subroutine streamingD2Q5_p

subroutine streamingD2Q5_con
use LBMset, ONLY: f_nega, f_posi
use grid, ONLY: nxp, nyp
! streaming
DO j=0,nyp
DO i=nxp,1,-1 !RIGHT TO LEFT
f_nega(1,i,j)=f_nega(1,i-1,j)
f_posi(1,i,j)=f_posi(1,i-1,j)
END DO
DO i=0,nxp-1 !LEFT TO RIGHT
f_nega(3,i,j)=f_nega(3,i+1,j)
f_posi(3,i,j)=f_posi(3,i+1,j)
END DO
END DO
DO j=nyp,1,-1 !TOP TO BOTTOM
DO i=0,nxp
f_nega(2,i,j)=f_nega(2,i,j-1)
f_posi(2,i,j)=f_posi(2,i,j-1)
END DO
END DO
DO j=0,nyp-1 !BOTTOM TO TOP
DO i=0,nxp
f_nega(4,i,j)=f_nega(4,i,j+1)
f_posi(4,i,j)=f_posi(4,i,j+1)
END DO
END DO
return
end subroutine streamingD2Q5_con

subroutine sfbound
use LBMset, ONLY: f
use grid, ONLY: nxp, nyp
implicit none
integer :: i,j
!$omp parallel private(i,j)
!$omp do
do j=0,nyp
! bounce back on west boundary
f(1,0,j)=f(1,nxp,j)
f(5,0,j)=f(5,nxp,j)
f(8,0,j)=f(8,nxp,j)
! bounce back on east boundary
f(3,nxp,j)=f(3,0,j)
f(7,nxp,j)=f(7,0,j)
f(6,nxp,j)=f(6,0,j)
end do
!$omp end do

!$omp do
! bounce back on south boundary u=v=0 NONSLIP BC
do i=0,nxp
f(2,i,0)=f(4,i,0)
f(5,i,0)=f(7,i,0)
f(6,i,0)=f(8,i,0)
end do
!$omp end do
!$omp do
! north boundary  NONSLIP BC
do i=0,nxp
!rhon=f(0,i,nyp)+f(1,i,nyp)+f(3,i,nyp)+2.0D0*(f(2,i,nyp)+f(6,i,nyp)+f(5,i,nyp))
f(4,i,nyp)=f(2,i,nyp)
f(8,i,nyp)=f(6,i,nyp)!+rhon*u0/6.0D0
f(7,i,nyp)=f(5,i,nyp)!-rhon*u0/6.0D0
end do
!$omp end do
!$omp end parallel 
return
end subroutine sfbound

subroutine f_potbound !execute before streaming!!!!
use LBMset, ONLY: f_potential,weight_p
use grid, ONLY: nxp, nyp
use physical_constant, ONLY: bottomBC, topBC
implicit none
integer::i,j,k
real*8 :: c_1_3=1.0D0/3.0D0
!$omp parallel private(i,j)
!$omp do
do j=0,nyp ! right B.C.    periodic
	f_potential(3,nxp,j)=f_potential(3,0,j)
end do
!$omp end do
!$omp do
do i=0,nxp ! bottom B.C.
	f_potential(2,i,0)=-f_potential(4,i,0)+c_1_3*bottomBC(i)
end do
!$omp end do
!$omp do
do j=0,nyp ! left boundary  periodic
	f_potential(1,0,j)=f_potential(1,nxp,j)
end do
!$omp end do
!$omp do
do i=0,nxp ! top B.C.  ! DIRICHLET B.C.
	f_potential(4,i,nyp)=-f_potential(2,i,nyp)+c_1_3*topBC(i)
end do
!$omp end do
!$omp end parallel 
return
end subroutine f_potbound

subroutine f_conbound !execute before streaming!!!!
use LBMset, ONLY: f_nega, f_posi, con_n, con_p, v, gradP_x,gradP_y
use grid, ONLY: nxp, nyp,dx
use physical_constant
implicit none
integer :: i,j,k
real*8 :: t2
!$omp parallel private(i,j)
!$omp do
do j=0,nyp ! right B.C.
	f_nega(3,nxp,j)=f_nega(3,0,j) ! 
	f_posi(3,nxp,j)=f_posi(3,0,j)   ! periodic B.C.
end do
!$omp end do
!$omp do
do j=0,nyp ! left boundary  ! periodic B.C.
	f_nega(1,0,j)=f_nega(1,nxp,j) ! periodic B.C.
	f_posi(1,0,j)=f_posi(1,nxp,j) ! periodic B.C.
end do
!$omp end do
!$omp do
!3rd order zero-flux B.C.
! for negative concentration
! do i=0,nxp
	! j=0
	! t2=Pe*v(i,j)-(-1.0D0)*gradP_y(i,j)*D_lb  ! effective velocity for negative ion at boundary
! f_nega(2,i,j)=D_lb*(18.0D0*con_n(i,j+1)-9.0D0*con_n(i,j+2)+2.0D0*con_n(i,j+3))/(6.0D0*dx*t2+11.0D0*D_lb)-f_nega(0,i,j)-f_nega(1,i,j)-f_nega(3,i,j)-f_nega(4,i,j)
	
	! j=nyp
	! t2=Pe*v(i,j)-(-1.0D0)*gradP_y(i,j)*D_lb  ! effective velocity for negative ion at boundary
! f_nega(4,i,j)=D_lb*(18.0D0*con_n(i,j-1)-9.0D0*con_n(i,j-2)+2.0D0*con_n(i,j-3))/(-6.0D0*dx*t2+11.0D0*D_lb)-f_nega(0,i,j)-f_nega(1,i,j)-f_nega(2,i,j)-f_nega(3,i,j)

! end do
! for positive concentration
! do i=0,nxp
	! j=0
	! t2=Pe*v(i,j)-(1.0D0)*gradP_y(i,j)*D_lb  ! effective velocity for positive ion at boundary

! f_posi(2,i,j)=D_lb*(18.0D0*con_p(i,j+1)-9.0D0*con_p(i,j+2)+2.0D0*con_p(i,j+3))/(6.0D0*dx*t2+11.0D0*D_lb)-f_posi(0,i,j)-f_posi(1,i,j)-f_posi(3,i,j)-f_posi(4,i,j)
	
	! j=nyp
	! t2=Pe*v(i,j)-(1.0D0)*gradP_y(i,j)*D_lb  ! effective velocity for negative ion at boundary
! f_posi(4,i,j)=D_lb*(18.0D0*con_p(i,j-1)-9.0D0*con_p(i,j-2)+2.0D0*con_p(i,j-3))/(-6.0D0*dx*t2+11.0D0*D_lb)-f_posi(0,i,j)-f_posi(1,i,j)-f_posi(2,i,j)-f_posi(3,i,j)
! end do
! DIRICHLET by poisson Boltzmann distribution
do i=0,nxp
	f_nega(2,i,0)=-f_nega(4,i,0)+con_n(i,0)/3.0D0
	f_posi(2,i,0)=-f_posi(4,i,0)+con_p(i,0)/3.0D0
	f_nega(4,i,nyp)=-f_nega(2,i,nyp)+con_n(i,nyp)/3.0D0
	f_posi(4,i,nyp)=-f_posi(2,i,nyp)+con_p(i,nyp)/3.0D0
end do
!$omp end do
!$omp end parallel
return
end subroutine f_conbound

subroutine potentialcalcu
use LBMset, ONLY: f_potential, psi
use grid, ONLY: nxp, nyp
use physical_constant, ONLY:zeta, topBC, bottomBC
implicit none
integer::i,j,k
!$omp parallel
!$omp do private(i,j)
do j=1,nyp-1
	do i=0,nxp
		psi(i,j)=f_potential(0,i,j)+f_potential(1,i,j)+f_potential(2,i,j)+f_potential(3,i,j)+f_potential(4,i,j)
	end do
end do
!$omp end do
!$omp end parallel
return
end subroutine potentialcalcu

subroutine concen_calcu
use LBMset, ONLY: f_nega, f_posi, con_n, con_p
use grid, ONLY: nxp, nyp
use physical_constant, ONLY:zeta, bulk_con
implicit none
integer::i,j,k
!$omp parallel private(i,j)
!$omp do 
do j=1,nyp-1
	do i=0,nxp
		! do k=0,4
			! ssum_con=ssum_con+f_posi(k,i,j)
		! end do
		! con_p(i,j)=ssum_con
		con_p(i,j)=f_posi(0,i,j)+f_posi(1,i,j)+f_posi(2,i,j)+f_posi(3,i,j)+f_posi(4,i,j)
	end do
end do
!$omp end do
!$omp do 
do j=1,nyp-1
	do i=0,nxp
		! do k=0,4
			! ssum_con=ssum_con+f_nega(k,i,j)
		! end do
		!con_n(i,j)=ssum_con
		con_n(i,j)=f_nega(0,i,j)+f_nega(1,i,j)+f_nega(2,i,j)+f_nega(3,i,j)+f_nega(4,i,j)
	end do
end do
!$omp end do
!$omp end parallel
return
end subroutine concen_calcu

subroutine rhouv
use LBMset , ONLY: f, u, v, ck, rho
use grid, ONLY: nxp, nyp
implicit none
integer::i,j,k
!$omp parallel private(i,j,k)
!$omp do 
do j=0,nyp
do i=0,nxp
rho(i,j)=f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
end do
end do
!$omp end do
!$omp do 
DO j=1,nyp-1
DO i=0,nxp
u(i,j)=(f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j))*ck/rho(i,j)
v(i,j)=(f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j))*ck/rho(i,j)
END DO
END DO
!$omp end do
!$omp end parallel
return
end subroutine rhouv
subroutine gradientP_x()  
use LBMset, ONLY: psi, gradP_x
use grid, ONLY: dx, nxp, nyp
implicit none
integer :: i,j,k
real*8 :: c_2_3_dx
real*8 :: c_12_dx
c_2_3_dx=2.0D0/3.0D0/dx
c_12_dx=1.0D0/12.0D0/dx
! 4th order FDM
do j=0,nyp
	i=0
	gradP_x(i,j)=-c_12_dx*(psi(i+2,j)-psi(nxp-2,j))+c_2_3_dx*(psi(i+1,j)-psi(nxp-1,j))
	i=1
	gradP_x(i,j)=-c_12_dx*(psi(i+2,j)-psi(nxp-1,j))+c_2_3_dx*(psi(i+1,j)-psi(i-1,j))
	do i=2,nxp-2
		gradP_x(i,j)=-c_12_dx*(psi(i+2,j)-psi(i-2,j))+c_2_3_dx*(psi(i+1,j)-psi(i-1,j))
	end do
	i=nxp-1
	gradP_x(i,j)=-c_12_dx*(psi(1,j)-psi(i-2,j))+c_2_3_dx*(psi(i+1,j)-psi(i-1,j))
	i=nxp
	gradP_x(i,j)=-c_12_dx*(psi(2,j)-psi(i-2,j))+c_2_3_dx*(psi(1,j)-psi(i-1,j))
end do
end subroutine gradientP_x
subroutine gradientP_y()
use LBMset, ONLY: psi, gradP_y
use grid, ONLY: dx, nxp, nyp
implicit none
integer :: i,j,k
real*8 :: c_12_0=11.0D0/6.0D0 
! way1 4th order explicit scheme
do i=0,nxp
	j=0
	gradP_y(i,j)=(psi(i,j+3)/3.0D0-1.50D0*psi(i,j+2)+3.0D0*psi(i,j+1)-c_12_0*psi(i,j))/dx
	j=1
	gradP_y(i,j)=(-psi(i,j+2)/6.0D0+psi(i,j+1)-0.50D0*psi(i,j)-psi(i,j-1)/3.0D0)/dx
	do j=2,nyp-2
		gradP_y(i,j)=-(psi(i,j+2)-psi(i,j-2))/12.0D0/dx+2.0D0*(psi(i,j+1)-psi(i,j-1))/3.0D0/dx
	end do
	j=nyp-1
	gradP_y(i,j)=(psi(i,j-2)/6.0D0-psi(i,j-1)+0.50D0*psi(i,j)+psi(i,j+1)/3.0D0)/dx
	j=nyp
	gradP_y(i,j)=(c_12_0*psi(i,j)-3.0D0*psi(i,j-1)+1.50D0*psi(i,j-2)-psi(i,j-3)/3.0D0)/dx
end do
end subroutine gradientP_y
