!This program is simulating the arrangement of clusters
!Adding one positive charge at the center in each cell, fixing the position
!The potential is given by -alpha*e^2/(2*epsilon^2*r^4)
!Calculating the RDF around the charge
!Using NVT ensemble
Program electroMD
implicit real(8)(a-h,o-z)
!epsilpn/kB=10.22K sigma=2.556 angstroms for helium
!mHe=6.65*10^-27 Kg.
!MD time unit=1.755 ps
!dt=3.5 fs
!MD pressure unit=84.4 bars
!dt=0.002 in MD units
real(8), parameter:: k_B=1.3806488d0 !Boltzmann constant in 10^-23 J/K
real(8), parameter:: m_He=6.65d0     !MD mass unit in 10^-24 gram
real(8), allocatable:: x(:),y(:),z(:)
real(8), allocatable:: fullx(:),fully(:),fullz(:)
real(8), allocatable:: vx(:), vy(:), vz(:)
real(8), allocatable:: fx(:), fy(:), fz(:)
real(8), allocatable:: g(:), rg(:)
real(8), allocatable:: disp(:,:)
real(8):: rcharge(3,27)   !only pick nearest ions.
real(8):: L
integer(4):: switch !controls the operation of RDF subroutine and lindemann_index subroutine
real(8), allocatable:: radSquared(:)

real(8):: initialx,initialy,initialz


!read input parameters***************************************
open(1,file='electroMD.inp',status='old')
read(1,*) Tunit !Epsilon/k_B
read(1,*) dunit !sigma
read(1,*) natm
read(1,*) rho   !initial density
read(1,*) T     !target temp.
read(1,*) num_equilibrium  
read(1,*) num_sample
read(1,*) dt
read(1,*) Rcut
close(1)
!*************************************************************
punit=k_B*Tunit*100.d0/(dunit**3)  !MD pressure unit
write(*,'("The pressure unit is:",F7.3)') punit
pressure=pressure/punit 
T=T/Tunit
rho=(dunit**3)/m_He*rho !convert to MD units
write(*,'("The initial density in MD units:",F7.4)') rho
L=(natm/rho)**(1.d0/3.d0)          !The dimension of box by initial number density
num=0
do i=1,3
    do j=1,3
        do k=1,3
            xc=-0.5d0*L+(i-1)*L
            yc=-0.5d0*L+(j-1)*L
            zc=-0.5d0*L+(k-1)*L
            num=num+1
            rcharge(1,num)=xc
            rcharge(2,num)=yc
            rcharge(3,num)=zc
        end do
    end do
end do

write(*,'("The initial size of the box:",F7.3)') L

T_avg=0.d0
nhis=500      !Num of bins used for RDF calculation
allocate(x(1:natm),y(1:natm),z(1:natm), &
fullx(1:natm),fully(1:natm),fullz(1:natm) , &
vx(1:natm),vy(1:natm),vz(1:natm),  &
fx(1:natm),fy(1:natm),fz(1:natm),g(1:nhis), rg(1:nhis),radSquared(1:natm))

allocate(disp(2,natm*(natm-1)/2))
call initialcondition_fcc(x,y,z,vx,vy,vz,natm,L,T)
call eforce(x,y,z,rcharge,fx,fy,fz,natm,L,Ep,Rcut)

open(2,file='fcc.dat',status='unknown')
write(2,'(3F7.3)') (x(i),y(i),z(i),i=1,natm)
close(2)

open(7,file='equilibrium.dat',status='unknown')
do i=1,num_equilibrium      !equilibrium time
   call  verlet(x,y,z,vx,vy,vz,fx,fy,fz,natm,dt,L,rcharge,Ek,Ep,Rcut)
   call  PBC(x,y,z,natm,L)
   call Berendsen_thermostat(vx,vy,vz,natm,dt,T,Ek)
   Tint=Ek/3.d0*2.0d0
   T_avg=T_avg+Tint*Tunit
   if(mod(i,2000)==0) write(*,'(I7,2F12.3)') i,T_avg/dfloat(i),Ep+Ek
   if(mod(i,1000)==0) write(7,'(I7,2F12.3)')  i,Tint,Ek+Ep 
end do
close(7)
T_avg=0.d0
switch=0
rho=natm/L**3
call RDF(x,y,z,natm,L,rho,rg,g,switch,nhis,ngr,delg)  !initialization of RDF
switch=1
print*, "The equilibrium time has finished."
 initialx = x(1)
 initialy = y(1)
 initialz = z(1)
do i=1,num_sample      !sample time
   call  verlet(x,y,z,vx,vy,vz,fx,fy,fz,natm,dt,L,rcharge,Ek,Ep,Rcut)
   call  PBC(x,y,z,natm,L)
   call  RDF(x,y,z,natm,L,rho,rg,g,switch,nhis,ngr,delg)
   call Berendsen_thermostat(vx,vy,vz,natm,dt,T,Ek)
   Tint=Ek/3.d0*2.0d0
   T_avg=T_avg+Tint*Tunit
  
  !radSquared(i)= (initialx-x(1))**2 + (initialy-y(1))**2 + (initialz-z(1))**2


   if(mod(i,5000) == 0) write(*,'(I7,2F12.3)') i,T_avg/dfloat(i),Ep+Ek
end do

  open(17,file='radSquared.dat',status='unknown')
  write(17,'(F7.3)') (radSquared(i),i=1,natm)
  close(17)
  
switch=2
call RDF(x,y,z,natm,L,rho,rg,g,switch,nhis,ngr,delg) 
T_avg=T_avg/num_sample

open(3,file='RDF.dat',status='unknown')
write(3,'(2F7.3)') (rg(i),g(i),i=1,nhis)
close(3)

open(4,file='geo.dat',status='unknown')
write(4,'(3F7.3)') (x(i),y(i),z(i),i=1,natm)
close(4)

open(8,file='Temp.dat',status='unknown')
write(8,'(2F9.3)') T_avg,rho*m_He/(dunit**3)
close(8) 
end Program electroMD


!Subroutine to calculate forces
subroutine force(x,y,z,fx,fy,fz,natm,L,Ep,virial,Rcut)
implicit real(8)(a-h,o-z)
real(8):: L
real(8):: x(natm),y(natm),z(natm)
real(8):: fx(natm), fy(natm), fz(natm)
pi=dacos(-1.d0)
fx=0.d0
fy=0.d0
fz=0.d0
epsilon_0=1.0d0; sigma=1.0d0  ! MD units
sigma2=sigma*sigma
Rcut2=Rcut*Rcut
Rcut_3=sigma2/Rcut2*sigma/Rcut
Rcut_6=Rcut_3*Rcut_3
Ep=0.d0
virial=0.d0
do i=1,natm-1
  do j=i+1,natm
     dx=x(i)-x(j)
	 dy=y(i)-y(j)
	 dz=z(i)-z(j)
	 !minimum image criterion
	 if(dabs(dx) > 0.5d0*L) dx=dx-sign(L,dx)
	 if(dabs(dy) > 0.5d0*L) dy=dy-sign(L,dy)
	 if(dabs(dz) > 0.5d0*L) dz=dz-sign(L,dz)
	 r2=dx*dx+dy*dy+dz*dz
	 if(r2 < Rcut2) then
	   fr2=sigma2/r2
       fr4=fr2*fr2
	   fr6=fr2*fr2*fr2
	   funit=48.d0*epsilon_0*fr6*(fr6-0.5d0)*fr2 
	   fxi=funit*dx
	   fyi=funit*dy
	   fzi=funit*dz

       fx(i)=fx(i)+fxi; fx(j)=fx(j)-fxi
	   fy(i)=fy(i)+fyi; fy(j)=fy(j)-fyi
	   fz(i)=fz(i)+fzi; fz(j)=fz(j)-fzi
	   Ep=Ep+4.d0*epsilon_0*fr6*(fr6-1.d0) 
	   virial=virial+funit*r2
	 end if
	end do
end do
Ep=Ep+8.d0*pi*natm**2/3.d0/L**3*(1.d0/3.d0*Rcut_6-1.d0)*Rcut_3 !tail correction for potential energy
Ep=Ep/natm
virial=virial/natm
end
!***************************************************************************************************

!Subroutine eforce including electrostatic forces.

subroutine eforce(x,y,z,rcharge,fx,fy,fz,natm,L,Ep,Rcut)
implicit real(8)(a-h,o-z)
real(8):: L
real(8):: x(natm),y(natm),z(natm)
real(8):: fx(natm), fy(natm), fz(natm)
real(8):: rcharge(3,27)
pi=dacos(-1.d0)
fx=0.d0
fy=0.d0
fz=0.d0
epsilon_0=1.0d0; sigma=1.0d0  ! MD units
E_cons=39.1921d0             !The coefficient in MD units
sigma2=sigma*sigma
Rcut2=Rcut*Rcut
Rcut_3=sigma2/Rcut2*sigma/Rcut
Rcut_6=Rcut_3*Rcut_3
Ep=0.d0
do i=1,natm-1
   do k=1,27              !calculating the electrostatic force
     dx=x(i)-rcharge(1,k)
     dy=y(i)-rcharge(2,k)
     dz=z(i)-rcharge(3,k)
     r2=dx*dx+dy*dy+dz*dz
     fr2=sigma2/r2
     fr4=fr2*fr2
     fr6=fr2*fr2*fr2
     fe=-4.d0*E_cons*fr6
     fx(i)=fx(i)+fe*dx
     fy(i)=fy(i)+fe*dy
     fz(i)=fz(i)+fe*dz
     Ep=Ep-E_cons*fr4
	 if(r2 < Rcut2) then  !the LJ interaction
	   funit=48.d0*epsilon_0*fr6*(fr6-0.5d0)*fr2 
       fxi=funit*dx
	   fyi=funit*dy
	   fzi=funit*dz
	   fx(i)=fx(i)+fxi
       fy(i)=fy(i)+fyi
       fz(i)=fz(i)+fzi
       Ep=Ep+4.d0*epsilon_0*fr6*(fr6-1.d0) 
	  end if
   end do
  do j=i+1,natm
     dx=x(i)-x(j)
	 dy=y(i)-y(j)
	 dz=z(i)-z(j)
	 !minimum image criterion
	 if(dabs(dx) > 0.5d0*L) dx=dx-sign(L,dx)
	 if(dabs(dy) > 0.5d0*L) dy=dy-sign(L,dy)
	 if(dabs(dz) > 0.5d0*L) dz=dz-sign(L,dz)
	 r2=dx*dx+dy*dy+dz*dz
	 if(r2 < Rcut2) then
	   fr2=sigma2/r2
       fr4=fr2*fr2
	   fr6=fr2*fr2*fr2
	   funit=48.d0*epsilon_0*fr6*(fr6-0.5d0)*fr2 
	   fxi=funit*dx
	   fyi=funit*dy
	   fzi=funit*dz

       fx(i)=fx(i)+fxi; fx(j)=fx(j)-fxi
	   fy(i)=fy(i)+fyi; fy(j)=fy(j)-fyi
	   fz(i)=fz(i)+fzi; fz(j)=fz(j)-fzi
	   Ep=Ep+4.d0*epsilon_0*fr6*(fr6-1.d0) 
	 end if
	end do
end do
do k=1,27              !calculating the electrostatic force for the last one
     dx=x(natm)-rcharge(1,k)
     dy=y(natm)-rcharge(2,k)
     dz=z(natm)-rcharge(3,k)
     r2=dx*dx+dy*dy+dz*dz
     fr2=sigma2/r2
     fr4=fr2*fr2
     fr6=fr2*fr2*fr2
     fe=-4.d0*E_cons*fr6
     fx(natm)=fx(natm)+fe*dx
     fy(natm)=fy(natm)+fe*dy
     fz(natm)=fz(natm)+fe*dz
     Ep=Ep-E_cons*fr4
	 if(r2 < Rcut2) then
	   funit=48.d0*epsilon_0*fr6*(fr6-0.5d0)*fr2 
	   fxi=funit*dx
	   fyi=funit*dy
	   fzi=funit*dz
	   fx(natm)=fx(natm)+fxi
       fy(natm)=fy(natm)+fyi
       fz(natm)=fz(natm)+fzi
       Ep=Ep+4.d0*epsilon_0*fr6*(fr6-1.d0) 
	 end if
end do
Ep=Ep+8.d0*pi*natm**2/3.d0/L**3*(1.d0/3.d0*Rcut_6-1.d0)*Rcut_3 !tail correction for potential energy
Ep=Ep/natm
end
!*************************************************************************************

!Subroutine to apply boundary conditions
subroutine PBC(x,y,z,natm,L)
implicit real(8)(a-h,o-z)
real(8):: L
real(8):: x(natm), y(natm), z(natm)
!real(8):: fullx(natm), fully(natm), fullz(natm)
do i=1,natm
  if(x(i) < 0.d0)  then
  x(i)=x(i)+L
  !fullx(i)= fullx(i)-L
  end if

  if(x(i) > L) then
  x(i)=x(i)-L
 !fullx(i)= fullx(i)+L
  end if

  if(y(i) < 0.d0) then
   y(i)=y(i)+L
  !fully(i)= fully(i)-L
  end if

  if(y(i) > L) then
   y(i)=y(i)-L
  !fully(i)= fully(i)+L
   end if

  if(z(i) < 0.d0) then
   z(i)=z(i)+L
   !fullz(i)= fullz(i)-L
  end if
  if(z(i) > L) then
  z(i)=z(i)-L
  !fullz(i)= fullz(i)+L
  end if
end do
end
!***********************************************************************

!Subroutine to do velocity-verlet step
subroutine verlet(x,y,z,vx,vy,vz,fx,fy,fz,natm,dt,L,rcharge,Ek,Ep,Rcut)
implicit real(8)(a-h,o-z)
real(8):: L
real(8):: x(natm), y(natm), z(natm)
!real(8):: fullx(natm), fully(natm), fullz(natm)
real(8):: vx(natm), vy(natm), vz(natm)
real(8):: fx(natm), fy(natm), fz(natm)
real(8):: rcharge(3,27)
pi=dacos(-1.d0)
half_dt=0.50*dt
!Rcut=2.5d0
sigma=1.d0
Rcut_3=sigma**3/Rcut**3
Rcut_6=Rcut_3*Rcut_3
do i=1,natm						!v(t+dt/2)
   vx(i)=vx(i)+fx(i)*half_dt
   vy(i)=vy(i)+fy(i)*half_dt
   vz(i)=vz(i)+fz(i)*half_dt
   x(i)=x(i)+vx(i)*dt	        !r(t+dt)
   !fullx(i) = fullx(i)+vx(i)*dt 
   y(i)=y(i)+vy(i)*dt
   !fully(i) = fully(i)+vy(i)*dt
   z(i)=z(i)+vz(i)*dt
   !fullz(i) = fullz(i)+vz(i)*dt
end do

call eforce(x,y,z,rcharge,fx,fy,fz,natm,L,Ep,Rcut)  !f(t+dt)
Ek=0.d0
do i=1,natm					   !v(t+dt)
   vx(i)=vx(i)+fx(i)*half_dt
   vy(i)=vy(i)+fy(i)*half_dt
   vz(i)=vz(i)+fz(i)*half_dt
   Ek=Ek+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
end do
Ek=Ek*0.5d0/natm
end
!************************************************************************


!Subroutine to do inialization simple cubic lattice
subroutine initialcondition_sc(x,y,z,vx,vy,vz,natm,L,T)
implicit real(8)(a-h,o-z)
real(8):: L
real(8):: x(natm), y(natm), z(natm)
real(8):: vx(natm), vy(natm), vz(natm)
sum_vx=0.d0
sum_vy=0.d0
sum_vz=0.d0
sum_v2=0.d0
n=idnint(natm**(1.d0/3.d0))
distance=L/dfloat(n-1)

num=0

do i=1,n						!cubic lattice position and random initial velocities
  do j=1,n
    do k=1,n
	  x1=(i-0.5d0)*distance
	  y1=(j-0.5d0)*distance
	  z1=(k-0.5d0)*distance
	  num=num+1
	  x(num)=x1
	  y(num)=y1
	  z(num)=z1
	  vx(num)=rand()-0.5d0
	  vy(num)=rand()-0.5d0
	  vz(num)=rand()-0.5d0
	  sum_vx=sum_vx+vx(num)
	  sum_vy=sum_vy+vy(num)
	  sum_vz=sum_vz+vz(num)
	  sum_v2=sum_v2+vx(num)*vx(num)+vy(num)*vy(num)+vz(num)*vz(num)
	end do
   end do
end do

vx_avg=sum_vx/natm
vy_avg=sum_vy/natm
vz_avg=sum_vz/natm
sum_v2=sum_v2/natm-vx_avg**2/natm-vy_avg**2/natm-vz_avg**2/natm
fs=dsqrt(3.d0*T/sum_v2)

do i=1,natm							!remove the velocity shift and scale by temperature
   vx(i)=fs*(vx(i)-vx_avg)
   vy(i)=fs*(vy(i)-vy_avg)
   vz(i)=fs*(vz(i)-vz_avg)
end do

end
!*********************************************************************************

!Subroutine to do inialization fcc lattice
subroutine initialcondition_fcc(x,y,z,vx,vy,vz,natm,L,T)
implicit real(8)(a-h,o-z)
real(8):: L
real(8):: x(natm), y(natm), z(natm)
real(8):: vx(natm), vy(natm), vz(natm)
real(8):: rcell(3,4)
sum_vx=0.d0
sum_vy=0.d0
sum_vz=0.d0
sum_v2=0.d0
n=idnint((natm/4.d0)**(1.d0/3.d0))
write(*,'("Num of cells along each direction:",I2)') n
dlat=L/(n+0.5d0) !size of unit cell of fcc lattice
rcell=reshape((/ 0.0d0 , 0.0d0 , 0.0d0 , &
                 0.5d0 , 0.5d0 , 0.0d0 , &
                 0.0d0 , 0.5d0 , 0.5d0 , &
                 0.5d0 , 0.0d0 , 0.5d0 /) , (/ 3,4 /) )
                 
num=0

do i=1,n						!cubic lattice position and random initial velocities
  do j=1,n
    do k=1,n
        do m=1,4
	      xcell=(i-1)*dlat
	      ycell=(j-1)*dlat
          zcell=(k-1)*dlat
          num=num+1
          x(num)=xcell+rcell(1,m)*dlat
          y(num)=ycell+rcell(2,m)*dlat
          z(num)=zcell+rcell(3,m)*dlat
          vx(num)=rand()-0.5d0
          vy(num)=rand()-0.5d0
	      vz(num)=rand()-0.5d0
	      sum_vx=sum_vx+vx(num)
	      sum_vy=sum_vy+vy(num)
	      sum_vz=sum_vz+vz(num)
	      sum_v2=sum_v2+vx(num)*vx(num)+vy(num)*vy(num)+vz(num)*vz(num)
        end do
	end do
   end do
end do

vx_avg=sum_vx/natm
vy_avg=sum_vy/natm
vz_avg=sum_vz/natm
sum_v2=sum_v2/natm-vx_avg**2/natm-vy_avg**2/natm-vz_avg**2/natm
fs=dsqrt(3.d0*T/sum_v2) !scaling factor by temperature

do i=1,natm			!remove the velocity shift and scale by temperature
   vx(i)=fs*(vx(i)-vx_avg)
   vy(i)=fs*(vy(i)-vy_avg)
   vz(i)=fs*(vz(i)-vz_avg)
end do

end

!****************************************************************************************

!Subroutine Berendsen(vx,vy,vz,natm,Ek)
subroutine Berendsen_thermostat(vx,vy,vz,natm,dt,T,Ek)
!tau constant is 1 ps
implicit real(8)(a-h,o-z)
real(8):: vx(natm),vy(natm),vz(natm)
real(8):: lambda
tau=0.6d0
Tint=Ek*2.d0/3.d0
lambda=dsqrt(dabs(1+dt/tau*(T/Tint-1.d0)))
do i=1,natm
  vx(i)=lambda*vx(i)
  vy(i)=lambda*vy(i)
  vz(i)=lambda*vz(i)
end do
end
!****************************************************************

!Subroutine for Berendsen barostat
subroutine Berendsen_barostat(x,y,z,pressure_int,pressure,natm,dt,L)
implicit real(8)(a-h,o-z)
real(8):: x(natm), y(natm), z(natm)
!real(8):: fullx(natm), fully(natm), fullz(natm)
real(8):: L
real(8):: kai
tau=0.3d0
kai=dabs((1-dt/tau*(pressure-pressure_int))**(1.d0/3.d0))
L=L*kai
do i=1,natm
    x(i)=x(i)*kai
	!fullx(i)=fullx(i)*kai
    y(i)=y(i)*kai
	!fully(i)=fully(i)*kai
    z(i)=z(i)*kai
	!fullz(i)=fullz(i)*kai

end do    
end    
!*****************************************************************************


!Subroutine to calculate radial distribution function
subroutine RDF(x,y,z,natm,L,rho,rg,g,switch,nhis,ngr,delg)
!switch indicates oerations: initialization, sample, determination
!ngis num of bins
!delg width of bin
!ngr sample numbers
!maximum interatomic distance L/2    
implicit real(8)(a-h,o-z)
integer(4):: switch
real(8):: L
real(8):: x(natm), y(natm), z(natm)
real(8):: g(nhis),rg(nhis)
pi=dacos(-1.d0)
half_L=0.5d0*L
if (switch == 0) then    !Initialization
    ngr=0
    delg=L/(2.d0*nhis)
    rg(1:nhis)=0.d0
    g(1:nhis)=0.d0
else if (switch == 1) then
    ngr=ngr+1
    do i=1,natm
          dx=x(i)-0.5d0*L   !RDF about the center of box
          dy=y(i)-0.5d0*L
          dz=z(i)-0.5d0*L
          r=dsqrt(dx*dx+dy*dy+dz*dz)
          if(r < half_L) then
              ig=idint(r/delg)+1 !notice the difference between nint() and int()
              g(ig)=g(ig)+1
          end if
    end do
else if (switch == 2) then
    do i=1,nhis
        rg(i)=(i-0.5d0)*delg
        vbin=(i**3-(i-1)**3)*delg**3
        fid=4.d0/3.d0*pi*vbin*rho   !normalization factor
        g(i)=g(i)/(ngr*fid)
    end do
end if
return
end
!********************************************************************************

!Subroutine to calculate Lindemann index
subroutine lindemann_index(x,y,z,natm,L,switch,ntime,disp,q_index)
implicit real(8)(a-h,o-z)
integer(4):: switch
real(8):: x(natm),y(natm),z(natm)
real(8):: L
real(8):: disp(2,natm*(natm-1)/2)

if (switch==0) then
   ntime=0
   disp=0.d0
   q_index=0.d0
end if 

if (switch==1) then
   ntime=ntime+1
   num=0
   do i=1,natm-1
     do j=i+1,natm
	    dx=x(i)-x(j)
        dy=y(i)-y(j)
        dz=z(i)-z(j)
	 !minimum image criterion
	 if(dabs(dx) > 0.5d0*L) dx=dx-sign(L,dx)
	 if(dabs(dy) > 0.5d0*L) dy=dy-sign(L,dy)
	 if(dabs(dz) > 0.5d0*L) dz=dz-sign(L,dz)
        r=dsqrt(dx*dx+dy*dy+dz*dz)
		num=num+1
		disp(1,num)=disp(1,num)+r
		disp(2,num)=disp(2,num)+r*r
	  end do
	end do
end if

if (switch==2) then
   do i=1,natm*(natm-1)/2
       disp(1,i)=disp(1,i)/ntime
       disp(2,i)=disp(2,i)/ntime
       q_index=q_index+dsqrt(disp(2,i)-disp(1,i)**2)/disp(1,i)
   end do
end if
end
!*************************************************************************************
