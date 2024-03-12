program velocidades

integer,parameter::nx=720,nz=160,ntime=1750,op=1
real,parameter::dx=50.0,dz=50.0,dt=0.004

real,dimension(nz,nx)::modevel
real,allocatable,dimension(:,:)::tti,velRMS

allocate(tti(nz,nx),velRMS(ntime,nx))

do iz=1,nz
do ix=1,nx
tti(iz,ix)=0.d0
end do
end do

do j=1,ntime
do i=1,nx
velRMS(j,i)=0.d0
end do
end do

OPEN(20,FILE='vel_model.dat',STATUS='old')
READ(20,*) ((modevel(iz,ix),iz=1,nz),ix=1,nx)
CLOSE(UNIT=20) 

call velo_tti_layers(modevel,nx,nz,dx,dz,tti,op)

call velo_rms_layers(modevel,tti,nx,nz,ntime,dx,dz,dt,velRMS,op)

print*,'Escrevendo saída para o arquivo!'

open(22,file='vel_tti.dat',status='unknown') 
write(22,'(f10.5)') ((tti(j,i),j=1,nz),i=1,nx)
close(22,status='keep')

open(24,file='vel_rms.dat',status='unknown') 
write(24,'(e10.2)') ((velRMS(j,i),j=1,ntime),i=1,nx)
close(24,status='keep')

deallocate(tti,velRMS)

end program velocidades