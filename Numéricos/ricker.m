function [ricker]=ricker(np,fr,dt);

npt=np*dt;
%t=(-npt/2):dt:(npt/2);
t=[1:np]*dt; %0:dt:npt;
t0=0.0; %(npt/2);
T=1/fr;
gama=4;
ni=0;

ricker=exp(-(pi*fr*(t-t0)).^2).*(1-2*pi^2*fr^2*(t-t0).^2);


