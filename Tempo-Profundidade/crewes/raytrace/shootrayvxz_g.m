function [t,r]=shootrayvxz_g(tstep,r0)
% SHOOTRAYVXZ_G more general raytracing in v(x,z).
%
% [t,r]=shootrayvxz_g(tstep,t0)
%
% General ray tracer for a 2D gridded velocity model. The model
% is built by defining a matrix of velocity values and calling
% RAYVELMOD. This establishes 3 global matrices containing v^2
% and the logarithmic derivatives of v with respect to x and z.
% (The velocity model matrix may be deleted after calling RAYVELMOD
% if memory is limited.) The raytracer implements an RK4 (4th order 
% Runge-Kutta) solution to the ray equations and, by default, uses
% nearest neighbor interpolation. Bilinear interpolation is available
% for more accuracy. To get bilinear interpolation, define a global
% variable called BILINEAR (all caps) and set its value to 1. This
% is quite a bit slower than nearest neighbor for the same grid size.
% If using nearest neighbor only, then SHOOTRAYVXZ is preferred
% because it is much faster. Rays will be terminated when they reach
% the maximum time or leave the bounds of the velocity model.
%
% tstep ... vector of time steps, (e.g [0:.01:tmax] steps from 0 to
%           tmax in 10 millisecond intervals.)
% r0 ... initial values of the ray vector. r0(1) and r0(2) are the
%        starting x and z coordinates and r0(3) and r0(4) are the
%        horizontal and vertical slownesses.
% t ... output time vector. If the ray stays within the bounds of
%        the model, this is the same as tstep, otherwise it may be
%        smaller.
% r ... output ray matrix. This is an N by 4 matrix, where N=length(t),
%       with each row being a ray vector for the corresponding time. The
%       ray vectors are as described in r0 above.
%
% G.F. Margrave and P.F. Daley, CREWES, June 2000
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

%

global RTV2 RTDLNVDX RTDLNVDZ RTDG BILINEAR RTX RTZ

if(isempty(RTDG))
	error('velocity model not defined. Use RAYVELMOD')
end


[m,n]=size(RTV2);
xmax=RTX(n-2);% set to abort within 2 dg of edge
zmax=RTZ(m-2);%
xmin=RTX(3);
zmin=RTZ(3);

r0=r0(:)';

r=zeros(length(tstep),length(r0));
r(1,:)=r0;

for k=2:length(tstep)
	tnow=tstep(k-1);
	rnow=r(k-1,:);
	h=tstep(k)-tstep(k-1);
	k1=h*drayvec(tnow,rnow)';
	k2=h*drayvec(tnow+.5*h,rnow+.5*k1)';
	k3=h*drayvec(tnow+.5*h,rnow+.5*k2)';
	k4=h*drayvec(tnow+h,rnow+k3)';
	r(k,:)=rnow+k1/6+k2/3+k3/3+k4/6;
	if(r(k,1)>xmax | r(k,2)>zmax | r(k,1)<xmin | r(k,2)<zmin)
		break
	end
end

if(k<length(tstep))
   r(k+1:length(tstep),:)=[];
   t=tstep(1:k);
else
   t=tstep;
end