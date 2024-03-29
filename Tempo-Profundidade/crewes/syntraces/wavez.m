function [w,tw,wd]=wavez(dt,fdom,tlength,m)
% WAVEZ: creates a zero phase w with a dominant frequency
%
% [w,tw]=wavez(dt,fdom,tlength,m)
% [w,tw]=wavez(dt,fdom,tlength)
% [w,tw]=wavez(dt,fdom) 
% [w,tw]=wavez(dt) 
% 
% WAVEZ returns a zero phase w with a realistic amplitude spectrum
%
% dt= desired temporal sample rate
% fdom= dominant frequency in Hz (default: 15 Hz)
% tlength= w length in seconds (default: 127*dt 
%                                     (ie a power of 2))
% m = exponent controlling spectral shape. See tntamp for a description
% ************ default 2 ************
% 
% 
% by G.F. Margrave, July 1991
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

if(nargin<4) m=2; end
 if(nargin<3)
   tlength=127.*dt;
 end
 if(nargin<2)
   fdom=15.; 
 end

% create a time vector
  nt=tlength/dt+1;
  nt=2.^nextpow2(nt);
  tmax=dt*(nt)/2;
  tw=-tmax+dt*(0:nt-1)';;
  fnyq=1./(2*(tw(2)-tw(1)));
  f=linspace(0.,fnyq,length(tw)/2+1)';
% create the amplitude spectrum spectrum
  aspec=tntamp(fdom,f,m);
% inverse transform
  w=ifftrl(aspec,f);
  w=fftshift(w);
  wd=ifftrl(i*f.*aspec,f);
  wd=fftshift(wd);
% make sure we don't return the Fourier pad if it was not asked for
  ind=between(-tlength/2,tlength/2,tw,2);
  w=w(ind);
  tw=tw(ind);
  wd=wd(ind);
% now normalize the w
  w=wavenorm(w,tw,2);