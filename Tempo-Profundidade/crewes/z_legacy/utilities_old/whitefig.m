function whitefig
% WHITEFIG Graphics figure default background color changed to white.
%
% WHITEFIG changes the current figure's default background color
%          to white as well as other properties so graphs look good.
%
%	G.F. Margrave
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
set(gcf,'color','white')
set(gcf,'defaulttextcolor','black')
set(gcf,'inverthardcopy','off')
%find all axes
hkids=get(gcf,'children');
axid=zeros(size(hkids));
for k=1:length(hkids)
	if(strcmp(get(hkids(k),'type'),'axes'))
		axid(k)=1;
	end
end
haxes=hkids(find(axid==1));
for kax=1:length(haxes)
	hax=haxes(kax);
	set(hax,'xcolor','black')
	set(hax,'ycolor','black')
	set(hax,'zcolor','black')
	htit=get(hax,'title');
	if( htit>0 )
		set(htit,'color','black');
	end
	h=get(hax,'children');
	for k=1:length(h)
		if( strcmp(get(h(k),'type'),'text') |  ...
			strcmp(get(h(k),'type'),'Text') )
			kol=get(h(k),'color');
			if( sum(kol)==3 )
				set(h(k),'color','black');
			end
		end
		if( strcmp(get(h(k),'type'),'line') | ...
			strcmp(get(h(k),'type'),'Line') )
			kol=get(h(k),'color');
			if( sum(kol)==3 )
				set(h(k),'color','black');
			end
		end
	
	end
end