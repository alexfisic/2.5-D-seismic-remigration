function fact=plotseis(smat,t,x,vaflag,fact,flipy,kolor,hax)
%
% fact=plotseis(smat,t,x,vaflag,fact,flipy,kolor,hax)
%
% PLOTSEIS does a quick plot of a seismic matrix in a figure window (made by plotseis).
% It plots the columns as wiggle traces centered at their x coordinates.
% Type plotseis with no arguments for a demo.
%
%       smat ... the seismic matrix to be plotted. Traces are assumed stored in
%       the columns smat.
%       t ... time coordinates of traces.
%       ****** default 0:.002:(nrows-1)*.002 where nrows = number of rows
%	in smat ****
%       x ... x coordinates of the traces
%       ****** default 1:ncols where ncols=number of columns *****
%       vaflag ... if 0, then traces are plotted as wt (wiggle traces)
%               if 1, then the traces are plotted wtva (shaded peaks)
%       if 2 then reverve polarity wtva (shaded troughs)
%               ************* default = 1 *********
%       fact ... scaling factors. Make fact(1) bigger for bigger wiggle
%               traces. fact(2) controls the overall plot scale. If not
%               provided, it is computed as max(abs(smat)). To scale two
%               plots the same, capture the return value from the first and
%               provide is as fact(2) for the second
%               ************* default 1.5 ***********
%       flipy ... if 1, then the y axis is reversed so that it goes from top
%               of window to the bottom
%          ************* default = 1 ***********
%  kolor ... color to plot the traces
%     ************* default = [1 0 0] (red) ************
%  hax ... handle of the axis to plot the seismic matrix in. If -1, then a
%       new figure is created.
%     ************* default = -1 **************
%
% NOTE: To scale two plots with respect to the maximum absolute value
% on the first plot, capture the return value from the first
% and provide it as the fifth argument for the second plot:
%       fact=plotseis(smat1,t1,x1);
%       plotseis(smat2,t2,x2,1,fact);
%
% G.F. Margrave, Aug 1995
% The CREWES Project
% Department of Geology and Geophysics
% The University of Calgary
% gary@geo.ucalgary.ca
%
if(nargin<1)
        % do a demo
                %Make a fake reflectivity
                t=0:.002:1.0;
                r=randn(size(t)).^5;
                %make a ricker wavelet
                tw=-.1:.002:.1;
                arg=(pi*15*tw).^2;
                w=(1-2.*arg).*exp(-arg);
                %convolve
                s=conv(r,w);
                s=s(51:length(t)+50)';
                s=s/max(s); %normalize

                smat=s*ones(1,20);
end
if(nargin < 8)
        hax=-1;
end
if(nargin<7)
        kolor=[0 0 0];
end
if(nargin<6)
        flipy=1;
end
if(nargin<5)
        fact=1.5;
end
if(nargin<4)
        vaflag=1;
end
if(nargin<3)
        ncols=size(smat,2);
        x=1:ncols;
end
if(nargin<2)
        nrows=size(smat,1);
        t=1:nrows;       
% t=0:.002:(nrows-1)*.002;
end

if(length(x)>1)
        bnds=(max(x)-min(x))/(length(x)+1);
else
        bnds=max(smat-min(smat))/2;
end

if( length(fact)<2 )
        ilive=find(~isnan(smat));
        trcmin=min(smat(ilive));
        trcmax=max(smat(ilive));
        s=max([abs(trcmax) abs(trcmin)]);
        fact=[fact s];
end

%if( hax==-1 )
%        figure;
%else
%        hfig=get(hax,'parent');
%        figure(hfig);
%        set(hfig,'currentaxes',hax);
%end

%subplot(313)

if(flipy)
        set(gca,'ydir','reverse');
end

%for k=1:length(x)
for k=length(x):-1:1
        trc=smat(:,k)/fact(2);
        ilive=find(~isnan(trc));
        m=mean(trc(ilive));
        trc=(trc-m)*bnds*fact(1) + x(k);
        if(~vaflag)
                line(trc(ilive),t(ilive),'color',kolor);
        elseif( vaflag==1)
                %replace NaN's with mean values
                if(length(ilive)< length(trc))
                        idead=find(isnan(trc));
                        trc(idead)=x(k)*ones(size(idead));
                end
                wtva(trc,t,kolor,x(k),1,1);                
        elseif(vaflag==2)
                wtva(trc(ilive),t(ilive),kolor,x(k),-1,1);
        end
end
axis('ij')
set(gca,'GridLineStyle','-')
set(gca,'YGrid','on')
set(gca,'XAxisLocation', 'top')
set(gca,'TickDir','out') 
set(gca,'FontName','times','FontSize',10)
%set(gca,'XTick',[.5 1. 1.5 2 2.4]) 
%set(gca,'XTickLabels',['0.5';'1.0';'1.5';'2.0';'2.4'])
%set(gca,'YTick',[2. 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9]) 
%set(gca,'YTickLabels',['2.0';'   ';'   ';'   ';'   ';'2.5';'   ';'   ';'   ';'2.9'])
 
