%% deconvolve noiseless and noisy traces with common parameters
makeQsynthetic

deconw_nonstat

names={'reflectivity','noiseless input','deconw',['deconw fmax=' int2str(fmax)]};
namesn={'reflectivity','noisey input','deconw',['deconw fmax=' int2str(fmaxn)]};

figure
fs=8;
x0=.02;y0=.1;wid=.4;sep=.08;ht=1-2*y0;
ya='n';
bc='w';
zt=10;
subplot('position',[x0,y0,wid,ht])
trplot(t,[r,sq,sd,sdb],'normalize',1,'order','d','names',names,'fontsize',fs,...
    'color',zeros(1,4),'tracespacing',1.25,'namesalign','left','nameshift',0,'yaxis',ya)
text(1,0,zt,str,'fontsize',fs,'backgroundcolor',bc)
text(1,-1.2,zt,strb,'fontsize',fs,'backgroundcolor',bc)
text(0.05,3.6,'A','fontsize',14,'fontweight','bold')
xlim([0 2.5])
subplot('position',[x0+wid+sep,y0,wid,ht])
trplot(t,[r,sqn,sdn,sdnb],'normalize',1,'order','d','names',namesn,'fontsize',fs,...
    'color',zeros(1,4),'tracespacing',1.25,'namesalign','left','nameshift',0,'yaxis',ya)
text(1,0,zt,strn,'fontsize',fs,'backgroundcolor',bc)
text(1,-1.2,zt,strnb,'fontsize',fs,'backgroundcolor',bc)
text(0.05,3.6,'B','fontsize',14,'fontweight','bold')
xlim([0 2.5])

prepfig

print -depsc decongraphics\deconwnonstat_time.eps

%spectral picture
figure
lw=[.5 1 1.2 1.5];
gl=[0,.25,.5,.7];
ls={':','-','-','-'};
subplot(1,2,1)
dbspec(t,[r sq sqn pad_trace(w,r)],'normoption',1,'graylevels',gl,'linewidths',lw,'linestyles',ls);
names={'reflectivity','trace noise-free',['trace s2n=' num2str(s2n)],...
    ['wavelet fdom=' num2str(fdom)]};
legend(names,'location','southwest')
title('Before deconw')
subplot(1,2,2)
dbspec(t,[r sd sdn],'normoption',1,'graylevels',gl(1:3),'linewidths',lw(1:3),'linestyles',ls(1:3));
names={'reflectivity','trace noise-free',['trace s2n=' num2str(s2n)]};
legend(names,'location','southwest')
title('After deconw')
prepfig

print -depsc decongraphics\deconwnonstat_freq.eps