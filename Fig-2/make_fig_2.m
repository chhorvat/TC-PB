clear

addpath('../Plot-Tools'); 
horvat_colors;
dataloc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/Vanuatu-TC-Blooms/';
load([dataloc '/IBTracs/stormmap_rad'],'tracker','LAT','LON','transr'); 

IB_tracker = tracker; 
IB_transr = transr; 

load([dataloc '/Pete-Cyclones/Pete-data.mat']); 
load([dataloc 'stormmap_rad'],'transfactor','hoverfactor','transr'); 
load([dataloc 'stormmap_rad'],'tracker'); 

hoverfactor(hoverfactor == 0) = nan;
hoverfactor(isinf(hoverfactor)) = nan; 

exp_start = 10; 
exp_x = 1:.5:max(hoverfactor); 

dat_fit = hoverfactor(hoverfactor > exp_start); 
dat_fit = dat_fit - min(dat_fit); 

muval = expfit(dat_fit); 

exp_mu = 1/muval; 
exp_vals = exp(-exp_mu*(exp_x - min(exp_x))); 

fprintf('Exponential fit has coefficient %2.3f \n',exp_mu); 

%%
close all
figure(1);

use_Pete = [1:11 13:15]; % Ignore Katrina


Ax{1} = subplot('position',[0.1 0.55 0.35 0.4]);

lat_lim_ALL = [-30 -10];
lon_lim_ALL = [150.4 180];

worldmap(lat_lim_ALL,lon_lim_ALL);
pcolorm(LAT,LON,IB_tracker); 
hold on
add_coastlines; 
scatterm(Pete_Lat,Pete_Lon,50*Pete_chl,'filled','markeredgecolor',clabs(3,:),'linewidth',2); 
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','w');
setm(gca,'PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','on','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0,'MLabelParallel','south');
title('Observed Storm \#','interpreter','latex');
cmap = cmocean('thermal');
% cmap(1:floor(0*256/max(get(gca,'clim'))),:) = 0;
load coastlines
plotm(coastlat, coastlon,'k','linewidth',2)
hold on
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[217,217,217]/256)


colormap(gca,cmap);
% colorbar
tightmap

Ax{2} = subplot('position',[0.55 0.55 0.375 0.4]);



worldmap(lat_lim_ALL,lon_lim_ALL);
pcolorm(LAT,LON,tracker); 
hold on
add_coastlines; 
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','w');
setm(gca,'PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','on','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0,'MLabelParallel','south');
tightmap

title('STORM Storm \#','interpreter','latex');
cmap = cmocean('thermal');
load coastlines
plotm(coastlat, coastlon,'k','linewidth',2)
hold on
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[217,217,217]/256)


colormap(gca,cmap);
% colorbar

Ax{3} = subplot('position',[0.1 0.05 0.8 0.4]); 


HPAR = Pete_hover(use_Pete);
CPAR = (Pete_chl(use_Pete)); % + VIIRS.MODIS_Aquachla05); 
CPAR = CPAR;
HBINS = [0:1:55 100];
mdl = fitlm(HPAR,CPAR); 
disp(['R-value is ' num2str(mdl.Rsquared.Adjusted^(1/2))]); 

fitpar = polyfit(HPAR,CPAR,1); 

[b,bint,r,rint,stats] = regress(CPAR,[ones(length(HPAR),1) HPAR]);


yyaxis right
set(gca,'ycolor','k'); 


% [a,b] = histcounts(hoverfactor,[0.001 1:25 1000]); 
[a,b] = histcounts(hoverfactor,HBINS); 
binm = b(1:end-1) + .5*diff(b);
binm(end) = binm(end-1) + .5*(binm(end-1) - binm(end-2)); 
a = log10(a) + 1; 
a(isinf(a)) = 0; 

bar(binm,a,2,'facecolor',clabs(3,:),'facealpha',0.25); 
set(gca,'ytick',[0 1 2 3 4 5 6],'yticklabel',{'0','1','10','100','1000','10000'}); 
ylim([0 5.5]);
set(gca,'ycolor',clabs(3,:)); 
ylabel('Storm Count','interpreter','latex','Color','k')
hold on

p = find(binm>min(2*exp_x),1); 
str = 10^(a(p)-1);

plot(exp_x,1+log10(str*exp_vals/max(exp_vals)),'-k'); 



yyaxis left
scatter(HPAR,CPAR,100,'filled','markeredgecolor',clabs(3,:)); 
hold on
plot(binm,polyval(fitpar,binm),'--k'); 
ylim([0 2.5]); 
xlim([1 max(binm)]); 
xlabel('Hover Parameter (hr$^2$/km)','interpreter','latex'); 
set(gca,'ycolor','k'); 
ylabel('Chl-a (mg/m$^3$','interpreter','latex','Color','k')
grid on; box on; 




pos = [6.5 4]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

tightfig;

for i = 1:length(Ax)
    
    if i~=3
        pos = get(Ax{i},'position');
        colorbar(Ax{i},'position',[pos(1)+pos(3) pos(2) .025 pos(4)]);
        set(Ax{i},'fontsize',9);
    end
    
end

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

%%
for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.02 posy(2)+posy(4)+.02 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');
    
end

saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/AGU Geophysical Research Letters_Extreme_TC-PB/Figures/Fig-2.pdf');
