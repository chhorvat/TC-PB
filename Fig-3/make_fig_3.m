clear

dataloc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/Vanuatu-TC-Blooms/';

load([dataloc 'MODIS_Vanu'],'area_Vanu');
load([dataloc 'stormmap_rad'],'hoverfactor','lathover','lonhover','LAT','LON','hoverr','tracker','omatrack','omatrackpos');

% transit_time_hover = hoverr./hoverv; 
% transit_time_all = trackr./trackv; 
% transit_time = transr./transv; 

hoverfactor(hoverfactor == 0) = nan;
hoverfactor(isinf(hoverfactor)) = nan; 

load([dataloc '/Pete-Cyclones/Pete-data.mat']); 

addpath('../Plot-Tools'); 
% Get locations where the storms hover

sublats = LAT(1:end,1:end); 
sublons = LON(1:end,1:end); 
earthellipsoid = referenceSphere('earth','km');
lldist = @(x,y) distance(x,y,earthellipsoid);


lat_lim_WVR = [-17.1 -14.1]; 
lon_lim_WVR = [164 165.4]; 

WVR_box = ones(size(LAT)); 
WVR_box(LAT < min(lat_lim_WVR)) = 0; 
WVR_box(LON < min(lon_lim_WVR)) = 0; 
WVR_box(LAT > max(lat_lim_WVR)) = 0; 
WVR_box(LON > max(lon_lim_WVR)) = 0; 

WVR_nanbox = WVR_box;
WVR_nanbox(WVR_nanbox == 0) = nan; 



lats = lathover(hoverfactor > 2); 
lons = lonhover(hoverfactor > 2); 
rads = hoverr(hoverfactor > 2); 

lats_Oma = lathover(hoverfactor > 10); 
lons_Oma = lonhover(hoverfactor > 10); 
rads_Oma = hoverr(hoverfactor > 10); 

M = createns([sublats(:) sublons(:)],'Distance',lldist);

locs_slow = 0*sublats; 
locs_oma = 0*sublons; 

for i = 1:length(rads)

    q = rangesearch(M,[lats(i) lons(i)],rads(i),'Distance',lldist); 
    locs_slow(q{1}) = locs_slow(q{1}) + 1;

end

for i = 1:length(rads_Oma)

    q = rangesearch(M,[lats_Oma(i) lons_Oma(i)],rads_Oma(i),'Distance',lldist); 
    locs_oma(q{1}) = locs_oma(q{1}) + 1;

end

N_TC = nanmean(nanmean(WVR_nanbox.*tracker)); 
N_PB = nanmean(nanmean(WVR_nanbox.*locs_slow)); 
N_OM = nanmean(nanmean(WVR_nanbox.*locs_oma)); 


fprintf('There were %d Storms. %d of them produce TC-PBs. %d were Oma-type \n',numel(hoverfactor),sum(hoverfactor>2),sum(hoverfactor>10)); 

fprintf('The WVR was reached by %d TCs, of which %d produced PBs and %d were Oma-type \n',round(N_TC),round(N_PB),round(N_OM)); 


%% Now show the cyclone frequency in this region
close all

horvat_colors; 
figure(1);
Ax{1} = subplot('position',[0.05 .55 .5 .4]);

cbars = [0 .5 1 2 3 4 5 6 7 8 9 10 20];

latter = linspace(min(lat_lim_WVR),max(lat_lim_WVR),100);
lonner = linspace(min(lon_lim_WVR),max(lon_lim_WVR),100);

latter = [latter 0*latter + max(latter) fliplr(latter) 0*latter + min(latter)];
lonner = [0*lonner + max(lonner) fliplr(lonner) 0*lonner + min(lonner) lonner];




lat_lim_Oma = [-22 -12];
lon_lim_Oma = [158 168]; 

lat_lim_ALL = [-30 -10]; 
lon_lim_ALL = [150.8 180]; 

worldmap(lat_lim_ALL,lon_lim_ALL)
% worldmap(lat_lim_Oma,lon_lim_Oma)

omatrack = reshape(omatrack,size(LAT)); 
omatrackpos = reshape(omatrackpos,size(LAT)); 

plotter = 10000./locs_slow; 
plotter(isinf(plotter)) = 20000; 
% 
% plotter(isinf(plotter)) = nan; 
% plotter(plotter >= 10000) = nan; 
pcolorm(sublats,sublons,locs_slow); 
shading flat
hold on
plotm(latter,lonner,'color',[.99 .99 .99],'linewidth',1);

conslon = [140:190]; 
conslat = -12.5 + 0*conslon; 
conslat2 = -20 + 0*conslon; 

% plotm(conslat,conslon,'color',clabs(3,:),'linewidth',2); 
% plotm(conslat2,conslon,'color',clabs(3,:),'linewidth',2); 

colorbar
colormap(gca,(cmocean('solar')));
title('Hovering TC','interpreter','latex');
setm(gca,'grid','on');
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','w');
setm(gca,'PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','off','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0,'MLabelParallel','south');
tightmap


% set(gca,'clim',[0 max(get(gca,'clim'))]); 
add_coastlines; 

Ax{2} = subplot('position',[.55 .575 .4 .35]);
plot(mean(sublats,1),mean(locs_slow,1),'k','linewidth',1); 
% xlabel('Latitude','interpreter','latex'); 
title('\# TC-PB','interpreter','latex'); 
set(gca,'xticklabel',''); 
grid on; box on; 
xline(Pete_Lat(1),'color',clabs(3,:),'linewidth',2,'Label','Oma','LabelHorizontalAlignment','Center','LabelVerticalAlignment','middle'); 

hold off

Ax{3} = subplot('position',[0.05 .1 .5 .4]);
worldmap(lat_lim_ALL,lon_lim_ALL)
hold on

pcolorm(sublats,sublons,locs_oma);%,[0:2:14])
shading flat
plotm(latter,lonner,'color',[.99 .99 .99],'linewidth',1);


% set(gca,'clim',[1 4]);
colorbar
% scatterm(,,transit_time(Oma_type),1:sum(Oma_type),'filled'); 
add_coastlines; 
colormap(gca,cmocean('solar')); 
title('Oma-type TC','interpreter','latex');
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','w');
setm(gca,'PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','on','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0,'MLabelParallel','south');
tightmap


Ax{4} = subplot('position',[.55 .15 .4 .35]);
hold on
grid on; box on; 
title('\# Oma-type TC-PB','interpreter','latex');
plot(mean(sublats,1),mean(locs_oma,1),'color','k','linewidth',1); 
xlabel('Latitude','interpreter','latex'); 
xline(Pete_Lat(1),'color',clabs(3,:),'linewidth',2,'Label','Oma','LabelHorizontalAlignment','Center','LabelVerticalAlignment','middle'); 



%
% subplot(223)
% [h,b] = histcounts(transit_time/24,[0:.5:12]); 
% mid = b(1:end-1) + diff(b); 
% 
% barc = cmocean('phase'); 
% fmult = ceil(256/length(h)); 
% bars = bar(mid,log10(h),'FaceColor','flat');
% bars.CData = barc(1:fmult:end,:);
% set(gca,'ytick',[0 1 2 3 4],'yticklabel',{'1','10','100','1000','10000'}); 
% xlim([0.5 10]);
% grid on
% box on
% drawnow
% 
% xline(100/24,'k','linewidth',2,'label','Oma'); 
% xlabel('Time (d)'); 
% ylabel('Number'); 
% title('Transit Times of Slow TCs','interpreter','latex'); 
% 
% subplot(224)
% [h,b] = histcounts(transit_time_hover/24,[0:.5:12]); 
% mid = b(1:end-1) + diff(b); 
% 
% barc = cmocean('phase'); 
% fmult = ceil(256/length(h)); 
% bars = bar(mid,log10(h),'FaceColor','flat');
% bars.CData = barc(1:fmult:end,:);
% set(gca,'ytick',[0 1 2 3 4],'yticklabel',{'1','10','100','1000','10000'}); 
% xlim([0.5 10]);
% grid on
% box on
% drawnow
% 
% xline(100/24,'k','linewidth',2,'label','Oma'); 
% xlabel('Time (d)'); 
% ylabel('Number'); 
% title('Transit Times of Hovering TCs','interpreter','latex'); 


   


pos = [6.5 3.5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1) posy(2)+posy(4)-.025 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');
    
end


saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/AGU Geophysical Research Letters_Extreme_TC-PB/Figures/Fig-3.pdf'); 

saveas(gcf,'Fig-3.pdf'); 