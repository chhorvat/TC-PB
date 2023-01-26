% Figure 1
clear

dataloc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/Vanuatu-TC-Blooms/';

%% Individual storm
% storm_index = 6068;
% stormdates = 6068:6081;
stormyr = 2019;
stormmo = 2;
stormdy = 21;
storm_duration = 1;


%% Load in Data - first Oma track

oma_sid = 1264;
OMA_lat = ncread([dataloc '/IBTracs/IBTrACS.SP.v04r00.nc'],'lat',[1 oma_sid],[Inf 1]);
OMA_lon = ncread([dataloc '/IBTracs/IBTrACS.SP.v04r00.nc'],'lon',[1 oma_sid],[Inf 1]);
OMA_cat = ncread([dataloc '/IBTracs/IBTrACS.SP.v04r00.nc'],'usa_sshs',[1 oma_sid],[Inf 1]);
Oma_time = ncread([dataloc '/IBTracs/IBTrACS.SP.v04r00.nc'],'iso_time',[1 1 oma_sid],[Inf Inf 1]);
OMA_timer = 85;
%% Second HIMAWARA data
SST_data = matfile([dataloc 'SST_Vanu']); 
HIMA_data = matfile([dataloc 'HIMAWARI_Vanu']);
SEAWIFS_data = matfile([dataloc 'SEAWIFS_Vanu']);
VIIRS_data = matfile([dataloc 'VIIRS_Vanu']);
MODIS_data = matfile([dataloc 'MODIS_Vanu']);

yr_HIMA = HIMA_data.yr;
mo_HIMA = HIMA_data.mo;
dy_HIMA = HIMA_data.dy;
storm_index_HIMA = find((yr_HIMA==stormyr)&(mo_HIMA == stormmo)&(dy_HIMA == stormdy));
stormdates_HIMA = storm_index_HIMA:(storm_index_HIMA + storm_duration-1);

HIMA_CHL = HIMA_data.chl_a(:,:,stormdates_HIMA);
HIMA_lat = HIMA_data.latvec;
HIMA_lon = HIMA_data.lonvec;

yr_VIIRS = VIIRS_data.yr;
mo_VIIRS = VIIRS_data.mo;
dy_VIIRS = VIIRS_data.dy;
storm_index_VIIRS = find((yr_VIIRS==stormyr)&(mo_VIIRS == stormmo)&(dy_VIIRS == stormdy));
stormdates_VIIRS = storm_index_VIIRS:(storm_index_VIIRS + storm_duration-1);

VIIRS_CHL = VIIRS_data.chl_a(:,:,stormdates_VIIRS);
VIIRS_lat = VIIRS_data.latvec;
VIIRS_lon = VIIRS_data.lonvec;

%%
% load([dataloc 'CHL_Vanu'],'chl_a','latvec','lonvec','mo','yr','dy','date_name');
load([dataloc 'MODIS_Vanu'],'latvec','lonvec','mo','yr','dy','date_name','area_Vanu');

storm_index_MODIS = find((yr==stormyr)&(mo == stormmo)&(dy == stormdy));
stormdates = storm_index_MODIS:(storm_index_MODIS + storm_duration-1);


%% Individual smoothed storm
chl_MODIS = (MODIS_data.chl_a(:,:,storm_index_MODIS));
chl_VIIRS =  (VIIRS_data.chl_a(:,:,storm_index_VIIRS));
chl_HIMA =  flipud(HIMA_data.chl_a(:,:,storm_index_HIMA));
%%

lat_lim_plot = [-30 -12];
lon_lim_plot = [150 172];


%% Compare to climatology for previous periods


%%

% vanu_peak_locs = blooming;
% vanu_peak_locs = CHL_stormper_rat > 6;
% vanu_peak_locs((latvec < -18) | (latvec > -12),:) = 0;
% vanu_peak_locs(:,(lonvec < 161) | (lonvec > 168)) = 0;
% vanu_peak_locs = double(vanu_peak_locs);
% vanu_peak_locs = imdilate(imerode(imfill(1*vanu_peak_locs),[1 1; 1 1; 1 1]),[1 1; 1 1 ; 1 1]);
% vanu_peak_locs(vanu_peak_locs == 0) = nan;

% WVR_box = ones(size(blooming)); 
% WVR_box(latvec < min(lat_lim_WVR),:) = 0; 
% WVR_box(:,lonvec < min(lon_lim_WVR)) = 0; 
% WVR_box(latvec > max(lat_lim_WVR),:) = 0; 
% WVR_box(:,lonvec > max(lon_lim_WVR)) = 0; 
% 
% vanu_peak_locs = double(blooming); 
% vanu_peak_locs(~WVR_box) = 0; 
% vanu_peak_locs(vanu_peak_locs == 0) = nan; 

% vanu_peak_locs = CHL_stormper_rat > 10;

% Reshape the infilled data to be nX by nt
% int_chl_vanu = reshape(chl_Oma_infilled,numel(vanu_peak_locs),[]);
% Then take the mean over all values that meat that criteria
% 
% disp('Calculating average CHL in the Vanu Box for different products'); 



% int_chl_vanu = squeeze(nanmean(nanmean(bsxfun(@times,area_Vanu.*WVR_box,chl_infilled),1),2));
% int_chl_vanu_SEAWIFS = squeeze(nanmean(nanmean(bsxfun(@times,WVR_box,chl_SEAWIFS_infilled),1),2));
% int_chl_vanu_VIIRS = squeeze(nanmean(nanmean(bsxfun(@times,WVR_box,chl_VIIRS_infilled),1),2));
% int_chl_vanu_HIMA = squeeze(nanmean(nanmean(bsxfun(@times,WVR_box,chl_HIMA_infilled),1),2))


%% 

%%
close
%
Ax{1} = subplot(131);

worldmap(lat_lim_plot,lon_lim_plot);
plotter = chl_MODIS;
pcolorm(latvec,lonvec,plotter);
hold on
% plotm(latter_WVR,lonner_WVR,'w','linewidth',2);

set(gca,'clim',[0 5]); %max(get(gca,'clim'))]);
setm(gca,'PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','off','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0);
setm(gca,'PlabelMeridian','West'); 
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','k');


title(['MODIS CHL-a'],'interpreter','latex');
cmap = cmocean('algae');
cmap(1:floor(0.5*256/max(get(gca,'clim'))),:) = 0;
colormap(gca,cmap);
add_coastlines; 


Ax{2} = subplot(132);
% hold on
% plotm(OMA_lat(OMA_cat > -2),OMA_lon(OMA_cat > -2),'-r','linewidth',1);
% 
% OMA_TClat = OMA_lat;
% OMA_TClat(OMA_cat <= 0) = nan;
% OMA_TClon = OMA_lon;
% OMA_TClon(OMA_cat <= 0) = nan;
% % plotm(OMA_TClat,OMA_TClon,'-b','linewidth',2);
% load coastlines
% plotm(coastlat, coastlon,'k','linewidth',4)
% hold on
% geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[200,200,200]/256)
worldmap(lat_lim_plot,lon_lim_plot);
plotter = chl_VIIRS; 
pcolorm(latvec,lonvec,plotter);
hold on
% plotm(latter_WVR,lonner_WVR,'w','linewidth',2);

set(gca,'clim',[0 5]); %max(get(gca,'clim'))]);
setm(gca,'ParallelLabel','off','PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','off','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0);
setm(gca,'PlabelMeridian','West'); 
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','k');
add_coastlines; 


title(['VIIRS CHL-a'],'interpreter','latex');
cmap = cmocean('algae');
cmap(1:floor(0.5*256/max(get(gca,'clim'))),:) = 0;
colormap(gca,cmap);

Ax{3} = subplot(133);

worldmap(lat_lim_plot,lon_lim_plot);
plotter = chl_HIMA; 
pcolorm(latvec,lonvec,plotter);
hold on
% plotm(latter_WVR,lonner_WVR,'w','linewidth',2);

set(gca,'clim',[0 5]); %max(get(gca,'clim'))]);
setm(gca,'ParallelLabel','off','PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','off','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0);
setm(gca,'PlabelMeridian','West'); 
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','k');


title(['HIMAWARI-8 CHL-a'],'interpreter','latex');
cmap = cmocean('algae');
cmap(1:floor(0.5*256/max(get(gca,'clim'))),:) = 0;
colormap(gca,cmap);
add_coastlines; 

colorbar('position',[.925 .25 .025 .5])

hold off
pos = [6.5 2];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

%%

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');

%     if i == 1
% 
    annotation('textbox',[posy(1)-.04 posy(2)+posy(4)-.01 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');
% 
%     else
% 
%      annotation('textbox',[posy(1) posy(2)+posy(4)-.01 .025 .025], ...
%         'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
%         'FontSize',8,'Tag','legtag');
%     end   
%     
end


saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/AGU Geophysical Research Letters_Extreme_TC-PB/Figures/Fig-S4.pdf');
saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/AGU Geophysical Research Letters_Extreme_TC-PB/Figures/Fig-S4.jpeg');


