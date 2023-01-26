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
load([dataloc 'MODIS_Vanu'],'chl_a','latvec','lonvec','mo','yr','dy','date_name','area_Vanu');

storm_index = find((yr==stormyr)&(mo == stormmo)&(dy == stormdy));
stormdates = storm_index:(storm_index + storm_duration-1);


%% Individual smoothed storm
chl_infilled = fillmissing(chl_a,'movmean',[storm_duration 0],3);
chl_SEAWIFS_infilled =  fillmissing(SEAWIFS_data.chl_a,'movmean',[storm_duration 0],3);
chl_VIIRS_infilled =  fillmissing(VIIRS_data.chl_a,'movmean',[storm_duration 0],3);
chl_HIMA_infilled =  flipud(fillmissing(HIMA_data.chl_a,'movmean',[storm_duration 0],3));
%%

fprintf('Begin date %s \n',date_name(stormdates(1),:));
fprintf('End date %s \n',date_name(stormdates(end),:));

disp(['Storm plotted on ' date_name(storm_index,:)]);

lat_lim_ALL = [-30 -10];
lon_lim_ALL = [150.4 180];

lat_lim_Oma = [-22 -12];
lon_lim_Oma = [158 172];

lat_lim_WVR = [-17.1 -14.1]; 
lon_lim_WVR = [164 165.4]; 

lat_lim_Win = [-18.3 -17.3];
lon_lim_Win = [171.6 172.6];

%% Compare to climatology for previous periods

% Ignore the storm in question
goodyrs = (yr~=stormyr + min(yr)-1);
goodmos = unique(mo(stormdates));

% Find the dates in each period that have the right days and mons
comp_dates = find( (mo == mo(stormdates(1)) & (dy == dy(stormdates(1)))));

yra = 1 + yr - min(yr);

CHL_stormper_clim = nan(size(chl_a,1),size(chl_a,2),max(yra));

for i = 1:length(comp_dates)
    CHL_stormper_clim(:,:,yra(comp_dates(i))) = nanmean(chl_a(:,:,comp_dates(i):comp_dates(i)+length(stormdates)-1),3);
end

CHL_stormper_clim(:,:,yra(stormdates(1))) = nan;
CHL_stormper_clim = nanmean(CHL_stormper_clim,3);

CHL_stormper = nanmean(chl_a(:,:,stormdates),3);
CHL_stormper_rat = CHL_stormper./CHL_stormper_clim;

%%
for i = 1:max(yra)
    
    CHL_ann(:,:,i) = nanmean(chl_a(:,:,yra==i),3);
    
end

yrlose = 1 + stormyr - min(yr);
CHL_ann(:,:,yrlose) = nan;

CHL_clim = nanmean(CHL_ann,3);

%% Get the maximum sum over each window

mov_amnt = nan(size(chl_infilled));



%%
for i = 1:size(chl_a,1)
    if mod(i,20) == 0
        disp(i)
    end
    %     for j = 1:size(chl_a,2)
    
    fld = chl_infilled(i,:,:);
    
    mov_amnt(i,:,:) = movmean(fld,storm_duration,3,'omitnan');
    
    %    end
end
%
% %%
%%
[a,maxloc] = nanmax(mov_amnt,[],3);

maxfrac = nanmax(mov_amnt,[],3)./CHL_clim;

% maxfrac = nanmax(mov_amnt,[],3)./nansum(chl_infilled,3);
maxnum = sum(~isnan(mov_amnt),3);

%%

blooming = chl_infilled(:,:,storm_index) > 0.5; 
blooming(latvec < min(lat_lim_Oma),:) = 0; 
blooming(:,lonvec < min(lon_lim_Oma)) = 0; 
blooming(latvec > max(lat_lim_Oma),:) = 0; 
blooming(:,lonvec > max(lon_lim_Oma)) = 0; 

area_blooming = nansum(nansum(blooming.*area_Vanu));
fprintf('Blooming area is %2.0f km2 \n',area_blooming);

% Find the locations where there is a lot of Chl on that day
% vanu_peak_locs = blooming;
% vanu_peak_locs = CHL_stormper_rat > 6;
% vanu_peak_locs((latvec < -18) | (latvec > -12),:) = 0;
% vanu_peak_locs(:,(lonvec < 161) | (lonvec > 168)) = 0;
% vanu_peak_locs = double(vanu_peak_locs);
% vanu_peak_locs = imdilate(imerode(imfill(1*vanu_peak_locs),[1 1; 1 1; 1 1]),[1 1; 1 1 ; 1 1]);
% vanu_peak_locs(vanu_peak_locs == 0) = nan;

WVR_box = ones(size(blooming)); 
WVR_box(latvec < min(lat_lim_WVR),:) = 0; 
WVR_box(:,lonvec < min(lon_lim_WVR)) = 0; 
WVR_box(latvec > max(lat_lim_WVR),:) = 0; 
WVR_box(:,lonvec > max(lon_lim_WVR)) = 0; 

vanu_peak_locs = double(blooming); 
vanu_peak_locs(~WVR_box) = 0; 
vanu_peak_locs(vanu_peak_locs == 0) = nan; 

% vanu_peak_locs = CHL_stormper_rat > 10;

% Reshape the infilled data to be nX by nt
% int_chl_vanu = reshape(chl_Oma_infilled,numel(vanu_peak_locs),[]);
% Then take the mean over all values that meat that criteria

disp('Calculating average CHL in the Vanu Box for different products'); 


int_chl_vanu = squeeze(nansum(nansum(area_Vanu.*WVR_box.*chl_infilled)) ...
    ./nansum(nansum(area_Vanu.*WVR_box)));

int_chl_vanu_SEAWIFS = squeeze(nansum(nansum(area_Vanu.*WVR_box.*chl_SEAWIFS_infilled)) ...
    ./nansum(nansum(area_Vanu.*WVR_box)));

int_chl_vanu_VIIRS = squeeze(nansum(nansum(area_Vanu.*WVR_box.*chl_VIIRS_infilled)) ...
    ./nansum(nansum(area_Vanu.*WVR_box)));

int_chl_vanu_HIMA = squeeze(nansum(nansum(area_Vanu.*WVR_box.*chl_HIMA_infilled)) ...
    ./nansum(nansum(area_Vanu.*WVR_box)));


% int_chl_vanu = squeeze(nanmean(nanmean(bsxfun(@times,area_Vanu.*WVR_box,chl_infilled),1),2));
% int_chl_vanu_SEAWIFS = squeeze(nanmean(nanmean(bsxfun(@times,WVR_box,chl_SEAWIFS_infilled),1),2));
% int_chl_vanu_VIIRS = squeeze(nanmean(nanmean(bsxfun(@times,WVR_box,chl_VIIRS_infilled),1),2));
% int_chl_vanu_HIMA = squeeze(nanmean(nanmean(bsxfun(@times,WVR_box,chl_HIMA_infilled),1),2))


%% 
fprintf('Area of WVR is %2.0f km2 \n',nansum(nansum(area_Vanu.*WVR_box))); 
area_usable = nansum(nansum(vanu_peak_locs.*area_Vanu));
fprintf('Area that bloomed on MODIS is %2.0f km2 \n',area_usable);

mean_VIIRS = nanmean(int_chl_vanu_VIIRS(1:2500)); 
std_VIIRS = nanstd(int_chl_vanu_VIIRS(1:2500)); 
peak_VIIRS = nanmax(int_chl_vanu_VIIRS); 

fprintf('VIIRS: Peak at %2.2f, Mean is %2.2f p/m %2.2f for z-score of %2.0f \n',peak_VIIRS,mean_VIIRS,std_VIIRS,(peak_VIIRS-mean_VIIRS)/std_VIIRS); 

mean_MODIS = nanmean(int_chl_vanu(1:6000)); 
std_MODIS = nanstd(int_chl_vanu(1:6000)); 
peak_MODIS = nanmax(int_chl_vanu); 

fprintf('MODIS: Peak at %2.2f, Mean is %2.2f p/m %2.2f for z-score of %2.0f \n',peak_MODIS,mean_MODIS,std_MODIS,(peak_MODIS-mean_MODIS)/std_MODIS); 

mean_HIMA = nanmean(int_chl_vanu_HIMA(1:1250)); 
std_HIMA = nanstd(int_chl_vanu_HIMA(1:1250)); 
peak_HIMA = nanmax(int_chl_vanu_HIMA); 

fprintf('HIMA: Peak at %2.2f, Mean is %2.2f p/m %2.2f for z-score of %2.0f \n',peak_HIMA,mean_HIMA,std_HIMA,(peak_HIMA-mean_HIMA)/std_HIMA); 

mean_SEAWIFS = nanmean(int_chl_vanu_SEAWIFS); 
std_SEAWIFS = nanstd(int_chl_vanu_SEAWIFS); 

fprintf('SEAWIFS: Mean is %2.2f p/m %2.2f \n',mean_SEAWIFS,std_SEAWIFS); 

fprintf('Mean Chl-a in the WVR is %2.2f \n',nansum(nansum(area_Vanu.*WVR_box.*chl_infilled(:,:,storm_index)))./nansum(nansum(area_Vanu.*WVR_box))); 
fprintf('Mean Clim Rate in the WVR is %2.2f \n',nansum(nansum(area_Vanu.*WVR_box.*CHL_stormper_rat))./nansum(nansum(area_Vanu.*WVR_box))); 

fprintf('Mean MEP in the S.P. is %2.2f \n',nanmean(nanmean(maxfrac(100:end,100:end)))); 
fprintf('Mean MEP in the WVR is %2.2f \n',nansum(nansum(area_Vanu.*WVR_box.*maxfrac))./nansum(nansum(area_Vanu.*WVR_box))); 

%%
addpath('../Plot-Tools');

close(figure(1))
figure(1)

Ax{1} = subplot('position',[0.06 0.55 0.25 0.4]);

cla
worldmap(lat_lim_ALL,lon_lim_ALL);

latter = linspace(min(lat_lim_Oma),max(lat_lim_Oma),100);
lonner = linspace(min(lon_lim_Oma),max(lon_lim_Oma),100);

latter = [latter 0*latter + max(latter) fliplr(latter) 0*latter + min(latter)];
lonner = [0*lonner + max(lonner) fliplr(lonner) 0*lonner + min(lonner) lonner];



latter_WVR = linspace(min(lat_lim_WVR),max(lat_lim_WVR),100);
lonner_WVR = linspace(min(lon_lim_WVR),max(lon_lim_WVR),100);

latter_WVR = [latter_WVR 0*latter_WVR + max(latter_WVR) fliplr(latter_WVR) 0*latter_WVR + min(latter_WVR)];
lonner_WVR = [0*lonner_WVR + max(lonner_WVR) fliplr(lonner_WVR) 0*lonner_WVR + min(lonner_WVR) lonner_WVR];


% plotter = ncread([dataloc '/MUR-SST/20190223090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'],'analysed_sst',[32000 5000 1],[Inf 4000 1]) ...
%     - ncread([dataloc '/MUR-SST/20190217090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'],'analysed_sst',[32000 5000 1],[Inf 4000 1]);
% 
% plot_SST_lat = double(ncread([dataloc '/MUR-SST/20190223090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'],'lat',5000,4000));% ,150,200)); 
% plot_SST_lon = double(ncread([dataloc '/MUR-SST/20190223090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'],'lon',32000,Inf));%,600,200));


%
plot_SST = (ncread([dataloc '/NOAA-SST/sst.day.mean.2019.nc'],'sst',[1 1 46],[Inf Inf 6]));%,[600 150 45],[200 200 14])); 
plot_SST_lat = double(ncread([dataloc '/NOAA-SST/sst.day.mean.2019.nc'],'lat'));% ,150,200)); 
plot_SST_lon = double(ncread([dataloc '/NOAA-SST/sst.day.mean.2019.nc'],'lon'));%,600,200));
plotter = plot_SST(:,:,end)' - plot_SST(:,:,1)';

pcolorm(plot_SST_lat,plot_SST_lon,plotter);
hold on
plotm(latter,lonner,'r','linewidth',2);
% plotm(latter_WVR,lonner_WVR,'w','linewidth',2);
set(gca,'clim',[-4.5 1]); %max(get(gca,'clim'))]);
setm(gca,'PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','on','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0);
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','k');

%
title(['$\Delta$ SST: 2/17/19-2/23/19 '],'interpreter','latex');
% cmap(1:floor(0.5*256/max(get(gca,'clim'))),:) = 0;
cmap = cmocean('curl','pivot',0);
colormap(gca,cmap);

OCscat = OMA_cat;
OCscat(OCscat <= 0) = 0.05;

scatterm(OMA_lat(OMA_timer),OMA_lon(OMA_timer),100,'markeredgecolor','r','linewidth',2);
hold on
plotm(OMA_lat(OMA_cat > -2),OMA_lon(OMA_cat > -2),'-r','linewidth',1);

OMA_TClat = OMA_lat;
OMA_TClat(OMA_cat <= 0) = nan;
OMA_TClon = OMA_lon;
OMA_TClon(OMA_cat <= 0) = nan;

plotm(OMA_TClat,OMA_TClon,'-b','linewidth',2);

load coastlines
plotm(coastlat, coastlon,'k','linewidth',4)
hold on
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[200,200,200]/256)

%
Ax{4} = subplot('position',[0.06 0.06 0.25 0.4]);
worldmap(lat_lim_Oma,lon_lim_Oma);

plotter = nanmax(chl_a(:,:,storm_index:storm_index),[],3);

pcolorm(latvec,lonvec,plotter);
hold on
% plotm(latter_WVR,lonner_WVR,'w','linewidth',2);

set(gca,'clim',[0 15]); %max(get(gca,'clim'))]);
setm(gca,'PlabelLocation',[-25 -15],'PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','on','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0);
setm(gca,'PlabelMeridian','West'); 
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','k');


title(['CHL-a ' date_name(storm_index(1),:)],'interpreter','latex');
cmap = cmocean('algae');
cmap(1:floor(0.5*256/max(get(gca,'clim'))),:) = 0;
colormap(gca,cmap);

hold on
% plotm(OMA_lat(OMA_cat > -2),OMA_lon(OMA_cat > -2),'-r','linewidth',1);

OMA_TClat = OMA_lat;
OMA_TClat(OMA_cat <= 0) = nan;
OMA_TClon = OMA_lon;
OMA_TClon(OMA_cat <= 0) = nan;
% plotm(OMA_TClat,OMA_TClon,'-b','linewidth',2);
load coastlines
plotm(coastlat, coastlon,'k','linewidth',4)
hold on
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[200,200,200]/256)

Ax{2} = subplot('position',[0.365 0.06 0.25 0.4]);
worldmap(lat_lim_Oma,lon_lim_Oma);

add_coastlines;

pcolorm(latvec,lonvec,CHL_stormper_rat);
hold on
plotm(latter_WVR,lonner_WVR,'color',[.99 .99 .99],'linewidth',1);

VPL = vanu_peak_locs; 
VPL(isnan(VPL)) = 0; 
% contourm(latvec,lonvec,VPL,[1 1],'w','linewidth',1);
% contourm(latvec,lonvec,chl_infilled(:,:,storm_index),[1 1],'w','linewidth',1);

set(gca,'clim',[0 30]);
setm(gca,'ParallelLabel','off','PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','on','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0);
setm(gca,'grid','on','GLineWidth',1,'GLineStyle','--','GColor','w');

title(['Mult. of Clim.'],'interpreter','latex');
% cmap(1:floor(1*256/max(get(gca,'clim'))),:) = 1;
load coastlines
plotm(coastlat, coastlon,'k','linewidth',4)
hold on
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[200,200,200]/256)
cmap = cmocean('thermal');
colormap(gca,cmap);


%
Ax{3} = subplot('position',[0.675 0.06 0.25 0.4]);
worldmap(lat_lim_Oma,lon_lim_Oma);
add_coastlines;

pcolorm(latvec,lonvec,100*maxfrac);
hold on
% contourm(latvec,lonvec,chl_a(:,:,storm_index),[5 5],'w','linewidth',1);
% contourm(latvec,lonvec,CHL_stormper_rat,[10 10],'w','linewidth',1);
plotm(latter_WVR,lonner_WVR,'color',[.99 .99 .99],'linewidth',1);

set(gca,'clim',[0 5000]);% max(get(gca,'clim'))]);
setm(gca,'ParallelLabel','off','PlineLocation',[-25 -15],'PLabelRound',0);
setm(gca,'MeridianLabel','on','MlabelLocation',[160 170],'MLineLocation',[160 170],'MLabelRound',0);
setm(gca,'grid','on','GLineWidth',1,'GLineStyle','--','GColor','w');

title('MEP (\%)','interpreter','latex');
cmap = cmocean('thermal');
% cmap(1:floor(10*256/max(get(gca,'clim'))),:) = 0;
colormap(gca,cmap);
% contourm(latvec,lonvec,VPL,[1 1],'w','linewidth',1);
load coastlines
plotm(coastlat, coastlon,'k','linewidth',4)
hold on
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[200,200,200]/256)
%

Ax{5} = subplot('position',[.375 .575 .55 .36]);

%
% int_chl_vanu = squeeze(nanmean(int_chl_vanu(vanu_peak_locs(:),:),1));
horvat_colors; 

% int_chl_vanu_yr = accumarray(yra',int_chl_vanu,[max(yra) 1],@nansum);
times_MODIS = datetime(date_name);
times_SEAWIFS = datetime(SEAWIFS_data.date_name);
times_VIIRS = datetime(VIIRS_data.date_name(:,1:10));
times_HIMA = datetime(HIMA_data.date_name,'InputFormat','yyyyMMdd');
plot(times_MODIS,int_chl_vanu,'k','linewidth',1);
hold on
plot(times_SEAWIFS,int_chl_vanu_SEAWIFS + 2,'color',clabs(2,:),'linewidth',1);
plot(times_VIIRS,int_chl_vanu_VIIRS + 4,'color',clabs(4,:),'linewidth',1);
plot(times_HIMA,int_chl_vanu_HIMA + 6,'color',clabs(5,:),'linewidth',1);
grid on; box on;
yline(0.5,'--k','label','Blooming','labelhorizontalalignment','left'); 
yline(2,'-k'); 
yline(2.5,'--k')
yline(4.5,'--k')
yline(6.5,'--k')
yline(2,'-k'); 
yline(4,'-k'); 
yline(6,'-k'); 
title('WVR Satellite CHL-a (mg/m$^3$)','interpreter','latex');

legend('MODIS','SEAWIFS + 2','VIIRS + 4','HIMAWARI-8 + 6','location','best','fontsize',6)
% ylabel('CHL$_a$','interpreter','latex');
hold off
pos = [6.5 3.75];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
%

for i = 1:length(Ax)
    
    if i~=5
        pos = get(Ax{i},'position');
        % colorbar(Ax{i},'position',[pos(1)+pos(3)+.01 pos(2)+.025 .015 .35],'orientation','horizontal');
        colorbar(Ax{i},'position',[pos(1)+ .01 pos(2)-.01 pos(3)-.02 .015],'orientation','horizontal');
        set(Ax{i},'fontsize',9);
    end
    
end

%%

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');

    if i == 1

    annotation('textbox',[posy(1)-.04 posy(2)+posy(4)-.01 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

    else

     annotation('textbox',[posy(1) posy(2)+posy(4)-.01 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');
    end   
    
end


saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/AGU Geophysical Research Letters_Extreme_TC-PB/Figures/Fig-1_1day.pdf');
saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/AGU Geophysical Research Letters_Extreme_TC-PB/Figures/Fig-1_1day.jpeg');


