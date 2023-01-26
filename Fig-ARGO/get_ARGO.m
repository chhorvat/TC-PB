
%%

clear

addpath('seawater_ver3_3.1/'); 
addpath('../Plot-Tools/'); 

chl_prof = ncread('6901656_Sprof.nc','CHLA_ADJUSTED'); % units mg/m^3
PRES = ncread('6901656_Sprof.nc','PRES_ADJUSTED'); % Pressure (decibar); 
PSAL = ncread('6901656_Sprof.nc','PSAL_ADJUSTED'); % salinity (psu)
TEMP = ncread('6901656_Sprof.nc','TEMP_ADJUSTED'); % Temp (deg C)
LAT = ncread('6901656_Sprof.nc','LATITUDE'); % 
LON = ncread('6901656_Sprof.nc','LONGITUDE'); % 
JULD = ncread('6901656_Sprof.nc','JULD'); % 
PAR = ncread('6901656_Sprof.nc','DOWN_IRRADIANCE490'); % Temp (deg C)


LAT = repmat(LAT',[size(PRES,1) 1]); 
LON = repmat(LON',[size(PRES,1) 1]); 

depth = sw_dpth(PRES,LAT); 


%%
diveuse = 1:25; 

date1 = datestr(datenum('1950/1/1') + JULD(diveuse(1)),'mm/yyyy'); 
date2 = datestr(datenum('1950/1/1') + JULD(diveuse(end)),'mm/yyyy'); 


zav = -nanmean(depth(:,diveuse),2); 
dz = -gradient(zav); 
cav = nanmean(chl_prof(:,diveuse),2);
tav = nanmean(TEMP(:,diveuse),2);
sav = nanmean(PSAL(:,diveuse),2);
pav = nanmean(PAR(:,:),2);
pav(isnan(pav)) = 0; 

fitter = fit(zav,pav,'exp1'); %  = expfit(pav./nanmax(pav)); 

sum_chl = nansum(dz.*cav); 

disp(['Dates are ' date1 ' to ' date2]); 
fprintf('Total column chl is %2.f mg \n',sum_chl); 
fprintf('Optical Depth is %2.f m \n',1./fitter.b); 

%% 


close all

Ax{5} = subplot(212);

worldmap([-30 -10],[120 200]); 
setm(gca,'parallellabel','off','meridianlabel','off');
add_coastlines; 

plotm(mean(LAT,1),mean(LON,1),'k','linewidth',.5); 
plotm(mean(LAT(:,diveuse),1),mean(LON(:,diveuse),1),'r','linewidth',1); 

Ax{1} = subplot(241);

plot(chl_prof(:,diveuse),-depth(:,diveuse),'color',[.5 .5 .5],'linewidth',.1); 
hold on
plot(cav,zav); 
ylim([-200 0]); 
grid on; box on; 
title('Chl-a (mg/m$^3$)','interpreter','latex')
ylabel('Depth (m)'); 

Ax{2} = subplot(242);

plot(TEMP(:,diveuse),-depth(:,diveuse),'color',[.5 .5 .5],'linewidth',.1); 
hold on
plot(tav,zav,'linewidth',2); 
ylim([-200 0]); 
grid on; box on; 
title('Temp ($^\circ$ C)','interpreter','latex')
set(gca,'yticklabel',''); 


Ax{3} = subplot(243);

plot(PSAL(:,diveuse),-depth(:,diveuse),'color',[.5 .5 .5],'linewidth',.1); 
hold on
plot(sav,zav,'linewidth',2); 
ylim([-200 0]); 
grid on; box on; 
title('Salinity (psu)','interpreter','latex')
set(gca,'yticklabel',''); 

Ax{4} = subplot(244);

plot(PAR(:,diveuse),-depth(:,diveuse),'color',[.5 .5 .5],'linewidth',.1); 
hold on
plot(pav,zav,'linewidth',2); 
ylim([-200 0]); 
grid on; box on; 
title('490 nm (W/m$^2$/nm)','interpreter','latex')
set(gca,'yticklabel',''); 

pos = [6.5 4]; 

set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.05 posy(2)+posy(4)+.015 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',10,'Tag','legtag');
    
end

saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Cyclones-Phytoplankton/Figures/Fig-ARGO.pdf');

