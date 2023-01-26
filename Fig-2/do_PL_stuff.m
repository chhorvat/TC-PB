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
HF = hoverfactor(~isnan(hoverfactor)); 
[alpha,xmin,L] = plfit(HF/min(HF)); 

%%
[p,gof] = plpva(HF/min(HF),xmin,'range',[4.001:0.001:6.001]); 