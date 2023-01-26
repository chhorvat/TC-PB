function create_local_map(lats,lons)

%subplot(22
axesm('mapprojection','pcarree','maplatlimit',lats,'maplonlimit',lons,...
    'Grid','on','MLineLocation',-135,'MeridianLabel','on',...
    'PLineLocation',[80],...
    'MLabelParallel',60,'ParallelLabel','On',...
    'labelformat','compass','labelrotation','on','glinestyle',':','glinewidth',0.5,...
    'gcolor',[0.500 0.500 0.500]);
set(gca,'fontname','helvetica','fontsize',12,'xminortick','on','yminortick','on')

hold on


% h1 = pcolorm(latback,lonback,ones(length(latback),length(lonback)));
% set(h1, 'facecolor', [198,219,239]/256);
box off

load coastlines
% plotm(coastlat, coastlon)
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[217,217,217]/256)
