load coastlines
plotm(coastlat, coastlon,'k','linewidth',2)
hold on
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor',[217,217,217]/256)
