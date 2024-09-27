%% Step 5: Sort based on PCA coefficients and redo rasters
%

% simple clustering 
thre1 = -0.02;
thre2 =  0.0;
unit1_i = find( vv(:,2) <  thre2 );
unit2_i = find( vv(:,2) >= thre2 );

% make unit raster
data_raster1 = zeros(size(data_raster));
data_raster2 = zeros(size(data_raster));

data_raster1( data_raster_ii(unit1_i) ) = 1;
data_raster2( data_raster_ii(unit2_i) ) = 1;

% plot the rasters
figure(21), clf,
subplot(311), plot( t, data_raster ), axis tight, grid on,
subplot(312), plot( t, data_raster1 ), axis tight, grid on,
subplot(313), plot( t, data_raster2 ), axis tight, grid on,


%% Step 6: Compute the spike-time-differences and autocorrelalogram
%

raster_timediff = (1e3/fs)*diff(find(data_raster));
raster1_timediff = (1e3/fs)*diff(find(data_raster1));
raster2_timediff = (1e3/fs)*diff(find(data_raster2));

raster_xcorr = xcorr( data_raster );
raster1_xcorr = xcorr( data_raster1 );
raster2_xcorr = xcorr( data_raster2 );

t_xc = [ -t(end:-1:2) t(1:end) ];

raster_xcorr( find( t_xc == t(1) ,1) ) = 0;
raster1_xcorr( find( t_xc == t(1) ,1) ) = 0;
raster2_xcorr( find( t_xc == t(1) ,1) ) = 0;


% make the unit time different histogram
figure(22), clf,
subplot(311), hist( raster_timediff, 100), xlim([0 300]),
subplot(312), hist( raster1_timediff, 100), xlim([0 300]),
subplot(313), hist( raster2_timediff ,100), xlim([0 300]),

% make the autocorrelalogram
figure(23), clf,
subplot(321), plot( t_xc*1e3, raster_xcorr ),
subplot(323), plot( t_xc*1e3, raster1_xcorr ),
subplot(325), plot( t_xc*1e3, raster2_xcorr ),
subplot(322), plot( t_xc*1e3, raster_xcorr ), xlim([-300 300]),
subplot(324), plot( t_xc*1e3, raster1_xcorr ), xlim([-300 300]),
subplot(326), plot( t_xc*1e3, raster2_xcorr ), xlim([-300 300]),


