%% Step 3: Analyze data - get snippets and do pile plots centered versus non-centered
%

snips_i = [-20:179];
t_snip = snips_i/fs;

snips = zeros( length(snips_i), n_raster_hits );
snips_minctr = snips;

for m = 1:n_raster_hits,
  snips(:,m) = data_filtered( data_raster_ii(m) + snips_i );
  [ ~, tmp_min_i ] = min( snips(:,m), [], 1 );
  snips_minctr(:,m) = data_filtered( data_raster_ii(m) + snips_i + tmp_min_i - 35 );
end


figure(5), clf,
subplot(211), pilePlot1(snips),
subplot(212), pilePlot1(snips_minctr),



%% Step 4: Analyze the snippets - PCA
%

% do pca here on the snippet matrix
[uu, ss, vv] = svd( snips );

% get the diagonal elements
ss_diag = diag( ss );

% plot the components
figure(11), clf,
plotmany( t_snip, uu(:,1:12) ),

% plot the sigma values and explained variance
figure(12), clf,
subplot(211), plot( ss_diag ),
subplot(212), plot( cumsum(ss_diag.^2)/sum(ss_diag.^2) )

% plot the components -- decide on clustering
figure(13), clf,
subplot(221), hist( vv(:,1), 100), 
axis tight, grid on, set(gca,'FontSize',12),
subplot(222), hist( vv(:,2), 100), 
axis tight, grid on, set(gca,'FontSize',12),
subplot(223), hist( vv(:,3), 100), 
axis tight, grid on, set(gca,'FontSize',12),
subplot(224), plot( vv(:,1), vv(:,2), 'x' ), 
axis tight, grid on, set(gca,'FontSize',12),