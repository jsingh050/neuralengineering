
% start clean
clear all
close all

load sim1_new_data2

fs = 24414;     % hz

t = [1:length(data)]/fs;

% filter the data for spike analysis
Wn = (300/fs) * 0.5;
[filt_b, filt_a] = butter( 2, Wn, 'high');

data_filtered = filter( filt_b, filt_a, data );

% plot the data and filtered data
figure(1), clf,
subplot(221), plot( t, data ),
axis tight, grid on, set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), title('Data','FontSize',14),
subplot(222), plot( t, data_filtered ),
axis tight, grid on, set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), title('Filtered Data','FontSize',14),
subplot(223), plot( t, data ),
axis tight, grid on, xlim([100e-3 1100e-3]), set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), 
subplot(224), plot( t, data_filtered ),
axis tight, grid on, xlim([100e-3 1100e-3]), set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), 

% find the negative deflections and make the raster plot
thr = mean(data_filtered) - 3.5*std( data_filtered);
thr = -40;
disp(sprintf('  threshold = %f', thr));

tmp_ii = find( data_filtered < thr );
tmp_di = find( diff( tmp_ii ) > 1);         % count the first crossing only
data_raster_ii = [ tmp_ii(1) tmp_ii(tmp_di+1) ];

data_raster = zeros(size(data));
data_raster( data_raster_ii ) = 1;
n_raster_hits = sum( data_raster );

% plot the raster
figure(2), clf,
subplot(311), plot(t, data_raster), 
axis tight, grid on, set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',14), ylabel('Raster Hits','FontSize',16), 
subplot(312), plot(t, data_raster), fatlines(1.5),
axis tight, grid on, xlim([100e-3 1100e-3]), set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',14), ylabel('Raster Hits','FontSize',16),
subplot(313), hist( 1e3*diff(data_raster_ii)/fs, 100),
axis tight, grid on, set(gca,'FontSize',12),
xlabel('Time Difference (ms)','FontSize',14), ylabel('Count','FontSize',16),


% get snippets and do pile plots centered versus non-centered
snips_i = [-20:179];
t_snip = snips_i/fs;

snips = zeros( length(snips_i), n_raster_hits );
snips_minctr = snips;

for m = 1:n_raster_hits,
  snips(:,m) = data_filtered( data_raster_ii(m) + snips_i );
end

% pile plot
figure(3), clf,
plot( 1e3*t_snip, snips, 'Color', [0.5 0.5 0.5]), axis tight, grid on, hold on,
plot( 1e3*t_snip, mean(snips ,2), 'LineWidth', 2, 'Color', [0 0 0.8]), hold off,
xlabel( 'Time (ms)' ), ylabel( 'Snip Amplitude' ),

% do pca here on the snippet matrix
[uu, ss, vv] = svd( snips );

% get the diagonal elements
ss_diag = diag( ss );

% plot components
figure(4), clf,
for m = 1:12, 
  subplot(3,4,m),
  plot( 1e3*t_snip, uu(:,m) ), axis tight, grid on,
  title( num2str( m ) ),
end

% plot the sigma values and explained variance
figure(5), clf,
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
figure(6), clf,
subplot(311), plot( t, data_raster ), axis tight, grid on,
subplot(312), plot( t, data_raster1 ), axis tight, grid on,
subplot(313), plot( t, data_raster2 ), axis tight, grid on,

figure(7), clf,
subplot(211),
plot( 1e3*t_snip, snips(:,unit1_i), 'Color', [0.5 0.5 0.5]), axis tight, grid on, hold on,
plot( 1e3*t_snip, mean(snips (:,unit1_i),2), 'LineWidth', 2, 'Color', [0 0 0.8]), hold off,
xlabel( 'Time (ms)' ), ylabel( 'Unit1 Amplitude' ),
subplot(212),
plot( 1e3*t_snip, snips(:,unit2_i), 'Color', [0.5 0.5 0.5]), axis tight, grid on, hold on,
plot( 1e3*t_snip, mean(snips(:,unit2_i) ,2), 'LineWidth', 2, 'Color', [0 0 0.8]), hold off,
xlabel( 'Time (ms)' ), ylabel( 'Unit2 Amplitude' ),


%
%
%

% filter the data for LFP analysis
Wn_lfp = (300/fs) * 0.5;
[filt_b, filt_a] = butter( 2, Wn_lfp, 'low');

data_filtered2 = filter( filt_b, filt_a, data );

figure(11), clf,
subplot(221), plot( t, data ),
axis tight, grid on, set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), title('Data','FontSize',14),
subplot(222), plot( t, data_filtered2 ),
axis tight, grid on, set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), title('Filtered Data','FontSize',14),
subplot(223), plot( t, data ),
axis tight, grid on, xlim([100e-3 1100e-3]), set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), 
subplot(224), plot( t, data_filtered2 ),
axis tight, grid on, xlim([100e-3 1100e-3]), set(gca,'FontSize',12),
xlabel('Time (sec)','FontSize',12), ylabel('Amplitude (uV)','FontSize',14), 


% make the frequency axis for the psd
f = [0:length(t)-1]*(fs/length(t));

% calculate the psd for the original data and filtered data
psd_data = abs(fft( data - mean(data) ));
psd_data = (psd_data/sum(psd_data)).^2;

psd_data_filt = abs(fft( data_filtered2 - mean(data_filtered2) ));
psd_data_filt = (psd_data_filt/sum(psd_data_filt)).^2;


% plot the psd in linear and log scale
figure(12),
subplot(221), plot(f, psd_data ), axis tight, grid on, xlim([0 100]),
ylabel('Magnitude'), xlabel('Frequency (Hz)'), 
subplot(222), plot(f, psd_data_filt ), axis tight, grid on, xlim([0 100]),
ylabel('Magnitude'), xlabel('Frequency (Hz)'), 
subplot(223), plot(f, 10*log10(psd_data) ), axis tight, grid on, xlim([0 100]),
ylabel('Log Magnitude'), xlabel('Frequency (Hz)'), 
subplot(224), plot(f, 10*log10(psd_data_filt) ), axis tight, grid on, xlim([0 100]),
ylabel('Log Magnitude'), xlabel('Frequency (Hz)'), 


% calculate the spectrogram

% setup the window length and skip
win_len = fs;
skp_len = fs/2;
win_n = floor((length(t)-win_len)/skp_len);
f_max = 100;

% setup the spectrogram time and frequency vectors
t_spec = [0:win_n-1]*(skp_len/fs);
f_win = [0:win_len-1]*(fs/win_len);
f_i = find( f_win <= f_max );

f_spec = f_win(f_i);
spectro = zeros( length(f_i), win_n );

% compute a simple spectrogram
for m = 1:win_n,
    tmpdata = data_filtered2( (m-1)*skp_len + [1:win_len] );
    tmpspec_psd = abs(fft(tmpdata));
    tmpspec_psd = (tmpspec_psd/sum(tmpspec_psd)).^2;
    spectro(:,m) = tmpspec_psd(f_i);
end

% show the result
figure(13), clf,
subplot(211),
  pcolor( t_spec, f_spec, spectro ),
  view(2), colormap jet,
  shading interp,
  axis tight,      
  xlabel('Time (sec)'),
  ylabel('Frequency (Hz)'),
  colorbar,
  set(gca,'CLim',[0 2e-3]),
subplot(212), 
  pcolor( t_spec, f_spec, 10*log10(spectro) ),
  view(2), colormap jet,
  shading interp,
  axis tight,      
  xlabel('Time (sec)'),
  ylabel('Frequency (Hz)'),
  colorbar,
  set(gca,'CLim',[-60 -30]),



