%% Load the data

load 0900702VisuallyEvoked.mat

%% Plot the raw data from one channel

channel_data = Wideband_data{1};
nsamples = numel(channel_data);
t = 0:(1/samprate):(nsamples - 1)/samprate; % time vector in seconds
xlabel('Time(mS)')
ylabel('Amplitude (mV)')
figure;
plot(t, channel_data);

%% Filter the data and plot from one channel

spk_cutoff = 300; % Hz
spk_filter_order = 2; %CHOSEN 2nd ORDER

[BSPK, ASPK] = butter(spk_filter_order, [spk_cutoff] / (samprate/2), 'low');

for n = 1:16
    DATAMATRIX(n,:) = Wideband_data{n, 1};
end

DATAMATRIX = double(DATAMATRIX);

% loop for filtering a channel at a time and then taking its respective mean and std
filtered_spk = {};
for i = 1:size(Wideband_data, 1)
    filtered_spk{i} = filter(BSPK, ASPK, DATAMATRIX(i,:));
    channelMeans(i) = mean(filtered_spk{i});
    channelSDs(i) = std(filtered_spk{i});
end

figure;
plot(t, filtered_spk{1});
% xlim([10, 11]);

%% Test the bandpass filter with white noise

n = randn(round(2*samprate), 20);
n_ms = n - mean(n, 1);
for m=1:20
    % noise
    y = filter(BSPK, ASPK, n(:,m));
    yf = fft(y);
    y_psd(:, m)= abs(yf.*yf);
    y_psd_shift(:, m) = abs(fftshift(yf.*yf, 1));
    % noise, mean-subtracted
    y_ms = filter(BSPK, ASPK, n_ms(:,m));
    yf_ms = fft(y_ms);
    y_psd_ms(:, m)= abs(yf_ms.*yf_ms);
    y_psd_shift_ms(:, m) = abs(fftshift(yf_ms.*yf_ms, 1));
end
x=[0:size(n, 1)-1]*samprate/size(n, 1);
x_shift=[-size(n, 1)/2:size(n, 1)/2-1]*samprate/size(n, 1);

figure; tiledlayout(2, 2);
nexttile; plot(x, mean(y_psd, 2)); xlabel("Frequency (Hz)"); ylabel("Power"); title("FFT");
nexttile; plot(x_shift, mean(y_psd_shift, 2)); xlabel("Frequency (Hz)"); ylabel("Power"); title("FFT (Shifted)"); xticks([-300, 0, 300]);
% Plot the mean-subtracted noise
nexttile; plot(x, mean(y_psd_ms, 2)); xlabel("Frequency (Hz)"); ylabel("Power");
nexttile; plot(x_shift, mean(y_psd_shift_ms, 2)); xlabel("Frequency (Hz)"); ylabel("Power");

%% Gather filtered signal into stimulation trials

n_channels = numel(filtered_spk);
n_stimulations = numel(TrigON);
stim_time = min(diff(TrigON)); % Use the minimum for setting the window size to avoid any overlaps
window_indices = (1 : round(stim_time*samprate)) - round(0.5*samprate);
stim_spk = {};
for ch = 1:n_channels
    for m = 1 : n_stimulations
        stimloc = find( t > TrigON(m), 1);
        stim_spk{ch}(:, m) = filtered_spk{ch}( stimloc + window_indices );
    end
end


%% Plot PSD for averaged stimuli and individual stimulus

ch = 1;
y_psd = zeros(size(stim_spk{1}, 1), n_stimulations);
y_psd_shift = zeros(size(stim_spk{1}, 1), n_stimulations);
for m= 1:n_stimulations
    y = stim_spk{ch}(:, m);
    yf = fft(y);
    y_psd(:, m)= abs(yf.*yf);
    y_psd_shift(:, m) = abs(fftshift(yf.*yf, 1));
end
x=[0:size(y, 1)-1]*samprate/size(y, 1);
x_shift=[-size(y, 1)/2:size(y, 1)/2-1]*samprate/size(y, 1);

example_stim = 10;
figure; tiledlayout(2, 2);
% Average PSD from all trials
nexttile;
plot(x, mean(y_psd, 2)); xlabel("Frequency (Hz)"); ylabel("Power"); title("FFT"); axis tight;
nexttile;
plot(x_shift, mean(y_psd_shift, 2)); xlabel("Frequency (Hz)"); ylabel("Power"); title("FFT"); xlim([-500, 500]);
% PSD from one trial
nexttile;
plot(x, y_psd(:, example_stim)); xlabel("Frequency (Hz)"); ylabel("Power"); title("FFT"); axis tight;
nexttile;
plot(x_shift, y_psd_shift(:, example_stim)); xlabel("Frequency (Hz)"); ylabel("Power"); title("FFT (Shifted)"); xlim([-500, 500]);

%% Compute the spectrogram over full signal

% calculate the spectrogram
% setup the window length and skip
win_len = samprate;
skp_len = samprate / 10;
win_n = floor((numel(t) - win_len) / skp_len);
f_max = 100;
% setup the spectrogram time and frequency vectors
t_spec = [0:win_n-1] * (skp_len/samprate);
f_win = [0:win_len-1] * (samprate/win_len);
f_i = find( f_win <= f_max );
f_spec = f_win(f_i);
spectro = zeros( numel(f_i), win_n );
% compute a simple spectrogram
for m = 1:win_n
    tmpdata = filtered_spk{4}( round((m-1) * skp_len + [1:win_len]) );
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
set(gca,'CLim',[0 5e-4]),
subplot(212),
pcolor( t_spec, f_spec, 10*log10(spectro) ),
view(2), colormap jet,
shading interp,
axis tight,
xlabel('Time (sec)'),
ylabel('Frequency (Hz)'),
colorbar,
set(gca,'CLim',[-60 -30]),

%% Compute the spectrogram over visual stimulation trials

% calculate the spectrogram
% setup the window length and skip
ch = 3;
win_len = samprate / 5;
skp_len = samprate / 10;
win_n = floor((size(stim_spk{ch}, 1) - win_len) / skp_len);
f_max = 100;
% setup the spectrogram time and frequency vectors
t_spec = [0:win_n-1] * (skp_len/samprate);
f_win = [0:win_len-1] * (samprate/win_len);
f_i = find( f_win <= f_max );
f_spec = f_win(f_i);
spectro = zeros( numel(f_i), win_n, 1 );
% compute a simple spectrogram
for s = 1:1
    for m = 1:win_n
        tmpdata = stim_spk{ch}( round((m-1) * skp_len + [1:win_len]), s);
        tmpspec_psd = abs(fft(tmpdata));
        tmpspec_psd = (tmpspec_psd/sum(tmpspec_psd)).^2;
        spectro(:, m, s) = tmpspec_psd(f_i);
    end
end
% show the result
figure(13), clf,
subplot(211),
pcolor( t_spec - (0.5 - win_len/samprate), f_spec, mean(spectro, 3) ),
view(2), colormap jet,
shading interp,
axis tight,
xlabel('Time from Stimulus Onset (sec)'),
ylabel('Frequency (Hz)'),
colorbar,
set(gca,'CLim',[0 2e-3]),
subplot(212),
pcolor( t_spec - (0.5 - win_len/samprate), f_spec, 10*log10(mean(spectro, 3)) ),
view(2), colormap jet,
shading interp,
axis tight,
xlabel('Time from Stimulus Onset (sec)'),
ylabel('Frequency (Hz)'),
colorbar,
% set(gca,'CLim',[-60 -30]),

%% Compute the spectrogram over all visual stimulation trials

% calculate the spectrogram
% setup the window length and skip
ch = 3;
win_len = samprate / 5;
skp_len = samprate / 10;
win_n = floor((size(stim_spk{ch}, 1) - win_len) / skp_len);
f_max = 100;
% setup the spectrogram time and frequency vectors
t_spec = [0:win_n-1] * (skp_len/samprate);
f_win = [0:win_len-1] * (samprate/win_len);
f_i = find( f_win <= f_max );
f_spec = f_win(f_i);
spectro = zeros( numel(f_i), win_n, n_stimulations );
% compute a simple spectrogram
for s = 1:n_stimulations
    for m = 1:win_n
        tmpdata = stim_spk{ch}( round((m-1) * skp_len + [1:win_len]), s);
        tmpspec_psd = abs(fft(tmpdata));
        tmpspec_psd = (tmpspec_psd/sum(tmpspec_psd)).^2;
        spectro(:, m, s) = tmpspec_psd(f_i);
    end
end
% show the result
figure(13), clf,
subplot(211),
pcolor( t_spec - (0.5 - win_len/samprate), f_spec, mean(spectro, 3) ),
view(2), colormap jet,
shading interp,
axis tight,
xlabel('Time from Stimulus Onset (Seconds)'),
ylabel('Frequency (Decibels)'),
title("Title 1", "Title 2", 'Color', 'b', 'FontWeight', 'bold'); 
colorbar,
set(gca,'CLim',[0 2e-3]),
subplot(212),
pcolor( t_spec - (0.5 - win_len/samprate), f_spec, 10*log10(mean(spectro, 3)) ),
view(2), colormap jet,
shading interp,
axis tight,
xlabel('Time from Stimulus Onset (seconds)'),
ylabel('Frequency (Decibels)'),
title("Title 1", "Title 2", 'Color', 'b', 'FontWeight', 'bold'); 
colorbar,
% set(gca,'CLim',[-60 -30]),

% Main title
sgt = sgtitle('Question 23: Average spectrogram for Stimuli from Trials', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial');
% annotation('textbox', [0.5, 0.03, 0.05, 0.05], 'String', 'Subtitle Text', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'center');

set(gcf,"color",      "w");

return;
%% Extra Credit: complete the process of affinity development 

%% Extra Credit IV: current source density of avg channel stimulation data 
% 
% for	m=1:size(stim_spk,1),
% 
%     csd(m,:)=-diff(stim_spk(m,:),2,2);
% 
% end;
% 
% ch = 5
% clf,
% pcolor(t,stim_spk{ch},csdF),
% view(2),
% shading('interp'),
% set(gca,'CLim',[-1	1]*0.5*min(csd(:))),

% Assuming you have your stim_spk data organized as stim_spk{channel#}(time)

% Assuming you have your stim_spk data organized as stim_spk{channel#}(time)
% Assuming you have your stim_spk data organized as stim_spk{channel#}(time)

% Define the channel number you want to analyze
ch = 5 ;

% Extract the data for the specified channel
data = stim_spk{ch};

% Compute the CSD for the specified channel
csd = -diff(data, 2, 2);

% Define the time vector based on the dimensions of csd
t = (1:size(csd, 2));

% Create a pcolor plot
clf;
pcolor(t, 1:size(csd, 1), csd');
view(2);
shading('interp');
% caxis([-1 1] * 0.5 * min(csd(:)); % Set color limits

xlabel('Time');
ylabel('Stimulus Onset (seconds)');
title(['Current Source Density (CSD) for Channel ' num2str(ch)]);

colorbar;





