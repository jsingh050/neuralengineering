%jaspreet singh
%assignment 4
%bioeng2615
%9/25/23

%% Load the data

load 0900702VisuallyEvoked.mat

%% Plot the raw data from one channel

channel_data = Wideband_data{1};
nsamples = numel(channel_data);
t = 0:(1/samprate):(nsamples - 1)/samprate; % time vector in seconds
figure;

plot(t, channel_data );

title("Raw Data from 16 Channels", FontSize=20)

xlabel("Time (mS)")
ylabel("Amplitude (mV)")

%set parameters to filter using spike data (same as past assignments)

spk_low_cutoff = 300; % Hz
spk_high_cutoff = 3000; % Hz
spk_filter_order = 2; %CHOSEN 2nd ORDER

[BSPK, ASPK] = butter(spk_filter_order/2, [spk_low_cutoff, spk_high_cutoff] / (samprate/2), 'bandpass');
%save data as doubles 
for n = 1:16
    DATAMATRIX(n,:) = Wideband_data{n, 1};
end

DATAMATRIX = double(DATAMATRIX);

% loop for filtering a channel at a time and then taking its respective
% mean and std (taken from last assignment)
filtered_spk = {};
for i = 1:size(Wideband_data, 1)
    filtered_spk{i} = filter(BSPK, ASPK, DATAMATRIX(i,:));
    channelMeans(i) = mean(filtered_spk{i});
    channelSDs(i) = std(filtered_spk{i});
end


%% Plot filtered data to compare to raw

figure;
title("Comparing Raw and Filtered Data for Visual Stimuli")
tiledlayout(3, 1);
x_limits = [20, 30];

nexttile;
plot(t, channel_data); % raw

title("Raw Data", FontSize=20);
xlabel("Time (mS)")
ylabel("Amplitude (mV)")
axis tight;

xlim(x_limits);
set(gca,"box",        "off",...
    "XColor",     "k",...
    "YColor",     "k",...
    "TickDir",    "out",...
    "TickLength", [0.01,0.01],...
    "FontSize",    13,...
    "LineWidth",   1);

nexttile;
plot(t, filtered_spk{1}); % filtered
title("Spike Filtered Data");
xlabel("Time (ms)")
ylabel("Amplitude (mV)")
axis tight;

xlim(x_limits);
set(gca,"box",        "off",...
    "XColor",     "k",...
    "YColor",     "k",...
    "TickDir",    "out",...
    "TickLength", [0.01,0.01],...
    "FontSize",    14,...
    "LineWidth",   1);

nexttile;
TrigONOFF = horzcat(TrigON', TrigOFF');
stim_binary = zeros(numel(channel_data), 1);
for stim = 1:size(TrigONOFF, 1)
    onoffsets = TrigONOFF(stim, :).*samprate;
    stim_binary(onoffsets(1):onoffsets(2)) = 1;
end
plot(t, stim_binary);
title("Visual Stimuli");
xlabel("Time (mS)")
ylabel("Amplitude (mV)")
axis tight;

xlim(x_limits);
set(gca,"box",        "off",...
    "XColor",     "k",...
    "YColor",     "k",...
    "TickDir",    "out",...
    "TickLength", [0.01,0.01],...
    "FontSize",    13,...
    "LineWidth",   1);

set(gcf,"color",      "w");

sgtitle("Filtered vs Raw with Visual Stimuli for Question 1", FontSize=20)
%% From the assignment ben sent this pseduocode for us to use in this assignment through part 3 

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


%% Spike detection and raster creation

detection_thr = 3; % detection threshold
raster_spk = cell(1, n_channels);

for ch = 1:n_channels
    negthresh = channelMeans(ch) - detection_thr * channelSDs(ch);
    raster_spk{ch} = zeros(numel(window_indices), n_stimulations);
    for m=1:n_stimulations
        ii = find(stim_spk{ch}(:, m) < negthresh);
        id = find(diff(ii) > 1);
        ii = [ii(1); ii(id + 1)];
        raster_spk{ch}(ii, m) = 1;
    end
end


titles = {'Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5', 'Channel 6', 'Channel 7', 'Channel 8', 'Channel 9', 'Channel 10', 'Channel 11', 'Channel 12', 'Channel 13', 'Channel 14', 'Channel 15', 'Channel 16'}
figure;
tld = tiledlayout(4, 4);
for ch = 1:n_channels
    nexttile;
    subtitle(titles(ch));
    
    
    plot_raster(window_indices./samprate, raster_spk{ch}(:, 1)); % Plot first stimulus as an example

    set(gca,"box",        "off",...
    "XColor",     "k",...
    "YColor",     "k",...
    "TickDir",    "out",...
    "TickLength", [0.01,0.01],...
    "FontSize",    10,...
    "LineWidth",   1);

set(gcf,"color",      "w");

end

title(tld, "Raster Plots for 16 Channels for Question #3", FontSize=20);



%% Plot PSTH for Question 4 

psth_dt = 50; % milliseconds
window_times = window_indices./samprate;
time_psth = min(window_times) : 50/1000 : max(window_times);
psth = cell(1, n_channels);
for ch = 1:n_channels
    psth{ch} = zeros(numel(time_psth), n_stimulations);
    for m = 1:numel(time_psth)
        hist_index = find((window_times >= time_psth(m)) & (window_times < time_psth(m) + psth_dt/1000));
        psth{ch}(m,:) = sum(raster_spk{ch}(hist_index, :), 1);
    end
end

figure;
tld = tiledlayout(4, 4);
for ch = 1:n_channels
    nexttile;
    
    plot(time_psth, psth{ch}(:, 1)); % Plot first stimulus as an example
    title(titles(ch));
    xlabel("Time")
    ylabel("Stimulus")
    axis tight;
    grid on;
end

sgtitle(tld, "PSTH for Each Channel for Question #4", FontSize=20);
set(gcf,"color",      "w");

% %% Plot average PSTH (!) fix using plot instead of title 
% n_channels = numel(filtered_spk);
% titles = {' 1', ' 2', ' 3', ' 4', ' 5', ' 6', ' 7', ' 8', ' 9', ' 10', ' 11', ' 12', ' 13', ' 14', ' 15', ' 16'}
% figure;
% title("Question #5");
% % tld = tiledlayout(4, 4);
% for ch = 1:n_channels
%     subplot(4, 4, ch);
%     title(titles{ch});
%      text(0.5, 0.5, titles{ch}, 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     % title(titles{ch});
%     % xlabel("x label")
%     % ylabel("y label")
%     plot(time_psth, mean(psth{ch}(:, :), 2)); % Take the mean along the  dimension for the stimulus (2nd)
%     axis tight;
%     grid on;
% end
% 
% 
% 
% 
%     set(gca,"box",        "off",...
%     "XColor",     "k",...
%     "YColor",     "k",...
%     "TickDir",    "out",...
%     "TickLength", [0.01,0.01],...
%     "FontSize",    12,...
%     "LineWidth",   1);
% 
% set(gcf,"color",      "w");
% 
% 
% return; 

%% Question 5 

n_channels = numel(filtered_spk);
titles = {' Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5', 'Channel 6', 'Channel 7', 'Channel 8', 'Channel 9', 'Channel 10', 'Channel 11', 'Channel 12', 'Channel 13', 'Channel 14', 'Channel 15', 'Channel 16'};
figure('Position', [100, 100, 800, 800]); % Set figure size as needed

for ch = 1:n_channels
    subplot(4, 4, ch);
    plot(time_psth, mean(psth{ch}(:, :), 2),'r'); % Take the mean along the dimension for the stimulus (2nd)
    axis tight;
    grid on;
    xlim([min(time_psth), max(time_psth)]); % Ensure consistent x-axis limits
    ylim([min(mean(psth{ch}(:, :), 2)), max(mean(psth{ch}(:, :), 2))]); % Ensure consistent y-axis limits
    
    % Add the label at the top with red color
    x = xlim;
    y = ylim;
    text(x(1) + 0.5*(x(2)-x(1)), y(2), titles{ch}, 'FontSize', 14, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'red');
    
    % Label the axes
    xlabel('Time');
    ylabel('Intensity');
end

set(gcf, 'color', 'w');

% Add the figure title
sgtitle('Question 5'); %%(!! edit)

%% Section 6  one channel 

%Define the time window for calculating average PSTH before and during stimulation
pre_stim_duration = 0.5;  % seconds
during_stim_duration = 1.0;  % seconds


% Initialize arrays to store average PSTH before and during stimulation for all channels
avg_psth_pre_stim = zeros(n_channels, 1);
avg_psth_during_stim = zeros(n_channels, 1);
avg_spike_rate_change = zeros(n_channels, 1);

% Calculate the average PSTH before and during stimulation for all channels
for ch = 1:n_channels
    % Extract the PSTH data for the current channel
    channel_psth = psth{ch}; % Use the correct data source, 'psth'

    % Find the indices corresponding to the pre-stimulation and during-stimulation periods
    pre_stim_indices = find(time_psth < (time_psth(1) + pre_stim_duration)); %we need to adjust time
    during_stim_indices = find(time_psth >= (time_psth(1) + pre_stim_duration) & time_psth < (pre_stim_duration + time_psth(1) + during_stim_duration));

    % Calculate the mean of the PSTH data for the pre-stimulation and during-stimulation indices
    avg_psth_pre_stim(ch) = mean(channel_psth(:, pre_stim_indices), 'all');
    avg_psth_during_stim(ch) = mean(channel_psth(:, during_stim_indices), 'all');

    % Calculate the change in spike rate (during vs. before) within the loop
    avg_spike_rate_change(ch) = avg_psth_during_stim(ch) - avg_psth_pre_stim(ch);



end
%
% Plot the change in spike rate as a function of depth (channel number)
figure;
plot(1:n_channels, avg_spike_rate_change, 'm');
xlabel('Channel Number (Depth)');
ylabel('Average Spike Rate Change (During vs. Before)');
% title('Average Spike Rate Change vs. Depth', FontSize=20);
grid on;
yline(0, "--k")

set(gcf, 'color', 'w');

% Add the figure title
sgtitle('Question 6');


% %%
% % Section 6
% 
% % Define the time window for calculating average PSTH before and during stimulation
% pre_stim_duration = 0.5;  % seconds (look at earlier plot to define .5 seconds)
% during_stim_duration = 1.0;  % seconds
% 
% % Initialize arrays to store average PSTH before and during stimulation for all channels
% avg_psth_pre_stim = zeros(n_channels, 5);
% avg_psth_during_stim = zeros(n_channels, 5);
% avg_spike_rate_change = zeros(n_channels, 5);
% 
% % Calculate the average PSTH before and during stimulation for all channels
% for ch = 1:n_channels
%     % Extract the PSTH data for the current channel
%     channel_psth = psth{ch}; % Use the correct data source, 'psth'
% 
%     % Find the indices corresponding to the pre-stimulation and during-stimulation periods
%     pre_stim_indices = find(time_psth < pre_stim_duration);
%     during_stim_indices = find(time_psth >= pre_stim_duration & time_psth < (pre_stim_duration + during_stim_duration));
% 
%     % Calculate the mean of the PSTH data for the pre-stimulation and during-stimulation indices
%     avg_psth_pre_stim(ch) = mean(channel_psth(:, pre_stim_indices), 'all');
%     avg_psth_during_stim(ch) = mean(channel_psth(:, during_stim_indices), 'all');
% 
%     % Calculate the change in spike rate (during vs. before) within the loop
%     avg_spike_rate_change(ch) = avg_psth_during_stim(ch) - avg_psth_pre_stim(ch);
% end
% 
% % Create one figure with 16 subplots
% figure;
% 
% for ch = 1:n_channels
%     subplot(4, 4, ch);
%     plot(time_psth, avg_spike_rate_change(ch), 'r');
%     xlabel('Time (s)');
%     ylabel('Average Spike Rate Change (During vs. Before)');
%     title(['Channel ' num2str(ch)], 'FontSize', 12);
%     grid on;
% end
% 
% % Add the figure title
% sgtitle('Question 6');


