%Jaspreet Singh
%BIOENG 2615 
%ASSIGNMENT V
%10/4/2023 

%% Load data

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

n_channels = numel(filtered_spk);
n_stimulations = numel(TrigON);
stim_time = min(diff(TrigON)); % Use the minimum for setting the window size to avoid any overlaps
window_indices = (1 : round(stim_time*samprate)) - round(0.5*samprate);
stim_spk = {};
idx_spk = {};
for ch = 1:n_channels
    for m = 1 : n_stimulations
        stimloc = find( t > TrigON(m), 1);
        stim_spk{ch}(:, m) = filtered_spk{ch}( stimloc + window_indices );
        idx_spk{ch}(:, m) = stimloc + window_indices;
    end
end

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

%% Extract snippets

snippet_window = [round(-4e-4 * samprate), round(1.1e-3 * samprate)]; % -0.4ms to + 1.1ms
snippets = [];
snippets_idx = [];
for ch = 7:7 % we'll use channel 7 for now
    negthresh = channelMeans(ch) - detection_thr * channelSDs(ch);
    for m=1:n_stimulations
        ii = find(stim_spk{ch}(:, m) < negthresh);
        id = find(diff(ii) > 1);
        ii = [ii(1); ii(id + 1)];
        for idx = 1:numel(ii)
            idx_snippet = ii(idx) + snippet_window;
            if idx_snippet(1) > 0 && idx_snippet(2) <= numel(window_indices)
                snippets = horzcat(snippets, stim_spk{ch}(idx_snippet(1):idx_snippet(2), m));
                snippets_idx = vertcat(snippets_idx, [ii(idx), m]);
            end
        end
    end
end

%% Plot all of the snippets (slow)

figure;
plot(1:38, snippets');

%% Step 10

% % Calculate the number of elements in the time vector
% num_time_points = length(snippets(:, snippet_num));
% 
% % Generate a time vector that matches the number of elements in the snippet data
% time_vector = linspace(snippet_window(1), snippet_window(2), num_time_points);
% % Define the snippet window
% snippet_window = [-0.4:1.1]; % in milliseconds
% 
% % Convert snippet_window to sample indices
% snippet_window_indices = round(snippet_window * (samprate / 1000));

% Plot the snippet with the updated time vector
time_vector = linspace(-0.4, 1.1, 38);
figure;
plot(time_vector, snippets(:, 2000));
axis tight;
title('Single Snippet from Channel 6');
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
grid on;

%%

% % Extract snippets for a specific channel (e.g., channel 1)
% channel_num = 1;
% snippets = zeros(length(snippet_window_indices), n_stimulations);
% 
% for m = 1:n_stimulations
%     stimloc = find(t > TrigON(m), 1);
%     snippets(:, m) = filtered_spk{channel_num}(stimloc + snippet_window_indices);
% end
% 
% % Plot a single snippet from the snippets matrix
% snippet_num = 1;
% figure;
% plot(snippet_window(1):1/samprate:snippet_window(2), snippets(:, snippet_num));
% title('Single Snippet from Channel 1');
% xlabel('Time (ms)');
% ylabel('Amplitude (mV)');
% grid on;

%% Step 11

% Calculate the minimum and maximum values for each snippet
min_values = min(snippets);
max_values = max(snippets);

% Plot min vs. max for each snippet
figure;
plot(min_values, max_values, 'o');
% plot(min_values, max_values, 'x')
title('Minimum vs. Maximum of Snippets','Question 11', FontSize=20);
xlabel('Minimum Value (mV)');
ylabel('Maximum Value (mV)');
grid on;

%% Step 12

% Perform PCA on the snippet data, SVD function discussed in class 
[A, e, W] = svd(snippets);
% Extract eigenvalues
eigenvalues = diag(e);
% Determine the explained variance
explained_variance = cumsum(eigenvalues.^2) / sum(eigenvalues.^2);
% Find the number of components that explain a desired variance (e.g., 95%)
desired_explained_variance = 0.95;
num_components = find(explained_variance >= desired_explained_variance, 1);
figure; plot(explained_variance); scatter(1:numel(explained_variance), explained_variance);
%Ben recommended fprintf for titles during OH
% Display the number of components that explain the desired variance
fprintf('Number of components to explain %.2f%% of variance: %d\n', desired_explained_variance * 100, num_components);

%% Step 14 (!!!)

% Plot 25 components just to have a look
figure; tiledlayout(5, 5, "TileSpacing", "tight");
for component_num = 1:25
    nexttile;
    plot(time_vector, A(:,component_num));
    title(['Component ' num2str(component_num)]);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    grid on;
    axis tight;
    
%     % Justification: You can visually inspect the components to decide which ones are significant.
%     % Typically, look for components with distinct waveforms or those that explain significant variance.
%     if component_num <= significant_components
%         disp(['Component ' num2str(component_num) ' is considered significant.']);
%     end
end

% Plot the eigen-vectors (components) and their magnitudes
significant_components = 3; % Choose the number of significant components to investigate
figure; hold on;
for component_num = 1:significant_components
    plot(time_vector, A(:,component_num));
end
title("Significant Components");
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;
axis tight;
legend(["1", "2", "3"])

%% Step 15 

% Select 2 or 3 significant components
selected_components = 3;

% Plot the amplitudes along a coordinated axis for the selected components
figure; tiledlayout(2, 3);
nexttile;
histogram(W(:, 1));
nexttile;
histogram(W(:, 2));
nexttile;
histogram(W(:, 3));
nexttile;
plot(W(:, 1), W(:, 2), 'x');
xlabel('Amplitude (Component 1)');
ylabel('Amplitude (Component 2)');
title('Amplitudes for Components 1 and 2');
nexttile;
plot(W(:, 2), W(:, 3), 'x');
xlabel('Amplitude (Component 2)');
ylabel('Amplitude (Component 3)');
title('Amplitudes for Components 2 and 3');
nexttile;
plot(W(:, 1), W(:, 3), 'x');
xlabel('Amplitude (Component 1)');
ylabel('Amplitude (Component 3)');
title('Amplitudes for Components 1 and 3');
% plot3(W(:, 1), W(:, 2), W(:, 3), 'x');
% xlabel('Amplitude (Component 1)');
% ylabel('Amplitude (Component 2)');
% zlabel('Amplitude (Component 3)');
% title('Amplitudes for Components 1, 2, and 3');
% if selected_components == 2
%     plot(W(:, 1), W(:, 2), 'x');
%     xlabel('Amplitude (Component 1)');
%     ylabel('Amplitude (Component 2)');
%     title('Amplitudes for Components 1 and 2');
% elseif selected_components == 3
%     plot3(W(:, 1), W(:, 2), W(:, 3), 'x');
%     xlabel('Amplitude (Component 1)');
%     ylabel('Amplitude (Component 2)');
%     zlabel('Amplitude (Component 3)');
%     title('Amplitudes for Components 1, 2, and 3');
% end
grid on;

% This plot shows a map of amplitudes for the selected components. Clusters may indicate similar spike shapes.

% To discriminate units based on these plots, you can use clustering algorithms (e.g., k-means) to group similar spikes together.


%% separate along component 3

idx_unit1 = (W(:, 3) > 0.005);
idx_unit2 = (W(:, 3) < 0.005);

% idx_unclustered = (idx_unit1 | idx_unit2 | idx_unit3) == 0;

% figure; tiledlayout(1, 2);
% nexttile;
% scatter(W(:, 1), W(:, 2), 50, idx, "Marker", "x"); colormap jet;
% xlabel('Amplitude (Component 1)');
% ylabel('Amplitude (Component 2)');
% title('Amplitudes for Components 1 and 2');

figure; tiledlayout(2,3); hold on;
colors = jet;
nexttile;
plot(time_vector, snippets(:, idx_unit1), "Color", [colors(1, :), 0.1]);
title("Cluster 1");
axis tight;
nexttile;
plot(time_vector, snippets(:, idx_unit2), "Color", [colors(end, :), 0.1]);
title("Cluster 2");
axis tight;
nexttile; hold on;
axis tight;
plot(time_vector, snippets(:, idx_unit1), "Color", [colors(1, :), 0.1]);
plot(time_vector, snippets(:, idx_unit2), "Color", [colors(end, :), 0.1]);
axis tight;

%% autocorrelelogram 

unit1_raster = zeros(size(raster_spk{7}));
unit2_raster = zeros(size(raster_spk{7}));

unit1_snippets_idx = sub2ind(size(raster_spk{7}), snippets_idx(idx_unit1, 1), snippets_idx(idx_unit1, 2));
unit1_raster(unit1_snippets_idx) = 1;

unit2_snippets_idx = sub2ind(size(raster_spk{7}), snippets_idx(idx_unit2, 1), snippets_idx(idx_unit2, 2));
unit2_raster(unit2_snippets_idx) = 1;

figure; tiledlayout(1, 3);
nexttile; imagesc(raster_spk{7}');
nexttile; imagesc(unit1_raster');
nexttile; imagesc(unit2_raster');

%% step 17 

% Get the time differences in full recording time
raster_timediff = (1e3/samprate)*diff(idx_spk{7}(raster_spk{7} == 1));
raster1_timediff = (1e3/samprate)*diff(idx_spk{7}(unit1_raster == 1));
raster2_timediff = (1e3/samprate)*diff(idx_spk{7}(unit2_raster == 1));

% Make spike vectors in full recording time
raster_full = zeros(1, size(DATAMATRIX, 2));
raster_full(idx_spk{7}(raster_spk{7} == 1)) = 1;
raster1_full = zeros(1, size(DATAMATRIX, 2));
raster1_full(idx_spk{7}(unit1_raster == 1)) = 1;
raster2_full = zeros(1, size(DATAMATRIX, 2));
raster2_full(idx_spk{7}(unit2_raster == 1)) = 1;

%%

raster_xcorr = xcorr(raster_full);
raster1_xcorr = xcorr(raster1_full);
raster2_xcorr = xcorr(raster2_full);

t_xc = [ -t(end:-1:2) t(1:end) ];

raster_xcorr( find( t_xc == t(1) ,1) ) = 0;
raster1_xcorr( find( t_xc == t(1) ,1) ) = 0;
raster2_xcorr( find( t_xc == t(1) ,1) ) = 0;

%%

figure(22), clf,
subplot(311), histogram( raster_timediff, 100) % , xlim([0 300]),
subplot(312), histogram( raster1_timediff, 100) % , xlim([0 300]),
subplot(313), histogram( raster2_timediff ,100) % , xlim([0 300]),

% make the autocorrelalogram
figure(23), clf,
subplot(321), plot( t_xc*1e3, raster_xcorr ),
subplot(323), plot( t_xc*1e3, raster1_xcorr ),
subplot(325), plot( t_xc*1e3, raster2_xcorr ),
subplot(322), plot( t_xc*1e3, raster_xcorr ), xlim([-300 300]),
subplot(324), plot( t_xc*1e3, raster1_xcorr ), xlim([-300 300]),
subplot(326), plot( t_xc*1e3, raster2_xcorr ), xlim([-300 300]),

%% trial one of  17 and 18
% % Step 17: Select snippets that meet specific criteria
% % You may need to adjust the criteria according to your dataset.
% unit1_ii = find((W(:, 1) > -0.1) & (W(:, 1) < 0.1) & (W(:, 2) > 0.1) & (W(:, 2) < 0.3));
% unit1_snips = snippets(:, unit1_ii);
% 
% % Make sure that ii and unit1_ii are of the same size
% ii = zeros(size(unit1_ii));
% 
% % Step 18: Compute and plot the autocorrelogram for the selected snippets
% unit1_raster = zeros(size(raster));
% unit1_raster(ii(unit1_ii)) = 1;
% unit1_xc = xcorr(unit1_raster);
% 
% % Compute the time vector for autocorrelation
% t_xc = [-t(end:-1:2), t(1:end)];
% 
% % Plot the autocorrelation for the selected snippets
% figure;
% plot(t_xc * 1e3, unit1_xc); % Multiply by 1e3 to convert to milliseconds
% xlabel('Time (ms)');
% ylabel('Autocorrelation');
% title('Autocorrelogram of Selected Snippets');
% 
% % Check if the selected snippets could represent individual units or clusters
% % You may need to analyze the resulting autocorrelation to make this determination.
return;

%% Extra Credit: 
% k-means

rng(1); % set random seed
[idx, C] = kmeans(W(:, 1:2), 2);

figure; tiledlayout(1, 2);
nexttile;
scatter(W(:, 1), W(:, 2), 50, idx, "Marker", "x"); colormap jet;
xlabel('Amplitude (Component 1)');
ylabel('Amplitude (Component 2)');
title('Amplitudes for Components 1 and 2');

nexttile; hold on;
colors = jet;
plot(time_vector, snippets(:, idx == 1), "Color", [colors(1, :), 0.1]);
plot(time_vector, snippets(:, idx == 2), "Color", [colors(end, :), 0.1]);





