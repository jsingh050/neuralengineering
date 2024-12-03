
%Jaspreet Singh
%BIOENG 2615
%9/14/2023
%Assignment 3

%load data
load 0900702VisuallyEvoked.mat


%variables 
%TRIGON= ; %put in the trig values (for this assignment, per the panapto
%video, none is needed)
%TRIGOFF= ; %look at Trig values 

%general information: raw data from channels is the wideband data

    
%% Outputting Histogram for 16 channels 

%mean, standard deviation video https://youtu.be/uM-ZdNe0YGA

spk_low_cutoff = 300; % Hz
spk_high_cutoff = 3000; % Hz
spk_filter_order = 2; %CHOSEN 2nd ORDER

[BSPK, ASPK] = butter(spk_filter_order/2, [spk_low_cutoff, spk_high_cutoff] / (samprate/2), 'bandpass');

for n = 1:16
    DATAMATRIX(n,:) = Wideband_data{n, 1};
end

DATAMATRIX = double(DATAMATRIX);

% loop for filtering a channel at a time and then taking its respective mean and std
for i = 1:size(Wideband_data, 1)
    filtered_spk(i,:) = filter(BSPK, ASPK, DATAMATRIX(i,:));
    channelMeans(i) = mean(filtered_spk(i,:));
    channelSDs(i) = std(filtered_spk(i,:));
end

%%
% 
% % 1) forloops are helpful for each channel
% for i = 1:16 %through channels 1 - 16
%     channelData = Wideband_data{i}; % Extract data for the current channel, each i can be diff, channel data variable
%     channelMeans(i) = mean(channelData); % Calculate the mean
%     channelStdDevs(i) = std(channelData); % Calculate the standard deviation
% end

% disp('Mean of 16 channels:');
% disp(channelMeans);
% 
% disp('Standard Deviation of 16 channels:');
% disp(channelStdDevs);


%2 plotting histogram



% Loop through each channel
for i = 1:16
    subplot(4, 4, i); 
    channelData = filtered_spk(i,:); 
    % Create a histogram with logarithmic Y-axis
    histogram(channelData);
    % set(gca,'YScale','log');
    title(sprintf('Channel #%d\nMean: %s\nStd: %s', i, channelMeans(i), channelSDs(i)));
end

%%
% 
% set(gca,'Yscale','log');
% han=axes(fig,'visible','off'); 
% 
% sgtitle('Histograms of Data for Each Channel');
% han.XLabel.Visible='on'
% han.YLabel.Visible='on'
% xlabel(han, 'Volts (V)');
% ylabel(han, 'Logarithmic  Frequency');

%% 

% % 3. Choose the threshold . 
% confidence_interval = [0.9990, 0.99994]; % Confidence interval range
% 
% % negative thresholds for each channel
% negative_thresholds = zeros(16, 1); %all channels negative threshold

threshold_value = 3; % how many STDs away

% Loop through 16 different channels for the negative threshold
for i = 1:16 %channels 1-16
    % Calculate the negative threshold for the current channel using
    % negative threshold func 
    negative_thresholds(i) = channelMeans(i) - threshold_value*channelSDs(i);
end

% Display or use the calculated negative thresholds as needed
disp('Negative Thresholds for Each Channel:');
disp(negative_thresholds);

%% 
%4/5: create index and find the time point and then make element in the
%index
% i need to store this information somewhereand find the indices of threshold crossings for each channel
%trying to make a cell array to store info

threshold_crossing_indices = cell(16, 1); %channels one to sixteen
crossed_time_points = cell(16,1)

% fourloop through each channel
for i = 1:16
    channelData = filtered_spk(i,:); % I want just data for the current channel from the larger raw dataset
    
    % use hte find function for the 16 channel to see where things were
    % negative 
    threshold_crossing_index = find(channelData < negative_thresholds(i));
    
    % Store the index in the cell array
    threshold_crossing_indices{i} = threshold_crossing_index;
    
    % crossed_time_points{i} = min(threshold_crossing_index);
end


%%

for i=1:16
    for j = 1:length(threshold_crossing_indices{i,1})
        if thresholcrossingindices{i,j+1}-thresholcrossingindices{i,j} >= 36 % 1.5ms window in terms of samples
            
            for n = 1:16
                    DATAMATRIX(n,:) = thresholcrossingindices {n, 1};
             end
              
                DATAMATRIX = double(DATAMATRIX);
               cross{} = threshold index crossing save to new data structure
        else
            continue;
        end
    end
end


%% as far as we got w Ben
%%

% 6. ignore threshold and adjust so we get just our appropriate snippet

% Calculate the Nyquist frequency: assuming the frequency needed is much
% higher than the max firing rate.
% 
% nyquist_frequency = max_firing_rate * 4 ;
% sampling_rate = 2 * nyquist_frequency;

% Calculate the minimum number of samples required to capture one cycle of the max frequency
%min_samples_per_cycle = round(sampling_rate / max_firing_rate);

% Define the size of the snippet window (you can adjust this as needed)
snippet_window_size = min_samples_per_cycle * 2; % the two woudl capture two cycles??



%7 create cell array and use index from 6 to get out the spike filter data
% 
% % spike_filtered_data = cell(16, 1); %filter for channel 1-16
% for i = 1:16
% channelData = spike_filtered_data{i};
% spike_filtered_data{i} = channelData(crossed_time_points{i} - snippet_window_size/2:crossed_time_points{i} + snippet_window_size/2);
% spike_filteredd_data{i} = channelData[spike_fliltered_data{i}]
% % end
% 

%8
minimums= [];
minimum_positions = [];
for i =1:16
    index_array = crossed_time_points{i} - snippet_window_size/2:crossed_time_points{i} + snippet_window_size/2;
    minimums = [minimums min(spike_filtered_data{i})];
    minimum_positions = [minimum_positions index_array(spike_filtered_data{i} ==min(spike_filtered_data{i}))];
    
end

%9
spike_filtered_data = cell(16, 1);
for i = 1:16
    channelData = spike_filtered_data{i}{i};
    spike_filtered_data{i} = channelData(minimum_positions(i) - snippet_window_size/2:minimum_positions(i) + snippet_window_size/2);
%     spike_filtered_data{i} = channelData[spike_filtered_data{i}]
end


%10

% Define the window size components for extracting data around the minimum value
windowComponentBefore = snippet_window_size/2;  

% Initialize an array to store the timestamps
timestamps = zeros(16, 1);

% Convert minimum position to timestamp for each spike
for i = 1:16
    channelData = Wideband_data{i};
    
    % Adjust for the 'window component before'
    adjustedPosition = crossed_time_points{i} - windowComponentBefore;
    
    % Ensure that the adjusted position is within bounds
    adjustedPosition = max(adjustedPosition, 1);
    
    % Grab the voltage value at the adjusted position
    voltageAtAdjustedPosition = channelData(adjustedPosition:minimum_positions(i)+1);

     if min(voltageAtAdjustedPosition) == min(spike_flitered_data{i})
         % Convert the adjusted position to time using the sampling rate
         timestamps(i) = adjustedPosition / samprate;
     else
         continue;
  %you can put something in the "else" or just enter "continue"
      
     end
end


%11 and 12

windowComponentAfter = snippet_window_size/2;
uniqueTimestamps = zeros(16, 1);
snippetMatrix = zeros(16, windowComponentBefore + windowComponentAfter + 1);

% Initialize counter variable
% counter = 1;

% Loop through each spike
for i = 1:16
    % Check if the current timestamp is different from the previous one
    if timestamps(i) ~= uniqueTimestamps(counter)
        % If different, store the timestamp and snippet in the matrix
        uniqueTimestamps(counter) = timestamps(i);
        snippetMatrix(counter, :) = spike_flitered_data{i};
        Matrixvariable('counter',:)
        
        % Increment the counter
        counter = counter + 1;
    end
end


% 14 Create pile plots for each channel
% %% STEP 14
% %%pile plots could be made using the stair function and hold on functions
% %%to overlay
% 
% % data for  sequences
% seq1 = ;
% seq2 = ;
% seq3 = ;
% 
% % Create a time vector (x-axis values)
% time = timestamps(i);
% 
% % Create a pile plot using the stairs function
% stairs(time, seq1, 'b'); % Blue
% hold on; % This allows you to overlay the other sequences
% stairs(time, seq2, 'g'); % Green
% stairs(time, seq3, 'r'); % Red
% 
% % Customize plot properties (labels, legend, etc.) as needed
% xlabel('Time (mS)');
% ylabel('Volts (mV)');
% title('Pile Plot for all 16 Channels');
% legend('Channel 1', 'Channel 2', 'Channel 3');
% 
% % Hold off to stop overlaying additional plots
% hold off;



% Loop through each channel and create and save separate pile plots
%%the general code for creating an "invisible figure" 
% % Create an invisible figure
% figure('visible','off');
% 
% % Plot some data on the invisible figure
% plot(x, y);
% 
% % Save the figure as an image without displaying it
% saveas(gcf, 'my_plot.png');
% 
% % Close the invisible figure
% close(gcf);

for i = 1:16
    % Create a new figure for each channel
    figure('visible','off');
    plot(snippetMatrix, (i, :));
    title(sprintf('Spike Channel %d', i));
    
    % Define the filename for saving the plot
    filename = sprintf('part14Channelpileplot.png', i);
    savePath = fullfile(saveFolder, filename);
    saveas(gcf, savePath);
    close(gcf);
end

% %% this is other work I tried -- ignore for now 
% %Spike plots
% sort(Wideband_data)
% mean_values = mean(Wideband_data,0,1 );
% Std_values = std(Wideband_data,0,1 );
% 
% % find the mean, standard deviation of each channel (16 times)
% %%
% %channel 1
% % mean_value1 = mean (Wideband_data, 1);
% % std1 = std(mean1)
% % mean_value1 = mean (channel,1,1)
% % std1 = std(mean channel1)
% % %%
%     % Plot each channel plot a histogram
%     for channel = 1:num_channels
%     % Get data for the current channel
%         channel_data = Wideband_data(channel, :);
% 
% for i = 1:16
% MEAN(i) = mean(filtered_spike(i,:));
% %gives mean of all 16 channels after filter for the channels and : means
% %all columns 
% SD(i) = std(filtered_spike (i, :));
% %try plotting the y axis in log 
% 
% %alternatively I could use set(gca, 'YScale', 'log')
% subplot (16,4,1)
% semilogy(x,y);
% %find filtered data less than negative threshold. We want the most neg
% %under threshold (look at negative threshold) 
% title(sprintf('Channel %d\nMean: %.2f\nStd: %.2f', i, mean_value, std_value));
% xlabel('x')
% ylabel('y')
% sgtitle('Histograms of Channel Data')
% %%
% %channel 16
% %histogram  x vs y is voltage vs count
% %should appear as a bell shaped curve 
% %chosen standard deviation 
% %% Creating Output 1: Histogram for 16 channels
% %frequncy range, filter the hertz 
% %figure out the chosen range 
% 
% % Plot a histogram of the data for each channel
% for channel = 1:size(Wideband_data, 1)
%     channel_data = Wideband_data(channel, :);
% end
% 
% 
%     % Plot histogram and add channel information to the title
% % According to HW3 video, you can use 'hist' or 'histogram' functions. The 
% %help hist function says "hist is not recommended. use histogram instead.
% %this helped me choose histogram. 
% %https://www.youtube.com/watch?v=4jjZ-G-eN_4
% histogram 
% 
% %chosen std threshold = 2.8 or something then -threshold = mean-std threshold selected *std of channel) mean (channel) -(*std(channel_data)
% % Adjust Y-axis to log scale if needed
% % Save/print the figure
% % Define the number of channels
% num_channels = size(Wideband_data, 1);
% 
% % Create a figure for the histograms
% figure;
% 
% 
% 
%     % Create a subplot for the current channel
%     subplot(4, 4, channel_1); % Change the layout as needed
% 
%     % Plot the histogram
%     histogram(channel_data, 'Normalization', 'probability');
% 
%     % Customize the plot
%     title(['Channel 1 '); % Add channel information to the title
%     xlabel('x label'); % Label the x-axis
%     ylabel('y label'); % Label the y-axis
% 
%     % Adjust histogram bin width or other thingss i need to 
% end
% 
% % Add a title to the entire figure
% suptitle('Histograms for 16 Channels');
% 
% % Adjust the layout of subplots for better visibility (optional)
% % Adjust margins and spacing as needed
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
% %% Indicies
% %
% 
% %%
% % Snippets: 
% %get small window so that only individual waveform is shown 
% %decide window 
% %take samples before a certain number (can it be arbitrary?) 
%     %how much waveform do you want to show for spike sorting? 
%     %sampling rate, time chosen, x many samples 
%     %samples on left then samples on the right 
% 
% %%group steps together 
% %%step 11:
% %after shifting to 0 point, check that we are not shifting window beyodn
% %the length of the data
% %Define LFP and Spike cutoff frequencies and filter orders, take from
% %assignment 2
% lfp_low_cutoff = 300; % Hz
% lfp_high_cutoff = 3000; % Hz
% lfp_filter_order = 8; %CHOSEN 10TH ORDER
% 
% spike_low_cutoff = 100; % Hz
% spike_high_cutoff = 5000; % Hz
% 
% spike_filter_order = 8; %CHOSEN SAME 10TH ORDER
% 
% %%
% %find negative threshold for every single channel 
% %find filterdata less than negative threshold 
% 
% 
% %25:52
% for i = 1:16
%     for j= 1:30000000
%     interval = channel(i,j+i)) = index(channel (i,j)); +
% 
% 
%     end
% 
