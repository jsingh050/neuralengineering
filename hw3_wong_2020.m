%% Ben Wong
% BioE 2615: Neural Engineering
% Assignment 3: Spike Sorting
% 9/9/20

%% Exercise 1

LFP_data = load('1000802VisuallyEvoked.mat'); % load my assigned data file
samprate = LFP_data.samprate; % extract the sampling rate

% for loop to extract the channel data into an accessible matrix
for i = 1:16
    channel{i} = double(LFP_data.Wideband_data{i}); % convert to double
    continue;
end

% cutoff freqs for Spike chosen 150Hz and 8500Hz
Wn_Spike = [150/(samprate/2),8500/(samprate/2)]; % normalized cutoff freqs in Wn_LFP vector
[b_spk,a_spk] = butter(1, Wn_Spike); % 1st order

% for loop to filter each channel's spike data and put into an accessible matrix
for k = 1:16 
    fil_spk{k} = filter(b_spk,a_spk,channel{k});
    continue;
end

% for loop to find mean and std of each channel and insert into matrices
for m = 1:16
    MN(m) = mean(fil_spk{m}); % mean
    STD(m) = std(fil_spk{m}); % standard deviation
    continue
end

%% Exercise 2

% for loop for each channel's histogram 
for n = 1:16
    subplot(4,4,n); % subplot each channel
    histogram(fil_spk{n}); % histogram
    set(gca,'Fontsize',6);
    title(sprintf('Ch %d: Mn=%dV Std=%dV',n,MN(n),STD(n)));
    xlabel('Relative Voltage (V)');
    ylabel('Count');
    continue
end

%% Exercise 3

% for loop for neg threshold of each channel
% detection threshold = 4 stds for 99.9%
for p = 1:16
    neg_thresh(p) = MN(p)-4*STD(p);
    continue
end

%% Exercise 4

% for loop to index where the spike crosses the threshold
for q = 1:16
    idx{q} = find(fil_spk{q} <= neg_thresh(q));
end

%% Exercise 5 & 6

% checking the APs are truly separate
for r = 1:16
  cross{r}(1) = idx{r}(1); % defining first index as 1st AP
  s = 2; % for count
  for u = 2:length(idx{r}) % starting after 1st AP
      if idx{r}(u)-idx{r}(u-1) > .001*samprate % greater than 1 ms time
          cross{r}(s) = idx{r}(u);
          s = s+1; % increase count
      else s = s;
      end
  end
end 

% window is symmetrical
window = .0015; % window size determined to be 1.5 ms
wdw_left = round((window/2)*samprate); % integer number of samples on left
wdw_right = round((window/2)*samprate); % integer number of samples on left

% indexing the AP crossing
for v = 1:16
    AP_av{v} = find(cross{v}-wdw_left >0 & cross{v}+wdw_right <= numel(fil_spk{v})); % checking all points if available
    for t = 1:length(AP_av{v})
        AP_idx{v}(t) = cross{v}(AP_av{v}(t)); %writing it to return index that cross point that would fit the available 
    end
end

%% Exercise 7

% for extracting APs within each window and index
for w = 1:16
    for x = 1:length(AP_idx{w}) % all of AP index
        AP{w,x} = fil_spk{w}(AP_idx{w}(x) - wdw_left:AP_idx{w}(x) + wdw_right);
    end
end

%% Exercise 8

% for loop for the min and min index
for y = 1:16
    for z = 1:length(AP{y})
        [minval{y,z} minidx{y,z}] = min(AP{y,z}); % return both min and index
        dif{y,z} = abs(minidx{y,z}-wdw_left-1); % measure the difference from min and center point
    end
end

%% Exercise 9

% loop to re-adjust the window around min
for aa = 1:16
    for bb = 1:length(AP{aa}) 
        shift{aa,bb} = AP_idx{aa}(bb)+dif{aa,bb}; % the shift total index and difference
        spk_min{aa,bb} = fil_spk{aa}(shift{aa,bb}-wdw_left:shift{aa,bb}+wdw_right); % shift around min on either window side
    end
end

%% Exercise 10-14

% time vector throughout the window
time = 0:1/samprate:.0015; % interval based on samprate

figure;
for dd = 1:16
    subplot(4,4,dd); % subplot for each channel
    for ee = 1:37
        plot(time,spk_min{dd,ee})
        title(sprintf('Channel %d',dd));
        set(gca,'XTickLabel',[0 .5 1 1.5]);
        xlabel('Time (ms)');
        ylabel('Amplitude (V)');
        hold on;
    end
end
