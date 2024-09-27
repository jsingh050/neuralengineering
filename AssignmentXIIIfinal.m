
%Jaspreet Singh
%November 21st, 2023
%%Assignment XIII
%BIOENG 2615

%%
% Notes: 
% red is blood vessels 
% and green is cCAMP (floresence) 
% - you want to observe the frame in the pdf of channel two nad the differnees in te channels so you want four things in one figure 

%% QUESTION 1
%take the literal differnce btween the frames 
% Load the data
load('tmp_tser8663_raw.mat');

% Display images
figure;
for mm = 1:size(data, 4)
    imagesc(squeeze(data(:,:,2,mm)), [0, 1]);
    title(['Motion Corrected Frame #', num2str(mm)]);
    t = 0 : 1/5 : size(data, 4)/5 - (1/5);
    xlabel(['Time= ', num2str(t(mm)), ' sec, Frame #', num2str(mm)]);
    axis image, colormap gray, colorbar;
    drawnow;
end

% STEP 1
hold on;
figure(1), clf,
subplot(1, 2, 1),
colormap gray,
imagesc(data(:,:,2,30)), axis image, colormap gray, colorbar,

title('Channel 2, Frame #30');
xlabel('Time (s), Frame #30');
ylabel('ylabel')


colormap gray,
subplot(1, 2, 2),
imagesc(data(:,:,2,30) - data(:,:,2,680)), axis image, colormap gray, colorbar,
title('Channel 2, Difference between Frame #680 and #30');
xlabel('Time (s), Difference Image')
hold off;

% Repeat for Channel 1
hold on

figure(2), clf
colormap gray,
subplot(1, 2, 1),
imagesc(data(:,:,1,30)), axis image, colormap gray, colorbar,
title('Channel 1, Frame #30');
xlabel('Time (s), Frame #30');

subplot(1, 2, 2),
colormap gray,
imagesc(data(:,:,1,30) - data(:,:,1,680)), axis image, colormap gray, colorbar,
title('Channel 1, Difference between Frame #680 and #30');
xlabel('Time (s), Difference Image')

hold off;

%%
%working in higher dimensions you can cpntorl wher ethe min and max of the
%colorbar adjusts how matlab sets the range 
% STEP 2
% Calculate average GCaMP image in channel 2
avgGCaMP = mean(data(:,:,2,:), 4);
figure(5), clf,
imagesc(avgGCaMP), axis image, colormap gray, colorbar,
title('Average GCaMP Image (Channel 2)')
xlabel('Vertical Posistion (microns)')
ylabel('Horizontal Posistion (microns)')

%% Didn't work for step 3 
% %% Step 3: Motion correction using imMotionCorrect
% [data_corrected, motion_est] = imMotionCorrect(data(:,:,2,:), 'reference', 1);
% figure(6), clf,
% plot(motion), xlabel('Frame'), ylabel('Motion (pixels)'),
% title('Estimated Motion Parameters')
% 
% % Extract green fluorescence data (channel 2)
% green_data = data(:,:,2,:);
% 
% % Perform motion correction using imMotionCorrect
% [data_corrected, motion_est] = imMotionCorrect(green_data, 'reference', 1);
% 
% % Display the estimated motion parameters
% figure,
% plot(motion_est)
% xlabel('Frame')
% ylabel('Motion (pixels)')
% title('Estimated Motion Parameters')
% 
% % Display the first two images side by side for channel 2 after motion correction
% figure,
% subplot(1, 2, 1),
% imagesc(mean(green_data, 3)), axis image, colormap gray, colorbar,
% title('Average GCaMP Image Before Motion Correction')
% 
% subplot(1, 2, 2),
% imagesc(mean(data_corrected, 3)), axis image, colormap gray, colorbar,
% title('Average GCaMP Image After Motion Correction')
%% STEP 3

%input your gcamp channel two data and then output the motion 
% Extract the green fluorescence data from channel 2
ch2_data = data(:,:,2,:);
ch2_data = squeeze(double(ch2_data)); % Convert double to perform operations


%reference frame (e.g., frame #1) for motion correction
reference_frame = 1;
reference_image = squeeze(ch2_data(:,:,reference_frame));

% Motion correction parameters
motion_correction_params = [0.01 0];

% Perform motion correction
% [data_corrected, motion_est] = imMotionCorrect(ch2_data, reference_image, motion_correction_params);
%[data_corrected, motion_est] = imMotionCorrect(ch2_data, reference_image, motion_correction_params(1), motion_correction_params(2));
[motion_est, data_corrected] = imMotionCorrect(ch2_data, reference_image, motion_correction_params);


% Images of estimated motion parameters
figure; hold on;
plot(motion_est(:, 1));  % Adjust this but it turns the 3D array to plot just 1D to correct error
plot(motion_est(:, 2));
legend("X Motion", "Y Motion");
ylabel('Motion (pixels)')
title('Estimated Motion Parameters')

% Images after the  motion correction
figure,
subplot(1, 2, 1),
imagesc(mean(ch2_data, 3)), axis image, colormap gray, colorbar,
title('Average GCaMP Image Before Motion Correction')

subplot(1, 2, 2),
imagesc(mean(data_corrected, 3)), axis image, colormap gray, colorbar,
title('Average GCaMP Image After Motion Correction')

%% STEP 4

% % Display estimated motion
% figure;
% plot(motion_est);
% title('Estimated Motion for each Frame');
% xlabel('Frame');
% ylabel('Motion (pixels)');

% Display motion-corrected images
figure;
for mm = 1:size(data_corrected, 3)
    imagesc(squeeze(data_corrected(:,:,mm)), [0, 1]);
    title(['Motion Corrected Frame #', num2str(mm)]);
    t = 0 : 1/5 : size(data_corrected, 3)/5 - (1/5);
    xlabel(['Time= ', num2str(t(mm)), ' sec, Frame #', num2str(mm)]);
    axis image, colormap gray, colorbar;
    drawnow;
end

% Calculate average GCaMP image using the motion-corrected data
avgGCaMP_corrected = mean(data_corrected, 3);

% Compare with the image from Question 2
figure;
subplot(1,2,1);
imagesc(avgGCaMP);
axis image, colormap gray, colorbar;
title('Average GCaMP Image (Original)');
xlabel('x')
ylabel('y')

subplot(1,2,2);
imagesc(avgGCaMP_corrected);
axis image, colormap gray, colorbar;
title('Average GCaMP Image (Motion-Corrected)');
xlabel('x')
ylabel('y')

%% STEP 4 (the correted image looks all wonky) (!!!!)

% % Compare the original and then the images that were corrected 
% figure;
% subplot(1,2,1);
% imshow(squeeze(ch2_data(:,:,reference_frame)), []);
% title('Original  Image');
% subplot(1,2,2);
% imshow(squeeze(data_corrected(:,:,reference_frame)), []);
% title('Corrected Image');
% sgtitle('Comparison between Original and Corrected Images');
% 
% % Calculate the average GCaMP image using the corrected data
% average_image_corrected = mean(data_corrected, 4);
% figure;
% imshow(average_image_corrected, []);
% title('Average GCaMP Image (Motion Corrected)');
% xlabel('x')
% ylabel('y')


%% STEP 5

% Using mask C
figure;
imagesc(maskC);
axis image;
colormap jet;
colorbar;
title('Cell Mask');

% Assuming t is defined as the time axis
% Extract time series for each segmented cell in maskC
num_cells = max(maskC(:));
data_C = zeros(size(data_corrected, 3), num_cells);

for mm = 1:size(data_corrected, 3)
    tmp_im = data_corrected(:,:,mm);
    for nn = 1:num_cells
        data_C(mm, nn) = mean(tmp_im(maskC == nn));
    end
end

% Plot time series from some cells
selected_cells = [1, 2, 3];  % Adjust as needed
figure;
plot(t, data_C(:, selected_cells));
xlabel('Time (seconds)');
ylabel('Fluorescence Intensity');
title('Time Series of Selected Cells');
legend('Cell 1', 'Cell 2', 'Cell 3');  % Adjust legend labels
axis tight;

%%
% %% STEP 6: MASK (redefine threshold or just use the dark spaces, maybe using greythresh and redefine 
% 
% % Calculate the average blood vessel image (channel 1 data)
% vessel_image = mean(data(:,:,1,:), 4);
% 
% % Generate a blood vessel mask by applying a threshold
% thr = .2   %(??)
% vessel_mask = vessel_image > thr;
% 
% % Generate a cell body mask using the variable maskC
% cell_mask = (maskC > 0);
% 
% % Generate the neuropil_mask using both the vessel_mask and cell_mask
% neuropil_mask = (~cell_mask) & (~vessel_mask);
% 
% % Display the neuropil mask
% figure;
% imshow(neuropil_mask, []);
% title('Neuropil Mask');
% 
% % Extract the neuropil time series
% gcamp_image = mean(data(:,:,2,:), 4);
% neuropil = zeros(size(data_corrected, 4), 1);
% 
% for mm = 1:size(data_corrected, 4)
%     tmp_im = data_corrected(:,:,mm);
%     neuropil(mm) = mean(tmp_im(neuropil_mask));
% end
% 
% % Plot the neuropil time series
% figure;
% plot(t, neuropil);
% xlabel('Time (seconds)');
% ylabel('Fluorescence Intensity');
% title('Neuropil Time Series');

%% STEP 6
% Generate a mask for a region-of-interest that contains the "neuropil"

% Calculate the average blood vessel image (channel 1 data)
vessel_image = mean(data(:,:,1,:), 4);

% Generate a blood vessel mask by applying a threshold to this image
thr = 0.08; % Adjust the threshold as needed
vessel_mask = vessel_image > thr;

% Generate a cell body mask using the variable maskC
cell_mask = (maskC > 0);

% Generate the neuropil_mask using both the vessel_mask and cell_mask
neuropil_mask = (~cell_mask) & (~vessel_mask);

% Generate a figure of your neuropil mask
figure;
subplot(1,3,1);
imagesc(vessel_mask);
axis image; colormap gray;
title('Vessel Mask');
subplot(1,3,2);
imagesc(cell_mask);
axis image; colormap gray;
title('Cell Mask');
subplot(1,3,3);
imagesc(neuropil_mask);
axis image; colormap gray;
title('Neuropil Mask');

% Extract the neuropil time series
gcamp_image = mean(data(:,:,2,:), 4);
for mm = 1:size(data_corrected, 3)
    tmpim = data_corrected(:,:,mm);
    neuropil(mm) = mean(tmpim(neuropil_mask));
end
figure;
plot(t, neuropil);
xlabel('Time (seconds)');
ylabel('Fluorescence Intensity');
title('Time Series of Neuropil');
axis tight;

%% STEP 7: K MEANS

% % Extract green fluorescence data (channel 2)
% ch2_data = data(:,:,2,:);
% 
% % Reshape the data for k-means clustering
% reshaped_data = reshape(ch2_data, [], size(ch2_data, 4));

% Number of clusters (this needs to be adjusted i think)
num_clusters = 6;

% Perform k-means clustering
kmeans_idx = kmeans(data_C.', num_clusters, 'distance', 'correlation');

kim = zeros(size(data,1), size(data,2)); 
for m = 1:max(maskC(:))
    kim(maskC==m) = kmeans_idx(m); 
end
figure;
imagesc(kim), axis image, colormap jet, colorbar;

tck = zeros(numel(t), num_clusters);
for m = 1:max(kmeans_idx(:))
    tck(:,m) = mean( data_C(:, kmeans_idx==m), 2); 
end
figure;
for m = 1:num_clusters
    subplot(3,3,m)
    plot( t, tck(:,m) )
    axis tight, grid on, drawnow, 
end

figure;
imagesc(corrcoef(tck)), axis image, colormap hot, colorbar;

% % Reshape the clustering results back to the original image size
% cluster_image = reshape(kmeans_idx, size(ch2_data, 1), size(ch2_data, 2));
% 
% % Display the results
% figure, clf,
% imagesc(cluster_image), axis image, colormap jet, colorbar,
% title('K-means Clustering on Original Data')

%% STEP 9: 
% Perform PCA
datasz = size(data_corrected);
newdata = reshape(data_corrected, [datasz(1)*datasz(2), datasz(3)]);
avgdata = mean(newdata, 2);
newdata = newdata - repmat(avgdata, [1, datasz(3)]);
[U, V, W] = svd(newdata', 'econ');
%%
% Plot mean intensity and the first 12 vectors of U
figure;
subplot(2,1,1);
plot(t, mean(newdata, 1));
title('Mean Intensity Over Time');
xlabel('Time (s)');
ylabel('Mean Intensity');
axis tight;

subplot(2,1,2);
for mm = 1:12
    plot(t, U(:,mm));
    hold on;
end
title('First 12 Vectors of U');
xlabel('Time (s)');
ylabel('Intensity');
axis tight;

%% STEP 11
% Display coefficient images for the first 12 components
figure;
for mm = 1:12
    subplot(3,4,mm);
    imagesc(reshape(W(:,mm), [datasz(1), datasz(2)]));
    axis image, colormap jet, colorbar;
end
%% Extra Credit: Step 14
% Assuming data is your imaging data variable
% Assuming channel 1 corresponds to blood vessels
vessel_data = squeeze(data(:,:,1,:));

% Calculate the average signal over time
average_vessel_signal = squeeze(mean(vessel_data, [1, 2]));

% Plot the intensity over time
figure;
plot(t, average_vessel_signal);
xlabel('Time (seconds)');
ylabel('Average Blood Vessel Intensity');
title('Average Blood Vessel Signal Over Time');

%% STEP 15
% Assuming average_vessel_signal and data_C are your average blood vessel signal and neuronal time series, respectively
correlation_coefficient = corr(average_vessel_signal, data_C);

% Display the correlation coefficient
disp(['Correlation Coefficient: ', num2str(correlation_coefficient)]);


