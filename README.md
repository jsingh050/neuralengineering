# neuralengineering

This repository is a collection of work I did for a neural engineering course during Fall 2023. I have only posted this here for reference to the type of neural analysis projects I have done before. I appreciate discretion. 

This homework introduces fundamental signal processing and image analysis techniques in MATLAB. The assignment begins by generating two sine waveforms: a basic sine wave A = sin(2π15t) and a distorted version B = sin(2π15t)^5. These waveforms are stored in a matrix along with their pointwise sum, product, and ratio. The resulting operations are visualized through subplots to highlight their differences in amplitude and shape. Next, three sine waves are plotted using different x-intervals (π spaced by 0.5, 0.1, and 0.01) to show how sampling resolution affects smoothness in waveform representation. Using hold on, the waves are overlaid with unique styles and then labeled appropriately.

In part (c), an image file (tire.mat) is loaded and visualized using imagesc, which displays intensity values using color mapping. A specific 3x2 sub-matrix is extracted and displayed to practice matrix indexing. In part (d), square waves are synthesized by thresholding sine waves: values above or below zero are replaced with ±1. This is based on Fourier theory, which shows that square waves can be approximated using sine components. Three plots compare 20 Hz and 30 Hz sine waves and the derived square wave.

Finally, part (e) revisits the tire image. The mean pixel intensity is computed, and a binary mask is created to highlight pixels above this mean. This mask is applied to the image to isolate and visualize regions of high intensity. Although the final section contains some syntax issues (like quoting variable names or misusing strings), the intention is to display both the original and masked images side-by-side for comparison, reinforcing the concept of thresholding in image processing.
Neural Signal Processing - Folder 3
This repository contains MATLAB code for neural spike detection and analysis based on wideband neural recordings from 16 channels. The script performs the following tasks:
 - Loads neural signal data from 0900702VisuallyEvoked.mat.
- Applies a bandpass filter (300–3000 Hz) to isolate spike activity.
        - bandpass filter removes Low-Frequency Artifacts and High-frequency components. The upper limit of 3000 Hz is chosen because most neural spikes occur between 300–3000 Hz.
  
Neuronal action potentials last about 1–2 milliseconds, and their energy is concentrated within this range.
Using a Butterworth filter, the bandpass ensures smooth but sharp filtering of spikes without distortion.


Filtering ensures we only analyze relevant spike frequencies.
Histograms help us set appropriate thresholds and understand data distributions.
- Computes mean and standard deviation for each channel to set spike detection thresholds.
- Detects threshold crossings, identifying spike events based on a 3-sigma threshold.
- Extracts spike snippets around detected events for further analysis.
- Generates histograms to visualize signal distributions across channels.
- Histograms help in visualizing the voltage distribution of the neural signals. Neural data contains both signal (spikes) and background noise. histograms, we see the distribution of voltage values, helping us separate spikes from noise. The threshold is typically set as mean ± (N × standard deviation) of the signal.
The histogram provides a visual guide to selecting this threshold optimally.
- Computes timestamps for detected spikes and removes redundant detections.
- Plots and saves pile plots of detected spikes.
The script is useful for spike sorting and neural signal processing.

Folder 5 
In this assignment, we analyze multi-channel electrophysiological data collected during visually evoked stimulation, focusing on spike detection, snippet extraction, and component-based classification. The script begins by loading raw neural data (Wideband_data) and plotting the unfiltered signal from one channel to inspect baseline activity. A bandpass Butterworth filter (300–3000 Hz) is then applied to isolate spike-relevant frequencies, removing low-frequency LFPs and high-frequency noise. The signal is segmented into stimulus-aligned windows using timestamps from TrigON, and the filtered data is stored for further analysis.

To detect spikes, the code defines a negative threshold based on the channel’s mean and standard deviation, marking voltage deflections that exceed this threshold. A binary raster matrix is constructed for each stimulation event, recording spike occurrences. From these rasters, snippets—brief segments of signal surrounding detected spikes—are extracted, focusing on channel 7. These snippets are then plotted to visualize waveform consistency.

Subsequent analysis includes plotting a single spike snippet with a real-time x-axis (in milliseconds), and calculating minimum vs. maximum voltage values across all snippets to assess their shape characteristics. To reduce dimensionality and identify dominant waveform features, Principal Component Analysis (PCA) is performed via SVD. The explained variance plot helps determine how many components account for most of the variability. The first 25 PCA components are visualized, and the most informative 3 are selected for further clustering analysis.

These components’ weights (W) are used to generate histograms and scatter plots, revealing clusters of spike waveforms, which may correspond to different neuronal units. Based on the third component, snippets are split into two clusters (unit 1 and unit 2), and their mean waveforms are compared. A raster plot visualizes spike activity for the full channel and each unit individually. To investigate temporal dynamics, the script calculates inter-spike intervals and generates autocorrelograms for each cluster, revealing refractory periods and firing patterns.

Finally, k-means clustering is applied to the PCA component space to validate manual separation, and the resulting clusters are overlaid with waveform plots. This pipeline demonstrates essential spike sorting techniques—filtering, thresholding, PCA, clustering, and autocorrelation—providing a comprehensive overview of how to identify and characterize single-unit activity from multi-channel recordings in response to sensory stimuli.
Neural Signal Processing - Folder 6
This fold has MATLAB code for analyzing visually evoked neural signals. 
Load neural data from 0900702VisuallyEvoked.mat.
Filter neural signals using a low-pass Butterworth filter to remove high-frequency noise.
Align neural signals with stimulus onset (TrigON) for event-related analysis.
Compute power spectral density (PSD) to examine the frequency response of neural signals.
Generate spectrograms to analyze how neural responses evolve over time.
Estimate Current Source Density (CSD) to determine where neural currents originate.
folder 7 
This script performs multistep image processing on high-bit-depth fluorescence microscopy images, particularly analyzing GFAP-labeled brain sections (IMPLANT_C004_GFAP.tif). The analysis begins by loading the original image and visualizing it with appropriate scaling to reflect its 12-bit depth. Spatial calibration is performed using a known scale of 0.621 microns per pixel, and the image is binned spatially into 2-micron bins to reduce resolution and allow easier quantification. This is done by averaging the intensity values within each bin, producing a lower-resolution but spatially meaningful image.

The script then thresholds the image based on intensity: pixels below a specified number of standard deviations below the mean are flagged, generating a binary mask that identifies regions with lower-than-average fluorescence (e.g., signal loss, background, or probe tract). A user-defined probe center is selected using ginput, from which a probe mask is built, followed by creating vertical bins along the tract. These bins are used to segment the image and calculate mean and standard deviation of intensity per bin, normalized to the threshold value. This provides a spatially resolved quantification of signal decay or variation along the depth of the probe.

In the next section, the script performs Singular Value Decomposition (SVD) across a full image stack (TSeries-10052015-1158-2053.tif) containing 1000 frames of 256x256 pixels. The 3D data is reshaped into a 2D matrix and decomposed using SVD. The singular values are plotted to visualize how much variance each component explains, and the image is reconstructed using a specified number of principal components (e.g., 65) to denoise the data while preserving spatial and temporal structure. The original and reconstructed frames are shown side-by-side to assess the effectiveness of dimensionality reduction.

Finally, an extra credit section introduces a cell-counting framework by defining a probe region and dividing it into multiple bins. It overlays bin boundaries on the image to visualize how cells or intensities could be quantified relative to the probe tract. This lays the foundation for future steps like automatic cell segmentation, quantification, or statistical comparison across treatment groups.


Neural Signal Processing - Folder 8 
This assignment involves analyzing time-series calcium imaging data to study neuronal and vascular activity using fluorescence signals. The data contains two channels: Channel 1 (red) captures blood vessel fluorescence, and Channel 2 (green) captures GCaMP expression, a proxy for intracellular calcium changes in neurons. The workflow begins by visualizing individual image frames to inspect differences over time and between channels. Specific frames (e.g., frame #30 and #680) are subtracted to observe localized activity changes, highlighting dynamic regions of interest. An average image of GCaMP fluorescence is computed to represent baseline neural activity. To correct for motion artifacts, a motion correction algorithm aligns frames to a reference, and motion parameters (X and Y displacement) are plotted. Before and after correction, average images are compared to assess the correction's effectiveness.

Next, time-series signals are extracted for each segmented cell region (using maskC) to analyze neural activity across time, and these traces are plotted for selected cells. In parallel, a neuropil mask is created by excluding both cell bodies and blood vessels, enabling the extraction of background calcium signals. These neuropil signals help distinguish true cellular activity from background fluctuations. Following that, k-means clustering is performed on the cell-wise activity data to uncover patterns of coordinated activity, and cluster identities are mapped back to spatial locations. Cluster-averaged signals are plotted to visualize distinct activity profiles, and their correlations are assessed.

The assignment continues with Principal Component Analysis (PCA), where the motion-corrected data is reshaped, mean-centered, and decomposed using singular value decomposition. This reduces dimensionality and highlights dominant spatiotemporal patterns. The mean intensity over time and the first 12 temporal components (U vectors) are plotted, followed by spatial coefficient maps (W vectors) for these components. Finally, the average blood vessel signal over time is plotted and correlated with neuronal time-series data to investigate neurovascular coupling, quantified using Pearson's correlation coefficient.

