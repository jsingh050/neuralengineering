# neuralengineering

This repository is a collection of work I did for a neural engineering course during Fall 2023. I have only posted this here for reference to the type of neural analysis projects I have done before. I appreciate discretion. 

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

Neural Signal Processing - Folder 6
This fold has MATLAB code for analyzing visually evoked neural signals. 
Load neural data from 0900702VisuallyEvoked.mat.
Filter neural signals using a low-pass Butterworth filter to remove high-frequency noise.
Align neural signals with stimulus onset (TrigON) for event-related analysis.
Compute power spectral density (PSD) to examine the frequency response of neural signals.
Generate spectrograms to analyze how neural responses evolve over time.
Estimate Current Source Density (CSD) to determine where neural currents originate.

Neural Signal Processing - Folder 8 
This assignment involves analyzing time-series calcium imaging data to study neuronal and vascular activity using fluorescence signals. The data contains two channels: Channel 1 (red) captures blood vessel fluorescence, and Channel 2 (green) captures GCaMP expression, a proxy for intracellular calcium changes in neurons. The workflow begins by visualizing individual image frames to inspect differences over time and between channels. Specific frames (e.g., frame #30 and #680) are subtracted to observe localized activity changes, highlighting dynamic regions of interest. An average image of GCaMP fluorescence is computed to represent baseline neural activity. To correct for motion artifacts, a motion correction algorithm aligns frames to a reference, and motion parameters (X and Y displacement) are plotted. Before and after correction, average images are compared to assess the correction's effectiveness.

Next, time-series signals are extracted for each segmented cell region (using maskC) to analyze neural activity across time, and these traces are plotted for selected cells. In parallel, a neuropil mask is created by excluding both cell bodies and blood vessels, enabling the extraction of background calcium signals. These neuropil signals help distinguish true cellular activity from background fluctuations. Following that, k-means clustering is performed on the cell-wise activity data to uncover patterns of coordinated activity, and cluster identities are mapped back to spatial locations. Cluster-averaged signals are plotted to visualize distinct activity profiles, and their correlations are assessed.

The assignment continues with Principal Component Analysis (PCA), where the motion-corrected data is reshaped, mean-centered, and decomposed using singular value decomposition. This reduces dimensionality and highlights dominant spatiotemporal patterns. The mean intensity over time and the first 12 temporal components (U vectors) are plotted, followed by spatial coefficient maps (W vectors) for these components. Finally, the average blood vessel signal over time is plotted and correlated with neuronal time-series data to investigate neurovascular coupling, quantified using Pearson's correlation coefficient.

