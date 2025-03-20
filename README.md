# neuralengineering

This repository is a collection of work I did for a neural engineering course during Fall 2024 at the University of Pittsburgh. 

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
