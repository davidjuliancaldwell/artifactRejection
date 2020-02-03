# Signal recovery from stimulation artifacts in intracranial recordings with dictionary learning
---
This is the readme file for the processing of stimulation related data in human electocorticography and deep brain stimulation datasets. The code is written in MATLAB. This is distributed as a MATLAB package, so simply download the package, and place it on the path. Only the root path is necessary, e.g. addpath '/path/to/artifactRejection', as all required functions are within the project. [MATLAB](https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html "MATLAB Packages") packages operate like namespaces (see link).

To use example data and to have all of the code required to perform the processing, clone the github repository, and download the data in the google drive link below. To run the code on the example data, place the data linked below in a folder named **+data**, in the same directory as the main script .

**https://drive.google.com/open?id=1yoOty7SPI9mcNtWc16nakC6kCDz_QtL4**

---
#### Description of folders and files

The script to call and run section by section is:
**artifact_rejection_script.m** - Running this script will generate figures in the style presenting in the manuscript (Figures 1, 2, 3, 4, 5 (panel a, b, c), 7, 8, and the supplemental figures)

 **haptic_touch_script.m** will create the natural haptic touch figures for Figure 5, panel d, e, and f, in the manuscript.

  **button_press_script.m** will create the figures for Figure 6 in the manuscript.

The **+analyFunc** folder contains the bulk of the analysis functions called for processing human stimulation electrocorticography data. Some of the functions in this folder are described below, and can be accessed by analyFunc.(____)

 **template_subtract.m** - this is a function which acts upon a samples x sensors x trial basis. Functions called by this include **get_artifact_indices**, **get_artifacts**, **template_equalize_length**, and either **template_average**, **template_trial**, or **template_dictionary**, depending on options selected.

**waveletWrapper** - this is a utility for performing wavelet decompositions of epoched data using Morlet wavelets, which calls **morletProcess** to perform the wavelet analysis.

The **+helpFunc** folder contains helper processing functions that are called by analysis scripts, and with modifications may prove useful to other researchers

The **+savitskyGolay** folder contains code to perform Savitszky-Golay filtering as part of the smoothing process for detecting artifact onset and offset, without additional MATLAB toolboxes.

The **+vizFunc** folder contains helper processing functions that prove useful in visualizing the data here, and with modifications may prove useful to other researchers. These include visualizing all epochs of data across different channels in the time domain, subselecting channels of interest, visualizing the power in different frequency bins, and visualizing time-frequency plots

The **+cbrewerHelper** folder contains the Colorbrewer colormaps for MATLAB.

The **+HDBSCAN** folder contains the code to perform HDBSCAN upon the data.

---

## Description of various parameters

We have found the parameters included for many of the variables which are used to control the artifact detection, clustering, scaling, and subtracting to work for a number of our datasets. In the manuscript, we include the values used for key parmaeters. Below, we will discuss a number of these parameters and provide guidelines for setting them. The code also contains comments throughout describing the purpose of parameters and reasonable values for them.

### artifact_rejection_script.m
**dataInt** = time x channels x epochs

**fsData** = sampling rate of the data Hz

**stimChans** = the channels used for stimulation . These should be noted and excluded from further analysis

**plotIt** = determines whether or not to plot the intermediate results of the functions.

**tEpoch** = epoched time window (s)

 **minDuration** = minimum duration of artifact in ms (0.5 ms as default for 200 &mu;s pulses, 0.25 ms for 60 &mu;s pulses )
#### parameters used for get_artifact_indices.m
 **pre** = default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (0.8 ms as default)

 **post** = default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (1 ms as default)

 **onsetThreshold** = This value is used as absolute valued z-score threshold to determine the onset of artifacts within a train. The differentiated smoothed signal is used to determine artifact onset. This is also used in determining how many stimulation pulses are within a train, by ensuring that beginning of each artifact is within a separate artifact pulse.   (Default, 1.5)

 **threshVoltageCut** = This is used to help determine the end of each individual artifact pulse dynamically. More specifically, this is a percentile value, against which the absolute valued, z-scored smoothed raw signal is compared to find the last value which exceeds the specified percentile voltage value. Higher values of this (i.e. closer to 100) result in a shorter duration artifact, while lower values result in a longer calculated duration of artifact. This parameter therefore should likely be set higher for more transient artifacts and lower for longer artifacts. (Default for first ECoG dataset, 75)

 **threshDiffCut** = This is used to help determine the end of each individual artifact pulse dynamically. More specifically, this is a percentile value, against which the absolute valued, z-scored differentiated smoothed raw signal is compared to find the last value which exceeds the specified percentile voltage value. Higher values of this (i.e. closer to 100) result in a shorter duration artifact, while lower values result in a longer calculated duration of artifact. This parameter therefore should likely be set higher for more transient artifacts and lower for longer artifacts. (Default for first ECoG dataset, 75)

 The longer artifact duration of the two values calculated by **threshVoltageCut** and **threshDiffCut** is used for determining the end of the artifact window.

**type** = this determines the type of template approach to use. The three options are **average**, **trial**, and **dictionary**.

#### parameters used in **analyFunc.template_subtract.m**

**average** is the simplest approach, and on a channel by channel basis it simply averages the artifact for each channel.

**trial** uses the stimulation pulse train within each trial.

**dictionary** is the most advanced, and uses a template matching algorithm with DBSCAN clustering to discover the family of artifacts. The manuscript primarily uses this approach. Below we highlight the options pertinent to this selection.

+ **eucl**, **cosine**, **corr**, for either euclidean distance, cosine similarity, or correlation for clustering and template matching. These are the metrics used if the dictionary method is selected for both clustering (**distanceMetricDbscan**) and subsequent template matching (**distanceMetricSigMatch**). We used euclidean distance for the clustering, and correlation for template matching for all of our data sets.

+ **distanceMetricDbscan** = (Default, 'eucl')
+ **distanceMetricSigMatch** = (Default, 'corr')
+ **amntPreAverage** = Number of samples over which to estimate the baseline of any given artifact window, starting from the first sample (Default, 3)
+ **normalize** = Method by which to normalize the data. Using **preAverage** uses the specified number of samples in **amntPreAverage**, while **firstSamp** just uses the first sample (Default, preAverage)

There are an additional set of HDBSCAN parameters and window selection for template matching.

**bracketRange** = This variable sets the number of samples around the maximum voltage deflection to use for template clustering and subsequent matching. The smaller this range, the lower the dimensionality used for clustering, and the fewer points used to calculate the best matching template. This value is used to try and ensure that non-informative points are not included in the clustering and template matching. This should be set to what looks like the approximate length of the artifact's largest part. (Default [-6:6])

**minPts** = This is a parameter that determines how many neighbors are used for core distance calculations. Increasing this parameter restricts clusters to increasingly dense areas. (Default, 2)

**minClustSize** = The minimum number of clustered artifact pulses for a cluster to be labelled as a true cluster. Increasing this number can reduce the number of clusters, and merges some clusters together that would have otherwise been labelled as individual clusters.   (Default, 3)

**outlierThresh** = Outlier parameter for labeling artifact pulses as noise in the HDBSCAN clustering. Any artifact pulse with an outlier score greater than this will be labelled as noise. Increasing this value results in fewer points being labelled as noise (Default, 0.95)



---

Note: For the recover_EP.m function, which is not used in the manuscript or for any figures, but could be helpful for processing signals with large exponential recovery component parts of the signals, the MATLAB curve fitting toolbox is required.

For this function, there are optional parameters for exponential recovery, not currently used. These could be helpful for signals with large exponential recoveries

recoverExp = This determines whether to try and do the exponential recovery on a pulse by pulse basis. (Default, 0)

expThreshVoltageCut = Parameter similar to the threshVoltageCut parameter outlined above (Default, 95)

expThreshDiffCut = Parameter similar to the threshDiffCut parameter outlined above (Default, 95)

---
### Questions and comments

Direct questions and comments to David Caldwell, at djcald@uw.edu

Authors: David J. Caldwell, Jeneva A. Cronin, Rajesh P.N. Rao, Kelly L. Collins, Kurt E. Weaver, Andrew L. Ko, Jeffrey G. Ojemann, Nathan J. Kutz, Bingni W. Brunton

___

BSD-3 License
