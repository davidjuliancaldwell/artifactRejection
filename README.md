# Signal recovery from stimulation artifacts in intracranial recordings with dictionary learning
---
This is the readme file for the processing of stimulation related data in human electocorticography and deep brain stimulation datasets. The code is written in MATLAB. This is distributed as a MATLAB package, so simply download the package, and place it on the path. Only the root path is necessary, e.g. addpath '/path/to/artifactRejection', as all required functions are within the project. [MATLAB](https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html "MATLAB Packages") packages operate like namespaces (see link).

To use example data and to have all of the code required to perform the processing, clone the github repository, and download the data in the google drive link below. To run the code on the example data, place the data linked below in a folder named **+data**, in the same directory as the main script .

**https://drive.google.com/open?id=1yoOty7SPI9mcNtWc16nakC6kCDz_QtL4**

---
#### Description of folders and files

The script to call and run section by section is:
**artifact_rejection_script.m** - Running this script will generate figures in the style presenting in the manuscript.

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

#### artifact_rejection_script.m
**dataInt** = time x channels x epochs

**fsData** = sampling rate of the data Hz

**stimChans** = the channels used for stimulation . These should be noted and excluded from further analysis

**plotIt** = determines whether or not to plot the intermediate results of the functions.

**tEpoch** = epoched time window (s)

 **minDuration** = minimum duration of artifact in ms (0.5 ms as default for 200 &mu;s pulses, 0.25 ms for 60 &mu;s pulses )

 **pre** = default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (0.8 ms as default)

 **post** = default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (1 ms as default)

**type** = this determines the type of template approach to use. The three options are **average**, **trial**, and **dictionary**.

This is used in the **analyFunc.template_subtract.m** function

**average** is the simplest approach, and on a channel by channel basis it simply averages the artifact for each channel.

**trial** uses the stimulation pulse train within each trial.

**dictionary** is the most advanced, and uses a template matching algorithm with DBSCAN clustering to discover the family of artifacts. The manuscript primarily uses this approach. Below we highlight the options pertinent to this selection.

+ **eucl**, **cosine**, **corr**, for either euclidean distance, cosine similarity, or correlation for clustering and template matching. These are the metrics used if the dictionary method is selected for both clustering (**distanceMetricDbscan**) and subsequent template matching (**distanceMetricSigMatch**). We used euclidean distance for the clustering, and correlation for template matching for all of our data sets. 

+ **distanceMetricDbscan** = (Default, 'eucl')
+ **distanceMetricSigMatch** = (Default, 'corr')
+ **amntPreAverage** = Number of samples over which to estimate the baseline of any given artifact window, starting from the first sample (Default, 3)
+ **normalize** = Method by which to normalize the data. Using **preAverage** uses the specified number of samples in **amntPreAverage**, while **firstSamp** just uses the first sample (Default, preAverage)



---

Note: For the recover_EP.m function, which is not used in the manuscript or for any figures, but could be helpful for processing signals with large exponential recovery component parts of the signals, the MATLAB curve fitting toolbox is required.

---
### Questions and comments

Direct questions and comments to David Caldwell, at djcald@uw.edu

Authors: David J. Caldwell, Jeneva A. Cronin, Rajesh P.N. Rao, Kelly L. Collins, Kurt E. Weaver, Andrew L. Ko, Jeffrey G. Ojemann, Nathan J. Kutz, Bingni W. Brunton

___

BSD-3 License
