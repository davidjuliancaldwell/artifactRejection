### Signal recovery from stimulation artifacts in intracranial recordings with dictionary learning

This is the readme file for the processing of stimulation related data in human electocorticography and deep brain stimulation datasets. The code is written in MATLAB. This is distributed as a MATLAB package, so simply download the package, and place it on the path. Only the root path is necessary, e.g. addpath '/path/to/artifactRejection', as all required functions are within the project. [MATLAB](https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html "MATLAB Packages") packages operate like namespaces (see link).

To use example data and to have all of the code required to perform the processing, clone the github repository, and download the data in the google drive link below. To run the code on the example data, place the data linked below in a folder named **+data**, in the same directory as the main script .

**https://drive.google.com/open?id=1yoOty7SPI9mcNtWc16nakC6kCDz_QtL4**

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

Note: For the recover_EP.m function, the MATLAB curve fitting toolbox is required.

Direct questions and comments to David Caldwell, at djcald@uw.edu

Authors: David J. Caldwell, Jeneva A. Cronin, Rajesh P.N. Rao, Kelly L. Collins, Kurt E. Weaver, Andrew L. Ko, Jeffrey G. Ojemann, Nathan J. Kutz, Bingni W. Brunton

___

BSD-3 License
