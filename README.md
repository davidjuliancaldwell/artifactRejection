### Processing of Electrical Stimulation Related Human Electrophysiologic Data

This is the readme file for the processing of stimulation related data in human electocorticography datasets. The data is written in MATLAB, and should require no special packages for use. This is distributed as a MATLAB package, so simply download the package, and place it on the path. Only the root path is necessary, e.g. addpath '/path/to/artifactRejection', as all required functions are within the project. [MATLAB](https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html "MATLAB Packages") packages operate like namespaces (see link).

To use example data and to have all of the code required to perform the processing, clone the github repository, and download the data in the google drive link below. To run the code on the example data, place the data linked below in a folder named **+data**

**https://drive.google.com/open?id=1yoOty7SPI9mcNtWc16nakC6kCDz_QtL4**

The script to call and run section by section is:
**artifact_rejection_script.m**

The **+helperFunc** folder contains helper processing functions that are called by analysis scripts, and with modifications may prove useful to other researchers.

The **+vizFunc** folder contains helper processing functions that prove useful in visualizing the data here, and with modifications may prove useful to other researchers.

The **+analyFunc** folder contains the bulk of the analysis functions called for processing human stimulation electrocorticography data.

Direct questions and comments to David Caldwell, at djcald@uw.edu

Authors: David Caldwell, Jeneva Cronin, Nathan Kutz, Bingni Brunton

___

BSD-3 License
