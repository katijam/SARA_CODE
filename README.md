# SARA_Code
Matlab code for SARA (Seismic Amplitude Ratio Analysis) used in the paper *"Assessing the potential to detect magma migration using seismicity through the analysis of the seismic network design"* 

The code in this repo uses an excerpt of the seismic network of Piton de la Fournaise (PdF). This can easily be changed to another network providing the user knows the location of the seismic stations as well as some basic parameters like the summit of the volcano.

## Details of the Files
* detect_positions_piton.m - runs the SARA detection capability analysis and plots the output
* detect_add_positions_piton.m - adds new stations to the PdF network, runs the the SARA detection capability analysis, and compares the results to results gained from detect_positions_piton.
* detect_remove_positions_piton.m - removes existing stations from the PdF network, runs the the SARA detection capability analysis, and compares the results to results gained from detect_positions_piton.

## How to Run
Running the code is pretty simple. The user can alter the parameters i.e. migration length, threshold, grid spacing if they wish.
Currently the code is set up to analyse vertical migrations of 1km using a threshold of 0.1

## External Functions
External functions are used in the code and the below explains where you can download each from. I would like to thank the authors of these codes.

* utm2deg.m and deg2utm.m (Rafael Palacios, 2006)
  * Rafael Palacios (2022). utm2deg (https://www.mathworks.com/matlabcentral/fileexchange/10914-utm2deg), MATLAB Central File Exchange. Retrieved November 24, 2022.
  * Rafael Palacios (2022). deg2utm (https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm), MATLAB Central File Exchange. Retrieved November 24, 2022.
  *
* GetSRTMData.m (Sebastian Hölz, 2004). Note for this function that the users will need NASA SRTM elevation files.
  * Sebastian Hölz (2022). GetSRTMData (https://www.mathworks.com/matlabcentral/fileexchange/5544-getsrtmdata), MATLAB Central File Exchange. Retrieved November 24, 2022.

* axesLabelsAlign3D.m (Matthew Arthington, 2010)
  *  Matthew Arthington (2022). Align axes labels in 3D plot (https://www.mathworks.com/matlabcentral/fileexchange/27450-align-axes-labels-in-3d-plot), MATLAB Central File Exchange. Retrieved November 24, 2022.

* distinguishable_colors.m (Timothy E. Holy, 2010-2011)
  * Tim Holy (2022). Generate maximally perceptually-distinct colors (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors), MATLAB Central File Exchange. Retrieved November 24, 2022.

* colorcet.m (Peter Kovesi, 2015)
  * Peter Kovesi. Good Colour Maps: How to Design Them. arXiv:1509.03700 [cs.GR] 2015 (https://colorcet.com/download/index.html)
