# SARA_Code
Matlab code for SARA (Seismic Amplitude Ratio Analysis) used in the paper *"Assessing the potential to detect magma migration using seismicity through the analysis of the seismic network design"* 
The code in this repo uses an excerpt of the seismic network of Piton de la Fournaise. This can easily be changed to another network providing the user knows the location of the seismic stations as well as some basic parameters like the summit of the volcano

## Details of the Files
* detect_positions_piton.m
* detect_add_positions_piton.m
* detect_remove_positions_piton.m

## How to Run
Running the code is pretty simple. The user can alter the parameters i.e. migration length, threshold, grid spacing if they wish.
Currently the code is set up to analyse vertical migrations of ... using a threshold of 

## External Functions
External functions are used in the code and the below explains where you can download each from. I would additionally like to thank the authors of these codes.

* "utm2deg" and "deg2utm" (Rafael Palacios, 2006)
  * Rafael Palacios (2022). utm2deg (https://www.mathworks.com/matlabcentral/fileexchange/10914-utm2deg), MATLAB Central File Exchange. Retrieved November 24, 2022.
  * Rafael Palacios (2022). deg2utm (https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm), MATLAB Central File Exchange. Retrieved November 24, 2022.
*
