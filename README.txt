********************
STARS - Slice-based Treatment-planning Algorithm for Radiosurgery Software- Version 0.2, "TomoAlpha"
********************
note: formerly known as 'Tomosurgery', prior to July 2012.

Author Info:
===================
Indraneel Gowdar, MD, MS
Department of Biomedical Engineering
Case Western Reserve University
ipg4@case.edu
707-646-9327

Background & Features:
===================
- Dynamic Slice-based Inverse Treatment-planning for Radiosurgical devices
- This software represents the latest iteration of 2+ years of on/off development. Originally written in MATLAB, this software

Instructions:
===================
To run a treatment plan optimization, a minimum of 2 pieces of information are needed:

1) DOSE KERNEL 
-----------------
A volumetric matrix representing the dose-rate at the isocenter of the desired radiosurgical device (i.e. the Leksell Gamma Knife). Choose 'File -> Load Dose Kernel...' and choose a *.txt file.

2) STRUCTURE SET
-----------------
A file containing the volume data for the lesion/tumor/ROI to be treated. Either this can be from a DICOM-RT file, such as 'RTSS.dcm', or it can be a separately processed text file containing the sequential values of a 3D volume (as long as a '*.h' header file is included). Choose 'File -> Load Structure...'.

Once these two inputs are loaded, a plan can be successfully run. However, DICOM image viewer functionality is included for visualization purposes.

TO RUN A PLAN:
-----------------
- 1) Load the two inputs described above.
- 2) Choose a desired slice thickness using the slider located at the bottom of the window.
- 3) Click the 'Plan' button (it should no longer be gray if both inputs are loaded correctly).
-- The structure will be divided into slices and shots will be planned according to the default parameters (Raster Width, Step Size, and Slice Thickness). 
- 4) Ensure that you are in the 'Plan' tab. You can now adjust the step-size and raster width of the plan for each slice individually using the parameter sliders. The red dots on the structure image will adjust in real-time. If you would like to see a dose overlay, click "Dose Preview" radio button.
- 5) Once you are satisfied with the plan, click 'Optimize' to run the weighting algorithm, which will choose the weights of each shot to optimize coverage and minimize dose external to the structure.
- 6) When finished, view the resulting dose and the associated coverage values listed in the list box. 