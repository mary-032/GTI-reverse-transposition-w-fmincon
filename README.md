# GTI-reverse-transposition-w-fmincon

This repository was created to provide complimetary material for the article cited below to the interested reader.

Insert proper citation later:
"Title", 2025
Rynoson, Marieke;
Ma Lu, Silvia;
Munkhammar, Joakim;
Campana, Pietro Elia;
doi:

It includes additional results which were not included in the article due to spatial constraints, 
as well as an example of the code used for the fmincon function used to estimate horizontal irradiance
components from global tilted irradiance.

FOLDER: HEATMAP PLOTS
Absolute and relative error for each model were plotted in relation to the corresponding solar azimuth and elevation.
This file contains both the individual plots and figures that show the performance of each model under the same circumstances.
The naming is as follows:
error type_estimated irradiance_starting GTI tilt angle_model

error type: Absolute (abs) and relative (rel)
est. irradiance: global horizontal (GHI) or diffuse horizontal (DHI)
starting GTI tilt angles: 30, 40, or 90 degrees
models: {[Erbs (E), Skartveit1 (S), Engerer2 (E2), Yang4 (Y4)] in combination with [Hay & Davies (HD), Perez1990 (P90)]}, GTI-DIRINT (D), Perez-Driesse (PD)

Also included: a heatmap plot of measured GHI. From this, the reader can estimate if a time period had mainly clear or cloudy days and compare that to when errors occured for the models.

FOLDER: STUDY RESULTS
Currently includes a boxplot which demonstrates the range of error (nRMSE) of each model for different tilt angles. 
We found that the accuracy generally decreases with increasing tilt angles. Naming of the plots similar to description above, but: model_starting GTI tilt angle
Not included in the paper due to purely spatial constraints.
