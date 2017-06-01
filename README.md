# OccupancyMS
An algorithm to calculate the site stoichiometry of a phospho-sites (or other modifications) measured under multiple conditions with quantitative mass spectrometry and provide confidence intervals. The approach for single sites is implemented in MATLAB for datasets of an arbitrary size. Similar code for multi-sites is forthcoming. 

Created by Marc Presler, Martin Wuehr, Allon Klein, Elizabeth Van Itallie. 

# To run:
-Download repository files and move to the same filepath of your choice. 

-Run script 'OccupancyCode_Superscript_v1' using Matlab, which calls the other necessary functions included in the download. 

-Sample data is provided in the file 'SampleFile_v1.xlsx.'

-Follow instructions in the comments. In short, code requires the user to provide the MS data for matched unmodified and modified forms of a peptide in an excel file. The variables for the data columns for the unmodified and modified forms are specified, as well as the Gene Symbol and Site Position for labeling purposes. Be sure to specify confidence interval parameters in Section 2 (see tips!). 

-Output is a .xlsx file of the estimated occupancy, the high and low bounds of the confidence intervals, as well as the original input data. Filename is records the parameters and input file use to generate the data along with the date.   

# Visualizations 
-Run with sample data to see the visualziations
-The unmodified vs. modified plot over the conditions. The values at each of these species at each conditions are used to perform the regression. 
-Regression plot with data, best fit, and confidence intervals (see tips!)
-Occupancy trend over conditions, as separate subplots or plotted together. 

# Tips:
-May required Matlab 2014 or later. Pseudocode for nonmatlab implementations are forthcoming. 

-For initial runs, limit the amount of bootstrapping iterations to 10-100 in Section 2, as the code is can be quite slow on most computers for larger datasets when using the recommended 1,000-10,000 iterations.

-For large datasets, set "plot_fitting_scatter' to 0 to avoid crashing matlab due to a graphics error. This functionality to most useful for visualizing the regression of a limited number of sites. 

-The data plotting functionality is compatable with large datasets, but may become impractical if the dataset exceeds 100 sites. Make sure to adjust plotting parameters to reflect your data. 

