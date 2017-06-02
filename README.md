# OccupancyMS
An algorithm to calculate the site stoichiometry of a phospho-sites (or other modifications) measured under multiple conditions with quantitative mass spectrometry and provide confidence intervals. The approach for single sites is implemented in MATLAB for datasets of an arbitrary size. The code for multi-sites will be posted shortly.  

Created by Marc Presler, Martin Wuehr, Allon Klein, Elizabeth Van Itallie. 

# To run:
1) Download repository files and move to the same filepath of your choice. 

2) Run script 'OccupancyCode_Superscript_v1' using Matlab, which calls the other necessary functions included in the download. 

  -Sample data is provided in the file 'SampleFile_v1.xlsx.'

3) Follow instructions in the comments. In brief, code requires the user to provide the MS data for matched unmodified and modified forms of a peptide in an excel file. The variables for the data columns for the unmodified and modified forms are specified, as well as the Gene Symbol and Site Position for labeling purposes. Be sure to specify confidence interval parameters in Section 2 (see tips!). 

4) Output is a .xlsx file of the estimated occupancy, the high and low bounds of the confidence intervals, as well as the original input data. Filename is records the parameters and input file use to generate the data along with the date.   

# Visualizations 
-Run code with 'SampleFile_v1.xlsx' to see the following visualziations:

1) The unmodified vs. modified plot over the conditions. The values at each of these species at each conditions are used to perform the regression. 

2) Regression plot with data, best fit, and confidence intervals (see tips!)

3) Occupancy trend over conditions, as separate subplots or plotted together. 

# Tips:
1) May require Matlab 2014 or later. Pseudocode for nonmatlab implementations is forthcoming. 

2) For initial runs, limit the amount of bootstrapping iterations to 10-100 in Section 2, as the code is can be quite slow on most computers for larger datasets when using the recommended 1,000-10,000 iterations.

3) For large datasets, set "plot_fitting_scatter' to 0 to avoid crashing matlab due to a graphics error. This functionality to most useful for visualizing the regression of a limited number of sites. 

4) The data plotting functionality is compatable with large datasets, but may become impractical if the dataset exceeds 100 sites. Make sure to adjust plotting parameters (i.e., x axis labels) to reflect your data. 

