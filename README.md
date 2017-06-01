# OccupancyMS
An algorithm to calculate the site stoichiometry of a phospho-sites or other modifications measured under multiple conditions with quantitative mass spectrometry. The approach for single sites is implemented in MATLAB for datasets of an arbitrary size. Code for multi-sites is forthcoming. 

Created by Marc Presler, Martin Wuehr, Allon Klein, Elizabeth Van Itallie. 

# To run:
-Download repository files and move to the same filepath of your choice. 
-Run script 'OccupancyCode_Superscript_v1', which calls the other necessary functions
-Sample data is provided in the file 'SampleFile_v1.xlsx'
-Follow instructions in the comments. In short, code requires the user to provide the MS data for matched unmodified and modified forms of a peptide in an excel file. The variables for the data columns for the unmodified and modified forms are specified, as well as the Gene Symbol and Site Position for labeling purposes. Be sure to specify confidence interval parameters in Section 2 (see tips!). 

# Tips:
-Pseudocode for nonmatlab implementations are forthcoming.
-For initial runs, limit the amount of bootstrapping iterations to 10-100 in Section 2, as the code is can be quite slow on most computers for larger datasets when using the recommended 1,000-10,000 bootstraps. 
-For large datasets, set "plot_fitting_scatter' to 0 to avoid crashing matlab due to a graphics error. This functionality to most useful for visualizing the regression of a limited number of sites. 
-The data plotting functionality is compatable with large datasets, but may become impractical if the dataset exceeds 100 sites. 

