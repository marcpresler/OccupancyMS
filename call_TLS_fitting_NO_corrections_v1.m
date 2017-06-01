function [occupancy_noncorrected_output] = call_TLS_fitting_NO_corrections_v1(input_matrix)


% call_TLS_fitting_NO_corrections_v1.m
% Marc Presler, December 16th, 2016
%
%Function calls Total Least Squares fitting function in 2D but allows for 
% impossible values greater than 0 and 100%, as there is no correction. 
%
% Inputs:
%   input matrix:   2 x number of conditions matrix of data. Row 1 is non-phos 
% and Row 2 is phos.
%
%
% Output:
%    ratio_output: 2 x 1 vector of occupancy for timepoint of interest.
%           Row 1 is phos occupancy. Row 2 is non-phos occupancy.



[~, Fit_slope_single_phos] = fit_2D_data_modified_v1(input_matrix(:,1),input_matrix(:,2),'no');


%solve for the slope of the perpendicular vector 
occupancy_ratio_single = -(1./Fit_slope_single_phos(1,1));


%Calculate percent occupancy 
occupancytrend_single = ((occupancy_ratio_single./(1+occupancy_ratio_single))*100);

occupancytrend_non = ((1./(1+occupancy_ratio_single))*100);

occupancy_noncorrected_output = [occupancytrend_single;occupancytrend_non];


end

