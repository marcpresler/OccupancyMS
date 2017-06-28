function [occupancy_trend_output] = call_fit_3D_data_SVD_v1(input_matrix_3D)

% call_fit_3D_data_SVD.m
% Marc Presler, December 28th, 2017
%
%Calls SVD-based orthogonal fitting in 3 dimensions for bootstrapping. 
%   Used to caclulates the normal vector to the plane then outputs the single and
%   double occupancy ratios, and is called when calculating confidence intervals.
%
%Function calls Singular Value Decomposition function in 3D and corrects for "impossible" values outside of 0 and 100%.  
%
% Inputs:
%   input matrix:   3 x "number of conditions" matrix of data. Row 1 is unmodified, 
% and Row 2 is first phos form, row 3 is the second phos form.
%
% Output:
%    ratio_output:  3 x 1 vector of occupancy for timepoint of interest.
%           Row 1 is unmodified occupancy. Row 2 is first phosform
%           occupancy. Row 3 is the second phosform occupancy. 


[normal_vector,~, ~,~] = fit_3D_data_modified_SVD_v1(input_matrix_3D(:,1),input_matrix_3D(:,2),input_matrix_3D(:,3));
 
%Calculates occupany from the orthogonal vector (i.e., normal) to the plane
%defined by the experimental values. We find the vector whose elements sum
%to 1 (similar to a stochastic or probability vector) to satisfy the conservation law.
%We multiply by 100 to convert to percent occupancy. 
occupancy_trend = normal_vector./repmat(sum(normal_vector),3,1)*100;

%Call each individual occupancy
occupancytrend_double = occupancy_trend(3);
occupancytrend_single = occupancy_trend(2);
occupancytrend_non = occupancy_trend(1);

%Correct for impossible values
if occupancytrend_non > 100
    occupancytrend_non = 100;
end

if occupancytrend_single > 100
    occupancytrend_single = 100;
end

if occupancytrend_double > 100
    occupancytrend_double = 100;
end


if occupancytrend_non < 0
    occupancytrend_non = 0;
end

if occupancytrend_single < 0
    occupancytrend_single = 0;
end

if occupancytrend_double < 0
    occupancytrend_double = 0;
end

%Store for output
occupancy_trend_output = [occupancytrend_non;occupancytrend_single;occupancytrend_double];
end

