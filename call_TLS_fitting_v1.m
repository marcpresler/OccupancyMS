function [occupancy_output] = call_TLS_fitting_v1(input_matrix, FLAG_repeated_points)

% call_TLS_fitting_v1.m
% Marc Presler, December 16th, 2016

%Function calls Total Least Squares fitting function in 2D.

% Inputs:
%   input matrix:   2 x number of conditions matrix of data. Row 1 is non-phos 
% and Row 2 is phos.
%
%   FLAG_repeated_points:   either 0 or 1 (only relevant if called during
%   bootstrapping
%       Default is 0 
%       1 if data is not unique

% Output:
%    ratio_output: 2 x 1 vector of occupancy for timepoint of interest.
%           Row 1 is phos occupancy. Row 2 is non-phos occupancy.


% If triggered, checks to make sure values are unique, if yes proceeds through function,
% if no, then returns 100% or 0% occupancy depending on input values
if FLAG_repeated_points
    
    %If line is vertical, give 0% occupancy 
    if size(unique(input_matrix(:,1)),1)==1

        occupancy_output = [1.e-7;99.9999999];
        return
    end


    % If line is horizantal, give 100% occupancy
    if size(unique(input_matrix(:,2)),1)==1
         
        occupancy_output = [99.9999999,1.e-7];
        return
    end
end


% Calculate slope of data using total least squares regression, nonphos vs phos
[Error, Fit_slope_single_phos] = fit_2D_data_modified_v1(input_matrix(:,1),input_matrix(:,2),'no');

%If error occurs because the same point was sampled 10 times during boostrapping, assigns a
%random occupancy value
if (Error)
    
    random_occ = randi(100,1);
    occupancy_output = [random_occ, 100-random_occ];
    return
    
end

%solve for the slope of the perpendicular vector 
occupancy_ratio_single = -(1./Fit_slope_single_phos(1,1));

        %Add correction for negative slopes
            %Set anything between 0 and -1 to ~0
            if occupancy_ratio_single < 0 && occupancy_ratio_single >=-1 
                occupancy_ratio_single = 1E-9;
            end
            
             %Set anything less than -1 to inifity-ish
            if occupancy_ratio_single < -1 
                occupancy_ratio_single = 1E9;
            end

%Calculate percent occupancy, which is the output 
occupancytrend_single = ((occupancy_ratio_single./(1+occupancy_ratio_single))*100);

occupancytrend_non = ((1./(1+occupancy_ratio_single))*100);

occupancy_output = [occupancytrend_single;occupancytrend_non];


end

