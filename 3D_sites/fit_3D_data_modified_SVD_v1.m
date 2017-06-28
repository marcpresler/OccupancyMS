function [normal_vector,MeanCenteredData,Vmatrix_3d, Smatrix_3d] = fit_3D_data_modified_SVD_v1(XData, YData, ZData)
%
% [normal_vector,MeanCenteredData,Vmatrix_3d, Smatrix_3d] = fit_3D_data_modified_SVD(XData, YData, ZData)
%
% Modified by Marc Presler from MathWorks File Exchance ID: #31109.
% June 28th, 2017
%
% Input parameters:
%  - XData: input data block -- x: axis
%  - YData: input data block -- y: axis
%  - ZData: input data block -- z: axis 

% Return parameters:
%  - normal_vector: the vector with the least variation determined by
%  Singulary value decompositon (i.e., the last "principal component").
%  This Vector us used to determine the occupancy of each species. 
%
% Other return parameters:
% - MeanCenteredData: Raw data with means substracted. 
% - Vmatrix_3d, Smatrix_3d: The V and S matrix of the SVD. 
%
% Original authors:
% Ivo Petras (ivo.petras@tuke.sk) 
% Dagmar Bednarova (dagmar.bednarova@tuke.sk)
% Date: 24/05/2006, 15/07/2009
%
% A modified version of the plotting technique is included in the main
% script.
%
% Note from original authors: Written for Matlab 7.0 (R14) with Statistics Toolbox
%
% We sincerely thank Peter Perkins, the author of the demo,
% and John D'Errico for their comments.
%

XData = XData - mean(XData);
YData = YData - mean(YData);
ZData = ZData - mean(ZData);

%
MeanCenteredData(:,1) = XData(:,1);
MeanCenteredData(:,2) = YData(:,1);
MeanCenteredData(:,3) = ZData(:,1);
%
[Umatrix_3d, Smatrix_3d, Vmatrix_3d] = svd(MeanCenteredData);

normal_vector = Vmatrix_3d(:,3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (c) 2011, Ivo Petras
% Copyright (c) 2016, Igor Podlubny
% All rights reserved.
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.