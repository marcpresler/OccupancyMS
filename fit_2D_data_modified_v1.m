function [Error, P, Yhat] = fit_2D_data_modified_v1(XData, YData, vizualization)
%
% fit_2D_data_modified_v1.m
% Modified by Marc Presler from MathWorks File Exchance ID: #31109.
% December 16, 2016

% Orthogonal linear regression method in 2D for model: y = a + bx   
%
% Input parameters:
%  - XData: input data block -- x: axis
%  - YData: input data block -- y: axis
%  - vizualization: figure ('yes','no')
% Return parameters:
%  - Err: error - sum of orthogonal distances 
%  - P: vector of model parameters [b-slope, a-offset] 
%
% Other return parameters:
% - Error: 0 or 1. O if default. 1 is return if a line cannot be fit
% because the data is 1 point
%
% Original authors:
% Ivo Petras (ivo.petras@tuke.sk) 
% Dagmar Bednarova (dagmar.bednarova@tuke.sk)
% Date: 24/05/2006, 15/07/2009

% Sets error flag to default of 0. 
Error = 0;

kx=length(XData);
ky=length(YData);
if kx ~= ky
   disp('Incompatible X and Y data.');
   close all;
end

sy=sum(YData)./ky;
sx=sum(XData)./kx;
sxy=sum(XData.*YData);
sy2=sum(YData.^2);
sx2=sum(XData.^2);
B=0.5.*(((sy2-ky.*sy.^2)-(sx2-kx.*sx.^2))./(ky.*sx.*sy-sxy));
b1=-B+(B.^2+1).^0.5;
b2=-B-(B.^2+1).^0.5;
a1=sy-b1.*sx;
a2=sy-b2.*sx;
R=corrcoef(XData,YData);


%This is triggered if the boostrapping calls the same point ten times, which will
%cause the corrcoef later in the code to fail. Flag this and return to calling script
if isnan(R(1,2))
    Error  = 1;
    %doesn't matter what p is
    P = [1 1];
    return;
    
end


if R(1,2) > 0 
    P=[b1 a1];
    Yhat = XData.*b1 + a1;
    Xhat = ((YData-a1)./b1);
end
if R(1,2) < 0
    P=[b2 a2];
    Yhat = XData.*b2 + a2;
    Xhat = ((YData-a2)./b2);
end   
alpha = atan(abs((Yhat-YData)./(Xhat-XData)));
d=abs(Xhat-XData).*sin(alpha);
Err=sum(abs(d));
%Err=sum(d.^2);

switch lower(vizualization)
     case {'yes'}
        plot(XData,YData,'blue*'); 
        hold on;
        plot(XData,Yhat,'black');
        hold off
     case {'no' }
%          disp('No vizualization.')
        return
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