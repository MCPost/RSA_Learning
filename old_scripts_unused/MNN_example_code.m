%% Multivariate Noise Normalization

% This code is basically the implementation of Guggenmos et al. (2018) as
% can be found (https://github.com/m-guggenmos/megmvpa). I adapted it for
% my analysis.
% 
% What it does is to estimate the noise variance within and the
% noise covariance between all channels. This covariance matrix is then
% used to weigh the data (normalize it). Since variances and covariances
% are estimated for each channel (here 40 channels) averaging trials (here
% 8 or less trials), the covariance matrix is likely to be ill-conditioned,
% rendering it non-invertible. To accomodate that, a shrikange parameter is
% estimated from the data and the covariance matrix is shrunk towards a
% simple diagonal matrix of variances until it is invertible.
%
% This is done seperately for each participant
%
%' Input Variables:
%       Data_EEG   -   Trial x Channel x Time points EEG data of one subject
%       TrialInfo  -   Cell array with descriptions of all (128) trials of
%                      one subject matching the order of Data_EEG (tells
%                      you whether trial 1 belongs to animate drawing or
%                      picture etc.)
%       curROI     -   index for a subset of channels (occipital region for
%                      example)

nonnan_ind = find(cell2mat(TrialInfo(:,1)) ~= 0);   % Indexes all NaN trials from the artifact correction
Sigma = zeros(16,length(curROI),length(curROI));    % Initialize Noise Covariance Matrix
for c = 1:16                                        % Loop over categories (the size of the final RSA Matrix)
    tmp = zeros(size(Data_EEG,3),length(curROI),length(curROI)); 
    for t = 1:size(Data_EEG,3)                                          % Loop over all timepoints
        dat = Data_EEG(intersect(nonnan_ind, (c-1)*8+(1:8)),curROI,t);  % Take all trials of one category (that are not NaN) per time point in one ROI
        tmp(t,:,:) = cov1para(dat(~sum(isnan(dat),2),:));               % This function estimates a covariance matrix of the channels (averaging over trials) for this time point
    end
    Sigma(c,:,:) = nanmean(tmp,1);                   % Average covariance matrices of all time points together
end
tmp = squeeze(mean(Sigma,1));   % Average over all categories
sigma_inv = tmp^-0.5;           % Take the square root of the inverse
Data_EEG_corr = zeros(size(Data_EEG,1), length(curROI), size(Data_EEG,3));  % Initialize the corrected Data set
for t = 1:size(Data_EEG,3)                                  % Loop over all time points 
    Data_EEG_corr(:,:,t) = Data_EEG(:,curROI,t)*sigma_inv;  % Correct the data of all trials and ROI channels per time point with the inverse square rooted covariance matrix
end



% I took this function from the Guggenmos code
% (https://github.com/m-guggenmos/megmvpa). Probably has to be copied into
% its own m-File to work.
function [sigma,shrinkage] = cov1para(x,shrink)

% function sigma=cov1para(x)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards one-parameter matrix:
%    all variances are the same
%    all covariances are zero
% if shrink is specified, then this value is used for shrinkage

% This version: 04/2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is released under the BSD 2-clause license.

% Copyright (c) 2014, Olivier Ledoit and Michael Wolf 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% de-mean returns
[t,n]=size(x);
meanx=mean(x);
x=x-meanx(ones(t,1),:);

% compute sample covariance matrix
sample=(1/t).*(x'*x);

% compute prior
meanvar=mean(diag(sample));
prior=meanvar*eye(n);

if (nargin < 2 | shrink == -1) % compute shrinkage parameters
  
  % what we call p 
  y=x.^2;
  phiMat=y'*y/t-sample.^2;
  phi=sum(sum(phiMat));
  
  % what we call r is not needed for this shrinkage target
  
  % what we call c
  gamma=norm(sample-prior,'fro')^2;

  % compute shrinkage constant
  kappa=phi/gamma;
  shrinkage=max(0,min(1,kappa/t));
    
else % use specified number
  shrinkage=shrink;
end

% compute shrinkage estimator
sigma=shrinkage*prior+(1-shrinkage)*sample;

end

