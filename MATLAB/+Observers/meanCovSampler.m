function X = meanCovSampler(estimate,numSamples)
% meanCovSampler  Sample generator given mean and covariance information.
%
%   X = meanCovSampler(estimate,numSamples) generates a set of samples X,
%   which is a matrix with the rows the variables and the columns the
%   samples, given an estimate structure with fields 'mean' and
%   'covariance'. numSamples samples are generated.
%
%   The estimate structure has two fields, 'mean' and 'covariance'. The
%   mean must be a column vector, and the covariance a square matrix of the
%   corresponding size.
% 
%   This function will work with state estimates from the Kalman, extended
%   Kalman, and unscented Kalman filters.
%
%   See also Observers.KalmanFilter, Observers.ExtendedKalmanFilter,
%   Observers.UnscentedKalmanFilter
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

X = repmat(estimate.mean,1,numSamples) + ...
    chol(estimate.covariance)'*randn(length(estimate.mean),numSamples);
