function Pxx = wcov(X,w,alpha,beta)
% wcov   Weighted covariance calculation
%
%   Pxx = wcov(X,w) computes the covariance matrix for samples X given
%   weights w.
%
%   Pxx = wcov(X,w,alpha,beta) computes the covariance matrix for samples X
%   given weights w and scaling parameters alpha and beta. This form is to
%   be used when the samples provided represent sigma points.
%
%   X must be formatted such that the rows are variables and the columns
%   are samples. If X is a vector, then both X and w can be either row or
%   column vectors.
%
%   See also Observers.wmean
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Ensure that if X is a vector, then both X and w are row vectors
if size(X,2)==1
	% Then X is a column vector, make it a row vector
	X = X';
end
if size(w,2)==1
	% Then w is a column vector, make it a row vector
	w = w';
end

% Compute weighted mean
mX = Observers.wmean(X,w);

% Initialize covariance
Pxx = 0;
for i=1:size(X,2)
	Pxx = Pxx + w(i)*(X(:,i)-mX)*(X(:,i)-mX)';
end

% Complete computation depending on whether samples are sigma points or not
if nargin<3
	Pxx = Pxx/(1-sum(w.^2));
else
	Pxx = Pxx + (1-alpha^2+beta)*(X(:,1)-mX)*(X(:,1)-mX)';
end
