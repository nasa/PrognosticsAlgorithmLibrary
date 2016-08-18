function mX = wmean(X,w)
% wmean   Weighted mean calculation
%
%   mX = wmean(X,w) computes the weighted mean given samples X and weights
%   w. 
%
%   X must be formatted such that the rows are variables and the columns
%   are samples. If X is a vector, then both X and w can be either row or
%   column vectors.
%
%   See also Observers.wcov
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Ensure that if X is a vector, then both X and w are column vectors
if size(X,2)==1
	% Then X is a column vector, make it a row vector (for the
	% multiplication)
	X = X';
end
if size(w,1)==1
	% Then w is a row vector, make it a column vector
	w = w';
end

% Compute weighted mean
mX = X*w;
