function [X,W,mX,PX] = computeSigmaPoints(mx,Pxx,method,parameter,alpha,sqrtFcn)
% computeSigmaPoints   Sigma point calculation
%
%   [X,W,mX,PX] = computeSigmaPoints(mx,Pxx,method,parameter,alpha,sqrtFcn)
%   computes a set of sigma points defined by samples X and weights W,
%   given mean mx, covariance Pxx, a sigma point method, free parameters, a
%   scaling parameter, and a matrix square root function.
%
%   X will be formatted such that the rows are variables and the columns
%   are samples. W will be a row vector. mx must be a column vector. Pxx
%   must be a square matrix.
%
%   Available methods include:
%
%     1. 'symmetric' = symmetric sigma points. Free parameters are kappa,
%        alpha, beta.
%
%     2. 'simplex' = minimal skew simplex sigma points. Free parameters are
%        W0, alpha, beta.
%
%     3. 'spherical' = spherical simplex sigma points. Free parameters are
%        W0, alpha, beta
%
%    alpha is optional scaling parameter. Default alpha=1 (no scaling).
%
%    sqrtFcn is optional handle to a matrix sqrt function to use (sqrtm or
%    chol). Default is chol.
%
%   The mean mX and covariance PX of the sigma points are also returned,
%   for comparison to the provided mean and covariance that the sigma
%   points should capture.
%
%   See also Observers.wmean, Observers.wcov,
%   Observers.UnscentedKalmanFilter
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Set default options
if nargin<5
	alpha = 1;
end
if nargin<6
	sqrtFcn = @chol;
end

% Ensure alpha is nonzero
if alpha==0
	error('alpha cannot be 0!');
end

% Compute number of samples
n = size(mx,1);

% Compute sigma points
switch method
	
    case 'symmetric'
		
		% Collect parameters
		kappa = parameter(1);
		
		% Compute unscaled sigma points
		X(:,1) = mx;
		matrixSq = sqrtFcn((n+kappa)*Pxx)';
		X(:,2:n+1) = repmat(mx,1,n) + matrixSq;
		X(:,n+2:2*n+1) = repmat(mx,1,n) - matrixSq;
		
        % Compute weights
		W(1) = kappa/(n+kappa);
		W(2:2*n+1) = .5/(n+kappa);
		
	case 'simplex'
		
		% Compute weights
		W(1) = parameter(1);
		W(2:3) = (1-W(1))/2^n;
		for i=4:n+2
			W(i) = 2^(i-3)*W(2);
		end
		
		% Initialize vector sequence
		X(1,1) = 0;
		X(1,2) = -1/sqrt(2*W(2));
		X(1,3) = 1/sqrt(2*W(2));
		
		% Construct vector sequence
		XOld = X;
		for j=2:n
			X = zeros(j,j+2);
			for i=0
				X(:,i+1) = [XOld(:,i+1); 0];
			end
			for i=1:j
				X(:,i+1) = [XOld(:,i+1); -1/sqrt(2*W(j+2))]; %j+2 here because 1-indexed
			end
			for i=j+1
				X(:,i+1) = [zeros(j-1,1); 1/sqrt(2*W(j+2))];
			end
			XOld = X;
		end
		
		% Transform to correct mean and covariance
		X = repmat(mx,1,n+2) + sqrtFcn(Pxx)'*X;
		
	case 'spherical'
		
		% Compute weights
		W(1) = parameter(1);
		W(2:n+2) = (1-W(1))/(n+1);
		
		% Initialize vector sequence
		X(1,1) = 0;
		X(1,2) = -1/sqrt(2*W(2));
		X(1,3) = 1/sqrt(2*W(2));
		
		% Construct vector sequence
		XOld = X;
		for j=2:n
			X = zeros(j,j+2);
			for i=0
				X(:,i+1) = [XOld(:,i+1); 0];
			end
			for i=1:j
				X(:,i+1) = [XOld(:,i+1); -1/sqrt(j*(j+1)*W(2))];
			end
			for i=j+1
				X(:,i+1) = [zeros(j-1,1); j/sqrt(j*(j+1)*W(2))];
			end
			XOld = X;
		end
		
		% Transform to correct mean and covariance
		X = repmat(mx,1,n+2) + sqrtFcn(Pxx)'*X;
		
end

% Scale the sigma points
p = size(X,2);
X = repmat(X(:,1),1,p) + alpha*(X-repmat(X(:,1),1,p));

% Compute weights
W(1) = W(1)/alpha^2 + (1-1/alpha^2);
W(2:end) = W(2:end)/alpha^2;

% Check for beta parameter
if length(parameter)==1
	beta = 0;
else
	beta = parameter(2);
end

% Compute weighted mean and covariance
mX = Observers.wmean(X,W);
PX = Observers.wcov(X,W,alpha,beta);
