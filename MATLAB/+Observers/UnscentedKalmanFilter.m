classdef UnscentedKalmanFilter < Observers.Observer
% UnscentedKalmanFilter   Class implementing unscented Kalman filter
%
% This class implements the unscented Kalman filter algorithm. It accepts a
% model of the explicit discrete time-variant form:
%   x(t+dt) = stateEqn(t,x(t),u(t),noise,dt)
%      y(t) = outputEqn(t,x(t),u(t),noise)
% where process and sensor noise are defined by covariance matrices Q and
% R.
%
% State and output equations must be defined in a vectorized form, i.e., so
% that they can take several samples of state, input, output, and noise.
% Matrices must be formed such that the row represents the variable and the
% column the sample.
%
% UnscentedKalmanFilter Methods:
%   initialize(UKF,t0,x0,u0) - Initialize filter given initial time, state,
%   and inputs
%   estimate(UKF,t,u,z) - Update the state and outpute estimates given new
%   input and output data.
%   getStateEstimate(UKF) - Return a state estimate structure with mean and
%   covariance. 
%
% See also Observers.Observer, Observers.KalmanFilter,
% Observers.ExtendedKalmanFilter, Observers.ParticleFilter,
% Observers.computeSigmaPoints
%
% Copyright (c) 2016 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% No copyright is claimed in the United States under Title 17, U.S.
% Code. All Other Rights Reserved.

    properties
        stateEqn;     % State equation function handle
        outputEqn;    % Output equation function handle
        x;            % State estimate
        z;            % Output estimate
        P;            % State covariance matrix estimate
        Q;            % Process noise covariance matrix
        R;            % Sensor noise covariance matrix
        method;       % Unscented transform method (string)
        parameters;   % Unscented transform free parameters
        alpha = 1;    % Unscented transform scaling parameter
        beta = 0;     % Unscented transform scaling parameter
        sigmaPoints;  % Sigma point struct array
    end
    
    methods
        
        function UKF = UnscentedKalmanFilter(stateEqn,outputEqn,modelVariance,sensorVariance,method,parameters,alpha,beta)
            % UnscentedKalmanFilter   Constructor
            %   Construct a UKF given a state equation, output equation,
            %   process and sensor noise covariance matrices, a sigma point
            %   selection method (a string), sigma point parameters (a
            %   vector), and (optional) scaling parameters alpha and beta.
            
            % Set state and output equations, noise covariances
            UKF.stateEqn = stateEqn;
            UKF.outputEqn = outputEqn;
            UKF.Q = modelVariance;
            UKF.R = sensorVariance;
            
            % Set unscented transform parameters
            UKF.method = method;
            UKF.parameters = parameters;
            if nargin>=7
                UKF.alpha = alpha;
            end
            if nargin==8
                UKF.beta = beta;
            end
            
            % Set initialization flag
            UKF.initialized = false;
        end
        
        
        function initialize(UKF,t0,x0,u0)
            % initialize   Initialize the UKF with initial time, state, and
            % inputs
            
            % Set initial values
            UKF.t = t0;
            UKF.x = x0;
            UKF.z = UKF.outputEqn(t0,x0,u0,0);
            UKF.u = u0;
            UKF.P = UKF.Q;
            
            % Intialize sigma points
            [X,W] = Observers.computeSigmaPoints(UKF.x,UKF.P,UKF.method,UKF.parameters,UKF.alpha);
            UKF.sigmaPoints.x = X;
            UKF.sigmaPoints.w = W;
            UKF.sigmaPoints.z = UKF.outputEqn(t0,X,u0,0);
            
            % Set filter to initialized
            UKF.initialized = true;
        end
        
        
        function estimate(UKF,t,u,z)
            % estimate   Update the state estimate given new time, inputs,
            % and outputs
            
            % Ensure UKF is initialized
            if ~UKF.initialized
                error('UnscentedKalmanFilter must be initialized first!');
            end
            
            dt = t-UKF.t;
            
            % 1. Predict
            xa = UKF.x;
            Pa = UKF.P;
            
            % Compute sigma points
            [X,W] = Observers.computeSigmaPoints(xa,Pa,UKF.method,UKF.parameters,UKF.alpha);
            
            % Propagate sigma points through transition function
            Xkk1 = X;
            Xkk1(1:length(UKF.x),:) = UKF.stateEqn(UKF.t,X(1:length(UKF.x),:),UKF.u,0,dt);
            
            % Recombine weighted sigma points to produce predicted state
            % and covariance
            xkk1 = Observers.wmean(Xkk1,W);
            Pkk1 = Observers.wcov(Xkk1,W,UKF.alpha,UKF.beta) + UKF.Q;
            
            % propagate sigma points through observation function
            Zkk1 = UKF.outputEqn(t,Xkk1(1:length(UKF.x),:),u,0);
            
            % Recombine weighted sigma points to produce predicted
            % measurement and covariance
            zkk1 = Observers.wmean(Zkk1,W);
            Pzz = Observers.wcov(Zkk1,W,UKF.alpha,UKF.beta) + UKF.R;
            
            % 2. Update
            
            % Compute state-measurement cross-covariance matrix
            Pxz = 0;
            for i=1:length(W)
                Pxz = Pxz + W(i)*(Xkk1(:,i)-xkk1)*(Zkk1(:,i)-zkk1)';
            end
            Pxz = Pxz + (1-UKF.alpha^2+UKF.beta)*(Xkk1(:,1)-xkk1)*(Zkk1(:,1)-zkk1)';
            
            % Compute kalman gain
            Kk = Pxz/Pzz;
            
            % Compute state and measurement estimates
            xk1 = xkk1 + Kk*(z-zkk1);
            UKF.x = xk1(1:length(UKF.x));
            UKF.z = UKF.outputEqn(t,UKF.x,u,0);
            Pk1 = Pkk1 - Kk*Pzz*Kk';
            UKF.P = Pk1(1:length(UKF.x),1:length(UKF.x));
            
            % Compute sigma points
            [X,W] = Observers.computeSigmaPoints(UKF.x,UKF.P,UKF.method,UKF.parameters,UKF.alpha);
            UKF.sigmaPoints.x = X;
            UKF.sigmaPoints.w = W;
            UKF.sigmaPoints.z = UKF.outputEqn(t,X,u,0);
            
            % Update time and inputs
            UKF.t = t;
            UKF.u = u;
        end
        
        
        function stateEstimate = getStateEstimate(UKF)
            % getStateEstimate  Return a structure with mean and covariance
            %   Returns the current state estimate as a structure with
            %   fields 'mean' and 'covariance'.
            stateEstimate.mean = UKF.x;
            stateEstimate.covariance = UKF.P;
        end
        
        
        function x = xMean(UKF)
            % xMean   Return mean state estimate
            x = UKF.x;
        end
        
        
        function x = xMax(UKF)
            % xMax   Return the maximum state in the sigma points
            x = max(UKF.sigmaPoints.x,[],2);
        end
        
        
        function x = xMin(UKF)
            % xMin   Return the minimum state in the sigma points
            x = min(UKF.sigmaPoints.x,[],2);
        end
        
        
        function z = zMean(UKF)
            % zMean   Return the mean output estimate
            z = UKF.z;
        end
        
        
        function z = zMax(UKF)
            % zMax   Return the maximum output in the sigma points
            z = max(UKF.sigmaPoints.z,[],2);
        end
        
        
        function z = zMin(UKF)
            % zMin   Return the minimum output in the sigma points
            z = min(UKF.sigmaPoints.z,[],2);
        end
        
        
        function v = xCov(UKF)
            % xCov   Return the weighted state covariance from the sigma
            % points
            v = Observers.wcov(UKF.sigmaPoints.x,UKF.sigmaPoints.w,UKF.alpha,UKF.beta);
        end
        
        
        function v = zCov(UKF)
            % zCov   Return the weighted output covariance from the sigma
            % points
            v = Observers.wcov(UKF.sigmaPoints.z,UKF.sigmaPoints.w,UKF.alpha,UKF.beta);
        end
        
    end
end
