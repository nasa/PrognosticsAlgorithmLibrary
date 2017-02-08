classdef ExtendedKalmanFilter < Observers.Observer
% ExtendedKalmanFilter   Class implementing the extended Kalman filter
% algorithm
%
% This class implements the extended Kalman filter algorithm. It accepts a
% model of the explicit discrete time-invariant form:
%   x(k+1) = stateEqn(x(k),u(k),dt) + N(0,Q)
%     y(k) = outputEqn(x(k),u(k)) + N(0,R)
% where N is the standard normal distribution with zero mean and covariance
% Q or R. This class implements the Observer interface.
%
% The state and output Jacobian matrices must also be provided. These are
% functions, both taking as input x and u.
%
% ExtendedKalmanFilter Methods:
%   initialize(EKF,t0,x0,u0) - Initialize filter given initial time, state,
%   and inputs
%   estimate(EKF,t,u,z) - Update the state and outpute estimates given new
%   input and output data.
%   getStateEstimate(EKF) - Return a state estimate structure with mean and
%   covariance.
%
% See also Observers.Observer, Observers.KalmanFilter,
% Observers.UnscentedKalmanFilter
%
% Copyright (c) 2016 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% No copyright is claimed in the United States under Title 17, U.S.
% Code. All Other Rights Reserved.

	properties
		stateEqn    % State update equation (function handle)
		outputEqn   % Output equation (function handle)
		xjacobian   % State Jacobian matrix
		zjacobian   % Output Jacobian matrix
		x           % State estimate
		z           % Output estimate
		P           % State covariance matrix estimate
		Q           % Process noise covariance matrix
		R           % Sensor noise covariance matrix
    end
    
	methods
		
		function EKF = ExtendedKalmanFilter(stateEqn,outputEqn,xjacobian,zjacobian,modelVariance,sensorVariance)
			% ExtendedKalmanFilter   Constructor
            %   Construct an EKF given the state equation, output equation,
            %   state Jacobian, output Jacobian, and the process and sensor
            %   noise covariance matrices.
            EKF.stateEqn = stateEqn;
			EKF.outputEqn = outputEqn;
			EKF.xjacobian = xjacobian;
			EKF.zjacobian = zjacobian;
			EKF.Q = modelVariance;
			EKF.R = sensorVariance;
			EKF.P = zeros(size(modelVariance));
			EKF.initialized = false;
        end
		
        
		function initialize(EKF,t0,x0,u0)
            % initialize   Initialize the EKF given initial time, states,
            % and inputs
            EKF.t = t0;
            EKF.x = x0;
            EKF.u = u0;
            EKF.z = EKF.outputEqn(x0,u0);
			EKF.initialized = true;
        end
		
        
		function estimate(EKF,t,u,z)
            % estimate   Update the state estimate for the current time
            % given inputs and outputs
            
            % Ensure EKF is initialized
            if ~EKF.initialized
                error('ExtendedKalmanFilter must be initialized first!');
            end
            
			deltaT = newT-EKF.t;
			
            % Predict
            Fk = EKF.xjacobian(EKF.x,EKF.u);
            xkk1 = EKF.stateEqn(EKF.x,EKF.u,deltaT);
            Pkk1 = Fk*EKF.P*Fk'+EKF.Q;
			
            % Update
			Hk = EKF.zjacobian(xkk1,u);
			yk = z-EKF.outputEqn(xkk1,u);
			Sk = Hk*Pkk1*Hk'+EKF.R;			
			Kk = Pkk1*Hk'/Sk;
			EKF.x = xkk1+Kk*yk;
			EKF.P = (eye(length(EKF.x))-Kk*Hk)*Pkk1;
			EKF.z = EKF.outputEqn(EKF.x,u);
			
			% Update time, inputs
            EKF.u = u;
            EKF.t = t;
        end
        
        
        function stateEstimate = getStateEstimate(EKF)
            % getStateEstimate  Return a structure with mean and covariance
            %   Returns the current state estimate as a structure with
            %   fields 'mean' and 'covariance'.
            stateEstimate.mean = EKF.x;
            stateEstimate.covariance = EKF.P;
        end
		
	end
end
		
