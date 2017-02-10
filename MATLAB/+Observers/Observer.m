classdef Observer < handle
% Observer   Interface class for state observer objects
% Abstract class for creating Observer objects that perform model-based
% state estimation. Subclass constructors should accept model equations
% and/or parameters.
%
% Observer Methods:
%   initialize(obj,t0,x0,u0) - Initialize filter given initial time, state,
%   and inputs
%   estimate(obj,t,u,z) - Update the state and outpute estimates given new
%   input and output data.
%   getStateEstimate(obj) - Return a state estimate structure with mean and
%   covariance. 
%
% See also Observers.KalmanFilter, Observers.ExtendedKalmanFilter,
% Observers.ParticleFilter, Observers.UnscentedKalmanFilter
%
% Copyright (c) 2016 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% No copyright is claimed in the United States under Title 17, U.S.
% Code. All Other Rights Reserved.

    properties
        t;            % Last updated time
        initialized;  % Flag describing whether observer is initialized
        u;            % Input at time t
    end
    
    methods (Abstract)
        
        % Given initial state and inputs, initialize observer
        initialize(obj,t,x,u)
        
        % Perform an estimation step given new time, inputs, and outputs
        estimate(obj,t,u,z)
        
        % Get state estimate
        getStateEstimate(obj)
    end
    
end
