classdef Predictor < handle
% Predictor   Class for implementing model-based prediction using sampling.
%
% This class defines an obect for performing model-based state prediction
% using sampling methods. A Predictor is constructed using a
% PrognosticsModel object, i.e., it must have a state equation, input
% equation, and threshold equation defined. The states are simulated until
% either a specified time horizon is met, or the threshold is reached, as
% defined by the threshold equation. The input equation is used to compute
% the inputs to the system at any given time point. The input equation
% should be defined to take in addition to the time, a set of "input
% parameters" that are used to determine what the input should be at a
% given time.
%
% The predict method takes the starting time of prediction, the horizon to
% predict to (i.e., it simulates to starting time + horizon), and functions
% to generate samples for the state at the time of prediction, the input
% parameters, and process noise. If horizon is inf, then prediction will
% continue until all samples hit the threshold. This is not recommended
% unless it is guaranteed that each sample will reach the threshold.
%
% The sample generator functions must be defined to each take only one
% argument, the number of samples to generate. This means that whatever
% information is needed to correctly generate samples should be built into
% these functions, such as distribution parameters.
%
% Predictor Properties:
%   predictions - Data structure containing results of the predict method.
%
% Predictor Methods:
%   predict - Predict the state evolution up to a threshold.
%
% See also Prognosis.Prognoser, Model.PrognosticsModel
%
% Copyright (c) 2016 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% No copyright is claimed in the United States under Title 17, U.S.
% Code. All Other Rights Reserved.
    
    properties
        model;                % PrognosticsModel object
        predictions;          % Prediction data structure
    end
    
    methods
        
        function P = Predictor(model)
            % Predictor   Constructor
            %   Construct a Predictor object given a PrognosticsModel
            %   object.
			P.model = model;
		end
        
        % Perform prediction given samples of states and input equations,
        % up to the threshold equation or a time horizon.
        function predict(P,t0,horizon,numSamples,stateSampler,...
                inputParameterSampler,processNoiseSampler)
            % predict   Predict the state evolution up to a threshold.
            %   The inputs to this function are as follows:
            %   - t0: starting time of prediction
            %   - horizon: time horizon of prediction (predicts to
            %     t0+horizon)
            %   - numSamples: number of samples for which to do prediciton
            %   - stateSampler: function that generates samples for the
            %     initial state at the start of prediction
            %   - inputParameterSampler: function that generates samples
            %     for the input parameters, which define possible future
            %     input trajectories for all future time points
            %   - processNoiseSampler: function that generates samples for
            %     the process noise. A default function for this purpose is
            %     available in the PrognosticsModel class
            % 
            % There are no outputs to this function, the prediction results
            % are stored in the class in the predictions property. The
            % following fields are computed:
            %    - t: the time of prediction
            %    - horizon: the time horizon used for prediction
            %    - numSamples: the number of samples used for prediction
            %    - thresholdReached: a logical row vector, with true/false
            %      for each sample, indicating whether the sample hit the
            %      threshold or not
            %    - thresholdTimes: a row vector, with for each sample the
            %      time that trajectory first reached the threshold. For
            %      any samples that did not, the value is inf.
            %    - states: a matrix with rows for the states and columns
            %      for the samples, capturing for each trajectory the value
            %      of the states at the time the threshold was first
            %      reached. If that trajectory did not reach the threshold,
            %      the values are all inf.
            %    - inputs: a matrix with rows for the inputs and columns
            %      for the samples, capturing for each trajectory the value
            %      of the inputs at the time the threshold was first
            %      reached. If that trajectory did not reach the threshold,
            %      the values are all inf.
            % 
            % Note that for a given sample trajectory and its threshold
            % time, states, and inputs, we can compute other values at that
            % time that are functions of these (e.g., the outputs).
            
            % Initialize predictions data structure with prediction
            % parameters
            P.predictions.t = t0;
            P.predictions.horizon = horizon;
            P.predictions.numSamples = numSamples;
            
            % Sample the initial states
            X = stateSampler(numSamples);
            
            % Sample the input parameters
            inputParameters = inputParameterSampler(numSamples);
            
            % Compute initial inputs
            U = P.model.inputEqn(t0,inputParameters);
            
            % Initialize prediction data
            P.predictions.thresholdReached = P.model.thresholdEqn(t0,X,U);
            P.predictions.thresholdTimes = inf(1,numSamples);
            P.predictions.thresholdTimes(P.predictions.thresholdReached) = t0;
            P.predictions.states = inf(size(X));
            P.predictions.inputs = inf(size(U));
            P.predictions.states(:,P.predictions.thresholdReached) = X(:,P.predictions.thresholdReached);
            P.predictions.inputs(:,P.predictions.thresholdReached) = U(:,P.predictions.thresholdReached);
            
            % Simulate all samples until threshold or horizon reached
            t = t0;
            dt = P.model.sampleTime;
            while sum(~P.predictions.thresholdReached) && t<t0+horizon
                % Sample process noise
                V = processNoiseSampler(numSamples);
                % Update state from time t to t+dt
                X = P.model.stateEqn(t,X,U,V,dt);
                % Increment time
                t = t + dt;
                % Get inputs for new time
                U = P.model.inputEqn(t,inputParameters);
                % Check threshold
                P.predictions.thresholdReached = P.predictions.thresholdReached | P.model.thresholdEqn(t,X,U);
                % Update threshold times, for those that have changed
                logicalIndices = P.predictions.thresholdTimes==inf & P.predictions.thresholdReached;
                P.predictions.thresholdTimes(logicalIndices) = t;
                % Update states and inputs at time of threshold crossing
                P.predictions.states(:,logicalIndices) = X(:,logicalIndices);
                P.predictions.inputs(:,logicalIndices) = U(:,logicalIndices);
            end
            
            % Compute probability of reaching threshold within horizon
            numReached = length(find(P.predictions.thresholdReached));
            P.predictions.probability = numReached/numSamples;
            
        end
        
    end
    
end
