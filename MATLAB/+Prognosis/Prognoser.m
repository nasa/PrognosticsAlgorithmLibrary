classdef Prognoser < handle
% Prognoser   Class for performing model-based prognosis.
%
% This class defines an object for performing model-based prognosis. It
% contains a system model, an Observer, a Predictor, and functions for
% generating samples for prediction. The system model must be defined using
% the PrognosticsModel class, i.e., it requires a state equation, output
% equation, input equation, and threshold equation.
%
% Prognoser Methods:
%   initialize - Initialize the prognoser with initial time, state, inputs
%   update - Update the prognoser and its state estimate with new data
%   predict - Predict the system evolution given the state estimate
%
% See also Model.PrognosticsModel, Observers.Observer, Prognosis.Predictor
%
% Copyright (c) 2016 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% No copyright is claimed in the United States under Title 17, U.S.
% Code. All Other Rights Reserved.

	properties
        t;                      % Last updated time
		model;                  % Model object for system under prognosis
        observer;               % Observer object for state estimation
        predictor;              % Predictor object for prediction
        numSamples;             % Number of samples for prediction
        horizon;                % Time horizon for prediction
        stateSampler;           % Sample generator for states given estimate
        inputParameterSampler;  % Sample generator for input parameters
        processNoiseSampler;    % Sample generator for process noise
        initialized = false;    % Initialization flag
	end
	
	methods
        
        function P = Prognoser(varargin)
            % Prognoser   Constructor
            % Construct a prognoser given various arguments. The arguments
            % are string, value pairs, e.g., Prognoser('model',model,...).
            % The following arguments are required:
            %   - model: the PrognosticsModel object
            %   - observer: the Observer object
            %   - numSamples: number of samples to use for prediction
            %   - horizon: time horizon to use for prediction
            %   - stateSampler: function to sample the states, given a
            %     state estimate structure and the number of samples to
            %     generate
            %   - inputParameterSampler: function to sample the input
            %     parameters
            %   - processNoiseSampler: function to sample the process noise
            %
            % The Predictor object is constructed automatically from the
            % given model, and does not need to be provided as input.
            
            % Create a structure from the string, value pairs
            args = struct(varargin{:});
            
            % Set object properties from the function arguments
            properties = fieldnames(args);
            for i=1:length(properties)
                P.(properties{i}) = args.(properties{i});
            end
            
            % Check for required properties
            requiredProperties = {'model' 'observer' 'numSamples'...
                'horizon' 'stateSampler' 'inputParameterSampler',...
                'processNoiseSampler'};
            for i=1:length(requiredProperties)
                if isempty(P.(requiredProperties{i}))
                    error('%s is a required property!',requiredProperties{i});
                end
            end
            
            % Create predictor from model
            P.predictor = Prognosis.Predictor(P.model);
        end
        
        
        function initialize(P,t,x,u)
            % initialize   Initialize the prognoser
            %    Initialize the prognoser given initial time, state, and
            %    inputs. This function will be called automatically the
            %    first time the update method is called, if the model
            %    contains an initialize equation, i.e., function that takes
            %    as input the current time, inputs, and outputs, and
            %    derives the state. If such a function is not available,
            %    then the initialize method should be directly called
            %    before update is called, otherwise an error will be
            %    thrown.
            
            % Initialize the time
            P.t = t;
            % Initialize the observer
            P.observer.initialize(t,x,u);
            % Set initialization flag
            P.initialized = true;
        end
        
        
        function update(P,t,u,z)
            % update   Update the prognoser and its state estimate
            %   Update the prognoser with new data. If not initialized,
            %   this function first attempts to initialize the prognoser
            %   using the given data through the model's initialize
            %   equation. If not available, an error is thrown. This
            %   function updates the internal prognoser time and calls the
            %   estimate method on the observer.
            
            % If not initialized, initialize based on initial data. This
            % requires an initialization equation.
            if ~P.initialized
                % Initialize using initialization equation
                if ~isempty(P.model.initializeEqnHandle)
                    x = P.model.initializeEqn(t,u,z);
                    P.initialize(t,x,u);
                    return
                else
                    error('Prognoser is not initialized, and initialize equation is not present!');
                end
            end
            
            % Update time
            P.t = t;
            
            % Estimate state
            P.observer.estimate(t,u,z);
        end
        
        
        function predict(P)
            % predict   Predict the state evolution up to a threshold.
            %   This method takes no arguments, as all the parameters for
            %   prediction are set up via the constructor. This method
            %   first constructs a sample generator for the initial state
            %   to use for prediction, based on the current state estimate
            %   from the observer. It then calls the predict method for the
            %   predictor object.
            
            % Here, the state sampler must be designed to take as input
            % state estimates from the observer and the number of samples
            % to generate.
            % For the predictor, only number of samples can be an argument,
            % so we construct the stateSampler as follows.
            stateSampler = @(N)P.stateSampler(P.observer.getStateEstimate(),N);
            
            % Now call the predict function for the predictor
            P.predictor.predict(P.t,P.horizon,P.numSamples,stateSampler,...
                P.inputParameterSampler,P.processNoiseSampler);
        end
        
	end
end
