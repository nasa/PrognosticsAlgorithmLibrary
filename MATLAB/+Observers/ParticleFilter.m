classdef ParticleFilter < Observers.Observer
% ParticleFilter   Class implementing the SIR particle filter
%
% This class implements the sampling-importance-resampling particle filter
% algorithm. It acccepts a model of the explicit discrete time-variant
% form:
%   x(t+dt) = stateEqn(t,x(t),u(t),noise,dt)
%     y(t) = outputEqn(t,x(t),u(t),noise)
%
% State and output equations must be defined in a vectorized form, i.e., so
% that they can take several samples of state, input, output, and noise.
% Matrices must be formed such that the row represents the variable and the
% column the sample.
%
% Assumptions:
% Currently assumes only Gaussian noise, specified using process and sensor
% noise variance vectors (so noise for each state/output is generated
% independently).
%
% Limitations:
% Calculation of likelihood is done in a way found to be most efficient,
% but uses a lot of memory, so this method will have to be changed if a
% very large number of particles are used (like 10,000).
%
% This class implements the Observer interface.
%
% ParticleFilter Methods:
%   initialize(PF,t0,x0,u0) - Initialize filter given initial time, state,
%   and inputs
%   estimate(PF,t,u,z) - Update the state and output estimates given new
%   input and output data.
%   getStateEstimate(PF) - Return a state estimate structure with samples
%   and weights fields
%
% Copyright (c) 2016 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% No copyright is claimed in the United States under Title 17, U.S.
% Code. All Other Rights Reserved.

	properties
		stateEqn;        % Handle to state equation
		outputEqn;       % Handle to output equation
		n;               % Sensor noise variance vector
		v;               % Process noise variance vector
		numParticles;    % Number of particles used by the filter
		particles = [];  % Particles structure
		minNEff;         % Effective particle number resampling threshold
        initGain = 1;    % Gain applied to sampling for particles o% initial states
	end
	
	methods
		
		function PF = ParticleFilter(stateEqn,outputEqn,processVariance,sensorVariance,numParticles,varargin)
			% ParticleFilter   Constructor
            %   Construct a particle filter given the state and output
            %   equations, process and sensor noise variance vectors, the
            %   number of particles, and additional optional arguments. The
            %   optional arguments are as follows:
            %   - minNEff: minimum number of effective particles before
            %     resampling. This defaults to numParticles/3
            %   - initGain: multiplcative gain to apply to initial state
            %     sampling. It samples using the process noise, which is
            %     then multiplied by this gain term.
            
            % Set up PF properties
            PF.stateEqn = stateEqn;
			PF.outputEqn = outputEqn;
			PF.v = processVariance;
			PF.n = sensorVariance;
			PF.numParticles = numParticles;
			
            % Set additional optional properties
            args = struct(varargin{:});
            PF.minNEff = numParticles/3;
            properties = fieldnames(args);
            for i=1:length(properties)
                PF.(properties{i}) = args.(properties{i});
            end
        end
		
        
		function noise = generateProcessNoise(PF)
            % generateProcessNoise   Generate process noise samples
            %   Sample from process noise distribution, which currently
            %   assumes independent Gaussian
			noise = sqrt(diag(PF.v))*randn(length(PF.v),PF.numParticles);
        end
		
        
		function noise = generateSensorNoise(PF)
			% generateSensorNoise   Generate sensor noise samples
            %   Sample from sensor noise distribution, which currently
            %   assumes independent Gaussian
			noise = sqrt(diag(PF.n))*randn(length(PF.n),PF.numParticles);
        end
		
        
		function initialize(PF,t0,x0,u0)
            % initialize   Initialize particle filter with initial time,
            % states, and inputs
            
			% Generate initial particle population
			x0 = repmat(x0,1,PF.numParticles);
			PF.particles.x = x0 + PF.initGain*PF.generateProcessNoise();
			PF.particles.z = PF.outputEqn(t0,PF.particles.x,u0,0);
			PF.particles.w = PF.likelihood(PF.outputEqn(t0,x0,u0,0),PF.particles.z);
			
            % Normalize weights
			PF.normalize();
            
            % Set time and inputs
            PF.t = t0;
            PF.u = u0;
            
            % Set initialization flag
            PF.initialized = true;
        end
		
        
		function normalize(PF)
            % normalize   Normalize the particle weights
			sumWeights = sum(PF.particles.w);
			PF.particles.w = PF.particles.w./sumWeights;
        end
		
        
		function l = likelihood(PF,zActual,zPredicted)
            % likelihood   Compute the likelihood of given observations for
            % predicted observations.
            
			% Note that zActual and zPredicted are matrices; diag is used
			% to extract correct multiplications. This is a more efficient
			% alternative to looping through each particle individually.
			l = 1/sqrt((2*pi)^size(zActual,1)*det(diag(PF.n))) * exp(-0.5*diag((zActual-zPredicted)'/diag(PF.n)*((zActual-zPredicted))));
        end
		
        
		function resample(PF)
            % resample   Resample the particles if the minimum effective
            % number of particles is below the threshold
            
			% Compute effective sample size
			nEff = PF.effectiveN();
			% Resample if nEff below threshold (systematic resampling)
			if nEff < PF.minNEff
				PF.systematicResample();
			end
        end
		
        
		function nEff = effectiveN(pf)
            % effectiveN   Compute the effective number of particles
			nEff = 1/sum(pf.particles.w.^2);
		end
			
		
		function estimate(PF,t,u,z)
            % estimate   Update the state estimate for the current time
            % given new inputs and outputs
            
			dt = t-PF.t;
			
			% Generate particles
            PF.particles.x = PF.stateEqn(PF.t,PF.particles.x,PF.u,...
                PF.generateProcessNoise(),dt);
            PF.particles.z = PF.outputEqn(t,PF.particles.x,u,0);
            
            % Compute likelihood
            lh = PF.likelihood(repmat(z,1,PF.numParticles),PF.particles.z);
            lh = lh + 1e-99;  % Enforce nonzero likelihood
            
            % Set particle weights, normalize and resample
            PF.particles.w = PF.particles.w.*lh;
            PF.normalize();
            PF.resample();
			
			% Update t and u
            PF.t = t;
            PF.u = u;
        end
		
        
		function systematicResample(pf)
            % systematicResample   Systematic resampling algorithm
            %   Resamples the particles to be distributed around the
            %   higher-weight particles, to increase the effective number
            %   of particles and reduce degeneracy.
            
			% Preallocate newParticles
			newParticles = pf.particles;
			
            % Construct CDF
			cs = cumsum(pf.particles.w);
			i = 1;
			
            % Draw starting point from U[0,1/Ns]
			u1 = rand/pf.numParticles;
			for j=1:pf.numParticles
				% Move along CDF
				u = u1+(j-1)/pf.numParticles;
				while u > cs(i)
					i = i+1;
				end
				% Reassign particle
				newParticles.x(:,j) = pf.particles.x(:,i);
				newParticles.z(:,j) = pf.particles.z(:,i);
			end
			newParticles.w = ones(size(pf.particles.w))/pf.numParticles;
			
            % Copy over resampled particles
			pf.particles = newParticles;
        end
        
        
        function stateEstimate = getStateEstimate(PF)
            % getStateEstimate  Return a structure with samples and weights
            %   Returns the current state estimate as a structure with
            %   fields 'samples' and 'weights'.
            stateEstimate.samples = PF.particles.x;
            stateEstimate.weights = PF.particles.w;
        end
		
        
		function m = xMean(pf)
			% xMean   Calculate mean through weighted sum of particles
			m = pf.particles.x*pf.particles.w;
        end
		
        
		function m = zMean(pf)
			% zMean   Calculate mean through weighted sum of particles
			m = pf.particles.z*pf.particles.w;
        end
		
        
		function x = xBest(pf)
			% xBest   Find particle with highest weight
			[~,i] = max(pf.particles.w);
			x = pf.particles.x(:,i);
        end
		
        
		function z = zBest(pf)
			% zBest   Find particle with highest weight
			[~,i] = max(pf.particles.w);
			z = pf.particles.z(:,i);
        end
		
        
		function x = xMax(pf)
			% xMax   Get max values for each xi in particles
			x = max(pf.particles.x,[],2);
        end
		
        
		function x = xMin(pf)
			% xMin   Get min values for each xi in particles
			x = min(pf.particles.x,[],2);
        end
		
        
		function z = zMax(pf)
			% zMax   Get max values for each zi in particles
			z = max(pf.particles.z,[],2);
        end
		
        
		function z = zMin(pf)
			% zMin   Get max values for each zi in particles
			z = min(pf.particles.z,[],2);
		end
		
	end
end
