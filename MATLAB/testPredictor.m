function testPredictor
% testPredictor   Test Predictor class for Battery model
%   
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Create battery model
battery = Battery.Create;

% Set variable input profile
loads = [8; 10*60; 4; 5*60; 12; 15*60; 5; 20*60; 10; 10*60];

% Get default initalization
x0 = battery.getDefaultInitialization();

% Create predictor
predictor = Prognosis.Predictor(battery);

% Create a sample generator for the initial state
% Return normally distributed samples with mean at the initial state and
% variance defined by the process noise variance vector.
stateSampler = @(N) repmat(x0,1,N) + repmat(sqrt(battery.V),1,N).*randn(length(x0),N);

% Create sample generator for input equation parameters
% For each of the 10 load segments, sample from a uniform distribution with
% the mean given in the loads vector and the range [-1,+1] W for load and
% [-60,+60] s for the durations.
gains = ones(10,1);
gains(2:2:end) = 60;
inputParameterSampler = @(N) repmat(loads,1,N) + repmat(gains,1,N).*(rand(10,N)-0.5);

% Create process noise generator
% Use default generator built into Model class
processNoiseSampler = @battery.generateProcessNoise;

% Predict
t = 0;
horizon = 5000;
numSamples = 100;
predictor.predict(t,horizon,numSamples,stateSampler,...
    inputParameterSampler,processNoiseSampler);

% Compute actual end of discharge time, giving the exact loading parameters
battery.inputEqnHandle = @(P,t)Battery.InputEqn(P,t,loads);
T = battery.simulateToThreshold();
trueEOD = T(end);

% Plot results: since the uncertainty was centered around the mean of the
% initial state and input parameters, the distribution should be spread
% around the mean, plotted as a red line.
figure;
dischargeTimes = predictor.predictions.thresholdTimes;
histogram(dischargeTimes,'BinWidth',60);
limits = axis;
line([trueEOD trueEOD],[limits(3) limits(4)],'Color','r','Linewidth',2);
fprintf('True EOD time: %g\n',trueEOD);
fprintf('Mean predicted EOD time: %g\n',mean(dischargeTimes));
