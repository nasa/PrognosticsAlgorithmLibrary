function testPrognoser
% testPrognoser   Test Prognoser class for Battery model
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Create battery model
battery = Battery.Create;

% Set variable input profile
loads = [8; 10*60; 4; 5*60; 12; 15*60; 5; 20*60; 10; 10*60];

% Set up noise covariance matrices
Q = diag(battery.V);
R = diag(battery.N);

% Create UKF
UKF = Observers.UnscentedKalmanFilter(@battery.stateEqn,@battery.outputEqn,...
    Q,R,'symmetric',3-8,1);

% Create sample generator for input equation parameters
% For each of the 5 load segments, sample from a uniform distribution with
% the mean given in the loads vector and the range [-1,+1] W for load and
% [-60,+60] s for the durations.
gains = ones(length(loads),1);
gains(2:2:end) = 60;
inputParameterSampler = @(N) repmat(loads,1,N) + repmat(gains,1,N).*(rand(length(loads),N)-0.5);

% Create Prognoser
horizon = 5000;
numSamples = 100;
prognoser = Prognosis.Prognoser('model',battery,'observer',UKF,...
    'horizon',horizon,'numSamples',numSamples,...
    'stateSampler',@Observers.meanCovSampler,...
    'inputParameterSampler',inputParameterSampler,...
    'processNoiseSampler',@battery.generateProcessNoise);

% Get initial state for battery simulation
t0 = 0;
[x0,u0,z0] = battery.getDefaultInitialization(t0,loads);

% Update/initialize prognoser based on initial data
prognoser.update(t0,u0,z0);

% Set up output data matrices
dt = 1;
T = t0:dt:3100;
X = zeros(length(x0),length(T));
Z = zeros(length(z0),length(T));
XEst = X;
ZEst = Z;
X(:,1) = x0;
Z(:,1) = z0;
XEst(:,1) = UKF.x;
ZEst(:,1) = UKF.z;
EODMean = [];
EODMax = [];
EODMin = [];
predictionTimes = [];

% Initialize simulation
x = x0;
u = u0;
z = z0;

% Simulate battery and run prognoser
for i=2:length(T)
    % Update state from T(i-1) to T(i)
    x = battery.stateEqn(T(i-1),x,u,battery.generateProcessNoise(),dt);
    % Get inputs for time T(i)
    u = battery.inputEqn(T(i),loads);
    % Compute outputs for time T(i)
    z = battery.outputEqn(T(i),x,u,battery.generateSensorNoise());
    
    % Update step for prognoser
    prognoser.update(T(i),u,z);
    
    % Predict once per minute
    if mod(i-1,60)==0
        prognoser.predict();
        
        % Print some status
        fprintf('Time: %g s\n',T(i));
        battery.printOutputs(z);
        fprintf('    EOD: %g s\n',mean(prognoser.predictor.predictions.thresholdTimes));
        
        % Save some prediction data
        predictionTimes(end+1) = T(i);
        EODMean(end+1) = mean(prognoser.predictor.predictions.thresholdTimes);
        EODMin(end+1) = min(prognoser.predictor.predictions.thresholdTimes);
        EODMax(end+1) = max(prognoser.predictor.predictions.thresholdTimes);
    end
    
    % Save data
    X(:,i) = x;
    Z(:,i) = z;
    XEst(:,i) = UKF.x;
    ZEst(:,i) = UKF.z;
end

% Plot output estimates
figure;
subplot(2,1,1);
plot(T,Z(1,:),'.',T,ZEst(1,:),'--');
title('Temperature Estimates');
xlabel('Time (s)')
ylabel('Temperature (deg C)');
subplot(2,1,2);
plot(T,Z(2,:),'.',T,ZEst(2,:),'--');
title('Voltage Estimates');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Measured','Estimated');

% Compute actual end of discharge time, giving the exact loading parameters
battery.inputEqnHandle = @(P,t)Battery.InputEqn(P,t,loads);
T = battery.simulateToThreshold();
trueEOD = T(end);

% Plot prediction results
figure;
plot(predictionTimes,EODMean,'o',predictionTimes,EODMin,'o',...
    predictionTimes,EODMax,'o',predictionTimes,...
    trueEOD*ones(size(predictionTimes)),'o');
legend('Predicted EOD Mean','Predicted EOD Max','Predicted EOD Min',...
    'True EOD');
axis tight;
xlabel('Time (s)')
ylabel('EOD (s)')
