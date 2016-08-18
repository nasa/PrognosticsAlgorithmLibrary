function testParticleFilter
% testParticleFilter   Test particle filter algorithm on Battery model
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Create battery model
battery = Battery.Create;
loads = [8; 10*60; 4; 5*60; 12; 15*60; 5; 20*60; 10; 10*60];
battery.inputEqnHandle = @(P,t)Battery.InputEqn(P,t,loads);

% Create PF
numParticles = 100;
PF = Observers.ParticleFilter(@battery.stateEqn,@battery.outputEqn,...
    battery.V,battery.N,numParticles);

% Get initial state for battery
t0 = 0;
[x0,u0,z0] = battery.getDefaultInitialization(t0);

% Initialize UKF
PF.initialize(t0,x0,u0);

% Set up output data matrices
dt = 1;
T = t0:dt:3100;
X = zeros(length(x0),length(T));
Z = zeros(length(z0),length(T));
XMean = X;
XMin = X;
XMax = X;
ZMean = Z;
ZMin = Z;
ZMax = Z;
X(:,1) = x0;
Z(:,1) = z0;
XMean(:,1) = PF.xMean();
XMax(:,1) = PF.xMax();
XMin(:,1) = PF.xMin();
ZMean(:,1) = PF.zMean();
ZMin(:,1) = PF.zMin();
ZMax(:,1) = PF.zMax();

% Initialize simulation
x = x0;
u = u0;
z = z0;

% Simulate battery and estimate state with UKF
for i=2:length(T)
    % Update state from T(i-1) to T(i)
    x = battery.stateEqn(T(i-1),x,u,0*battery.generateProcessNoise(),dt);
    % Get inputs for time T(i)
    u = battery.inputEqn(T(i));
    % Compute outputs for time T(i)
    z = battery.outputEqn(T(i),x,u,battery.generateSensorNoise());
    
    % Estimation step for UKF
    PF.estimate(T(i),u,z);
    
    % Save data
    X(:,i) = x;
    Z(:,i) = z;
    XMean(:,i) = PF.xMean();
    XMin(:,i) = PF.xMin();
    XMax(:,i) = PF.xMax();
    ZMean(:,i) = PF.zMean();
    ZMin(:,i) = PF.zMin();
    ZMax(:,i) = PF.zMax();
end

% Plot output estimates
figure;
subplot(2,1,1);
plot(T,Z(1,:),'.',T,ZMean(1,:),'--',T,ZMin(1,:),'k:',T,ZMax(1,:),'k:','LineWidth',2);
title('Temperature Estimates');
xlabel('Time (s)')
ylabel('Temperature (deg C)');
axis tight;
subplot(2,1,2);
plot(T,Z(2,:),'.',T,ZMean(2,:),'--',T,ZMin(2,:),'k:',T,ZMax(2,:),'k:','LineWidth',2);
title('Voltage Estimates');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Measured','Mean Estimated','Min Estimated','Max Estimated');
axis tight;
