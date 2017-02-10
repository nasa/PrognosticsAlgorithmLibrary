function testKalmanFilter
% tesKalmanFilter   Test Kalman filter class with a simple model.
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Create state-space model with two states, two inputs, and one output
t0 = 0;
x0 = [5; 5];
u0 = [0; 0];
A = [-1 2; 2 -4]/10;
B = [1 0; 0 0];
C = [3 1];
D = [0 0];
Q = [1e-4 0; 0 2e-4];
R = [2e-2];

% Create Kalman filter
KF = Observers.KalmanFilter(A,B,C,D,Q,R);

% Initialize the filter - provide incorrect initial state to ensure that
% the estimate converges
KF.initialize(t0,[0; 0],u0);

% Initialize the system simulation
x = x0;
u = u0;
z = C*x + D*u;

% Set up some output data to assess filter performance
T = t0:1:30;
X = zeros(length(x),length(T));
Z = zeros(length(z),length(T));
XEst = X;
ZEst = Z;

% Set initial data
X(:,1) = x;
Z(:,1) = z;
XEst(:,1) = [0; 0];
ZEst(:,1) = C*[0; 0] + D*u0;

% Simulate the system, and generate state and output estimates using the
% Kalman filter
for i=2:length(T)
    % Simulate data for the new time step
    x = A*x + B*u + chol(Q)*randn(length(x),1);
    u = u + 0.01;
    z = C*x + D*u + chol(R)*randn(length(z),1);
    t = T(i);
    
    % Run an estimation step for the KF
    KF.estimate(t,u,z);
    
    % Save the data
    X(:,i) = x;
    Z(:,i) = z;
    XEst(:,i) = KF.x;
    ZEst(:,i) = KF.z;
end

% Plot the results
figure;
plot(T,X,'.',T,XEst,'--');
title('State Estimates');
xlabel('Time');
legend('Simulated State 1', 'Simulated State 2',...
    'Estimated State 1','Estimated State 2');
axis tight;

figure;
plot(T,Z,'.',T,ZEst,'--');
title('Output Estimates');
xlabel('Time');
legend('Simulated Output','Estimated Output');
axis tight;
