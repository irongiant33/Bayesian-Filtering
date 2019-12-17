clc; clear variables; close all; dbstop if error

%generate random sequence consistently.
%comment out for more "unpredictable" scenario
rng(1);

parameters.numSteps = 33;   % how long the test lasts
parameters.scanTime = 1;   % duration between time steps
parameters.limits = [-1000 1000;-1000 1000; -10 10; -10 10];   % limits of surveillance area
parameters.sigmaDrivingNoise = .1;   % uncertainty in target motion (maneuverability of the target)        

parameters.detectionProbability = [.9 .5];
parameters.meanClutter = [10 7];
parameters.numDevices = 3; % number of devices in the network
parameters.networkMode = 1; % 0 if networking not enabled. 1,2,3 if networking enabled. 1 if device 1 transmits first, 2 if 2 transmits first, etc.
parameters.measurementRange = 200;
parameters.velocityRange = 20;
parameters.sigmaMeasurementNoiseRange = 1;
parameters.sigmaMeasurementNoiseVelocity = 1;
parameters.sigmaMeasurementNoiseBearing = 1;

parameters.startState = [0;0;0;0];   % xPos, yPos, xVel, yVel
parameters.sensorPosition = [0 50;25 25;-25 25];   % xPos, yPos of sensor
parameters.priorCovariance = diag([100;100;5;5]);   % initial uncertainty of xPos,yPos,xVel,yVel of where target is
parameters.priorMean = parameters.startState + sqrt(parameters.priorCovariance)*randn(4,1);   % initial belief of target state

parameters.birthProb = 0.01;   %probability that the target appears at an individual time step
parameters.survivalProb = 0.999;   %depents on how long a target is expected to stay in the surveillance area
parameters = getBirthComponents(parameters,100,10);

parameters.detectionThreshold = 0.5;
parameters.maxComponents = 900;

% generate track
apperanceFromTo = [1;70];
trueTracks = getTrueTrack(parameters,apperanceFromTo);

% generate observations
correctObservations = getObservations(trueTracks,parameters);
observations = getClutteredObservations(correctObservations,parameters);

% run Bernoulli filter
[estimatedTracks,existences] = bernoulliFilter(observations,parameters);

% calculate error
rmse = getError(trueTracks,estimatedTracks);
meanRmse = mean(rmse(~isnan(rmse)));

% show results
showDynamic = false;
plotOverTime(trueTracks,observations,estimatedTracks,showDynamic);

