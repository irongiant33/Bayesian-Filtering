clc; clear variables; close all;

%generate random sequence consistently.
%comment out for more "unpredictable" scenario
%rng(50);

%how long the test lasts
parameters.numSteps = 100;
%how frequently you make observations
parameters.scanTime = 1; 

%maximum surveillance area
parameters.limits = [-1000 1000;-1000 1000; -10 10; -10 10];

%maneuverability of the target (uncertainty in target position)
parameters.sigmaDrivingNoise = .1;         
%uncertainty in measurement
parameters.sigmaMeasurementNoise = 1;

parameters.detectionProbability = .8 ;
parameters.meanClutter = 3;
parameters.measurementRange = 200;


%xPos, yPos, xVel, yVel
parameters.startState = [0;0;0;0];
%xPos, yPos of sensor
parameters.sensorPosition = [0;50];
%initial uncertainty of xPos,yPos,xVel,yVel of where target is
parameters.priorCovariance = diag([100;100;5;5]);
%initial belief of target state
parameters.priorMean = parameters.startState + sqrt(parameters.priorCovariance)*randn(4,1);


trueTracks = getTrueTrack(parameters);
correctObservations = getObservations(trueTracks,parameters);
observations = getClutteredObservations(correctObservations,parameters);

estimatedTracks = pdaf(observations,parameters);
rmse = getError(trueTracks,estimatedTracks);

figure(1)
%plot XY of true tracks
plot(trueTracks(1,:),trueTracks(2,:))
axis([-100 100 -100 100])
hold on
%plot XY of estimated tracks to compare
plot(estimatedTracks(1,:),estimatedTracks(2,:))


estRmse = mean(rmse)