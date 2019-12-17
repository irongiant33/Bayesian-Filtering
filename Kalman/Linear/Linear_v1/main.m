clc; clear variables; close all; 

%generate random sequence consistently.
%comment out for more "unpredictable" scenario
%rng(5);

parameters.numSteps = 100;
parameters.scanTime = 1;             

parameters.sigmaDrivingNoise = .1;         
parameters.sigmaMeasurementNoise = 1;

parameters.startState = [0;0;0;0];
parameters.priorCovariance = diag([100;100;1;1]);
parameters.priorMean = parameters.startState + sqrt(parameters.priorCovariance)*randn(4,1);

trueTracks = getTrueTrack(parameters);
observations = getObservations(trueTracks,parameters);

%estimatedTracks = kalmanFilterNew(observations,parameters);
%estimated tracks with only two measurements:
estimatedTracks = kalmanFilter(observations,parameters);
rmse1 = getError(trueTracks,estimatedTracks);
rmse2 = getError(trueTracks,observations);

figure(1)
%plot XY of true tracks
plot(trueTracks(1,:),trueTracks(2,:))
axis([-100 100 -100 100])
hold on
%plot XY of estimated tracks to compare
plot(estimatedTracks(1,:),estimatedTracks(2,:))
estRmse=mean(rmse1)
obsRmse=mean(rmse2)


