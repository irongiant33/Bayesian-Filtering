clc; clear variables; close all; 
rng(5)

parameters.numSteps = 100;
parameters.scanTime = 1;             

parameters.sigmaDrivingNoise = .1;         
parameters.sigmaMeasurementNoiseRange = 5;
parameters.sigmaMeasurementNoiseBearing = .5;
parameters.numMeasurements = 2;

parameters.startState = [0;0;0;0];
parameters.sensorPosition = [0;50];
parameters.priorCovariance = diag([100;100;1;1]);
parameters.priorMean = parameters.startState + sqrt(parameters.priorCovariance)*randn(4,1);
parameters.dimensionality = 2;
parameters.alpha = 1e-3;
parameters.beta = 2;
parameters.kappa = 0;

trueTracks = getTrueTrack(parameters);
observations = getObservations(trueTracks,parameters);

estimatedTracks = zeros(4,parameters.numSteps+1);
estimatedTracks(:,1) = parameters.priorMean;
estCov = parameters.priorCovariance;
noiseGain = [0.5*parameters.scanTime^2 0;0 0.5*parameters.scanTime^2;parameters.scanTime 0;0 parameters.scanTime];
Q = noiseGain*diag([parameters.sigmaDrivingNoise;parameters.sigmaDrivingNoise])*noiseGain';
R = diag([parameters.sigmaMeasurementNoiseRange;parameters.sigmaMeasurementNoiseBearing]);
f = @(x)[x(1)+parameters.scanTime*x(3);x(2)+parameters.scanTime*x(4);x(3);x(4)];
h = @(x)[sqrt((x(1)-parameters.sensorPosition(1))^2+(x(2)-parameters.sensorPosition(2))^2);atan2d(x(1)-parameters.sensorPosition(1),x(2)-parameters.sensorPosition(2))];
for i=2:parameters.numSteps+1
    obs = observations(:,i-1);
    [estimatedTracks(:,i),estCov]=ukf(f,estimatedTracks(:,i-1),estCov,h,obs,Q,R);
end
rmse = getError(trueTracks,estimatedTracks);

figure(1)
plot(trueTracks(1,:),trueTracks(2,:))
%axis([-100 100 -100 100])
hold on
plot(estimatedTracks(1,:),estimatedTracks(2,:))
estRmse=mean(rmse)
%obsRmse - how to measure?


