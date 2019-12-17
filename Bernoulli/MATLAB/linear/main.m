clc; clear variables; close all;

%generate random sequence consistently.
%comment out for more "unpredictable" scenario
rng(50);

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

parameters.detectionProbability = .99;
parameters.meanClutter = 0.1;
parameters.measurementRange = 200;

%xPos, yPos, xVel, yVel
parameters.startState = [0;0;0;0];
%xPos, yPos of sensor
parameters.sensorPosition = [0;50];
%initial uncertainty of xPos,yPos,xVel,yVel of where target is
parameters.priorCovariance = diag([100;100;5;5]);
%initial belief of target state
parameters.priorMean = parameters.startState + sqrt(parameters.priorCovariance)*randn(4,1);

parameters.birthProb = 0.01; %how often you expect a target to show up. 
parameters.survivalProb = 0.999; %how long you expect a target to stay in the field of view
parameters.priorExistence = 0; %probability that target exists when you start scanning

[parameters.birthMean,parameters.birthCov,parameters.birthGaussians] = getBirthStats(parameters,50,10);
parameters.birthWeights = (1/parameters.birthGaussians)*ones(1,parameters.birthGaussians);


trueTracks = getTrueTrack(parameters);
correctObservations = getObservations(trueTracks,parameters);
observations = getClutteredObservations(correctObservations,parameters);
[estimatedTracks,existence,deltaK] = bernoulliFilter(observations,parameters);

%estimatedTracks = pdaf(observations,parameters);
rmse = getError(trueTracks,estimatedTracks);

figure(1)
plot=1;
if(plot==0)
    %plot XY of true tracks
    plot(trueTracks(1,:),trueTracks(2,:))
    axis([-100 100 -100 100])
    hold on
    %plot XY of estimated tracks to compare
    plot(estimatedTracks(1,:),estimatedTracks(2,:))
else
    plotOverTime(trueTracks,observations,estimatedTracks);
end


estRmse = mean(rmse(~isnan(rmse)))

function [birthMean, birthCov, birthGaussians] = getBirthStats(parameters,rangeStd,velStd)
    %make this a grid of births such that 1 standard deviation of
    %each birth mean covers the entire area.
    limits = parameters.limits;
    [numMeasurements,~] = size(parameters.startState);
    rangeStd; %range standard deviation
    velStd; %velocity standard deviation
    birthGaussians = round(((limits(1,2)-limits(1,1))/(2*rangeStd))*((limits(2,2)-limits(2,1))/(2*rangeStd))*((limits(3,2)-limits(3,1))/(2*velStd))*((limits(4,2)-limits(4,1))/(2*velStd)),0);
    birthMean = zeros(numMeasurements,birthGaussians);
    birthCov = zeros(numMeasurements,numMeasurements,birthGaussians);
    index = 1;
    for x=(limits(1,1)+rangeStd):2*rangeStd:limits(1,2)
        for y = (limits(2,1)+rangeStd):2*rangeStd:limits(2,2)
            for xVel=(limits(3,1)+velStd):2*velStd:limits(3,2)
                for yVel=(limits(4,1)+velStd):2*velStd:limits(4,2)
                    birthMean(:,index) = [x;y;xVel;yVel];
                    birthCov(:,:,index) = diag([rangeStd^2;rangeStd^2;velStd^2;velStd^2]);
                    index = index + 1;
                end
            end
        end
    end
end