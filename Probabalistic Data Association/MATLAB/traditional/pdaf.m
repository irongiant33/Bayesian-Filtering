% Implementation of the Unscented Probabilistic Data Association Filter
% (PDAF). By Sam Miller, 6/26/19

function estimatedState = pdaf(observations, parameters)
    numSteps = parameters.numSteps;
    priorMean = parameters.priorMean;
    scanTime = parameters.scanTime;
    estimatedCov = parameters.priorCovariance;
    dVar = parameters.sigmaDrivingNoise^2;
    mVar = parameters.sigmaMeasurementNoise^2;
    stateTransform = diag([1;1;1;1]);
    measurementCov = diag([mVar;mVar;mVar;mVar]);
    limits = parameters.limits;
    meanClutter = parameters.meanClutter; %average number of clutter points
    lambda = meanClutter / ( (limits(1,2)-limits(1,1)) * (limits(2,2)-limits(2,1)) * (limits(3,2)-limits(3,1)) * (limits(4,2)-limits(4,1)));
    pd = parameters.detectionProbability;
    [stateTransition,noiseMatrix] = getTransitionMatrices(scanTime);
    drivingCov = noiseMatrix*diag([dVar;dVar])*noiseMatrix';
    
    estimatedState = zeros(length(priorMean),numSteps+1);
    estimatedState(:,1) = priorMean;
    
    for i=2:(numSteps+1)
        currentObs = cell2mat(observations(i-1));
        [numMeasurements, numObs] = size(currentObs);
        
        %% prediction
        priorMean = stateTransition * estimatedState(:,i-1);
        priorCov = stateTransition * estimatedCov * stateTransition' + drivingCov;
        measurementPrediction = stateTransform * priorMean;
        measurementPredictionCov = stateTransform * priorCov*stateTransform' + measurementCov; %aka innovation cov
        
        %% data association
        associationProbabilities = zeros(1,numObs+1); %last index corresponds to no detections
        likelihood = zeros(1,numObs);
        %evaluate likelihood for each measurement (SOMETHING IS WRONG HERE)
        for k=1:numObs
            likelihood(k) = gaussian(measurementPrediction,measurementPredictionCov,currentObs(:,k));
        end 
        likelihood = (likelihood * pd)/lambda;
        %calculate association probabilities for all association events
        totalLikelihood = sum(likelihood);
        for k=1:numObs
            associationProbabilities(k) = likelihood(k)/(1-pd+totalLikelihood);
        end
        associationProbabilities(numObs+1) = (1-pd)/(1-pd+totalLikelihood);
        
        %% pre-update calculations
        combinedInnovation = zeros(numMeasurements,1);
        innovation = zeros(numMeasurements,numObs);
        for k=1:numObs
            innovation(:,k) = currentObs(:,k)-measurementPrediction;
            combinedInnovation = combinedInnovation + associationProbabilities(k) * innovation(:,k);
        end
        combinedSum = zeros(numMeasurements,numMeasurements);
        for k=1:numObs
            combinedSum = combinedSum + associationProbabilities(k)*innovation(:,k)*innovation(:,k)';
        end
        combinedSum = combinedSum - combinedInnovation*combinedInnovation';
        
        %% update
        kalmanGain = priorCov * stateTransform' * (measurementPredictionCov^-1);
        estimatedState(:,i) = priorMean + kalmanGain * combinedInnovation;
        updatedCov = priorCov - kalmanGain * measurementPredictionCov * kalmanGain';
        originUncertaintyOffset = kalmanGain * (combinedSum)*kalmanGain'; %effect of measurement origin uncertainty
        estimatedCov = associationProbabilities(numObs+1) * priorCov + (1-associationProbabilities(numObs+1))*updatedCov+originUncertaintyOffset;
    end
    estimatedState = estimatedState(:,2:length(observations)+1);
end

function [A,W] = getTransitionMatrices(scanTime)
    A = diag([1;1;1;1]);
    A(1,3) = scanTime;
    A(2,4) = scanTime;
    
    %W is noise gain
    W = zeros(4,2);
    W(1,1) = 0.5*scanTime^2;
    W(2,2) = 0.5*scanTime^2;
    W(3,1) = scanTime;
    W(4,2) = scanTime;
end

function likelihood=gaussian(mean,covariance,x)
    [dim,~]=size(covariance);
    detInverse = (det(covariance))^-0.5;
    coeff = 1/((2*pi)^(dim/2)) * detInverse;
    diff = x - mean;
    likelihood = coeff * exp(-0.5 * diff' * covariance^-1 * diff);
end
