% Implementation of the Unscented Probabilistic Data Association Filter
% (PDAF). By Sam Miller, 6/26/19

function estimatedTracks = unscentedPdaf(observations, parameters,correctObservations)
    numSteps = parameters.numSteps;
    priorMean = parameters.priorMean;
    scanTime = parameters.scanTime;
    priorCov = parameters.priorCovariance;
    drivingVar = parameters.sigmaDrivingNoise^2;
    rangeVar = parameters.sigmaMeasurementNoiseRange^2;
    angleVar = parameters.sigmaMeasurementNoiseBearing^2;
    measurementCov = diag([rangeVar;angleVar]);
    sensorPosition = parameters.sensorPosition;
    meanClutter = parameters.meanClutter;
    [stateTransition,noiseMatrix] = getTransitionMatrices(scanTime);
    drivingCov = noiseMatrix*diag([drivingVar;drivingVar])*noiseMatrix';
    pd = parameters.detectionProbability;
    lambda = meanClutter / ((2*parameters.measurementRange)^2);
    
    estimatedState = zeros(length(priorMean),numSteps+1);
    estimatedState(:,1) = priorMean;
    estimatedCov = priorCov;
    
    for i=2:(numSteps+1)
        currentObs = cell2mat(observations(i-1));
        [numMeasurements,numObs] = size(currentObs);
        
        %% prediction
        predictedState = stateTransition * estimatedState(:,i-1);
        predictedCov = stateTransition * estimatedCov * stateTransition' + drivingCov;
        
        %% fetch sigma points and propogate through nonlinear model
        [sigmas, weights] = getSigmas(predictedState,predictedCov);
        [stateSize,numPoints] = size(sigmas);
        
        %propogate sigma points through nonlinear model
        predictedMeasurements = zeros(numMeasurements,length(sigmas));
        predictedMeasurements(1,:) = sqrt((sigmas(1,:)-sensorPosition(1)).^2+(sigmas(2,:)-sensorPosition(2)).^2);
        predictedMeasurements(2,:) = atan2d(sigmas(1,:)-sensorPosition(1),sigmas(2,:)-sensorPosition(2));
        
        %% calculate statistics of predictedMeasurements
        meanMeasurements = zeros(numMeasurements,1);
        for k=1:numPoints
            meanMeasurements(1) = meanMeasurements(1) + weights(k) * predictedMeasurements(1,k);
        end
        meanMeasurements(2) = angleFromPoints(predictedMeasurements(2,:)',weights);
        
        covMeasurements = zeros(numMeasurements,numMeasurements);
        diffMeasurements = zeros(numMeasurements,length(sigmas));
        for k=1:numPoints
            diffMeasurements(:,k) = (predictedMeasurements(:,k) - meanMeasurements);
            diffMeasurements(2,k) = wrapTo180(diffMeasurements(2,k));
            covMeasurements = covMeasurements + weights(k) * (diffMeasurements(:,k) * diffMeasurements(:,k)');
        end
        covMeasurements = covMeasurements + measurementCov;
        
        %calculate cross covariance matrix
        covCross = zeros(stateSize,2);
        for k=1:numPoints
            covCross = covCross + weights(k) * (sigmas(:,k) - predictedState) * diffMeasurements(:,k)';
        end
        
        %% data association
        associationProbabilities = zeros(1,numObs+1); %last index corresponds to no detections
        likelihood = zeros(1,numObs);
        %evaluate likelihood for each measurement
        for k=1:numObs
            likelihood(k) = gaussian(meanMeasurements,covMeasurements,currentObs(:,k));
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
            innovation(:,k) = currentObs(:,k)-meanMeasurements;
            innovation(2,k) = wrapTo180(innovation(2,k));
            combinedInnovation = combinedInnovation + associationProbabilities(k) * innovation(:,k);
        end
        combinedSum = zeros(numMeasurements,numMeasurements);
        for k=1:numObs
            combinedSum = combinedSum + associationProbabilities(k)*innovation(:,k)*innovation(:,k)';
        end
        combinedSum = combinedSum - combinedInnovation*combinedInnovation';
        
        %% update
        kalmanGain = covCross * (covMeasurements^-1);
        estimatedState(:,i) = predictedState + kalmanGain * combinedInnovation;
        updatedCov = predictedCov - kalmanGain * covMeasurements * kalmanGain';
        originUncertaintyOffset = kalmanGain * (combinedSum)*kalmanGain'; %effect of measurement origin uncertainty
        estimatedCov = associationProbabilities(numObs+1) * predictedCov + (1-associationProbabilities(numObs+1))*updatedCov+originUncertaintyOffset;
    end
    estimatedTracks = estimatedState(:,2:length(estimatedState));
end

function [sigmas,weights] = getSigmas(predictedState,predictedCov)
    [L,~] = size(predictedState);
    numSigmas = 2*L;
    weights = (1/numSigmas)*ones(numSigmas,1);
    sigmas = zeros(L,numSigmas);
    sigmas(:,1:L) = repmat(predictedState,1,L) + sqrt(L) * sqrtm(predictedCov);
    sigmas(:,L+1:numSigmas) = repmat(predictedState,1,L) - sqrt(L) * sqrtm(predictedCov);
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

function [angle] = angleFromPoints(points, weights)
    [numPoints,~] = size(points);

    signs = sign(points);
    if(sum(signs*signs(1)) ~= numPoints)
        totalWeightPositive = sum(weights(signs>=0));
        positiveMean = sum(points(signs>=0).*weights(signs>=0)/totalWeightPositive);
        totalWeightNegative = sum(weights(signs<0));
        negativeMean = sum(points(signs<0).*weights(signs<0)/totalWeightNegative);

        if((positiveMean - negativeMean) >= 180)
            negativeMean = negativeMean+360;
        end
        angle = positiveMean*totalWeightPositive+negativeMean*totalWeightNegative;
    else
        angle = sum(points.*weights);
    end
    angle = wrapTo180(angle);
end

function likelihood=gaussian(mean,covariance,x)
    [dim,~]=size(covariance);
    detInverse = (det(covariance))^-0.5;
    coeff = 1/((2*pi)^(dim/2)) * detInverse;
    diff = x - mean;
    diff(2) = wrapTo180(diff(2));
    likelihood = coeff * exp(-0.5 * diff' * covariance^-1 * diff);
end