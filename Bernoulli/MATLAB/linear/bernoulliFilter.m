function [estimatedState,existence,dk] = bernoulliFilter(observations,parameters)
    numSteps = parameters.numSteps;
    priorMean = parameters.priorMean;
    scanTime = parameters.scanTime;
    priorExistence = parameters.priorExistence;
    ps = parameters.survivalProb;
    pb = parameters.birthProb;
    estimatedCov = parameters.priorCovariance;
    dVar = parameters.sigmaDrivingNoise^2;
    mVar = parameters.sigmaMeasurementNoise^2;
    stateTransform = diag([1;1;1;1]);
    measurementCov = diag([mVar;mVar;mVar;mVar]);
    limits = parameters.limits;
    meanClutter = parameters.meanClutter; %average number of clutter points
    lambdaC = meanClutter / ((limits(1,2)-limits(1,1)) * (limits(2,2)-limits(2,1)) * (limits(3,2)-limits(3,1)) * (limits(4,2)-limits(4,1)));
    pd = parameters.detectionProbability;
    [stateTransition,noiseMatrix] = getTransitionMatrices(scanTime);
    drivingCov = noiseMatrix*diag([dVar;dVar])*noiseMatrix';
    birthGaussians = parameters.birthGaussians;
    percentageGaussians = birthGaussians;
    numGaussians = 0;
    birthMean = parameters.birthMean;
    birthCov = parameters.birthCov;
    spatialMean = birthMean;
    spatialCov = birthCov;
    birthWeights = parameters.birthWeights;
    spatialWeights = birthWeights;
    
    estimatedState = zeros(length(priorMean),numSteps+1);
    estimatedState(:,1) = priorMean;
    existence = zeros(1,numSteps+1);
    dk = zeros(1,numSteps + 1);
    dk(1) = 0;
    existence(1) = priorExistence;
    
    for i=2:(numSteps+1)  
        i
        currentObs = cell2mat(observations(i-1));
        [numMeasurements, numObs] = size(currentObs);
        
        %% declare variables  
        % Make room for more birth gaussians in the new spatial pdf.
        predictedSpatialMean = zeros(numMeasurements,numGaussians + birthGaussians);
        predictedSpatialCov = zeros(numMeasurements,numMeasurements,numGaussians + birthGaussians);
        % Make room for more birth weights in the new spatial weight array. 
        weights = zeros(1,numGaussians + birthGaussians);
        weights(1,numGaussians+1:numGaussians+birthGaussians) = birthWeights;
        
        % If there were prior gaussians in the spatial PDF, assign them to
        % the spatial PDF
        if(numGaussians > 0)
            predictedSpatialMean(:,1:numGaussians) = spatialMean;
            predictedSpatialCov(:,:,1:numGaussians) = spatialCov;
            weights(1,1:numGaussians) = spatialWeights;
        end
        
        %% predict
        % Predict existence probability
        predictedExistence = pb * (1 - priorExistence) + ps * priorExistence;
        
        % If there were prior gaussians in the spatial PDF, predict their
        % mean and covariance
        for n=1:numGaussians
            priorMean = stateTransition * predictedSpatialMean(:,n);
            priorCov = stateTransition * predictedSpatialCov(:,:,n) * stateTransition' + drivingCov;
            
            survivalComp = (ps * priorExistence)/predictedExistence;
            weights(n) = survivalComp * spatialWeights(n);
            predictedSpatialMean(:,n) = priorMean;
            predictedSpatialCov(:,:,n) = priorCov;
        end
        % update birth weights and fill in birth mean/covariance
        birthComp = (pb * (1 - priorExistence))/predictedExistence;
        weights(1+numGaussians:birthGaussians+numGaussians) = birthComp * birthWeights;
        predictedSpatialMean(:,1+numGaussians:birthGaussians+numGaussians) = birthMean;
        predictedSpatialCov(:,:,1+numGaussians:birthGaussians+numGaussians) = birthCov;
        
        % update total number of gaussians in the new spatial pdf and
        % remember the number of gaussians in the prior spatial pdf
        priorNumGaussians = numGaussians;
        numGaussians = numGaussians + birthGaussians;
        
        %% pre-update calculations   
        %calculate deltaK
        deltaK = 0;
        addedTerm = zeros(1,priorNumGaussians * numObs);
        for obs=1:numObs 
            index = 1; 
            for n=(1 + (obs - 1) * priorNumGaussians):(priorNumGaussians + (obs - 1) * priorNumGaussians)
                predictedMean = stateTransform * predictedSpatialMean(:,index);
                predictedCov = stateTransform * predictedSpatialCov(:,:,index) * stateTransform' + measurementCov; 
                likelihood = gaussian(currentObs(:,obs),predictedMean,predictedCov);
                addedTerm(n) = (spatialWeights(index) * likelihood)/lambdaC;
                deltaK = deltaK + addedTerm(n);
                index = index + 1;
            end
        end
        deltaK = pd - pd * deltaK;
        dk(i) = deltaK;
        
        % calculate measurement likelihood for each observation according to
        % the new spatial pdf.
        predictedMeasurementMean = zeros(numMeasurements,numGaussians * numObs);
        predictedMeasurementCov = zeros(numMeasurements,numMeasurements,numGaussians * numObs);
        measurementLikelihood = zeros(1,numGaussians * numObs);
        for obs=1:numObs            
            index = 1;
            for n=(1 + numGaussians * (obs - 1)):(numGaussians + numGaussians * (obs - 1))
                predictedMeasurementMean(:,n) = stateTransform * predictedSpatialMean(:,index);
                predictedMeasurementCov(:,:,n) = stateTransform * predictedSpatialCov(:,:,index) * stateTransform' + measurementCov;
                measurementLikelihood(n) = gaussian(currentObs(:,obs),predictedMeasurementMean(:,n),predictedMeasurementCov(:,:,n));
                index = index + 1;
            end
        end
        
        %% update
        spatialMean = zeros(numMeasurements,numGaussians + numObs * numGaussians);
        spatialCov = zeros(numMeasurements,numMeasurements,numGaussians + numObs * numGaussians);
        spatialWeights = zeros(1,numGaussians + numObs * numGaussians);
        
        % update existence
        priorExistence = ((1 - deltaK)/(1 - deltaK * predictedExistence)) * predictedExistence;
        existence(i) = priorExistence;
        
        % assign all predicted spatial statistics to the end of the updated
        % array for ease of parsing in the next step. This is the missed
        % detection component of the updated spatial pdf
        spatialComp = (1 - pd)/(1 - deltaK);
        spatialMean(:,(1 + numObs * numGaussians):(numGaussians + numObs * numGaussians)) = predictedSpatialMean;
        spatialCov(:,:,(1 + numObs * numGaussians):(numGaussians + numObs * numGaussians)) = predictedSpatialCov;
        spatialWeights(1,(1 + numObs * numGaussians):(numGaussians + numObs * numGaussians)) = spatialComp * weights;
        
        % this is the detection component of the updated spatial pdf
        updateComp = pd/(1 - deltaK);
        for obs=1:numObs
            index = 1;
            for n=(1 + (obs - 1) * numGaussians):(numGaussians + (obs - 1) * numGaussians)
                spatialWeights(1,n) = updateComp * ((weights(1,index) * measurementLikelihood(n))/lambdaC);
                kalmanGain = (predictedSpatialCov(:,:,index) * stateTransform')/predictedMeasurementCov(:,:,n);
                spatialMean(:,n) = (predictedSpatialMean(:,index) + kalmanGain * (currentObs(:,obs) - predictedMeasurementMean(:,n)));
                spatialCov(:,:,n) = (predictedSpatialCov(:,:,index) - kalmanGain * stateTransform * predictedSpatialCov(:,:,index));
                index = index + 1;
            end
        end
        %update total number of gaussians in the updated spatial pdf
        numGaussians = numGaussians + numObs * numGaussians;
        
        %% Calculate estimate
        % estimate is a weighted average over all gaussians in the updated
        % spatial pdf. The prior existence also serves as a weight to
        % determine how much it should rely on the previous known state (if
        % the target doesn't exist) and how much it should rely on the
        % spatial pdf (if the target exists)
        for n=1:numGaussians
            estimatedState(:,i) = (estimatedState(:,i)+spatialWeights(n) * spatialMean(:,n)) *priorExistence + estimatedState(:,i-1) * (1 - priorExistence);
        end
        
        %% pruning
        % keep only top N weights where N is the number of Gaussians
        % added at each new birth
        threshold = min(maxk(spatialWeights,percentageGaussians));
        tempMean = spatialMean;
        tempCov = spatialCov;
        tempWeights = spatialWeights;
        index = 1;
        for k=1:numGaussians
            if(tempWeights(k)>=threshold)
                spatialMean(:,index) = tempMean(:,k);
                spatialCov(:,:,index) = tempCov(:,:,k);
                spatialWeights(index) = tempWeights(k);
                index = index + 1;
            end
        end
        spatialMean = spatialMean(:,1:index-1);
        spatialCov = spatialCov(:,:,1:index-1);
        % Renormalize Weights
        spatialWeights = spatialWeights(1:index-1)/sum(spatialWeights(1:index-1));
        numGaussians = index - 1; 
    end
    estimatedState = estimatedState(:,2:length(observations)+1);
    existence = existence(1,2:length(observations)+1);
    dk = dk(1,2:length(observations)+1);
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

function likelihood=gaussian(x,mean,covariance)
    [dim,~]=size(covariance);
    detInverse = (det(covariance))^-0.5;
    coeff = 1/((2*pi)^(dim/2)) * detInverse;
    diff = x - mean;
    likelihood = coeff * exp(-0.5 * diff' * covariance^-1 * diff);
end
