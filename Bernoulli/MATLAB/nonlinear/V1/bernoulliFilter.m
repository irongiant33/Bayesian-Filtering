function [estimatedState,existence,dk] = bernoulliFilter(observations,parameters)
    %% Initialize parameters
    numSteps = parameters.numSteps;
    priorMean = parameters.priorMean;
    scanTime = parameters.scanTime;
    dVar = parameters.sigmaDrivingNoise^2;
    rangeVar = parameters.sigmaMeasurementNoiseRange^2;
    angleVar = parameters.sigmaMeasurementNoiseBearing^2;
    measurementRange = parameters.measurementRange;
    meanClutter = parameters.meanClutter; %average number of clutter points
    pd = parameters.detectionProbability;
    birthGaussians = parameters.birthGaussians;
    percentageGaussians = birthGaussians;
    numGaussians = 0;
    birthMean = parameters.birthMean;
    birthCov = parameters.birthCov;
    spatialMean = birthMean;
    spatialCov = birthCov;
    birthWeights = parameters.birthWeights;
    spatialWeights = birthWeights;
    sensorPosition = parameters.sensorPosition;
    priorExistence = parameters.priorExistence;
    ps = parameters.survivalProb;
    pb = parameters.birthProb;
    
    %% Initialize variables
    [stateTransition,noiseMatrix] = getTransitionMatrices(scanTime);
    drivingCov = noiseMatrix*diag([dVar;dVar])*noiseMatrix';
    lambdaC = meanClutter /  (measurementRange * 360);
    measurementCov = diag([rangeVar;angleVar]);
    
    estimatedState = zeros(length(priorMean),numSteps+1);
    estimatedState(:,1) = priorMean;
    existence = zeros(1,numSteps+1);
    dk = zeros(1,numSteps + 1);
    dk(1) = 0;
    existence(1) = priorExistence;
    
    for i=2:(numSteps+1)  
        i
        % Get observations from current time step
        currentObs = cell2mat(observations(i-1));
        [numMeasurements, numObs] = size(currentObs);
        [numStates,~] = size(spatialMean);
        
        %% Initialize temporary variables  
        % Make room for more birth gaussians in the new spatial pdf.
        predictedSpatialMean = zeros(numStates,numGaussians + birthGaussians);
        predictedSpatialCov = zeros(numStates,numStates,numGaussians + birthGaussians);
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
        % initialize temporary variables
        deltaK = 0;
        addedTerm = zeros(1,priorNumGaussians * numObs);
        predictedMean = zeros(numMeasurements,priorNumGaussians * numObs);
        predictedCov = zeros(numMeasurements,numMeasurements,priorNumGaussians * numObs);
        likelihood = zeros(1,priorNumGaussians * numObs);
        
        % Calculate deltaK
        for obs=1:numObs 
            index = 1; 
            for n=(1 + (obs - 1) * priorNumGaussians):(priorNumGaussians + (obs - 1) * priorNumGaussians)
                % Generate sigma points
                [sigmas,sigmaWeights] = getSigmaPoints(predictedSpatialMean(:,index),predictedSpatialCov(:,:,index));
                numPoints = size(sigmas,2);

                % Propagate sigma points through observation model
                pointsMeasurements = zeros(2,numPoints); %2x5 matrix
                pointsMeasurements(1,:) = sqrt((sigmas(1,:) - sensorPosition(1)).^2+(sigmas(2,:) - sensorPosition(2)).^2);
                pointsMeasurements(2,:) =  atan2d(sigmas(1,:) - sensorPosition(1),sigmas(2,:) - sensorPosition(2));

                %calculate statistics of transformed state and prediction state
                for k=1:numPoints
                    predictedMean(1,n) = predictedMean(1,n) + sigmaWeights(k) * pointsMeasurements(1,k);
                end
                predictedMean(2,n) = angleFromPoints(pointsMeasurements(2,:)', sigmaWeights);
                
                diffMeasurements = zeros(2,numPoints);
                for k=1:numPoints
                    diffMeasurements(:,k) = (pointsMeasurements(:,k) - predictedMean(:,n));
                    diffMeasurements(2,k) = wrapTo180(diffMeasurements(2,k));
                    predictedCov(:,:,n) = predictedCov(:,:,n) + sigmaWeights(k) * (diffMeasurements(:,k) * diffMeasurements(:,k)');
                end
                predictedCov(:,:,n) = predictedCov(:,:,n) + measurementCov;
                 
                %calculate likelihood of observation and add to deltaK
                likelihood(n) = gaussian(currentObs(:,obs),predictedMean(:,n),predictedCov(:,:,n));
                addedTerm(n) = (spatialWeights(index) * likelihood(n))/lambdaC;
                deltaK = deltaK + addedTerm(n);
                index = index + 1;
            end
        end
        deltaK = pd - pd * deltaK;
        dk(i) = deltaK;
        
        % Initialize temporary variables
        predictedMeasurementMean = zeros(numMeasurements,numGaussians * numObs);
        predictedMeasurementCov = zeros(numMeasurements,numMeasurements,numGaussians * numObs);
        covCross = zeros(numStates,numMeasurements,numGaussians * numObs);
        innovation = zeros(numMeasurements,numGaussians * numObs);
        measurementLikelihood = zeros(1,numGaussians * numObs);
        
        % calculate measurement likelihood for each observation according to
        % the new spatial pdf.
        for obs=1:numObs            
            index = 1;
            for n=(1 + numGaussians * (obs - 1)):(numGaussians + numGaussians * (obs - 1))
                [sigmas,sigmaWeights] = getSigmaPoints(predictedSpatialMean(:,index),predictedSpatialCov(:,:,index));
                numPoints = size(sigmas,2);

                %propagate sigma points through observation model
                pointsMeasurements = zeros(2,numPoints); %2x5 matrix
                pointsMeasurements(1,:) = sqrt((sigmas(1,:) - sensorPosition(1)).^2+(sigmas(2,:) - sensorPosition(2)).^2);
                pointsMeasurements(2,:) =  atan2d(sigmas(1,:) - sensorPosition(1),sigmas(2,:) - sensorPosition(2));

                %calculate statistics of transformed state and prediction state
                for k=1:numPoints
                    predictedMeasurementMean(1,n) = predictedMeasurementMean(1,n) + sigmaWeights(k) * pointsMeasurements(1,k);
                end
                predictedMeasurementMean(2,n) = angleFromPoints(pointsMeasurements(2,:)', sigmaWeights);
                
                diffMeasurements = zeros(2,numPoints);
                for k=1:numPoints
                    diffMeasurements(:,k) = (pointsMeasurements(:,k) - predictedMeasurementMean(:,n));
                    diffMeasurements(2,k) = wrapTo180(diffMeasurements(2,k));
                    predictedMeasurementCov(:,:,n) = predictedMeasurementCov(:,:,n) + sigmaWeights(k) * (diffMeasurements(:,k) * diffMeasurements(:,k)');
                end
                predictedMeasurementCov(:,:,n) = predictedMeasurementCov(:,:,n) + measurementCov; 
                
                %calculate cross covariance matrix
                for k=1:numPoints
                    covCross(:,:,n) = covCross(:,:,n) + sigmaWeights(k) * (sigmas(:,k) - predictedSpatialMean(:,index)) * diffMeasurements(:,k)';
                end
                
                innovation(:,n) = currentObs(:,obs)-predictedMeasurementMean(:,n);
                innovation(2,n) = wrapTo180(innovation(2,n));
                
                %calculate measurement likelihood for later use in the
                %update step. 
                measurementLikelihood(n) = gaussian(currentObs(:,obs),predictedMeasurementMean(:,n),predictedMeasurementCov(:,:,n));
                index = index + 1;
            end
        end
        
        
        %% update
        % Initialize update variables
        spatialMean = zeros(numStates,numGaussians + numObs * numGaussians);
        spatialCov = zeros(numStates,numStates,numGaussians + numObs * numGaussians);
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
                kalmanGain = covCross(:,:,n)/predictedMeasurementCov(:,:,n);
                spatialMean(:,n) = predictedSpatialMean(:,index) + kalmanGain * innovation(:,n);
                spatialCov(:,:,n) = predictedSpatialCov(:,:,index) - kalmanGain  * predictedMeasurementCov(:,:,n) * kalmanGain';
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
    % save final variables for output
    estimatedState = estimatedState(:,2:length(observations)+1);
    existence = existence(1,2:length(observations)+1);
    dk = dk(1,2:length(observations)+1);
end

function [A,W] = getTransitionMatrices(scanTime)
    % Next state is previous state plus scanTime * previous state velocity.
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

% Important to note this gaussian wraps around for angle
function likelihood=gaussian(x,mean,covariance)
    [dim,~]=size(covariance);
    detInverse = (det(covariance))^-0.5;
    coeff = 1/((2*pi)^(dim/2)) * detInverse;
    diff = x - mean;
    diff(2) = wrapTo180(diff(2));
    likelihood = coeff * exp(-0.5 * diff' * covariance^-1 * diff);
end

% Produce 8 sigma points with the same mean and covariance as the parameters.
function [points, weights] = getSigmaPoints(mean, covariance)
    [lengthState, ~] =  size(mean);
    numSigmaPoints = 2*lengthState;

    %determine sigma point weights
    weights = 1/numSigmaPoints * ones(numSigmaPoints,1);

    %calculate sigma points
    sqrtmCovariance = sqrtm(covariance); %matrix square root
    points = zeros(lengthState,2*lengthState);
    points(:,1:lengthState) = repmat(mean,1,lengthState) + sqrt(lengthState)*sqrtmCovariance;
    points(:,lengthState+1:2*lengthState) = repmat(mean,1,lengthState) - sqrt(lengthState)*sqrtmCovariance;
end

% this function calculates the mean angle from a set of weighted points
function [angle] = angleFromPoints(points, weights)
    numPoints = size(points,2);

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
