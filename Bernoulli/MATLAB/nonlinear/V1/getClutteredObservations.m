function [clutteredObservations] = getClutteredObservations(correctObservations, parameters)
detectionProbability = parameters.detectionProbability;
meanClutter = parameters.meanClutter;
measurementRange = parameters.measurementRange;
[~, numSteps] = size(correctObservations);

clutteredObservations = cell(numSteps,1);
for step = 1:numSteps
    numFalseAlarms = poissrnd(meanClutter);
    falseAlarms = zeros(2,numFalseAlarms);
    falseAlarms(1,:) = measurementRange * rand(1,numFalseAlarms); %false alarms in range
    falseAlarms(2,:) = 360 * rand(1,numFalseAlarms) - 180; %false alarms in angle
    
    allMeasurements = falseAlarms;
    if(all(~isnan(correctObservations(:,step))))
        if(rand < detectionProbability)
            allMeasurements = [allMeasurements, correctObservations(:,step)]; %include the correct observation amidst the false alarms
        end
    end
    allMeasurements = allMeasurements(:,randperm(size(allMeasurements,2))); %randomly rearranges everything so the correct observation isn't always last.
    clutteredObservations{step} = allMeasurements;
end

end

