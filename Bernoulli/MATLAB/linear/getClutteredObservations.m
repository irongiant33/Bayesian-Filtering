function [clutteredObservations] = getClutteredObservations(correctObservations, parameters)
detectionProbability = parameters.detectionProbability;
meanClutter = parameters.meanClutter;
limits = parameters.limits;
[~, numSteps] = size(correctObservations);

clutteredObservations = cell(numSteps,1);
for step = 1:numSteps
    
    numFalseAlarms = poissrnd(meanClutter);
    falseAlarms = zeros(2,numFalseAlarms);
    falseAlarms(1,:) = (limits(1,2) - limits(1,1)) * rand(1,numFalseAlarms) + limits(1,1); %false alarms in range
    falseAlarms(2,:) = (limits(2,2) - limits(2,1)) * rand(1,numFalseAlarms) + limits(2,1); %false alarms in angle
    falseAlarms(3,:) = (limits(3,2) - limits(3,1)) * rand(1,numFalseAlarms) + limits(3,1);
    falseAlarms(4,:) = (limits(4,2) - limits(4,1)) * rand(1,numFalseAlarms) + limits(4,1);
    
    allMeasurements = falseAlarms;
    
    if(all(~isnan(correctObservations(:,step))))
        if(rand < detectionProbability)
            allMeasurements = [allMeasurements, correctObservations(:,step)]; %include the correct observation amidst the false alarms
        end
        allMeasurements = allMeasurements(:,randperm(size(allMeasurements,2))); %randomly rearranges everything so the correct observation isn't always last.
    end
    
    clutteredObservations{step} = allMeasurements;
end

end

