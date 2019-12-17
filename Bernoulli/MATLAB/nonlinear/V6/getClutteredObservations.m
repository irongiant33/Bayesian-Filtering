function [clutteredObservations] = getClutteredObservations(correctObservations, parameters)
detectionProbability = parameters.detectionProbability;
meanClutter = parameters.meanClutter;
measurementRange = parameters.measurementRange;
velocityRange = parameters.velocityRange;
[~, numSteps, numDevices] = size(correctObservations);

clutteredObservations = cell(numSteps,numDevices);
mode = 1;
for step = 1:numSteps
    for device = 1:numDevices
        if(mod(step-1,3)==(device-1))
            mode = 1;
            numFalseAlarms = poissrnd(meanClutter(1));
            falseAlarms = zeros(3,numFalseAlarms);
            falseAlarms(1,:) = measurementRange * rand(1,numFalseAlarms); %false alarms in range
            falseAlarms(3,:) = 2*velocityRange * rand(1,numFalseAlarms) - velocityRange; %false alarms in velocity
            falseAlarms(2,:) = 360 * rand(1,numFalseAlarms) - 180; %false alarms in angle
        else
            mode =2;
            numFalseAlarms = poissrnd(meanClutter(2));
            falseAlarms = zeros(3,numFalseAlarms);
            falseAlarms(1,:) = measurementRange * rand(1,numFalseAlarms); %false alarms in range
            falseAlarms(2,:) = 360 * rand(1,numFalseAlarms) - 180; %false alarms in angle
            falseAlarms(3,:) = zeros(1,size(falseAlarms(1,:),2));
        end    
        allMeasurements = falseAlarms;
        if(all(~isnan(correctObservations(:,step,device))))
            if(rand < detectionProbability(mode))
                allMeasurements = [allMeasurements, correctObservations(:,step,device)]; %include the correct observation amidst the false alarms
            end
        end
        allMeasurements = allMeasurements(:,randperm(size(allMeasurements,2))); %randomly rearranges everything so the correct observation isn't always last.
        clutteredObservations{step,device} = allMeasurements;
    end
end

end

