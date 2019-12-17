function plotOverTime(trueTracks,observations,estimatedTracks,isMoving)
[~,numSteps]=size(trueTracks);

figure(1)
if(isMoving)
    for step=1:numSteps
        currentObs = cell2mat(observations(step));
        [~,numObs] = size(currentObs);
        q = plot(trueTracks(1,1:step),trueTracks(2,1:step),'-k');        
        hold on
        axis([-100 100 -100 100])
        p = plot(estimatedTracks(1,1:step),estimatedTracks(2,1:step),'-r');
        estColor = get(p,'color');
        trueColor = get(q,'color');
        scatter(estimatedTracks(1,step),estimatedTracks(2,step),50,estColor,'x');
        scatter(trueTracks(1,step),trueTracks(2,step),50,trueColor,'x');
        
        for obs=1:numObs
            scatter(currentObs(1,obs),currentObs(2,obs),30,'b');
        end
        hold off
        
        pause(0.1)
    end
else
    plot(trueTracks(1,:),trueTracks(2,:),'-k')
    axis([-100 100 -100 100])
    hold on
    plot(estimatedTracks(1,:),estimatedTracks(2,:),'-r')
    hold off
end

end

%maybe improve this so that the estimated tracks' size of the scatter point
%is related to the estimated covariance at that step to get an idea of
%how the estimate improves over time. Would need to save covariance over
%time in the filter code though.