function [ effDeathAsso ] = h_GetDeathFeat( imgLab, imgDeath, nFrame )

if( isempty( imgLab ) );effDeathAsso = [];return;end
currIds = unique(imgLab(:));
currIds = currIds(currIds>0);
currIds = sort(currIds);
effDeathAsso = -1000*ones(nFrame,1,length(currIds)); % if zero labels, it will be empty
for jj = 1:length(currIds)
    for i=1:nFrame
        currLabImage = imgLab(:,:,i);
        binImage = (currLabImage==currIds(jj));
        if(~any(binImage));continue;end

        effDeathAsso(i,1,jj) = sum(sum(binImage.*imgDeath(:,:,i)))/sum(binImage(:));

    end

end

end