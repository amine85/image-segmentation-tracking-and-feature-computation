function asso = h_GetAssocFeat( tarImg, tarLab, effImg, nFrame )

% global nFrame

if( isempty( tarLab ) );asso = [];return;end
currIds = unique(tarLab(:));
currIds = currIds(currIds>0);
currIds = sort(currIds);
asso = -1000*ones(nFrame,1,length(currIds)); % if zero labels, it will be empty
for jj = 1:length(currIds)
    for i=1:nFrame
        currLabImage = tarLab(:,:,i);
        binImage = (currLabImage==currIds(jj));
        if(~any(binImage));continue;end

        asso(i,1,jj) = GetAssociation( binImage, effImg(:,:,i) );
    end

end

end


function intensity = GetAssociation( binImage, effImg )

     D = bwdist(binImage);
     intensity = 0;
     for r = 1:10
         maskImage = double(D>0 & D<=r);
         intensity = intensity + (sum(sum(maskImage.*effImg))/(r^2));
     end
     intensity = intensity/(sum(effImg(:)+eps));

end
