function [gfeat,numCell] = h_GetGeomFeats(labImage,intImage,nFrame)
% gfeat, size [nFrame,9 (number of features), # cells ]

if( isempty( labImage ) );
    gfeat = [];
    numCell = 0;
    return;
end
currIds = unique(labImage(:));
currIds = currIds(currIds>0);
currIds = sort(currIds);
% gfeat = -1000*ones(nFrame,9,length(currIds)); % if zero labels, it will be empty
gfeat = -1000*ones(nFrame,4,length(currIds)); % if zero labels, it will be empty
for jj = 1:length(currIds)
    for i=1:nFrame
        currLabImage = labImage(:,:,i);
        binImage = (currLabImage==currIds(jj));
        if(~any(binImage));continue;end

        props = regionprops(binImage,'Centroid', 'Area','Orientation','EquivDiameter','MinorAxisLength','MajorAxisLength','Perimeter');
        eccentr = props.MinorAxisLength/props.MajorAxisLength;
        gfeat(i,1:3,jj) = [
            props.Centroid(1); % fet_1_2
            props.Centroid(2);
            eccentr; % fet_3
%             props.Area; % fet_4
%             props.Orientation; % fet_5
%             props.EquivDiameter; % fet_6
%             props.Perimeter; % fet_7
        ];
%         gfeat(i,8,jj) = sum( sum( binImage.*intImage(:,:,i) ) ); % fet_8
    end
   % compute the displacement
   xy = gfeat(:,1:2,jj);
   dist = circshift(xy,[-1 0])-xy;
   dist = sqrt(sum(dist.^2,2));
   dist(end,1) = 0;

   % FIXME, what happens when there is a missing time point ? 
   distMiss = sum(abs(circshift(xy,[-1 0]))+abs(xy),2);

   dist( distMiss > 1000,1 ) = -1000;
   gfeat(:,4,jj) = dist;  % fet_4

end


numCell = length(currIds);
end
