function  CleanLabelImage = h_CleanLabel(LabelImage,nFrame)
% global nFrame

CleanLabelImage = zeros(size(LabelImage));
for tt = 1:nFrame
    currImg = uint16(LabelImage(:,:,tt));
    currIds = unique(currImg(:));
    currIds = currIds(currIds>0);
    
    labelImage = zeros(size(currImg));
    for ii = 1:length(currIds)
        curr_id = currIds(ii);
        binImg = (currImg == curr_id);   
        
        labelImg = bwlabel(binImg,4);
        
        props = regionprops(labelImg,'Area','PixelIdxList');
        
        if(length(props)>1)
            Areas = zeros(length(props),1);
            for jj = 1:length(props)
                Areas(jj) = props(jj).Area;    
            end
            [~,idx] = sort(Areas,'descend');
            for kk=2:length(idx)
                labelImg(props(idx(kk)).PixelIdxList)=0;
            end
%             fprintf('cleaned: %d labels\n',length(idx)-1);

        end
        
        labelImage = labelImage + double(curr_id)*double(labelImg>0);
    end

    CleanLabelImage(:,:,tt) = labelImage;
   
    
end
