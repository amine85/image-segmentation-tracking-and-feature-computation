function img = h_read3DImage(fName,num_slices)

    if(~exist(fName,'file'))
        fName
        error('file does not exist');
    end
    temp = double(imread(fName,1));
    [m,n] = size(temp);
    img = zeros(m,n,num_slices);
    
    img(:,:,1) = temp;

    for zz = 2:num_slices
        img(:,:,zz) = double(imread(fName,zz));
    end
% InfoImage=imfinfo(fName);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% % num_slices=length(InfoImage);
% img=zeros(nImage,mImage,num_slices,'double');
% 
% TifLink = Tiff(fName, 'r');
% for i=1:num_slices
%    TifLink.setDirectory(i);
%    img(:,:,i)=double(TifLink.read());
% end
% TifLink.close();

end









% fName='ImageStack.tif';
% InfoImage=imfinfo(fName);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% num_slices=length(InfoImage);
% img=zeros(nImage,mImage,num_slices,'uint16');
% 
% TifLink = Tiff(fName, 'r');
% for i=1:num_slices
%    TifLink.setDirectory(i);
%    img(:,:,i)=TifLink.read();
% end
% TifLink.close();









