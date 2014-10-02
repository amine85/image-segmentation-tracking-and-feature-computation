function h_save3DImage(fName,imdata,bittype)

switch (bittype)
    case '16'
        for zz = 1:size(imdata,3)
              imwrite( uint16(imdata(:,:,zz)),fName,'Compression','none','WriteMode','append');              
        end
    case '8'
        for zz = 1:size(imdata,3)
              imwrite( uint8(imdata(:,:,zz)),fName,'Compression','none','WriteMode','append');              
        end
end

end