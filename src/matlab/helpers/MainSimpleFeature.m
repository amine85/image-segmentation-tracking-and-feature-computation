function MainSimpleFeature( dd, dirList, N, w_hsize, timeToCountSeeds, dataPath, process_channels_simple_fea, nFrame, num_well )

block_number = dd - 2;
block_name = dirList(dd).name;
if( isempty(strfind(block_name,'B')) ),return,end
%     disp(['In Block ',block_name])

savePath = [dataPath,block_name,'/crops/'];

% Read text file of the well detection
shift = dlmread([dataPath,block_name,'/crops/','a_shifts.txt']); % time,xx_shift,yy_shift
center  = dlmread([dataPath,block_name,'/crops/','a_centers.txt']); % Row, col, x_coor, y_coord

% Load the well data and apply the shift
numWells = size(center,1);
well_location_local = cell(size(center,1),1);
for ww = 1:numWells
    if( center(ww,3) == -1000 || center(ww,4) == -1000 )
        well_location_local{ww,1} = [-1000,1000]; % Not a valid well
        continue;
    end
    temp = [];
    shiftAcc = [0,0];
    for time = 1:nFrame
        shiftAcc = shiftAcc+shift(time,2:3);
        temp = [temp;center(ww,3:4) - shiftAcc]; % FIXME, the shift is with respect to the previous one
    end
    well_location_local{ww,1} = temp;
end

%% Output files
% 2 rows per channel
numRowsPerWell = 2*5;
outFeatures = -1000*ones( num_well*num_well*numRowsPerWell, nFrame+3 );
outFeatures( :, 1 ) = block_number;


%     parfor ww = 1:numWells%1:length(well_info_local)         % process this well


% now read channels, segment and track %    
for ch = 1:length(process_channels_simple_fea)
    for ww = 1:numWells
        well_location = well_location_local{ww,1};
        if( well_location(1) == -1000 || well_location(2) == -1000 )
            continue;
        end

%         tableRange = 1+(numRowsPerWell)*(ww-1):(numRowsPerWell)*(ww);
        tableRange = 1+(numRowsPerWell)*(ww-1):(numRowsPerWell)*(ww);
    %     ch1Raw = 1+(num_well*num_well)*(w-1);
    %     ch2Raw = 1+(num_well*num_well)*(w-1)+2;
    %     ch3Raw = 1+(num_well*num_well)*(w-1)+4;
    %     ch4Raw = 1+(num_well*num_well)*(w-1)+6;
        chRaw = 1+(numRowsPerWell)*(ww-1)+2*(ch-1);
        chBin = 1+(numRowsPerWell)*(ww-1)+2*(ch-1)+1;

        % Put the labels of wells (only valid ones)
        % The centers are stored in base 0 (c++)
        row = center(ww,1) + 1; % The name is 1 based
        col = center(ww,2) + 1;
        outFeatures( tableRange, 2 ) = row;
        outFeatures( tableRange, 3 ) = col;
        
        % Read binary image
        if( ch == 1 || ch == 2 ) % The preprocess channels, add the bg
            img_fname = [savePath,'imgR',num2str(row),'C',num2str(col),process_channels_simple_fea{ch},'bg.tif'];
        else
            img_fname = [savePath,'imgR',num2str(row),'C',num2str(col),process_channels_simple_fea{ch},'.tif'];
        end
        bin_fname = [savePath,'imgR',num2str(row),'C',num2str(col),process_channels_simple_fea{ch},'bin.tif'];
        
        % Raw image
        if(~exist(img_fname,'file')),continue;end
        img_data = double(h_read3DImage(img_fname,nFrame));
        img_sum = squeeze(sum(sum(img_data)))';
        outFeatures( chRaw, 3+1:end ) = img_sum/(size(img_data,1)*size(img_data,2));
            
        % Bin image
        if(~exist(bin_fname,'file')),continue;end
        bin_data = double(h_read3DImage(bin_fname,nFrame));
        bin_sum = squeeze(sum(sum(bin_data)))';
        outFeatures( chBin, 3+1:end ) = bin_sum/(size(img_data,1)*size(img_data,2));
        
    end
end

% Save features
feName = [dataPath,block_name,'/a_PixelFeatures.txt'];
if( ~isempty( outFeatures ));dlmwrite(feName,outFeatures,'delimiter','\t');end













