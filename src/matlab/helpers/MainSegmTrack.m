function MainSegmTrack( dd, dirList, N, w_hsize, timeToCountSeeds, dataPath, process_channels_seg_trc, nFrame )
    
block_number = dd - 2;
block_name = dirList(dd).name;
if( isempty(strfind(block_name,'B')) ),return,end
%     disp(['In Block ',block_name])

savePath = [dataPath,block_name,'/crops/'];

% Read text file of the well detection
shift = dlmread([dataPath,block_name,'/crops/','a_shifts.txt']); % time,xx_shift,yy_shift (Shift have to be applied negative)
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

%% now read channels, segment and track %    
for ch = 1:length(process_channels_seg_trc)
%         tt2 = clock;

    % Read seed information
    seed_fname = [dataPath,block_name,'/bg_',block_name,process_channels_seg_trc{ch},'.txt'];  % Stored like col (x-axis), row (y-axis)
    seed_data = dlmread(seed_fname);
    
    for ww = 1:numWells
        well_location = well_location_local{ww,1};
        if( well_location(1) == -1000 || well_location(2) == -1000 )
            continue;
        end

        % The centers are stored in base 0 (c++)
        row = center(ww,1) + 1; % The name is 1 based
        col = center(ww,2) + 1;

        % Test if binary exists
        fName = [savePath,'imgR',num2str(row),'C',num2str(col),process_channels_seg_trc{ch},'bin.tif'];
        if(~exist(fName)),continue,end


        cell_count = -ones(nFrame,1);
        valid = true;

        for time = 1:nFrame
            timeBaseZero = time - 1;
            seeds = seed_data(seed_data(:,1)==timeBaseZero,2:3);

            % FIXME, what was this thing for
%                 if(time==1)
%                     t = 1;
%                 else
%                     t = time-1;
%                 end

%             % FIXME, CHECK THAT THE ROWS AND COLS ARE CORRECT
%             rowc = well_location(time,1);
%             colc = well_location(time,2);
% 
%             rowctemp = max(w_hsize+1,rowc);
%             rowctemp = min(N-(w_hsize+1),rowctemp);
% 
%             rmin = rowctemp - w_hsize;
%             rmax = rowctemp + w_hsize;
%             if(rmin<1 || rmax>N)
%                 valid = false;
% %                     valid
%                 break;
%             end
% 
%             colctemp = max(w_hsize+1,colc);
%             colctemp = min(N-(w_hsize+1),colctemp);
% 
%             cmin = colctemp - w_hsize;
%             cmax = colctemp + w_hsize;
%             if(cmin<1 || cmax>N)
%                 valid = false;
% %                     valid
%                 break;
%             end
            
            xCenter = well_location(time,1);
            yCenter = well_location(time,2);

            xmin = xCenter - w_hsize;
            xmax = xCenter + w_hsize;
            
            ymin = yCenter - w_hsize;
            ymax = yCenter + w_hsize;

            condition = logical(seeds(:,1)>xmin & seeds(:,1)<xmax &...
                                seeds(:,2)>ymin & seeds(:,2)<ymax);
            cell_count(time,1) = sum(condition);
        end

% %             % find out the number of cells: 
% %             cell_count = [cell_count(end);cell_count;cell_count(1)];
% %             cell_count = medfilt1(cell_count,3);
% %             cell_count = cell_count(2:end-1);
        cell_count = cell_count(1:timeToCountSeeds,1);

        cell_count = cell_count(cell_count>=0);

        counts = sort(unique(cell_count));
        count_hist = zeros(length(counts),1);
        for ii = 1:length(counts)
            count_hist(ii) = sum(cell_count==counts(ii));                
        end
        count_hist = count_hist/sum(count_hist);
        [likelihood,idx] = max(count_hist);
        num_cells = counts(idx);
        if(num_cells==0)
            continue;
        end

        % Read binary image
        binaryImage = double(h_read3DImage(fName,nFrame));

        % segment the binary image %
        if(num_cells>0 && valid)
            labelImage = h_SpectralClusteringSegmentation(binaryImage,num_cells,nFrame);
            labelImage = h_CleanLabel(labelImage,nFrame);
            labelImage = h_GetTracks(labelImage,nFrame);
        else
            labelImage = binaryImage;
        end

        % Save the label image
        fName = [savePath,'imgR',num2str(row),'C',num2str(col),process_channels_seg_trc{ch},'label.tif'];
        if(exist(fName,'file')),delete(fName);end
        h_save3DImage(fName,labelImage,'16');    
    end
%         fprintf('\tCH: %d, t: %d\n',ch,etime(clock,tt2))
end
%     fprintf('Block: %d, t: %d\n',dd-2,etime(clock,tt))
