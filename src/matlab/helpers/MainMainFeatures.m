function MainMainFeatures( dd, mainDirList,nFrame,num_well,dataPath,process_channels)

    block_number = dd - 2;
    block_name = mainDirList(dd).name;
    if( isempty(strfind(block_name,'B')) ),return,end

    Table_E1_T1 = [];
    Table_E1_T2 = [];
    Table_E1_T3 = [];
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chEfe = find(strcmp(process_channels,'CH1'));
    chTar = find(strcmp(process_channels,'CH2'));
    chDea = find(strcmp(process_channels,'CH3'));

    chBead1 = find(strcmp(process_channels,'CH4'));
    chBead2 = find(strcmp(process_channels,'CH5'));

    % All the features to compute
    numRows = 6+8*5+1+2;
    % 6 Effector
    % 8x5 Target x5
    % 2 Beads
    % 1 Time points
    
    % To store the features of the well
    Table_Well = [];
    
    curr_directory_crops = [dataPath,block_name,filesep,'crops',filesep];
%     curr_directory_featu = [dataPath,block_name,filesep,'features',filesep];
%     fprintf('In dir: %s\n',block_name);
    
    for row =  1:num_well
        for col = 1:num_well
            imgOrg = cell(1,length(process_channels));
            imgLab = cell(1,length(process_channels));
            numEfe = [0];
            numTar = [0];
            for ch = 1:length(process_channels)
            
                fNameOrg = [curr_directory_crops,'imgR',num2str(row),'C',num2str(col),process_channels{ch},'*.tif'];
                dirList = dir(fNameOrg);
                if isempty( dirList )
                    continue;
                end
                if( ch == 1 || ch == 2 ) % The preprocess channels, add the bg
                    fNameOrg = [curr_directory_crops,'imgR',num2str(row),'C',num2str(col),process_channels{ch},'bg.tif'];
                else
                    fNameOrg = [curr_directory_crops,'imgR',num2str(row),'C',num2str(col),process_channels{ch},'.tif'];
                end
                if( ~exist(fNameOrg,'file'))
                    continue;
                end
                imgOrg{1,ch} = h_read3DImage(fNameOrg,nFrame);
%                 fprintf('Row: %d, Col: %d, Ch: %d\n',[col,row,ch]);

                fNameLab = [curr_directory_crops,'imgR',num2str(row),'C',num2str(col),process_channels{ch},'label.tif'];
                if( ~exist(fNameLab,'file'))
                    continue;
                end
                imgLab{1,ch} = h_read3DImage(fNameLab,nFrame);
            end
            
            % Effector features ( DEBUG: for now separated )
            effFeat = [];
            if( chEfe ~= 0 )
                [effFeat,numEfe] = h_GetGeomFeats(imgLab{1,chEfe},imgOrg{1,chEfe},nFrame);
            end
            
            effDeathAsso = [];
            if( ~isempty( effFeat) && ~isempty(imgOrg{1,chDea}) ) %chDea ~=0
                [ effDeathAsso ] = h_GetDeathFeat( imgLab{1,chEfe},imgOrg{1,chDea},nFrame);
            end
            
            % Target features ( DEBUG: for now separated )
            tarFeat = [];
            if( chTar ~= 0 )
                [tarFeat,numTar] = h_GetGeomFeats(imgLab{1,chTar},imgOrg{1,chTar},nFrame);
            end
            
            tarDeathAsso = [];
            if( ~isempty( tarFeat) && ~isempty(imgOrg{1,chDea}) ) %chDea ~=0
                [ tarDeathAsso ] = h_GetDeathFeat( imgLab{1,chTar},imgOrg{1,chDea},nFrame);
            end
            
            % Association with a different channel
            % In general target is the base channel, and effector is the
            % asociated channel
            tarFeatAsso = [];
            if( ~isempty( effFeat) && ~isempty( tarFeat ) )
                [ tarFeatAsso ] = h_GetAssocFeat(imgOrg{1,chTar},imgLab{1,chTar},imgOrg{1,chEfe},nFrame);
            end
            
            % Beads
            % FIXME, only works with one bead, otherwise just chooses the
            % first one to store the intensity
            bead1Asso = [];
            if( chBead1 ~=0 )
                [ bead1Asso ] = h_GetDeathFeat( imgLab{1,chBead1},imgOrg{1,chBead1},nFrame);
            end
            
            bead2Asso = [];
            if( chBead2 ~=0 )
                [ bead2Asso ] = h_GetDeathFeat( imgLab{1,chBead2},imgOrg{1,chBead2},nFrame);
            end
            
            %% Put results together
            
            if numEfe > 1 || numTar>5
                continue
            end
            % FIXME, allocate the space before putting results together
%             FinalTable = -1000*ones(numRows,3+nFrame);

            % Empty part effector
            if numEfe == 0 
                effFeat = -1000*ones(4,nFrame)';
            end
            % Empty part deadmarker for effector
            if isempty(effDeathAsso) 
                 effDeathAsso = -1000*ones(1,nFrame)';
            end
            
            % 6 rows
            FinalTable = [effFeat';
                effDeathAsso';
                -1*ones(1,nFrame)];
            
            % Empty part deadmarker for target
            if isempty(tarDeathAsso)  
                for cnt = 1:numTar
                    tarDeathAsso(:,:,cnt) = -1000*ones(1,nFrame)';
                end
            end
            
            % Empty part targetfeatureassociation
            if isempty(tarFeatAsso)  
                for cnt = 1:numTar
                   tarFeatAsso(:,:,cnt) = -1000*ones(1,nFrame)';
                end
            end
            
            % 8 per target, for 5 targets
            for cnt = 1:numTar
                FinalTable = [FinalTable;
                    tarFeat(:,:,cnt)';
                    tarDeathAsso(:,:,cnt)';
                    -1*ones(1,nFrame);
                    tarFeatAsso(:,:,cnt)';
                    -2*ones(1,nFrame);];
            end
            
            %empty part target
            for cnt = 1:5-numTar
                FinalTable = [FinalTable;-1000*ones(8,nFrame);];
            end
            
            % Beads, two rows 
            if( isempty( bead1Asso ) )
                FinalTable = [FinalTable;-1000*ones(1,nFrame);];
            else
                FinalTable = [FinalTable;bead1Asso(:,:,1)'];
            end
                    
            if( isempty( bead2Asso ) )
                FinalTable = [FinalTable;-1000*ones(1,nFrame);];
            else
                FinalTable = [FinalTable;bead2Asso(:,:,1)';];
            end
            
            % Separate by the number of time points
            FinalTable = [FinalTable;(1:nFrame);];
            % Add the block, row, col in the left
            FinalTable = [ones(numRows,1)*str2double(block_name(2:end)),ones(numRows,1)*row,ones(numRows,1)*col,FinalTable];
            % Write localy
%             outName = [curr_directory_featu,'imgR',num2str(col),'R',num2str(row),'features.txt'];
%             if( exist(outName,'file'));delete( outName );end
%             dlmwrite(outName,FinalTable,'delimiter','\t');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % All the wells
            Table_Well = [Table_Well;FinalTable];
            % Specific classification
            if( (numTar == 1) && (numEfe == 1 ) )
                Table_E1_T1 = [Table_E1_T1;FinalTable];
            end
            if( (numTar == 2) && (numEfe == 1 ) )
                Table_E1_T2 = [Table_E1_T2;FinalTable];
            end
            if( (numTar == 3) && (numEfe == 1 ) )
                Table_E1_T3 = [Table_E1_T3;FinalTable];
            end

        end
    end
    
    outName = [dataPath,block_name,filesep,'a_table_block.txt'];
    if( exist(outName,'file'));delete( outName );end
    if( ~isempty( Table_Well ));dlmwrite(outName,Table_Well,'delimiter','\t');end
    
    outName = [dataPath,block_name,filesep,'a_table_e1_t1.txt'];
    if( exist(outName,'file'));delete( outName );end
    if( ~isempty( Table_E1_T1 ));dlmwrite(outName,Table_E1_T1,'delimiter','\t');end
    
    outName = [dataPath,block_name,filesep,'a_table_e1_t2.txt'];
    if( exist(outName,'file'));delete( outName );end
    if( ~isempty( Table_E1_T2 ));dlmwrite(outName,Table_E1_T2,'delimiter','\t');end
    
    outName = [dataPath,block_name,filesep,'a_table_e1_t3.txt'];
    if( exist(outName,'file'));delete( outName );end
    if( ~isempty( Table_E1_T3 ));dlmwrite(outName,Table_E1_T3,'delimiter','\t');end
    
end