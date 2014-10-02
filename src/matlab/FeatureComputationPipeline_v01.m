clear; close all;clc
%% Parameters to be changed

nFrame = 9;
num_well = 8; % Num of well in each col, and each row
dataPath = '/home/melquiades/aLab/00-nano-team/results/20140425ILZENKA_5_RES/';
process_channels = {'CH1','CH2','CH3','CH4','CH5'};
preprocess_channels = {'CH1','CH2'}; % Not implemented, by default 1 and 2, use the bg sufix, in the crop image

%% Open pool (TO open pool, use once)
% matlabpool
% parpool
if (~matlabpool('size') > 0 )
    error('ERROR: No pool open')
end
addpath('helpers')

%% From here, DO NOT CHANGE

DEBUG = 1;
mainDirList = dir(dataPath);

%% Main loop, simple features
disp('Computing features')
parfor dd = 3:length(mainDirList)
    tt = clock;
    MainMainFeatures( dd, mainDirList,nFrame,num_well,dataPath,process_channels);
    fprintf('Block: %d, t: %d\n',dd-2,etime(clock,tt))
end

%% Put things together
disp('Putting results together')
Table_Exp = [];
Table_E1_T1 = [];
Table_E1_T2 = [];
Table_E1_T3 = [];
% Put results together
for dd = 3:length(mainDirList)
    tt = clock;
    
    block_name = mainDirList(dd).name;
    if( isempty(strfind(block_name,'B')) ),continue,end
    
    outName = [dataPath,block_name,filesep,'a_table_block.txt'];
    if( exist(outName,'file'));
        Table_Exp = [Table_Exp;dlmread( outName )];
    end
    
    outName = [dataPath,block_name,filesep,'a_table_e1_t1.txt'];
    if( exist(outName,'file'));
        Table_E1_T1 = [Table_E1_T1;dlmread( outName )];
    end
    
    outName = [dataPath,block_name,filesep,'a_table_e1_t2.txt'];
    if( exist(outName,'file'));
        Table_E1_T2 = [Table_E1_T2;dlmread( outName, '\t' )];
    end
    
    outName = [dataPath,block_name,filesep,'a_table_e1_t3.txt'];
    if( exist(outName,'file'));
        Table_E1_T3 = [Table_E1_T3;dlmread( outName )];
    end
    fprintf('Block: %d, t: %d\n',dd-2,etime(clock,tt))
end
 
 %% Save the global results
mkdir( dataPath,'features')

outName = [dataPath,'features',filesep,'Table_Exp.txt'];
if( exist(outName,'file'));delete( outName );end
dlmwrite(outName,Table_Exp,'delimiter','\t');

outName = [dataPath,'features',filesep,'Table_E1_T1.txt'];
if( exist(outName,'file'));delete( outName );end
dlmwrite(outName,Table_E1_T1,'delimiter','\t');

outName = [dataPath,'features',filesep,'Table_E1_T2.txt'];
if( exist(outName,'file'));delete( outName );end
dlmwrite(outName,Table_E1_T2,'delimiter','\t');

outName = [dataPath,'features',filesep,'Table_E1_T3.txt'];
if( exist(outName,'file'));delete( outName );end
dlmwrite(outName,Table_E1_T3,'delimiter','\t');





























