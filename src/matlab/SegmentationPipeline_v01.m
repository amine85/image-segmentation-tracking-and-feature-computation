clear;close all;
%% Parameters to be changed

% nFrame = 65;
% num_well = 8; % Num of well in each col, and each row
% dataPath = '/data/pdata/nrey/00-nano-team/results/20140425ILZENKA_RES2/';
nFrame = 70;
num_well = 7; % Num of well in each col, and each row
dataPath = '/data/pdata/nrey/00-nano-team/results/20140323GR_RES3/';
process_channels_seg_trc = {'CH1','CH2'};
preprocess_channels = {'CH1','CH2'}; % Not implemented, by default 1 and 2, use the bg sufix, in the crop image for the pixel features
timeToCountSeeds = 7;

%% Open pool (TO open pool, use once)
% matlabpool
% parpool
if (~matlabpool('size') > 0 )
    error('ERROR: No pool open')
end
addpath('helpers')

%% From here, DO NOT CHANGE
process_channels_simple_fea = {'CH1','CH2','CH3','CH4','CH5'};

dirList = dir(dataPath);
N = 2048;
if( num_well == 8 )
    w_hsize = 115; % (231-1)/2
elseif( num_well == 7 )
    w_hsize = 125; % (231-1)/2
end

if(~exist(dataPath,'dir'))
    error('directory does not exist');
end

%% Main Loop Segmentation Tracking
disp('Segmentation and Tracking')
parfor dd = 3:length(dirList)%ini:fin%3:24+2%length(dirList)%23:length(dirList)%3:length(dirList)
    tt = clock;
    MainSegmTrack( dd, dirList, N, w_hsize, timeToCountSeeds, dataPath, process_channels_seg_trc, nFrame );
    fprintf('Block: %d, t: %d\n',dd-2,etime(clock,tt))
end
    
%% Main loop, simple features
disp('Pixel features')
parfor dd = 3:length(dirList)
    tt = clock;
    MainSimpleFeature( dd, dirList, N, w_hsize, timeToCountSeeds, dataPath, process_channels_simple_fea, nFrame, num_well );
    fprintf('Block: %d, t: %d\n',dd-2,etime(clock,tt))
end

%% Put things together
Table_Pixel = [];
% Put results together
for dd = 3:length(dirList)
    
    block_name = dirList(dd).name;
    if( isempty(strfind(block_name,'B')) ),continue,end
    
    outName = [dataPath,block_name,filesep,'/a_PixelFeatures.txt'];
    if( exist(outName,'file'));
        Table_Pixel = [Table_Pixel;dlmread( outName )];
    end
end
 
 %% Save the global results
mkdir( dataPath,'features')

outName = [dataPath,'features',filesep,'Table_Pixel.txt'];
if( exist(outName,'file'));delete( outName );end
dlmwrite(outName,Table_Pixel,'delimiter','\t');










