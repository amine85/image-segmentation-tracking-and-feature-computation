
################################################################### This Part is for Users ##############################################

# Adquisition conditions
numb_well_per_row = 8
num_of_time_decimals = 2
num_of_block_decimals = 3

# Directories and Path:

# -------- Experiments
#microscope = 'leica'
#data_id = '20140530CARbeads'
#root_path = '/home/melquiades/aLab/00-nano-team/data/leica/'

microscope = 'zeiss'
data_id = '20140425ILZENKA'
root_path = '/data/shared/iliadi/'

#microscope = 'zeiss'
#data_id = '20140326GRCARN6calcium'
#root_path = '/data/shared/gromain/'

# -------- paths
data_path = root_path + data_id +'/'

save_root_path = '/data/pdata/nrey/00-nano-team/results/'
out_path = save_root_path + data_id +'_RES3/'

exe_path= '/data/pdata/nrey/00-nano-team/exes/'
param_path = '/data/pdata/nrey/00-nano-team/parameters/'

## Run Flags ##
## This has to be 1 at least once for each data set (subfolders)
saveFileNames = 1
## save filenames and run image to stack should be enabled/disabled at the same time
runImageToStack = 1
##
runUnmixing = 1
##
runBackgroundSubstraction = 1
## If this is 1, make sure the runBackgroundSubstraction has also been 1 at least once
runSegmentation = 1
##
runWellCropping = 1

## Channel Dictionary: ##
## change the channel name according to the actual fluophore used in the experiment ##

## Leica Dictionary
#channel_dict = {"bright_field":"BF",
    #"effectors":"GFP",
    #"targets":"RFP",
    #"death":"",
    #"beads_1":"CY5",
    #"beads_2":""}

##Zeiss Dictionary
channel_dict = {"bright_field":"c1_ORG",
                "effectors":"c4_ORG",
                "targets":"c3_ORG",
                "death":"c2_ORG",
                "beads_1":"--",
                "beads_2":"--"}

## Indicate which Channels to Process: ##
stack_tuple = ("bright_field","targets","effectors","death") ## eg in a time lapse, we wantto stack multiple snap shots together to create a video

preprocess_tuple = ("effectors","targets")## remiving background

segment_tuple = ("effectors","targets") ## for segmenting

unmix_tuple = ("targets","effectors") ## assume the first channel leaks into the second one,  removing spectral overlap

crop_tuple = ("bright_field","death") # Preproces and segment tuple will be segmented, plus this one

# Analysis
range_blocks = range(1,8) ## it can be a range or a list


