#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, re, operator, pdb, subprocess, multiprocessing
import sys, os, time
import numpy
import helper
import time
from final_pipeline_params import *
#from helper import SaveFiles, RunParallel

if __name__ == '__main__':

    ################################################################### This Part is for Code Readers ##############################################

    numprocessors = 7
    numprocessorsForCropping = 20

    channel_naming_dict = {"bright_field":"0",
        "effectors":"1",
        "targets":"2",
        "death":"3",
        "beads_1":"4",
        "beads_2":"5"}

    ################################################################### This Part is for Code Readers ##############################################
    if runImageToStack:
       saveFileNames = 1

    # Do a path check
    if not os.path.isdir(exe_path):
        print "executable directory:\n"+ exe_path +"\ndoes not exist"
        sys.exit()
    if not os.path.isdir(root_path):
        print "root directory:\n"+ root_path +"\ndoes not exist"
        sys.exit()
    if saveFileNames:
        if not os.path.isdir(data_path):
            print "data directory:\n"+ data_path +"\ndoes not exist"
            sys.exit()
    if not os.path.isdir(param_path):
        print "parameters directory:\n"+ param_path +"\ndoes not exist"
        sys.exit()
    if not os.path.isdir(save_root_path):
        print "parameters directory:\n"+ save_root_path +"\ndoes not exist"
        sys.exit()

    if not os.path.isdir(out_path):
        os.makedirs(os.path.join(out_path))

    if saveFileNames:
        helper.SaveFiles( data_path, save_root_path, data_id, range_blocks, out_path, stack_tuple, channel_dict, channel_naming_dict, microscope, num_of_time_decimals, num_of_block_decimals )
        file_list_path = os.path.join(save_root_path+data_id+'_FileList/')

########################################## Main Processing Loop ##################################

    block_list = os.listdir(out_path)
    block_list = sorted(block_list)
    #for block_dir in block_list[1]:
    t1 = time.time()
    commandsGlobal = []
    for block_dir in block_list:
        print 'In Block '+str(block_dir)
        if not os.path.isdir(os.path.join(out_path,block_dir)):
            os.makedirs(os.path.join(out_path,block_dir))

        filename_dict = {}
        for ch in stack_tuple:
            filename_dict[ch] = block_dir +'CH'+channel_naming_dict[ch]+'.tif'

        #### Stack Image Data ####################################
        ompNumThreads = 1 #FIXME is giving an error when reading in parallel
        commands = []
        if runImageToStack:
            for ch in stack_tuple:
                temp = []
                temp.append(os.path.join(exe_path,'image_to_stack'))
                temp.append(os.path.join(file_list_path,'inputfnames_'+block_dir+ '_' +channel_naming_dict[ch]+'.txt'))
                temp.append(os.path.join(out_path,block_dir))
                temp.append(filename_dict[ch])
                temp.append( str(ompNumThreads) )
                print '\timgToStack, cmd: '+" ".join(temp)

                #subprocess.call(temp)
                temp = " ".join(temp)
                commandsGlobal.append(temp)
                #commands.append(temp)
            #helper.RunParallel(commands,1)
    helper.RunParallel(commandsGlobal,numprocessors)
    print "1-STACK TIME: " +str(time.time() - t1)

    t1 = time.time()
    commandsGlobal = []
    for block_dir in block_list:
        filename_dict = {}
        for ch in stack_tuple:
            filename_dict[ch] = block_dir +'CH'+channel_naming_dict[ch]+'.tif'
      ###### spectral unmixing ############################################################
        ompNumThreads = 11
        umx_prefix = 'umx_'
        commands = []
        if runUnmixing:
            temp = []
            temp.append( os.path.join(exe_path,'unmix16') )
            temp.append( os.path.join(out_path,block_dir,filename_dict[unmix_tuple[0]]) )
            temp.append( os.path.join(out_path,block_dir,filename_dict[unmix_tuple[1]]) )
            temp.append( os.path.join(out_path,block_dir,umx_prefix + filename_dict[unmix_tuple[0]]) )
            temp.append( os.path.join(out_path,block_dir,umx_prefix + filename_dict[unmix_tuple[1]]) )
            temp.append( os.path.join(param_path,str(data_id)+'_mixing_matrix.txt' ) )
            temp.append( os.path.join(out_path,block_dir,'a_newMixing.tif' ) )
            temp.append( str(ompNumThreads) )

            print '\tunmix, cmd: '+" ".join(temp)

            #subprocess.call(temp)
            temp = " ".join(temp)
            commandsGlobal.append(temp)
            #commands.append(temp)
        #helper.RunParallel(commands,2)
    helper.RunParallel(commandsGlobal,numprocessors)
    print "2-UNMIX TIME: " +str(time.time() - t1)

    t1 = time.time()
    commandsGlobal = []
    for block_dir in block_list:
        filename_dict = {}
        for ch in stack_tuple:
            filename_dict[ch] = block_dir +'CH'+channel_naming_dict[ch]+'.tif'
      ##### background subtraction ############################################################
        ompNumThreads = 11
        bg_param = 40
        bg_prefix = 'bg_'
        commands = []
        if runBackgroundSubstraction:
            for ch in preprocess_tuple:
                temp = []
                temp.append(os.path.join(exe_path,'background_subtraction'))
                if ch in unmix_tuple:
                    temp.append( os.path.join(out_path,block_dir,umx_prefix + filename_dict[ch]) )
                else:
                    temp.append( os.path.join(out_path,block_dir,filename_dict[ch]) )
                temp.append( os.path.join(out_path,block_dir,bg_prefix + filename_dict[ch]) )
                temp.append( str(bg_param) )
                temp.append( str(ompNumThreads) )
                print '\tbacksubs, cmd: '+" ".join(temp)

                #subprocess.call(temp)
                temp = " ".join(temp)
                commandsGlobal.append(temp)
                #commands.append(temp)
        #helper.RunParallel(commands,2)
    helper.RunParallel(commandsGlobal,numprocessors)
    print "3-BS TIME: " +str(time.time() - t1)

    t1 = time.time()
    commandsGlobal = []
    for block_dir in block_list:
        filename_dict = {}
        for ch in stack_tuple:
            filename_dict[ch] = block_dir +'CH'+channel_naming_dict[ch]+'.tif'
      ####### mixture segmentation #####################################################################################
        ompNumThreads = 11
        clean_prefix = 'bin_'
        commands = []
        if runSegmentation:
            for ch in segment_tuple:
                temp = []
                temp.append(os.path.join(exe_path,'mixture_segment'))
                temp.append(os.path.join(out_path,block_dir,bg_prefix + filename_dict[ch]))
                temp.append(os.path.join(out_path,block_dir,clean_prefix + filename_dict[ch]))
                temp.append(os.path.join(param_path,str(data_id)+'_segmentation_paramters.txt'))
                temp.append( str(ompNumThreads) )
                print '\tsegment, cmd: '+" ".join(temp)

                #subprocess.call(temp)
                temp = " ".join(temp)
                commandsGlobal.append(temp)
                #commands.append(temp)
        #helper.RunParallel(commands,2)
    helper.RunParallel(commandsGlobal,numprocessors)
    print "4-SEGM TIME: " +str(time.time() - t1)

    t1 = time.time()
    commandsGlobal = []
    for block_dir in block_list:
        filename_dict = {}
        for ch in stack_tuple:
            filename_dict[ch] = block_dir +'CH'+channel_naming_dict[ch]+'.tif'
      ####### well cropping #####################################################################################
        ompNumThreads = 3
        commands = []
        if runWellCropping:
            if not os.path.isdir(os.path.join(out_path,block_dir,'crops')):
                os.makedirs( os.path.join(out_path,block_dir,'crops') )
            if not os.path.isdir(os.path.join(out_path,block_dir,'features')):
                os.makedirs( os.path.join(out_path,block_dir,'features') )
            temp = []
            temp.append(os.path.join(exe_path,'wellDetection'))
            temp.append(os.path.join(param_path,str(data_id)+'_P_rom.tif'))
            temp.append(os.path.join(param_path,str(data_id)+'_P_squ.tif'))
            temp.append(os.path.join(out_path,block_dir,'crops' ))
            temp.append('0' )  #starting slice
            temp.append( str(numb_well_per_row) ) #number of wells in each column/row
            temp.append( str(ompNumThreads) )
            for ch in crop_tuple:
                temp.append(os.path.join(out_path,block_dir, block_dir+'CH'+channel_naming_dict[ch]+'.tif' ))
            for ch in segment_tuple:
                temp.append(os.path.join(out_path,block_dir, 'bin_'+ block_dir+'CH'+channel_naming_dict[ch]+'.tif' ))
            for ch in preprocess_tuple:
                temp.append(os.path.join(out_path,block_dir, 'bg_'+ block_dir+'CH'+channel_naming_dict[ch]+'.tif' ))

            print temp
            #subprocess.call(temp)
            temp = " ".join(temp)
            commandsGlobal.append(temp)
            #commands.append(temp)
        #helper.RunParallel(commands,2)
    helper.RunParallel(commandsGlobal,numprocessorsForCropping)
    print "5-CROP TIME: " +str(time.time() - t1)








