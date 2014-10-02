# -*- coding: utf-8 -*-
import shutil, fnmatch, os, subprocess, os.path, time, glob, sys, inspect, filecmp, datetime, getpass, multiprocessing
import os, re, operator, pdb, subprocess


def SaveFiles( data_path, save_root_path, data_id, range_blocks, out_path, stack_tuple, channel_dict, channel_naming_dict, microscope, num_of_time_decimals, num_of_block_decimals ):

    if 'leica' in microscope:
        #num_of_time_decimals = 3
        #num_of_block_decimals = 3
        blocks = range_blocks
        # read files #
        file_list = os.listdir(data_path)
        if not file_list:
            print "empty directory:\n"+ data_path 
            sys.exit()
            
        file_list_path = os.path.join(save_root_path+data_id+'_FileList/')
        if not os.path.isdir( file_list_path ):
          os.makedirs( file_list_path )
        new_blocks = []
        for b in blocks:
            for ch in stack_tuple:
                if num_of_block_decimals == 2:
                    if b<10:                       
                        tempbb = '0'+str(b)          
                    elif b<100:
                        tempbb = ''+str(b)
                    else:
                        tempbb = str(b)
                if num_of_block_decimals == 3:
                    if b<10:                       
                        tempbb = '00'+str(b)          
                    elif b<100:
                        tempbb = '0'+str(b)
                    else:
                        tempbb = str(b)
                block_dir = 'B'+tempbb

                channel_list = [filename for filename in file_list \
                    if channel_dict[ch] in filename and "_t"+str(b)+".TIF" in filename and channel_dict[ch] ] # and channel_dict[ch], added for the cases of empy dictiory
                if channel_list:
                    new_blocks.append(b)
                    print "In Block (1.stack): "+str(b)+", ch: "+ch
                    time_list = []
                    #print len(channel_list)
                    for filename in channel_list:
                        temp = filename.split("_w") 
                        temp = temp[0]
                        temp = temp[-num_of_time_decimals:]                # get the 3 characters before _w
                        if re.search(r'\d+', temp):
                            time_list.append(int(re.search(r'\d+', temp).group()))
                        else :
                            time_list.append(-1)

                    # sort the files and get the index:
                    index = [time_list.index(x) for x in sorted(time_list)]
                    index = [x for x in index if x != time_list.index(-1)]
                    # sort the channel list:
                    channel_list = [channel_list[i] for i in index]
                    fp = open(os.path.join(file_list_path,'inputfnames_'+block_dir+ '_' +channel_naming_dict[ch]+'.txt'),'w')
                    for filename in channel_list:
                        fp.write(os.path.join(data_path,filename))
                        #print filename_list_txt[(ch,b)]
                        fp.write('\n')
                    fp.close()

        for b in new_blocks:
            if num_of_block_decimals == 2:
                if b<10:                       
                    tempbb = '0'+str(b)          
                elif b<100:
                    tempbb = ''+str(b)
                else:
                    tempbb = str(b)
            if num_of_block_decimals == 3:
                if b<10:                       
                    tempbb = '00'+str(b)          
                elif b<100:
                    tempbb = '0'+str(b)
                else:
                    tempbb = str(b)
            block_dir = 'B'+tempbb+'/'
            if not os.path.isdir(os.path.join(out_path,block_dir)):
                os.makedirs(os.path.join(out_path,block_dir))
                
                
    if 'zeiss' in microscope:
        #num_of_time_decimals = 2
        #num_of_block_decimals = 2
        blocks = range_blocks
        # read files #
        file_list = os.listdir(data_path)
        if not file_list:
          print "empty directory:\n"+ data_path 
          sys.exit()

        file_list_path = os.path.join(save_root_path+data_id+'_FileList/')
        if not os.path.isdir(file_list_path):
          os.makedirs(file_list_path )
        new_blocks = []
        for b in blocks:
            for ch in stack_tuple:
                if num_of_block_decimals == 2:
                    if b<10:                       
                        tempbb = '0'+str(b)          
                    elif b<100:
                        tempbb = ''+str(b)
                    else:
                        tempbb = str(b)
                if num_of_block_decimals == 3:
                    if b<10:                       
                        tempbb = '00'+str(b)          
                    elif b<100:
                        tempbb = '0'+str(b)
                    else:
                        tempbb = str(b)
                block_dir = 'B'+tempbb

                tempbb = "_s"+tempbb
                channel_list = [filename for filename in file_list
                                if channel_dict[ch] in filename and tempbb in filename and ".tif" in filename and channel_dict[ch] ] # and channel_dict[ch], added for the cases of empy dictiory
                if channel_list:
                    print channel_list
                    new_blocks.append(b)
                    print "In Block (1.stack): "+str(b)+", ch: "+ch
                    time_list = []
                    for filename in channel_list:
                        temp = filename.split("c") 
                        temp = temp[0]
                        temp = temp[-num_of_time_decimals:]             # get the 3 characters before _w
                        print temp
                        if re.search(r'\d+', temp):
                          time_list.append(int(re.search(r'\d+', temp).group()))
                          #print int(re.search(r'\d+', temp).group())
                        else :
                          time_list.append(0)

                    # sort the files and get the index:
                    index = [time_list.index(x) for x in sorted(time_list)]
                    # sort the channel list:
                    channel_list = [channel_list[i] for i in index]
                    fp = open(os.path.join(file_list_path,'inputfnames_'+block_dir+ '_' +channel_naming_dict[ch]+'.txt'),'w')
                    for filename in channel_list:
                        fp.write(os.path.join(data_path,filename))
                        #print filename_list_txt[(ch,b)]
                        #print os.path.join(data_path,filename)
                        fp.write('\n')
                    fp.close()    

        for b in new_blocks:
            if num_of_block_decimals == 2:
                if b<10:                       
                    tempbb = '0'+str(b)          
                elif b<100:
                    tempbb = ''+str(b)
                else:
                    tempbb = str(b)
            if num_of_block_decimals == 3:
                if b<10:                       
                    tempbb = '00'+str(b)          
                elif b<100:
                    tempbb = '0'+str(b)
                else:
                    tempbb = str(b)
            block_dir = 'B'+tempbb+'/'
            if not os.path.isdir(os.path.join(out_path,block_dir)):
              os.makedirs(os.path.join(out_path,block_dir))



def RunParallel(subprocess_command_list, numThreads):
  #print subprocess_command_list
  
  #spawn and add subprocesses (each representing an execution of register_pair.exe) into a list
  subprocess_list=[]
  file_handle_list=[]
  launched_processes_list=[]
  launched_subp_list=[]
  still_launched_subp_list=[]

  idx = 0
  threads_launched = 0
  contador = 0
  flagAllStart = 0
  for subprocess_command in subprocess_command_list:
    #while loop to check active processes to see if we can launch more threads
    process = subprocess.Popen('free -k|grep Mem:', shell=True, stdout=subprocess.PIPE)
    lista = process.communicate()[0].split()

    memoryAvilable = 1
    if int(lista[3]) < 150000000:
      memoryAvilable = 0
      memoryAvilable = 1
    print lista[3]+str(memoryAvilable)  

    while (threads_launched >= numThreads ) | (memoryAvilable == 0): # at least one thread, but no more than 3/4 of the machine
    #while (threads_launched >= max(1,multiprocessing.cpu_count())-5 ) | (memoryAvilable == 0): # at least one thread, but no more than 3/4 of the machine
    #while (threads_launched >= max(1,1) ) | (memoryAvilable == 0): # at least one thread, but no more than 3/4 of the machine
      if flagAllStart == 0:
        flagAllStart = 1
        print "\tALL THE PROCESS ARE STARTED"
      threads_launched = 0
      still_launched_subp_list = []
      for launched_subp in launched_processes_list:                  
        if launched_subp.poll() == None:
          threads_launched = threads_launched + 1
          still_launched_subp_list.append(launched_subp)
      launched_processes_list = still_launched_subp_list[:]
      #print "Threads launched: " + str(threads_launched)
      time.sleep(0.1)
      process = subprocess.Popen('free -k|grep Mem:', shell=True, stdout=subprocess.PIPE)
      lista = process.communicate()[0].split()
      memoryAvilable = 1
      if int(lista[3]) < 150000000:
        memoryAvilable = 0
        memoryAvilable = 1
    #fh = open(os.getcwd() + '/debug/debug_from_' + from_image_list[idx] + '_to_' + to_image_list[idx] + '.txt', 'w')
    contador = contador + 1
    #print "Launching " + subprocess_command + ", "+str(contador)+" of "+str(len(subprocess_command_list))
    subp = subprocess.Popen(subprocess_command, shell=True)
    #subp = subprocess.Popen(subprocess_command)
    
    subprocess_list.append(subp)
    launched_processes_list.append(subp)
    #file_handle_list.append(fh)
    threads_launched = threads_launched + 1;
    idx = idx + 1
    time.sleep(0.1)
    
  #loop through entire list to check if process are done
  processes_done = False
  while not processes_done:
    processes_done = True
    time.sleep(1) #sleep 1 second before polling all processes
    for subp in subprocess_list:
      if subp.poll() == None:
        processes_done = False

