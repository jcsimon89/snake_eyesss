"""
For large images such as a anat images, each step can be super slow with not
too many timepoints.
It would be better to have multiprocessing that just creates and index such
as [0,1,..n] and feeds them to the motion correction function.
There's much more overhead but because each process takes several minutes that
should be managable.



Benchmark:

Newest: Yandan's func data (2x10.87Gb files)
threads: 32, 8 cores: 00:55:45 Memory Utilized: 21.58 GB

threads: 32, 8 cores:  00:09:51

Using Yandan's data (Previously: 05:42:11, Memory Efficiency: 5.94% of 70.93 GB)
threads: 32, 8 cores: 00:55:36 Memory Efficiency: 37.53% of 111.46 GB (before optimizing with del!)
threads: 32, 16 cores:  00:54:51, Memory Efficiency: Memory Efficiency: 19.38% of 111.46 GB (after optimizing)

Large Yandan's data (previously 09:22:01, Memory Efficiency: 85.64% of 171.34 GB)
threads: 32, 8 cores 02:48:45, Memory Efficiency: 51.83% of 250.00 GB
threads: 32, 16 cores -> failed SILENTLY with memory error
--> I need to add a control mechanism, for loop checking if index is missing and call moco_func
before combining
threads: 32, 16 cores, more memory: 02:41:33, Memory Efficiency: 78.12% of 156.04 GB

-> It seems that ant.registration is taking advantage of multiple available cores! No decrease
when going from 8 to 16 parallel processes! But this is still much faster and memory requirements much
lower than before!
"""
import os
import nibabel as nib
import pathlib
import ants
import numpy as np
import time
import multiprocessing
import natsort
import shutil
import itertools
import argparse
import sys
from pystackreg import StackReg 
from para_stack_reg import ParaReg

# To import files (or 'modules') from the brainsss folder, define path to scripts!
# path of workflow i.e. /Users/dtadres/snake_brainsss/workflow
#scripts_path = pathlib.Path(#
#    __file__
#).parent.resolve()
#sys.path.insert(0, pathlib.Path(scripts_path, "workflow").as_posix())
parent_path = str(pathlib.Path(pathlib.Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, parent_path)
#print(sys.path)
# This just imports '*.py' files from the folder 'brainsss'.
from brainsss import moco_utils
from brainsss import utils
###
# Global variable
###
WIDTH = 120 # Todo: make a global parameter class somewhere to keep track of this variable!
#TESTING = False

def moco_slice(
    moving_path,
    functional_channel_paths,
    par_output,
    moco_settings
):
    """
    Loop doing the motion correction for each slice.
    This is the function that is doing the heavy lifting of the multiprocessing.
    Saves the transformed files as nii (float32).
    Saves transformation matrices as tmats_struct or func.npy (float32)
    """
    # Unpack functional paths
    if functional_channel_paths is None:
        functional_path_one = None
        functional_path_two = None
    elif len(functional_channel_paths) == 1:
        functional_path_one = functional_channel_paths[0]
        functional_one_proxy = nib.load(functional_path_one)
        functional_data_one_final = np.empty(functional_one_proxy.shape)
        functional_path_two = None
    elif len(functional_channel_paths) == 2:
        functional_path_one = functional_channel_paths[0]
        functional_one_proxy = nib.load(functional_path_one)
        functional_data_one_final = np.empty(functional_one_proxy.shape)
        functional_path_two = functional_channel_paths[1]
        functional_two_proxy = nib.load(functional_path_two)
        functional_data_two_final = np.empty(functional_two_proxy.shape)
    else:
        # Fix this, should be identical to if!
        functional_path_one = None
        functional_path_two = None

    # Load moving proxy in this process
    moving_proxy = nib.load(moving_path) #xyzt
    n_timepoints = moving_proxy.shape[-1]
    # Determine number of slices
    if len(moving_proxy.shape)==4: # xyzt
        n_slices = moving_proxy.shape[2]
    elif len(moving_proxy.shape)==3: # xyt
        n_slices = 1
    else:
        print('Error: Could not determine number of slices in moving data.')
    
    moving_data_final = np.empty(moving_proxy.shape)

    tmats_final = np.empty([3,3,n_slices,n_timepoints]) # for rigid, 3x3 tmat for each slice and timepoint
    # for a rigid transformation, tmat is a 3x3 matrix [[cos(r) −sin(r) tx],[sin(r) cos(r) ty],[0 0 1]] where r is angle, t is translation
    for slice in range(n_slices):
        # Keeping track of time
        t_function_start = time.time()
        if n_slices==1: # data is single plane
            moving_data = np.squeeze(np.asarray(moving_proxy.dataobj[:,:,:],dtype='float32')) #source data xyt
        else: # data is volume
            moving_data = np.squeeze(np.asarray(moving_proxy.dataobj[:,:,slice,:],dtype='float32')) #source data xyzt into xyt
            moving_data = np.moveaxis(moving_data, -1, 0) #rearrange moving axes to t,x,y
        
        pr = ParaReg(reg_mode=moco_settings['reg_mode'],
                     smooth=moco_settings['smooth'],
                     avg_wid=moco_settings['avg_wid'],
                     n_proc=moco_settings['n_proc'],
                     mean_frames=moco_settings['moco_mean_frames']
                     )
        pr.register(moving_data)

        # apply transform, reorder axes back to xyt
        moving_data = pr.transform(moving_data)
        moving_data = np.moveaxis(moving_data, 0, -1)#rearrange moving axes back to x,y,t

        moving_data_final[:,:,slice,:] = moving_data
        

        if functional_path_one is not None:
            functional_data_one = np.squeeze(np.asarray(functional_one_proxy.dataobj[:,:,slice,:],dtype='float32')) #source data xyzt into xyt
            functional_data_one = np.moveaxis(functional_data_one, -1, 0) #rearrange moving axes to t,x,y

            # apply transform, reorder axes back to xyt
            functional_data_one = pr.transform(functional_data_one)
            functional_data_one = np.moveaxis(functional_data_one, 0, -1) #rearrange moving axes back to x,y,t
            functional_data_one_final[:,:,slice,:] = functional_data_one

            if functional_path_two is not None:
                functional_data_two = np.squeeze(np.asarray(functional_two_proxy.dataobj[:,:,slice,:],dtype='float32')) #source data xyzt into xyt
                functional_data_two = np.moveaxis(functional_data_two, -1, 0) #rearrange moving axes to t,x,y

                # apply transform, reorder axes back to xyt
                functional_data_two = pr.transform(functional_data_two)
                functional_data_two = np.moveaxis(functional_data_two, 0, -1) #rearrange moving axes back to x,y,t
                functional_data_two_final[:,:,slice,:] = functional_data_two

        # save transform params for slice to tmats_final
        for t_ind in range(n_timepoints):
            tmats_final[:,:,slice,t_ind] = np.asarray(pr._tmats[t_ind],dtype='float32')
            # for a rigid transformation, tmat is a 3x3 matrix [[cos(r) −sin(r) tx],[sin(r) cos(r) ty],[0 0 1]] where r is angle, t is translation

        print('Motion correction for ' + moving_path.as_posix()
            + 'at slice ' + str(slice) + ' took : '
            + repr(round(time.time() - t_function_start, 1))
            + 's\n')
    
    # save tmats_final to tmats.npy file
    np.save(pathlib.Path(par_output), np.squeeze(tmats_final))

    # save nifti files
    if n_slices==1:
        aff = np.eye(3)
    else:
        aff = np.eye(4)

    nib.Nifti1Image(np.squeeze(moving_data_final), aff).to_filename(moving_output_path)
    
    if functional_path_one is not None:
        nib.Nifti1Image(np.squeeze(functional_data_one_final), aff).to_filename(functional_channel_output_paths[0])
        
        if functional_path_two is not None:
            nib.Nifti1Image(np.squeeze(functional_data_two_final), aff).to_filename(functional_channel_output_paths[1])

    #derive recording metadata path for moco plot - from moco directory

    metadata_dir = os.path.join(moving_output_path.parent.parent, 'imaging')
    moco_dir = os.path.dirname(moving_output_path)

    moco_utils.save_moco_figure_stackreg_rigid(
    transform_matrix=tmats_final,
    metadata_dir=metadata_dir,
    moco_dir=moco_dir,
    printlog=printlog,
    n_slices=n_slices
    )
    

if __name__ == '__main__':
    ############################
    ### Organize shell input ###
    ############################
    parser = argparse.ArgumentParser()
    parser.add_argument("--fly_directory", help="Folder of fly to save log")
    parser.add_argument("--dataset_path", nargs="?", help="Folder pointing 'preprocessed'")

    parser.add_argument("--brain_paths_ch1", nargs="?", help="Path to ch1 file, if it exists")
    parser.add_argument("--brain_paths_ch2", nargs="?", help="Path to ch2 file, if it exists")
    parser.add_argument("--brain_paths_ch3", nargs="?", help="Path to ch3 file, if it exists")

    parser.add_argument("--STRUCTURAL_CHANNEL", nargs="?", help="variable with string containing the structural channel")
    parser.add_argument("--FUNCTIONAL_CHANNELS", nargs="?", help="list with strings containing the functional channel")

    parser.add_argument("--moco_path_ch1", nargs="?", help="Path to ch1 moco corrected file, if Ch1 exists")
    parser.add_argument("--moco_path_ch2", nargs="?", help="Path to ch2 moco corrected file, if Ch2 exists")
    parser.add_argument("--moco_path_ch3", nargs="?", help="Path to ch3 moco corrected file, if Ch3 exists")

    parser.add_argument("--par_output", nargs="?", help="Path to parameter output")

    parser.add_argument("--moco_temp_folder", nargs="?", help="Where to save the temp file")

    # moco settings
    parser.add_argument("--moco_transform_type", nargs="?", help="Type of transformation to use") # default is rigid
    parser.add_argument("--moco_smooth", nargs="?", help="Whether to smooth the data") # default is True for func, False for anat
    parser.add_argument("--moco_avg_wid", nargs="?", help="Width of the average") # if smooth=True, default is 3 for func, 1 for anat
    parser.add_argument("--moco_mean_frames", nargs="?", help="Number of frames to average for moco fixed brain") # default is 40 for func, 40 for anat
    parser.add_argument("--cores", nargs="?", help="Number of cores to use") # default is 40

    args = parser.parse_args()

    print('args: ' + repr(args))

    par_output = args.par_output

    #####################
    ### SETUP LOGGING ###
    #####################
    fly_directory = pathlib.Path(args.fly_directory)
    logfile = utils.create_logfile(fly_directory, function_name="moco_parallel")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, "moco_parallel")
    ####################################
    ### Identify the structural channel ###
    ####################################
    # Normally, we would have one structural channel
    if args.STRUCTURAL_CHANNEL is not None:
        if 'channel_1' == args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.brain_paths_ch1)
            moving_output_path = pathlib.Path(args.moco_path_ch1)
        elif 'channel_2' ==  args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.brain_paths_ch2)
            moving_output_path = pathlib.Path(args.moco_path_ch2)
        elif 'channel_3' ==  args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.brain_paths_ch3)
            moving_output_path = pathlib.Path(args.moco_path_ch3)
        # Convert the string represenation of a list to a list - it's either
        # ['channel_1',] or ['channel_1','channel_2'] or similar
        functional_channel_paths = []
        functional_channel_output_paths = []
        if 'channel_1' in args.FUNCTIONAL_CHANNELS and 'channel_1' not in args.STRUCTURAL_CHANNEL:
            # It possible to record i.e. anat scan with ONLY the structural
            # marker. Hence, check here if the functional channel exists.
            if args.brain_paths_ch1 is not None:
                functional_channel_paths.append(pathlib.Path(args.brain_paths_ch1))
                functional_channel_output_paths.append(pathlib.Path(args.moco_path_ch1))
            else:
                print('Information: Channel 1 in func channel, but this channel does not exist')
        if 'channel_2' in args.FUNCTIONAL_CHANNELS and 'channel_2' not in args.STRUCTURAL_CHANNEL:
            if args.brain_paths_ch2 is not None:
                functional_channel_paths.append(pathlib.Path(args.brain_paths_ch2))
                functional_channel_output_paths.append(pathlib.Path(args.moco_path_ch2))
            else:
                print('Information: Channel 2 in func channel, but this channel does not exist')
        if 'channel_3' in args.FUNCTIONAL_CHANNELS and 'channel_3' not in args.STRUCTURAL_CHANNEL:
            if args.brain_paths_ch3 is not None:
                functional_channel_paths.append(pathlib.Path(args.brain_paths_ch3))
                functional_channel_output_paths.append(pathlib.Path(args.moco_path_ch3))
            else:
                print('Information: Channel 3 in func channel, but this channel does not exist: ')
        #else:
        # If no functional path has been found, length of list == 0 and set both to None
        # to explictly state that no functional path has been defined
        if len(functional_channel_paths) == 0:
            functional_channel_paths = None
            functional_channel_output_paths = None
            print('No functional channel defined!')
    else:
        printlog('"structural_channel" NOT DEFINED!!! \n'
                 'You must define the "structural_channel" in the "fly.json" file!')
    """
    # With the 'STRUCTURAL_CHANNEL' we are now enforcing that the user must define
    # the channel to be used as the structural channel (previously 'anatomy_channel').
    # Even if you don't have a dedicated structural marker/channel, you must define
    # one of the functional channels as the structural channel.  Moco will be done on the structural
    # channel and the resulting transforms will be applied to any other functional channels
    """

    print("moving_path" + repr(moving_path))
    print("moving_output_path" + repr(moving_output_path))

    ########################################
    ### Multiprocessing code starts here ###
    ########################################

    #cores = 40 (set in moco settings in fly.json now)

    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    n_timepoints = brain_shape[-1]
    if len(brain_shape)==4: # xyzt
        n_slices = brain_shape[2]
    elif len(brain_shape)==3: # xyt
        n_slices = 1


    # moco settings
        # unpack transform type from string
    if args.moco_transform_type == 'StackReg.RIGID_BODY':
        moco_transform_type = StackReg.RIGID_BODY
    moco_settings = {}
    moco_settings['reg_mode'] = moco_transform_type
    if args.moco_smooth == 'True':
        moco_settings['smooth'] = True
    elif args.moco_smooth == 'False':
        moco_settings['smooth'] = False
    else:
        print('Error: could not interpret moco_smooth, should be True or False, datatype=string')
    moco_settings['avg_wid'] = int(args.moco_avg_wid)
    moco_settings['n_proc'] = int(args.cores)
    moco_settings['moco_mean_frames'] = int(args.moco_mean_frames)

    print('moco_settings: ' + repr(moco_settings))

    print('Will perform motion correction on a total of ' + repr(n_timepoints) + ' timepoints and ' + repr(n_slices) + ' slice(s).')

    print('Starting MOCO')
    time_start = time.time()

    # DO MOCO
    moco_slice(
             moving_path=moving_path,
             functional_channel_paths=functional_channel_paths,
             par_output = par_output,
             moco_settings = moco_settings,
             )
    
    print('Took: ' + repr(round(time.time() - time_start,1)) + 's\n to motion correct files')
    print('Motion correction done.')