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
TESTING = False

def moco_slice(
    n_proc,
    fixed_path,
    moving_path,
    functional_channel_paths,
):
    """
    Loop doing the motion correction for a given set of index.
    This is the function that is doing the heavy lifting of the multiprocessing
    Saves the chunks as npy files in the 'temp_save_path' folder with the index contained
    encoded in the filename.
    :param index:
    :return:
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
    
    # Load meanbrain (fixed)
    fixed_proxy = nib.load(fixed_path)
    n_slices = fixed_proxy.shape[-1] #xyz
    # Load moving proxy in this process
    moving_proxy = nib.load(moving_path)
    
    moving_data_final = np.empty(moving_proxy.shape)

    for slice in range(n_slices):
        # Keeping track of time
        t_function_start = time.time()
        fixed_data = np.squeeze(np.asarray(fixed_proxy.dataobj[:,:,slice],dtype='float32')) #source data xyz into xy

        moving_data = np.squeeze(np.asarray(moving_proxy.dataobj[:,:,slice,:],dtype='float32')) #source data xyzt into xyt
        moving_data = np.moveaxis(moving_data, -1, 0) #rearrange moving axes to t,x,y

        pr = ParaReg(reg_mode=StackReg.RIGID_BODY, smooth=smooth, avg_wid=avg_wid, n_proc=n_proc)
        print('size of moving_data: ' + str(moving_data.shape))
        print('size of fixed_data: ' + str(fixed_data.shape))
        pr.register(moving_data,fixed_data)

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

        #TODO: save transform params

        print('Motion correction for ' + moving_path.as_posix()
            + 'at slice ' + str(slice) + ' took : '
            + repr(round(time.time() - t_function_start, 1))
            + 's\n')
    
    # save nifti files
    aff = np.eye(4)
    nib.Nifti1Image(moving_data_final, aff).to_filename(moving_output_path)
    
    if functional_path_one is not None:
        nib.Nifti1Image(functional_data_one_final, aff).to_filename(functional_channel_output_paths[0])
        
        if functional_path_two is not None:
            nib.Nifti1Image(functional_data_two_final, aff).to_filename(functional_channel_output_paths[1])


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

    parser.add_argument("--mean_brain_paths_ch1", nargs="?", help="Path to ch1 meanbrain file, if it exists")
    parser.add_argument("--mean_brain_paths_ch2", nargs="?", help="Path to ch2 meanbrain file, if it exists")
    parser.add_argument("--mean_brain_paths_ch3", nargs="?", help="Path to ch3 meanbrain file, if it exists")

    parser.add_argument("--STRUCTURAL_CHANNEL", nargs="?", help="variable with string containing the structural channel")
    parser.add_argument("--FUNCTIONAL_CHANNELS", nargs="?", help="list with strings containing the functional channel")

    parser.add_argument("--moco_path_ch1", nargs="?", help="Path to ch1 moco corrected file, if Ch1 exists")
    parser.add_argument("--moco_path_ch2", nargs="?", help="Path to ch2 moco corrected file, if Ch2 exists")
    parser.add_argument("--moco_path_ch3", nargs="?", help="Path to ch3 moco corrected file, if Ch3 exists")

    parser.add_argument("--par_output", nargs="?", help="Path to parameter output")

    parser.add_argument("--moco_temp_folder", nargs="?", help="Where to save the temp file")

    args = parser.parse_args()

    #####################
    ### SETUP LOGGING ###
    #####################
    fly_directory = pathlib.Path(args.fly_directory)
    logfile = utils.create_logfile(fly_directory, function_name="moco_parallel")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, "moco_parallel")

    ####################################
    ### Identify the anatomy channel ###
    ####################################
    # Normally, we would have one anatomy channel
    if args.STRUCTURAL_CHANNEL is not None:
        if 'channel_1' == args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.brain_paths_ch1)
            fixed_path = pathlib.Path(args.mean_brain_paths_ch1)
            moving_output_path = pathlib.Path(args.moco_path_ch1)
        elif 'channel_2' ==  args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.brain_paths_ch2)
            fixed_path = pathlib.Path(args.mean_brain_paths_ch2)
            moving_output_path = pathlib.Path(args.moco_path_ch2)
        elif 'channel_3' ==  args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.brain_paths_ch3)
            fixed_path = pathlib.Path(args.mean_brain_paths_ch3)
            moving_output_path = pathlib.Path(args.moco_path_ch3)
        # If we have an antomy channel, we can use the parameters to
        # motion correct the (generally) less clearly visible functional channel(s)
        #if args.FUNCTIONAL_CHANNELS != ['']:
        # Convert the string represenation of a list to a list - it's either
        # ['channel_1',] or ['channel_1','channel_2'] or similar
        # if
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
    # the channel to be used as the structural channel (previously 'anatomy_channel')
    else:
        # However, sometimes, we don't have an anatomy channel (e.g.
        # when we only have GCAMP expressed and not anatomical marker)
        # Note that it'll just take the first channel as the channel to
        # perform moco and will ignore the rest. I.e. if your
        # 'functional_channel: "['channel_1', 'channel_2']" we currently
        # only do moco on channel_1
        if 'channel_1' in args.FUNCTIONAL_CHANNELS:
            moving_path = pathlib.Path(args.brain_paths_ch1)
            fixed_path = pathlib.Path(args.mean_brain_paths_ch1)
            moving_output_path = pathlib.Path(args.moco_path_ch1)
        elif 'channel_2' == args.FUNCTIONAL_CHANNELS:
            moving_path = pathlib.Path(args.brain_paths_ch2)
            fixed_path = pathlib.Path(args.mean_brain_paths_ch2)
            moving_output_path = pathlib.Path(args.moco_path_ch2)
        elif 'channel_3' == args.FUNCTIONAL_CHANNELS:
            moving_path = pathlib.Path(args.brain_paths_ch3)
            fixed_path = pathlib.Path(args.mean_brain_paths_ch3)
            moving_output_path = pathlib.Path(args.moco_path_ch2)
        functional_channel_paths = None
        functional_channel_output_paths = None
    """

    print("moving_path" + repr(moving_path))
    print("fixed_path" + repr(fixed_path))
    print("moving_output_path" + repr(moving_output_path))

    param_output_path = args.par_output

    ##########################
    ### ORGANIZE TEMP PATH ###
    ##########################
    # Path were the intermediate npy files are saved.
    # It is important that it's different for each run.
    # We can just put in on scratch
    # This will only work if we have a folder called trc and data is in /data, of course
    #relevant_temp_save_path_part = moving_path.as_posix().split('trc/data/')[-1]
    dataset_path = pathlib.Path(args.dataset_path)
    relevant_temp_save_path_part = moving_path.as_posix().split(dataset_path.as_posix())[-1]
    relevant_temp_save_path_part=relevant_temp_save_path_part[1::] # remove first forward slash
    ###################
    # DON'T CHANGE THIS-if this points to your actual experimental folder, the shutil.rmtree
    # below will DELETE YOUR DATA. THIS MUST BE A TEMPORARY PATH
    #temp_save_path = pathlib.Path('/scratch/groups/trc', relevant_temp_save_path_part).parent
    temp_save_path = pathlib.Path(args.moco_temp_folder, relevant_temp_save_path_part).parent
    ##################

    if TESTING:
        temp_save_path = pathlib.Path('C:/Users/jcsimon/.snakemake/temp')
        if temp_save_path.is_dir():
            shutil.rmtree(temp_save_path)

    if temp_save_path.is_dir():
        shutil.rmtree(temp_save_path)
    #else:
    #    print('WARNING: Did not remove files in ' + str(temp_save_path))
    #    print('Only remove folders that are on scratch to avoid accidentally deleting source data')

    # create dir
    # No need for exist_ok=True because the dir should have been deleted just now
    temp_save_path.mkdir(parents=True)

    ########################################
    ### Multiprocessing code starts here ###
    ########################################
    # Note on cores: I benchmarked using a 30 minute functional dataset (256,128,49,~3300)
    # with each channel about 10Gb and a large anatomical dataset (1024, 512, 100, 100)
    # with each channel about 25Gb. I didn't find a difference using 8 vs 16 cores and
    # the 16 cores try with anatomical actually threw a memory error (with 256Gb of RAM).
    # There seems to be some multiprocessing going on in ants.registration which might explain
    # why we don't see improved times during benchmarking.
    # Hence, I fixed the core count at 8.
    cores = 40

    # create an index going from [0,1,...,n]
    time_index = moco_utils.prepare_time_index(moving_path)
    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    nslices = brain_shape[2] #assumes data is x,y,z,t

    if TESTING:
        cores = 40
        time_index = list(range(0,1000,1))#[0,1,2,3,4,5,6,7]

    smooth = False
    avg_wid = 10


    print('Will perform motion correction on a total of ' + repr(len(time_index)) + ' timepoints and ' + repr(nslices) + ' slices.')

    time_start = time.time()
    print('Starting MOCO')

    # DO MOCO
    moco_slice(
             n_proc = cores,
             fixed_path=fixed_path,
             moving_path=moving_path,
             functional_channel_paths=functional_channel_paths
             )

    print('Motion correction done.')
    print('Took: ' + repr(round(time.time() - time_start,1)) + 's\n to motion correct files')