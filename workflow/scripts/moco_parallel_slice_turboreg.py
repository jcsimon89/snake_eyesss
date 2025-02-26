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
    slice,
    fixed_path,
    moving_path,
    functional_channel_paths,
    temp_save_path,
    fly_directory
):
    """
    Loop doing the motion correction for a given set of index.
    This is the function that is doing the heavy lifting of the multiprocessing
    Saves the chunks as npy files in the 'temp_save_path' folder with the index contained
    encoded in the filename.
    :param index:
    :return:
    """

    # Keeping track of time
    t_function_start = time.time()
    # Load meanbrain (fixed)
    fixed_proxy = nib.load(fixed_path)
    fixed_data = np.squeeze(np.asarray(fixed_proxy.dataobj[:,:,slice],dtype='float32')) #source data xyz into xy


    # Load moving proxy in this process
    moving_proxy = nib.load(moving_path)
    moving_data = np.squeeze(np.asarray(moving_proxy.dataobj[:,:,slice,:],dtype='float32')) #source data xyzt into xyt
    moving_data = np.moveaxis(moving_data,[0, 2], [2, 0]) #rearrange moving axes to t,x,y

    pr = ParaReg(eg_mode=StackReg.RIGID_BODY, smooth=smooth, avg_wid=avg_wid, n_proc=cores)
    pr.register(moving_data,fixed_data)

    # apply transform, reorder axes back to xyt
    moving_data = pr.transform(moving_data)
    moving_data = np.moveaxis(moving_data,[0, 2], [2, 0]) #rearrange moving axes back to x,y,t

    # save warped slice in temp folder
    np.save(pathlib.Path(temp_save_path, moving_path.name + '_slice{}_'.format(slice)),
            moving_data)

    # Unpack functional paths
    if functional_channel_paths is None:
        functional_path_one = None
        functional_path_two = None
    elif len(functional_channel_paths) == 1:
        functional_path_one = functional_channel_paths[0]
        functional_path_two = None
    elif len(functional_channel_paths) == 2:
        functional_path_one = functional_channel_paths[0]
        functional_path_two = functional_channel_paths[1]
    else:
        # Fix this, should be identical to if!
        functional_path_one = None
        functional_path_two = None
    if functional_path_one is not None:
        # Load functional one proxy in this process
        functional_one_proxy = nib.load(functional_path_one)
        if functional_path_two is not None:
            functional_two_proxy = nib.load(functional_path_two)

    if functional_path_one is not None:
        functional_data_one = np.squeeze(np.asarray(functional_one_proxy.dataobj[:,:,slice,:],dtype='float32')) #source data xyzt into xyt
        functional_data_one = np.moveaxis(functional_data_one,[0, 2], [2, 0]) #rearrange moving axes to t,x,y

        # apply transform, reorder axes back to xyt
        functional_data_one = pr.transform(functional_data_one)
        functional_data_one = np.moveaxis(functional_data_one,[0, 2], [2, 0]) #rearrange moving axes back to x,y,t

        # save warped slice in temp folder
        np.save(pathlib.Path(temp_save_path, functional_path_one.name + '_slice{}_'.format(slice)),
            functional_data_one)

        if functional_path_two is not None:
            functional_data_two = np.squeeze(np.asarray(functional_two_proxy.dataobj[:,:,slice,:],dtype='float32')) #source data xyzt into xyt
            functional_data_two = np.moveaxis(functional_data_two,[0, 2], [2, 0]) #rearrange moving axes to t,x,y

            # apply transform, reorder axes back to xyt
            functional_data_two = pr.transform(functional_data_two)
            functional_data_two = np.moveaxis(functional_data_two,[0, 2], [2, 0]) #rearrange moving axes back to x,y,t
            
            # save warped slice in temp folder
            np.save(pathlib.Path(temp_save_path, functional_path_two.name + '_slice{}_'.format(slice)),
            functional_data_two)

    #TODO: save transform params

    print('Motion correction for ' + moving_path.as_posix()
          + 'at slice ' + str(slice) + ' took : '
          + repr(round(time.time() - t_function_start, 1))
          + 's\n')


def find_missing_temp_files(fixed_path,
                            moving_path,
                            functional_channel_paths,
                            temp_save_path,
                            fly_directory
                            ):
    """
    It can happen that registration of one or a few images fails.
    This function takes the temp folder and runs through it and
    tries to identify missing files (e.g. if 'index_0', 'index_1' and 'index_3'
    exist and 'index_2' is missing.
    If it find a missing frame, call the motion_correction function on that
    index, wait and try again.
    :param temp_save_path:
    :return:
    """

    start_time = time.time()
    # A list to keep track of missing files!
    slice_of_missing_files = []

    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    nslices = brain_shape[2] #assumes data is x,y,z,t

    for slice in range(nslices):
        slice_tracker = 0
        for current_file in natsort.natsorted(temp_save_path.iterdir()):
            #print('Finding missing files: current_file ' + current_file.name)
            # Check if moving_path.name, for example channel_1.nii is in filename
            # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
            if '.npy' in current_file.name and moving_path.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                # Extract index number and slice number
                if slice == slice_tracker:
                    # Great!
                    pass
                else:
                    # in case more than one file (e.g. 1 & 2) are missing!
                    while slice > slice_tracker:
                        print('Missing slice: ' + repr(slice_tracker))
                        slice_of_missing_files.append(slice_tracker)
                        slice_tracker+=1 #
                # add 1 to be prepared for the next loop!
                slice_tracker+=1
        # it's possible that we are missing only functional files but not anatomical files.
        # Also collect those
        if functional_channel_paths is None:
            functional_path_one = None
            functional_path_two = None
        elif len(functional_channel_paths) == 1:
            functional_path_one = functional_channel_paths[0]
            functional_path_two = None
        elif len(functional_channel_paths) == 2:
            functional_path_one = functional_channel_paths[0]
            functional_path_two = functional_channel_paths[1]
        else:
            functional_path_one = None
            functional_path_two = None
        if functional_path_one is not None:
            slice_tracker = 0
            for current_file in natsort.natsorted(temp_save_path.iterdir()):
                if '.npy' in current_file.name and functional_path_one.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                    if slice == slice_tracker:
                        # Great!
                        pass
                    else:
                        # in case more than one file (e.g. 1 & 2) are missing!
                        while slice > slice_tracker:
                            print('Missing slice: ' + repr(slice_tracker))
                            slice_of_missing_files.append(slice_tracker) 
                            slice_tracker += 1  #TODO: change for multiple time indices
                    #add 1 to be prepared for the next loop!
                    slice_tracker += 1
        if functional_path_two is not None:
            slice_tracker = 0
            for current_file in natsort.natsorted(temp_save_path.iterdir()):
                if '.npy' in current_file.name and functional_path_two.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                    if slice == slice_tracker:
                        # Great!
                        pass
                    else:
                        # in case more than one file (e.g. 1 & 2) are missing!
                        while slice > slice_tracker:
                            slice_of_missing_files.append(slice_tracker)
                            slice_tracker += 1  ##TODO: change for multiple time indices
                    # Once slice == slice_tracker, add 1 to be prepared for the next loop!
                    slice_tracker += 1
                    # remove duplicate entries
    slice_of_missing_files = dict.fromkeys(slice_of_missing_files)
    slice_of_missing_files = list(slice_of_missing_files.keys())

    # loop through index_of_missing_files. If it's an empty list, don't loop and skip
    if len(slice_of_missing_files) > 0:
        print('==========================================================================================')
        print('WARNING: Not all files that should have been created in the temp folder have been create')
        print('This might be due to a memory (RAM) error in some registration calls.')
        print('Will now try to run the missing files but might run out of allocated time as this is done serially.')
        print('SUGGESTION: Increase available RAM for the motion correction call!')
        print('==========================================================================================')
        for i in range(len(slice_of_missing_files)):
            # Call motion_correction function on index of missing files
            # THIS IS SLOW AS IT'S NOT PARALLELIZED. Hopefully this only is used
            # in very rare circumstances.
            current_slice = slice_of_missing_files[i]
            print('Missing slice, currently working on: slice ' + str(current_slice))
            moco_slice(
                            slice=current_slice,
                            fixed_path=fixed_path,
                            moving_path=moving_path,
                            functional_channel_paths=functional_channel_paths,
                            temp_save_path=temp_save_path,
                            fly_directory=fly_directory)
    print('Checking for missing files took: ' + repr(round(time.time() - start_time,2)))

def combine_temp_files(moving_path,
                       functional_channel_paths,
                       temp_save_path,
                       moving_output_path,
                       functional_channel_output_paths,
                       param_output_path):
    """
    This function crawls through the temp_save_path and stitched the individual frames
    together.

    In order to save RAM, it does one imaging file after the other.

    :param moving_path:
    :param functional_channel_paths:
    :param temp_save_path:
    :param moving_output_path:
    :param functional_channel_output_paths:
    :param param_output_path:
    :return:
    """
    time_start = time.time()
    ####
    # STITCH STRUCTURAL_CHANNEL
    ####
    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    nslices = brain_shape[2]
    # Preallocate array for anatomy...
    stitched_anatomy_brain = np.zeros((brain_shape[0],brain_shape[1],
                                       brain_shape[2], brain_shape[3]),
                                      dtype=np.float32)
    #TODO: transform matrix
    # Loop through all files. Because it's sorted we don't have to worry about
    # the index!
    printlog('Start combining ' + moving_output_path.name)
    
    for slice in range(nslices):
        slice_tracker = 0
        for current_file in natsort.natsorted(temp_save_path.iterdir()):
            # Check if moving_path.name, for example channel_1.nii is in filename
            # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
            if '.npy' in current_file.name and moving_path.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                stitched_anatomy_brain[:,:,slice,:] = np.load(current_file)
                slice_tracker += 1
            
            # and collect motcorr_params, this is tiny so no worries about space here
            # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
            #TODO: transform parameters
            # elif 'motcorr_params' in current_file.name and '_slice{}_'.format(slice) in current_file.name:
            #     index = moco_utils.index_from_filename(current_file)
            #     transform_matrix[slice,index,:] = np.load(current_file)


    if slice_tracker != nslices:
            print('There seems to be a problem with the temp files for: ' + repr(current_file))
            print('only {} out of {} slices found'.format(slice_tracker,slice))

    # SAVE
    # we create a new subfolder called 'moco' where the file is saved
    # Not a great solution - the definition of the output file should happen in the snakefile, not hidden
    # in here!
    savepath_root = pathlib.Path(moving_path.parents[1], 'moco')
    savepath_root.mkdir(parents=True, exist_ok=True)
    # Prepare nifti file
    aff = np.eye(4)
    # Save anatomy channel: - unfortunately this is super memory intensive! We
    # make a copy of the huge array before we are able so save it.
    #stitched_anatomy_brain_nifty = nib.Nifti1Image(stitched_anatomy_brain, aff)
    #stitched_anatomy_brain_nifty.to_filename(moving_output_path)
    # Try this and check if we still need 4 times memory
    nib.Nifti1Image(stitched_anatomy_brain, aff).to_filename(moving_output_path)


    del stitched_anatomy_brain # Explicitly release the memory (it might not be
    # immediately released, though. Test this.
    #del stitched_anatomy_brain_nifty

    printlog('Done combining and saving ' + moving_output_path.name)

    #####
    # Work on functional paths
    #####
    # Unpack functional paths
    if functional_channel_paths is None:
        functional_path_one = None
        functional_path_two = None
    elif len(functional_channel_paths) == 1:
        functional_path_one = functional_channel_paths[0]
        functional_path_two = None
    elif len(functional_channel_paths) == 2:
        functional_path_one = functional_channel_paths[0]
        functional_path_two = functional_channel_paths[1]
    else: # todo fix, should be same as the 'if' above.
        functional_path_one = None
        functional_path_two = None


    ####
    # STITCH FUNCTIONAL 1
    ####

    if functional_path_one is not None:
        printlog('Start combining ' + functional_channel_output_paths[0].name)
        # Essentially identical to the code above
        stitched_functional_one = np.zeros((brain_shape[0],brain_shape[1],
                                       brain_shape[2], brain_shape[3]),
                                      dtype=np.float32)
        nslices = brain_shape[2]
        for slice in range(nslices):
                slice_tracker = 0
                for current_file in natsort.natsorted(temp_save_path.iterdir()):
                    # Check if moving_path.name, for example channel_1.nii is in filename
                    # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
                    if '.npy' in current_file.name and moving_path.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                        stitched_anatomy_brain[:,:,slice,:] = np.load(current_file)
                        slice_tracker += 1

        if slice_tracker != nslices:
                print('There seems to be a problem with the temp files for: ' + repr(current_file))
                print('only {} out of {} slices found'.format(slice_tracker,slice))


        #savepath_func_one = pathlib.Path(savepath_root, functional_path_one.stem + '_moco.nii')
        #stitched_functional_one_nifty = nib.Nifti1Image(stitched_functional_one, aff)
        #stitched_functional_one_nifty.to_filename(functional_channel_output_paths[0])
        nib.Nifti1Image(stitched_functional_one, aff).to_filename(functional_channel_output_paths[0])

        del stitched_functional_one
        #del stitched_functional_one_nifty

        printlog('Done combining and saving ' + functional_channel_output_paths[0].name)
        ####
        # STITCH FUNCTIONAL 2
        ####
        if functional_path_two is not None:
            stitched_functional_two = np.zeros((brain_shape[0],brain_shape[1],
                                           brain_shape[2], brain_shape[3]),
                                          dtype=np.float32)
            # Collect data for second functional channel
            printlog('Start combining ' + functional_channel_output_paths[1].name)
            nslices = brain_shape[2]
            for slice in range(nslices):
                    slice_tracker = 0
                    for current_file in natsort.natsorted(temp_save_path.iterdir()):
                        # Check if moving_path.name, for example channel_1.nii is in filename
                        # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
                        if '.npy' in current_file.name and moving_path.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                            stitched_anatomy_brain[:,:,slice,:] = np.load(current_file)
                            slice_tracker += 1

            if slice_tracker != nslices:
                    print('There seems to be a problem with the temp files for: ' + repr(current_file))
                    print('only {} out of {} slices found'.format(slice_tracker,slice))
            # Save the nifty file
            # stitched_functional_two_nifty = nib.Nifti1Image(stitched_functional_two, aff)
            # stitched_functional_two_nifty.to_filename(functional_channel_output_paths[1])
            nib.Nifti1Image(stitched_functional_two, aff).to_filename(functional_channel_output_paths[1])

            del stitched_functional_two
            #del stitched_functional_two_nifty

            printlog('Done combining and saving ' + functional_channel_output_paths[1].name)

    # After saving the stitched file, delete the temporary files
    shutil.rmtree(temp_save_path)

    # Save transform matrix:
    #param_savename = savepath_root, 'motcorr_params.npy'
    # np.save(param_output_path, transform_matrix) #TODO: save transforms

    print('Took: ' + repr(time.time() - time_start) + 's to combine files')

    t0 = time.time()
    # ### MAKE MOCO PLOT ###
    # moco_utils.save_moco_figure( #TODO: change for slice moco
    #     transform_matrix=transform_matrix,
    #     parent_path=moving_path.parent,
    #     moco_dir=moving_output_path.parent,
    #     printlog=printlog,
    # )

    # print('took: ' + repr(round(time.time() - t0,2)) + ' s to plot moco')


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
    if TESTING:
        cores = 40

    smooth = False
    avg_wid = 10

    # create an index going from [0,1,...,n]
    time_index = moco_utils.prepare_time_index(moving_path)
    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    nslices = brain_shape[2] #assumes data is x,y,z,t

    print('Will perform motion correction on a total of ' + repr(len(time_index)) + ' timepoints and ' + repr(nslices) + ' slices.')
    if TESTING:
        time_index = list(range(0,100,1))#[0,1,2,3,4,5,6,7]

    time_start = time.time()
    print('Starting MOCO')

    # DO MOCO
    for current_slice in range(nslices):
        moco_slice(slice=current_slice,
             fixed_path=fixed_path,
             moving_path=moving_path,
             functional_channel_paths=functional_channel_paths,
             temp_save_path=temp_save_path,
             fly_directory=fly_directory)


    print('Took: ' + repr(time.time() - time_start) + 's to motion correct files')
    print('Motion correction done. Checking for missing files now')

    find_missing_temp_files(fixed_path,
                            moving_path,
                            functional_channel_paths,
                            temp_save_path,
                            fly_directory
                            )

    print('Checked for missing files. See above if something was missing. Combining files now')
    combine_temp_files(moving_path, functional_channel_paths, temp_save_path,
                       moving_output_path, functional_channel_output_paths, param_output_path)
    print('files combined')