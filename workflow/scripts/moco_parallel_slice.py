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
type_of_transform = "SyN"
flow_sigma = 3
total_sigma = 0
aff_metric = 'mattes'

TESTING = False

def moco_slice(
    index,
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
    # Load meanbrain (fixed) and put into ants format
    fixed_proxy = nib.load(fixed_path)
    # Load data to memory. Dtypes is a bit confusing here: Meanbrain comes as uint16...
    fixed_data = fixed_proxy.dataobj[:,:,slice]
    # However, ants seems to require float32 (I think)
    fixed_ants = ants.from_numpy(np.asarray(fixed_data, dtype=np.float32))

    # Load moving proxy in this process
    moving_proxy = nib.load(moving_path)

    # Load data in a given process
    current_moving = moving_proxy.dataobj[:,:,slice,index]
    # Convert to ants images
    moving_ants = ants.from_numpy(np.asarray(current_moving, dtype=np.float32))

    #t0 = time.time()
    # Perform the registration
    moco = ants.registration(fixed_ants, moving_ants,
                             type_of_transform=type_of_transform,
                             flow_sigma=flow_sigma,
                             total_sigma=total_sigma,
                             aff_metric=aff_metric)
    #print('Registration took ' + repr(time.time() - t0) + 's')

    # Save warped image in temp_save_path with index in filename.
    np.save(pathlib.Path(temp_save_path, moving_path.name + '_slice' + str(slice) + '_index'
                         + str(index)),
            moco["warpedmovout"].numpy())

    #t0 = time.time()
    # Next, use the transform info for the functional image
    transformlist = moco["fwdtransforms"]

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
        current_functional_one = functional_one_proxy.dataobj[:,:,slice,index]
        moving_frame_one_ants = ants.from_numpy(np.asarray(current_functional_one, dtype=np.float32))
        # to motion correction for functional image
        moving_frame_one_ants = ants.apply_transforms(fixed_ants, moving_frame_one_ants, transformlist)
        # put moco functional image into preallocated array
        #moco_functional_one[:, :, :, counter] = moving_frame_one_ants.numpy()
        #print('apply transforms took ' + repr(time.time() - t0) + 's')
        np.save(pathlib.Path(temp_save_path, functional_path_one.name + '_slice' + str(slice)+ '_index'
                             + str(index)),
                moving_frame_one_ants.numpy())

        if functional_path_two is not None:
            current_functional_two = functional_two_proxy.dataobj[:,:,slice, index]
            moving_frame_two_ants = ants.from_numpy(np.asarray(current_functional_two, dtype=np.float32))
            moco_functional_two = ants.apply_transforms(fixed_ants, moving_frame_two_ants, transformlist)
            #moco_functional_two[:,:,:, counter] = moco_functional_two.numpy()
            np.save(pathlib.Path(temp_save_path, functional_path_two.name + '_slice' + str(slice) + '_index'
                                 + str(index)),
                    moco_functional_two.numpy())

    #t0=time.time()
    # delete writen files:
    # Delete transform info - might be worth keeping instead of huge resulting file? TBD
    for x in transformlist:
        if ".mat" in x:
            # Keep transform_matrix, I think this is used to make the plot
            # called 'motion_correction.png'
            temp = ants.read_transform(x)
            #transform_matrix[counter, :] = temp.parameters
            param_savename = pathlib.Path(temp_save_path, "motcorr_params" + '_slice' + str(slice) + '_index'
                                          + str(index))
            np.save(param_savename, temp.parameters) # that's the transform_matrix in brainsss

        # lets' delete all files created by ants - else we quickly create thousands of files!
        pathlib.Path(x).unlink()
    print('Motion correction for ' + moving_path.as_posix()
          + 'at slice ' + str(slice) + ', at index ' + str(index) + ' took : '
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
    slice_and_index_of_missing_files = []

    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    nslices = brain_shape[2] #assumes data is x,y,z,t

    for slice in range(nslices):
        index_tracker = 0
        for current_file in natsort.natsorted(temp_save_path.iterdir()):
            #print('Finding missing files: current_file ' + current_file.name)
            # Check if moving_path.name, for example channel_1.nii is in filename
            # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
            if '.npy' in current_file.name and moving_path.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                # Extract index number and slice number
                index = moco_utils.index_from_filename(current_file)
                if index == index_tracker:
                    # Great!
                    pass
                else:
                    # in case more than one file (e.g. 1 & 2) are missing!
                    while index > index_tracker:
                        print('Missing files: ' + repr(index_tracker))
                        slice_and_index_of_missing_files.append(tuple([slice,index_tracker]))
                        index_tracker+=1 #
                # Once index == index_tracker, add 1 to be prepared for the next loop!
                index_tracker+=1
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
            index_tracker = 0
            for current_file in natsort.natsorted(temp_save_path.iterdir()):
                if '.npy' in current_file.name and functional_path_one.name in current_file.name:
                    # Extract index number
                    index = moco_utils.index_from_filename(current_file)
                    if index == index_tracker:
                        # Great!
                        pass
                    else:
                        # in case more than one file (e.g. 1 & 2) are missing!
                        while index > index_tracker:
                            slice_and_index_of_missing_files.append(tuple([slice,index_tracker]))
                            index_tracker += 1  #
                    # Once index == index_tracker, add 1 to be prepared for the next loop!
                    index_tracker += 1
        if functional_path_two is not None:
            index_tracker = 0
            for current_file in natsort.natsorted(temp_save_path.iterdir()):
                if '.npy' in current_file.name and functional_path_two.name in current_file.name:
                    # Extract index number
                    index = moco_utils.index_from_filename(current_file)
                    if index == index_tracker:
                        # Great!
                        pass
                    else:
                        # in case more than one file (e.g. 1 & 2) are missing!
                        while index > index_tracker:
                            slice_and_index_of_missing_files.append(tuple([slice,index_tracker]))
                            index_tracker += 1  #
                    # Once index == index_tracker, add 1 to be prepared for the next loop!
                    index_tracker += 1
                    # remove duplicate entries
        slice_and_index_of_missing_files = np.unique(np.asarray(slice_and_index_of_missing_files))

    # loop through index_of_missing_files. If it's an empty list, don't loop and skip
    if len(slice_and_index_of_missing_files) > 0:
        print('==========================================================================================')
        print('WARNING: Not all files that should have been created in the temp folder have been create')
        print('This might be due to a memory (RAM) error in some ants.registration calls.')
        print('Will now try to run the missing files but might run out of allocated time as this is done serially.')
        print('SUGGESTION: Increase available RAM for the motion correction call!')
        print('==========================================================================================')
        for i in range(len(slice_and_index_of_missing_files)):
            # Call motion_correction function on index of missing files
            # THIS IS SLOW AS IT'S NOT PARALLELIZED. Hopefully this only is used
            # in very rare circumstances.
            current_slice = slice_and_index_of_missing_files[i][0]
            current_index = slice_and_index_of_missing_files[i][1]
            print('Missing (slice,index) currently working on: (' + str(current_slice) + ',' + str(current_index) + ')')
            moco_slice(current_index,
                                current_slice,
                                fixed_path,
                                moving_path,
                                functional_channel_paths,
                                temp_save_path,
                                fly_directory)
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
    # ...and transform matrix.
    transform_matrix = np.zeros((nslices,brain_shape[3],6))
    # Loop through all files. Because it's sorted we don't have to worry about
    # the index!
    printlog('Start combining ' + moving_output_path.name)
    
    for slice in range(nslices):
        index_tracker = 0
        for current_file in natsort.natsorted(temp_save_path.iterdir()):
            # Check if moving_path.name, for example channel_1.nii is in filename
            # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
            if '.npy' in current_file.name and moving_path.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                # Extract index number
                index = moco_utils.index_from_filename(current_file)
                stitched_anatomy_brain[:,:,slice,index] = np.load(current_file)
                # Just a sanity check! E.g. for first image we expect '0'
                if index_tracker != index:
                    print('There seems to be a problem with the temp files for: ' + repr(current_file))
                    print('Previous index was ' + str(index_tracker - 1))
                    print('Current slice is {}'.format(slice))
                    print('Next index (based on filename) was ' + str(index))

                index_tracker += 1

            # and collect motcorr_params, this is tiny so no worries about space here
            # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
            elif 'motcorr_params' in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                index = moco_utils.index_from_filename(current_file)
                transform_matrix[slice,index,:] = np.load(current_file)
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
            index_tracker = 0
            for current_file in natsort.natsorted(temp_save_path.iterdir()):
                    # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
                    if '.npy' in current_file.name and functional_path_one.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                        index = moco_utils.index_from_filename(current_file)
                        stitched_functional_one[:,:,slice,index] = np.load(current_file)
                        # Just a sanity check! E.g. for first image we expect '0'
                        if index_tracker != index:
                            print('There seems to be a problem with the temp files for: ' + str(current_file))
                            print('Previous index was ' + str(index_tracker - 1))
                            print('Current slice is {}'.format(slice))
                            print('Next index (based on filename) was ' + str(index))
                        index_tracker += 1

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
                index_tracker = 0
                for current_file in natsort.natsorted(temp_save_path.iterdir()):
                    # NOTE: '_slice{}_' is essential, otherwis searching for slice1 returns slice1 and slice10
                    if '.npy' in current_file.name and functional_path_two.name in current_file.name and '_slice{}_'.format(slice) in current_file.name:
                        index = moco_utils.index_from_filename(current_file)
                        stitched_functional_two[:, :, slice, index] = np.load(current_file)
                        # Just a sanity check! E.g. for first image we expect '0'
                        if index_tracker != index:
                            print('There seems to be a problem with the temp files for: ' + str(current_file))
                            print('Previous index was ' + str(index_tracker - 1))
                            print('Current slice is {}'.format(slice))
                            print('Next index (based on filename) was ' + str(index))
                        index_tracker += 1
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
    np.save(param_output_path, transform_matrix)

    print('Took: ' + repr(time.time() - time_start) + 's to combine files')

    t0 = time.time()
    ### MAKE MOCO PLOT ###
    moco_utils.save_moco_figure(
        transform_matrix=transform_matrix,
        parent_path=moving_path.parent,
        moco_dir=moving_output_path.parent,
        printlog=printlog,
    )

    print('took: ' + repr(round(time.time() - t0,2)) + ' s to plot moco')


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
        cores = 4

    # create an index going from [0,1,...,n]
    time_index = moco_utils.prepare_time_index(moving_path)
    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    nslices = brain_shape[2] #assumes data is x,y,z,t

    print('Will perform motion correction on a total of ' + repr(len(time_index)) + ' volumes.')
    if TESTING:
        time_index = [0,1,2,3,4,5,6,7]

    # Manual multiprocessing, essentially copy-paste from answer here:
    # https://stackoverflow.com/questions/23119382/how-can-i-multithread-a-function-that-reads-a-list-of-objects-in-python-astroph/23436094#23436094
    pool = multiprocessing.Pool(cores)
    # Keep track of number of processes
    running_processes = 0
    counter = 0
    # Make a list to keep track of the spawned child processes
    child_processes = []
    # Define what the max number of processes is
    max_processes = cores

    # Loop thorugh slices
    for current_slice in range(nslices):
    # Loop through index, yield 0, 1 etc.
        for current_index in time_index:
            # Run until break
            while True:
                # Only fork a new process is there are less processes running than max_processes
                if running_processes < max_processes:
                    # Define process: Target is a function, args are the arguments
                    p = multiprocessing.Process(target=moco_slice,
                                                args=(current_index,
                                                    current_slice,
                                                    fixed_path,
                                                    moving_path,
                                                    functional_channel_paths,
                                                    temp_save_path,
                                                    fly_directory))
                    # start process
                    p.start()
                    # To keep track of child processes, add to list
                    child_processes.append(p)
                    # to keep track of running_processes
                    running_processes += 1
                    counter += 1
                    # get out of the current 'while' loop and go back to the for loop
                    break
                # Processor wait loop if we don't have running_processes < max_processes
                else:
                    # Stay here until break is called
                    while True:
                        # loop through the child_processes
                        for current_child_process in range(len(child_processes)):
                            # Check if process is still running
                            if child_processes[current_child_process].is_alive():
                                # Continue for loop (i.e. check next child_process)
                                continue
                            else:
                                # If it's found that a child process isn't running anymore,
                                # remove the item at the current index
                                child_processes.pop(current_child_process)
                                # Subtract running processes by one
                                running_processes -= 1
                                # and break the for loop
                                break
                        # We are here either because the for loop finished or becuse
                        # it was found that a process is not running anymore.
                        # Check if we have less running processes than max processes
                        if running_processes < max_processes:
                            # If yes, break this inner while loop and go back to the
                            # outer while loop that keeps to start a new child process
                            break
                            # Else stay in this while loop and check again for processes
                            # that are finished.
    print('Submitted all indeces. Waiting for remaining processes to complete')
    # wait for remaining processes to complete --> this is the same code as the
    # processor wait loop above
    while len(child_processes) > 0:
        for next in range(len(child_processes)):
            if child_processes[next].is_alive():
                continue
            else:
                child_processes.pop(next)
                running_processes -= 1
                break  # need to break out of the for-loop,
                # because the child_list index is changed by pop

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