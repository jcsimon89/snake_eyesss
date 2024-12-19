import numpy as np
import pathlib
import sys
import time
import traceback
import matplotlib as mpl
mpl.use("agg")  # Agg, is a non-interactive backend that can only write to files.
# Without this I had the following error: Starting a Matplotlib GUI outside of the main thread will likely fail.
import matplotlib.pyplot as plt
import nibabel as nib
import shutil
import ants
import h5py
from scipy import ndimage
import sklearn
import skimage.filters
import sklearn.feature_extraction
import sklearn.cluster
import os
import multiprocessing
import natsort
import itertools
import argparse

# To import files (or 'modules') from the visanalysis folder, define path to scripts!
# path of workflow i.e. /Users/jcsimon89/snake_visanalysis/workflow
#scripts_path = pathlib.Path(
#    __file__
#).parent.resolve()
#sys.path.insert(0, pathlib.Path(scripts_path, "workflow"))
parent_path = str(pathlib.Path(pathlib.Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, parent_path)
# This just imports '*.py' files from the folder 'brainsss'.

#from visanalysis import ____ 

from brainsss import utils
#from brainsss import moco_utils

####################
# GLOBAL VARIABLES #
####################
WIDTH = 120  # This is used in all logging files
# Bruker gives us data going from 0-8191 so we load it as uint16.
# However, as soon as we do the ants.registration step (e.g in
# function motion_correction() we get floats back.
# The original brainsss usually saved everything in float32 but some
# calculations were done with float64 (e.g. np.mean w/o defining dtype).
# To be more explicit, I define here two global variables.
# Dtype is what was mostly called when saving and loading data
DTYPE = np.float32
# Dtype_calculation is what I explicity call e.g. during np.mean
DTYPE_CACLULATIONS = np.float32

TESTING = False

def make_mean_brain(fly_directory,
                    meanbrain_n_frames,
                    path_to_read, path_to_save,
                    rule_name
):
    """
    Function to calculate meanbrain.
    This is based on Bella's meanbrain script.
    :param fly_directory: pathlib.Path object to a 'fly' folder such as '/oak/stanford/groups/trc/data/David/Bruker/preprocessed/fly_001'
    :param meanbrain_n_frames: First n frames to average over when computing mean/fixed brain | Default None (average over all frames).
    :param path_to_read: Full path as a list of pathlib.Path objects to the nii to be read
    :param path_to_save: Full path as a list of pathlib.Path objects to the nii to be saved
    :param rule_name: a string used to save the log file
    """

    ####
    # LOGGING
    ####
    logfile = utils.create_logfile(fly_directory, function_name=rule_name)
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, rule_name)


    #####
    # CONVERT PATHS TO PATHLIB.PATH OBJECTS
    #####
    path_to_read = utils.convert_list_of_string_to_posix_path(path_to_read)
    path_to_save = utils.convert_list_of_string_to_posix_path(path_to_save)

    for current_path_to_read, current_path_to_save in zip(path_to_read, path_to_save):
        brain_data = None  # make sure the array is empty before starting another iteration!
        ###
        # Read imaging file
        ###
        printlog("Currently looking at: " + repr(current_path_to_read.name))
        if current_path_to_read.suffix == ".nii":
            # Doesn't load anything, just points to a given location
            brain_proxy = nib.load(current_path_to_read)
            # Load data, it's np.uint16 at this point, no point changing it.
            brain_data = np.asarray(brain_proxy.dataobj, dtype=np.uint16)
        elif current_path_to_read.suffix == ".h5":
            # Original: Not great because moco brains are saved as float32
            #with h5py.File(current_path_to_read, "r") as hf:
            #    brain_data = np.asarray(hf["data"][:], dtype="uint16")
            with h5py.File(current_path_to_read, "r") as hf:
                brain_data = np.asarray(hf["data"][:], dtype=DTYPE)
        else:
            printlog("Current file has suffix " + current_path_to_read.suffix)
            printlog("Can currently only handle .nii and .h5 files!")
        # brain_data = np.asarray(nib.load(path_to_read).get_fdata(), dtype='uint16')
        # get_fdata() loads data into memory and sometimes doesn't release it.

        ###
        # CREATE MEANBRAIN
        ###
        if meanbrain_n_frames is not None:
            # average over first meanbrain_n_frames frames
            meanbrain = np.mean(brain_data[..., : int(meanbrain_n_frames)], axis=-1, dtype=DTYPE_CACLULATIONS)
        else:  # average over all frames
            meanbrain = np.mean(brain_data, axis=-1, dtype=DTYPE_CACLULATIONS)

        printlog("Datatype of meanbrain: " + repr(meanbrain.dtype))
        ###
        # SAVE MEANBRAIN
        ###
        aff = np.eye(4)
        meanbrain_nifty = nib.Nifti1Image(
            meanbrain, aff
        )
        meanbrain_nifty.to_filename(current_path_to_save)

        ###
        # LOG SUCCESS
        ###
        fly_print = pathlib.Path(fly_directory).name
        func_print = str(current_path_to_read).split("/imaging")[0].split("/")[-1]
        # func_print = current_path_to_read.name.split('/')[-2]
        printlog(
            f"meanbrn | COMPLETED | {fly_print} | {func_print} | {brain_data.shape} ===> {meanbrain.shape}"
        )
            
def bleaching_qc(
    fly_directory,
    path_to_read,
    path_to_save,
):
    """
    Perform bleaching quality control
    This is based on Bella's 'bleaching_qc.py' script
    Input are all nii files per folder (e.g. channel_1.nii and channel_2.nii) and output is single 'bleaching.png' file.
    Bleaching is defined as overall decrease of fluorescence in a given nii file.

    :param fly_directory: a pathlib.Path object to a 'fly' folder such as '/oak/stanford/groups/trc/data/David/Bruker/preprocessed/fly_001'
    :param path_to_read: list of paths to images to read, can be more than one
    :param path_to_save: list of paths to the 'bleaching.png' file
    :param functional_channel_list: list with channels marked as functional channels by experimenter
    :param anatomical_channel: the channel marked as the anatomy channel by the experimenter
    """
    ####
    # LOGGING
    ####
    logfile = utils.create_logfile(fly_directory, function_name="bleaching_qc")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, "bleaching_qc")

    #####
    # CONVERT PATHS TO PATHLIB.PATH OBJECTS
    #####
    path_to_read = utils.convert_list_of_string_to_posix_path(path_to_read)
    path_to_save = utils.convert_list_of_string_to_posix_path(path_to_save)
    # There is only one path to save! It comes as a list
    path_to_save = path_to_save[0]

    #####
    # CALCULATE MEAN FLUORESCENCE PER TIMEPOINT
    #####
    data_mean = {}
    # For each path in the list
    for current_path_to_read in path_to_read:
        printlog(f"Currently reading: {current_path_to_read.name:.>{WIDTH - 20}}")
        # Doesn't load anything to memory, just a pointer
        brain_proxy = nib.load(current_path_to_read)
        # Load data into memory, brain at this point is half of uint14, no point doing float
        brain = np.asarray(brain_proxy.dataobj, dtype=np.uint16)
        utils.check_for_nan_and_inf_func(brain)
        # calculate mean over time
        data_mean[current_path_to_read.name] = np.mean(brain, axis=(0, 1, 2))

    ##############################
    ### OUTPUT BLEACHING CURVE ###
    ##############################
    # plotting params
    plt.rcParams.update({"font.size": 24})
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    signal_loss = {}

    for filename in data_mean:
        xs = np.arange(len(data_mean[filename]))
        color = "k"
        # I slightly changed the colornames to make it obvious that data was analyzed
        # with a pipeline different from the normal one
        if "channel_1" in filename:
            color = "tomato"  # on our scope this is red
        elif "channel_2" in filename:
            color = "lime"  # on our scope this is green
        elif "channel_3" in filename:
            color = "darkred"  # on our scope this should be an IR channel
        ax.plot(data_mean[filename], color=color, label=filename)
        try:
            # Fit polynomial to mean fluorescence.
            linear_fit = np.polyfit(xs, data_mean[filename], 1)
            # and plot it
            ax.plot(np.poly1d(linear_fit)(xs), color="k", linewidth=3, linestyle="--")
            # take the linear fit to calculate how much signal is lost and report it as the title
            signal_loss[filename] = (
                linear_fit[0] * len(data_mean[filename]) / linear_fit[1] * -100
            )
        except:
            print("unable to perform fit because of error: ")
            print(traceback.format_exc())
            print(
                "\n Checking for nan and inf as possible cause - if no output, no nan and inf found"
            )
            print(utils.check_for_nan_and_inf_func(data_mean[filename]))

    ax.set_xlabel("Frame Num")
    ax.set_ylabel("Avg signal")
    try:
        loss_string = ""
        for filename in data_mean:
            loss_string = (
                loss_string
                + filename
                + " lost"
                + f"{int(signal_loss[filename])}"
                + "%\n"
            )
    except:  # This happens when unable to peform fit (see try..except above).
        pass
    ax.set_title(loss_string, ha="center", va="bottom")

    ###
    # SAVE PLOT
    ###
    save_file = pathlib.Path(path_to_save)
    fig.savefig(save_file, dpi=300, bbox_inches="tight")

    ###
    # LOG SUCCESS
    ###
    printlog(f"Prepared plot and saved as: {str(save_file):.>{WIDTH - 20}}")

def fictrac_qc(
        fly_directory,
        fictrac_file_path,
        fictrac_fps
):
    """
    Perform fictrac quality control.
    This is based on Bella's fictrac_qc.py  script.
    :param fly_directory: a pathlib.Path object to a 'fly' folder such as '/oak/stanford/groups/trc/data/David/Bruker/preprocessed/fly_001'
    :param fictrac_file_paths: a list of paths as pathlib.Path objects
    :param fictrac_fps: frames per second of the videocamera recording the fictrac data, an integer
    """
    ####
    # LOGGING
    ####
    logfile = utils.create_logfile(fly_directory, function_name="fictrac_qc")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, "fictrac_qc")

    #####
    # CONVERT PATHS TO PATHLIB.PATH OBJECTS
    #####
    fictrac_file_path = utils.convert_list_of_string_to_posix_path(fictrac_file_path)
    # Organize path - there is only one, but it comes as a list
    fictrac_file_path = fictrac_file_path[0]

    ####
    # QUALITY CONTROL
    ####
    #for current_file in fictrac_file_paths:
    printlog("Currently looking at: " + repr(fictrac_file_path))
    fictrac_raw = fictrac_utils.load_fictrac(fictrac_file_path)
    # This should yield something like 'fly_001/func0/fictrac
    full_id = ", ".join(fictrac_file_path.parts[-4:-2])

    resolution = 10  # desired resolution in ms # Comes from Bella!
    expt_len = fictrac_raw.shape[0] / fictrac_fps * 1000
    behaviors = ["dRotLabY", "dRotLabZ"]
    fictrac = {}
    for behavior in behaviors:
        if behavior == "dRotLabY":
            short = "Y"
        elif behavior == "dRotLabZ":
            short = "Z"
        fictrac[short] = fictrac_utils.smooth_and_interp_fictrac(
            fictrac_raw, fictrac_fps, resolution, expt_len, behavior
        )
    time_for_plotting = np.arange(0, expt_len, resolution) # comes in ms

    # Call these helper functions for plotting
    fictrac_utils.make_2d_hist(
        fictrac, fictrac_file_path, full_id,  fixed_crop=True
    )
    fictrac_utils.make_2d_hist(
        fictrac, fictrac_file_path, full_id, fixed_crop=False
    )
    fictrac_utils.make_velocity_trace(
        fictrac, fictrac_file_path, full_id, time_for_plotting,
    )

    ###
    # LOG SUCCESS
    ###
    printlog(f"Prepared fictrac QC plot and saved in: {str(fictrac_file_path.parent):.>{WIDTH - 20}}")

def moco_slice(
    index,
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
    fixed_data = fixed_proxy.dataobj
    # However, ants seems to require float32 (I think)
    fixed_ants = ants.from_numpy(np.asarray(fixed_data, dtype=np.float32))

    # Load moving proxy in this process
    moving_proxy = nib.load(moving_path)

    # Load data in a given process
    current_moving = moving_proxy.dataobj[:,:,:,index]
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
    np.save(pathlib.Path(temp_save_path, moving_path.name + 'index_'
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
        current_functional_one = functional_one_proxy.dataobj[:,:,:,index]
        moving_frame_one_ants = ants.from_numpy(np.asarray(current_functional_one, dtype=np.float32))
        # to motion correction for functional image
        moving_frame_one_ants = ants.apply_transforms(fixed_ants, moving_frame_one_ants, transformlist)
        # put moco functional image into preallocated array
        #moco_functional_one[:, :, :, counter] = moving_frame_one_ants.numpy()
        #print('apply transforms took ' + repr(time.time() - t0) + 's')
        np.save(pathlib.Path(temp_save_path, functional_path_one.name + 'index_'
                             + str(index)),
                moving_frame_one_ants.numpy())

        if functional_path_two is not None:
            current_functional_two = functional_two_proxy.dataobj[:,:,:, index]
            moving_frame_two_ants = ants.from_numpy(np.asarray(current_functional_two, dtype=np.float32))
            moco_functional_two = ants.apply_transforms(fixed_ants, moving_frame_two_ants, transformlist)
            #moco_functional_two[:,:,:, counter] = moco_functional_two.numpy()
            np.save(pathlib.Path(temp_save_path, functional_path_two.name + 'index_'
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
            param_savename = pathlib.Path(temp_save_path, "motcorr_params" + 'index_'
                                          + str(index))
            np.save(param_savename, temp.parameters) # that's the transform_matrix in brainsss

        # lets' delete all files created by ants - else we quickly create thousands of files!
        pathlib.Path(x).unlink()
    print('Motion correction for ' + moving_path.as_posix()
          + ' at index ' + str(index) + ' took : '
          + repr(round(time.time() - t_function_start, 1))
          + 's\n')
