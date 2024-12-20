import numpy as np
import pathlib
import sys
import time
import traceback
import matplotlib as mpl
print('initial backend for preprocessing is: ' + repr(mpl.get_backend()))
mpl.use("agg")
print('setting backend for preprocessing to: ' + repr(mpl.get_backend()))
 # Agg, is a non-interactive backend that can only write to files.
# Without this I had the following error: Starting a Matplotlib GUI outside of the main thread will likely fail.
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
    print('test printing before log definition')
    logfile = utils.create_logfile(fly_directory, function_name="bleaching_qc")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, "bleaching_qc")
    print('test printing after log definition')
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
    #mpl.pyplot.rcParams.update({"font.size": 24})
    print('reseting backend for bleaching qc, current backend is: ' + repr(mpl.get_backend()))
    mpl.pyplot.switch_backend('agg')
    print('backend for bleaching qc set to: ' + repr(mpl.get_backend()))
    fig = mpl.pyplot.figure(figsize=(10, 10))
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
