import numpy as np
import pathlib
import sys
import time
import traceback
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("agg")  # Agg, is a non-interactive backend that can only write to files.
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

def moco_slice(
    fly_directory,
    path_to_read,
    brain_master,
    brain_mirror = 'none'
):
    """
    :param fly_directory:
    :param path_to_read:
    """
    #####################
    ### SETUP LOGGING ###
    #####################

    logfile = utils.create_logfile(fly_directory, function_name="moco_slice")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    utils.print_function_start(logfile, "moco_slice")

    ##########
    ### Convert list of (sometimes empty) strings to pathlib.Path objects
    ##########
    path_to_read = utils.convert_list_of_string_to_posix_path(path_to_read)
    #save_path_cluster_labels = utils.convert_list_of_string_to_posix_path(
    #    save_path_cluster_labels
    #)
    #save_path_cluster_signals = utils.convert_list_of_string_to_posix_path(
    #    save_path_cluster_signals
    #)

    # Can have more than one functional channel, hence loop!
    for (
        current_path_to_read,
        current_path_to_save_labels,
        current_path_to_save_signals,
    ) in zip(path_to_read, save_path_cluster_labels, save_path_cluster_signals):
        printlog("Loading " + repr(current_path_to_read))
        ### LOAD BRAIN ###

        # brain_path = os.path.join(func_path, 'functional_channel_2_moco_zscore_highpass.h5')
        t0 = time.time()
        # with h5py.File(brain_path, 'r+') as h5_file:
        #with h5py.File(current_path_to_read, "r+") as file:
        #    # Load everything into memory, cast as float 32
        #    brain = file["data"][:].astype(np.float32)
        #    # Convert nan to num, ideally in place to avoid duplication of data
        #    brain = np.nan_to_num(brain, copy=False)
        #    # brain = np.nan_to_num(h5_file.get("data")[:].astype('float32'))

        # Everything is only nifty in this pipeline! Define proxy
        brain_master_data_proxy = nib.load(current_path_to_read)
        # Load everything into memory, cast DTYPE
        brain_master = np.asarray(brain_master_data_proxy.dataobj, dtype=DTYPE)
        # Convert nan to num, ideally in place to avoid duplication of data
        brain_master = np.nan_to_num(brain_master, copy=False)

        printlog("brain shape: {}".format(brain_master.shape))

        # load mirror brain (if it exists)
        if brain_mirror not 'none':
            # Everything is only nifty in this pipeline! Define proxy
            brain_mirror_data_proxy = nib.load(current_path_to_read)
            # Load everything into memory, cast DTYPE
            brain_mirror = np.asarray(brain_mirror_data_proxy.dataobj, dtype=DTYPE)
            # Convert nan to num, ideally in place to avoid duplication of data
            brain_mirror = np.nan_to_num(brain_mirror, copy=False)
            printlog("brain mirror shape: {}".format(brain_mirror.shape))


        printlog("load duration: {} sec".format(time.time() - t0))

        ind_slice = brain_master.shape[2]

        ### MAKE MOCO DIRECTORY ###

        current_path_to_save_labels.parent.mkdir(exist_ok=True, parents=True)

        nifti1_limit = (2**16 / 2)

        ### RUN MOCO ###

        for slice in range(ind_slice):

            printlog("Running Moco on slice {}".format(slice))
            
            t0 = time.time()

            brain_master_slice = brain_master[:,:,slice,:]

            if np.any(np.array(brain_slice.shape) >= nifti1_limit):  # Need to save as nifti2
                nib.save(nib.Nifti2Image(slice.astype('float32'), np.eye(4)), os.path.join(dir,savefile_master))
            else:  # Nifti1 is OK
                nib.save(nib.Nifti1Image(slice.astype('float32'), np.eye(4)), os.path.join(dir,savefile_master))
            
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

def background_subtract(
    input
):
    """
    :param input:

    """
    #####################
    ### SETUP LOGGING ###
    #####################

    logfile = utils.create_logfile(fly_directory, function_name="background_subtract")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    utils.print_function_start(logfile, "background_subtract")

    ##########
    ### Convert list of (sometimes empty) strings to pathlib.Path objects
    ##########
    path_to_read = utils.convert_list_of_string_to_posix_path(path_to_read)