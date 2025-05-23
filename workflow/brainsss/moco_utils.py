import h5py
import numpy as np
import matplotlib as mpl
mpl.use("agg")
# Agg, is a non-interactive backend that can only write to files.
# Without this I had the following error: Starting a Matplotlib GUI outside of the main thread will likely fail.
import pathlib
import time
import nibabel as nib
import math

# To import brainsss, define path to scripts!
import sys

scripts_path = pathlib.Path(
    __file__
).parent.resolve()  # path of workflow i.e. /Users/dtadres/snake_brainsss/workflow
sys.path.insert(0, pathlib.Path(scripts_path, "workflow"))

from brainsss import utils

def prepare_time_index(moving_path):
    """
    Returns a list with length (time dimension) of dataset. List starts at 0 and
    goes to data.shape[-1] like this: [0,1,..n]
    :param moving_path: path to data used to extract time
    :return:
    """
    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    # last dimension is time, indicating the amount of volumes in the dataset
    experiment_total_frames = brain_shape[-1]
    # make a list that represents the index of the total_frames
    time_index = list(np.arange(experiment_total_frames))

    return(time_index)

def prepare_time_index_chunks(moving_path,chunk):
    """
    Returns a list with length (time dimension) of dataset. List starts at 0 and
    goes to data.shape[-1] like this: [0,chunk-1,2*chunk-1,...n].
    It includes the last element even if the chunk size changes (if total frames
    divided by chunk size is not an integer)
    :param moving_path: path to data used to extract time
    :return: time_index (list of tuples of length eperiment_total_frames, where each tuple is pair of start and end inds for that chunk)
    """
    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    # get total number of frames
    experiment_total_frames = brain_shape[3] #assumes data is x,y,z,t
    # make a list that represents the index of the total_frames
    n_chunks = math.ceil(experiment_total_frames/chunk)
    start_ind = list(range(0,experiment_total_frames-1,chunk))
    end_ind = list(range(0+chunk-1,experiment_total_frames,chunk))
    if len(start_ind)!=len(end_ind): #experiment_total_frames is not a multiple of chunk (missing last end_ind)
        end_ind.append(experiment_total_frames-1)
    time_index = []
    for i in range(0,n_chunks):
        current_ind_pair = (start_ind[i],end_ind[i])
        time_index.append(current_ind_pair)

    return(time_index)

def index_from_filename(filename):
    """
    Temporary files are saved with a characteristic file name indicating their index
    in the original array such as: 'channel_1.niiindex_0.npy'
    This function just extracts the index (in this case 0) from the filename.
    :param filename:
    :return:
    """
    index = int(filename.name.split('index')[-1].split('.npy')[0])

    return(index)

def index_range_from_filename(filename):
    """
    Temporary files are saved with a characteristic file name indicating their index
    in the original array such as: '...index0-99.npy'
    This function just extracts the index range (in this case (0,99)) from the filename.
    :param filename:
    :return: index range for file as tuple (start_index,end_index)
    """
    start_index = int(filename.name.split('index')[-1].split('.npy')[0].split('-')[0])
    end_index = int(filename.name.split('index')[-1].split('.npy')[0].split('-')[1])
    index_range = (start_index,end_index)

    return(index_range)


'''
def make_empty_h5(directory, file, brain_dims): 
    moco_dir = pathlib.Path(directory, "moco")
    moco_dir.mkdir(exist_ok=True, parents=True)

    # savefile = os.path.join(moco_dir, file)
    savefile = pathlib.Path(moco_dir, file)
    with h5py.File(savefile, "w") as f:
        dset = f.create_dataset("data", brain_dims, dtype="float32", chunks=True)
    return moco_dir, savefile'''


def save_moco_figure(transform_matrix, parent_path, moco_dir, printlog):
    """

    :param transform_matrix:
    :param parent_path:
    :param moco_dir:
    :param printlog:
    :return:
    """

    xml_path = None
    # Get voxel resolution for figure
    for current_file in parent_path.iterdir():
        if "recording_metadata.xml" in current_file.name:
            xml_path = current_file

            # if xml_path == None:
            # 	printlog('Could not find xml file for scan dimensions. Skipping plot.')
            # 	return
            # elif not xml_path.is_file():
            # 	printlog('Could not find xml file for scan dimensions. Skipping plot.')
            # 	return

            printlog(f"Found xml file.")
            x_res, y_res, z_res = utils.get_resolution(xml_path)

            # Save figure of motion over time
            # save_file = os.path.join(moco_dir, 'motion_correction.png')
            save_file = pathlib.Path(moco_dir, "motion_correction.png")
            fig = mpl.pyplot.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)
            ax.plot(
                transform_matrix[:, 9] * x_res, label="y"
            )  # note, resolutions are switched since axes are switched
            ax.plot(transform_matrix[:, 10] * y_res, label="x")
            ax.plot(transform_matrix[:, 11] * z_res, label="z")
            ax.set_ylabel("Motion Correction, um")
            ax.set_xlabel("Time")
            ax.set_title(moco_dir)
            mpl.pyplot.legend()
            fig.savefig(save_file, bbox_inches="tight", dpi=300)

            return
    printlog("Could not find xml file for scan dimensions. Skipping plot.")

def save_moco_figure_stackreg_rigid(transform_matrix, metadata_dir, moco_dir, printlog, n_slices):
    """

    :param transform_matrix:
    :param parent_path:
    :param moco_dir:
    :param printlog:
    :return:
    """

    xml_path = None
    # Get voxel resolution for figure
    for current_file in pathlib.Path(metadata_dir).iterdir():
        if "recording_metadata.xml" in current_file.name:
            xml_path = current_file

            # if xml_path == None:
            # 	printlog('Could not find xml file for scan dimensions. Skipping plot.')
            # 	return
            # elif not xml_path.is_file():
            # 	printlog('Could not find xml file for scan dimensions. Skipping plot.')
            # 	return

            # transform_matrix dimensions: 3,3,slice,timepoint
            printlog(f"Found xml file.")
            x_res, y_res, z_res = utils.get_resolution(xml_path)

            for slice in range(n_slices):
                # Save figure of motion over time for each slice
                # save_file = os.path.join(moco_dir, 'motion_correction.png')
                save_file_translation = pathlib.Path(moco_dir, "motion_correction_slice_{}_translation.png".format(slice))
                fig1 = mpl.pyplot.figure(figsize=(10, 10))
                ax1 = fig1.add_subplot(111)
                ax1.plot(transform_matrix[0][2][slice][:] * x_res, label="x")
                ax1.plot(transform_matrix[1][2][slice][:] * y_res, label="y")
                ax1.set_ylabel("Motion Correction, um")
                ax1.set_xlabel("Frame")
                ax1.set_title(moco_dir + "moco_slice_{}_translation".format(slice))
                mpl.pyplot.legend()
                fig1.savefig(save_file_translation, bbox_inches="tight", dpi=300)
                
                save_file_angle = pathlib.Path(moco_dir, "motion_correction_slice_{}_angle.png".format(slice))
                fig2 = mpl.pyplot.figure(figsize=(10, 10))
                ax2 = fig2.add_subplot(111)
                transform_angle = []
                for time_ind in range(transform_matrix.shape[-1]):
                    transform_angle.append(math.acos(transform_matrix[0][0][slice][time_ind]))
                ax2.plot(transform_angle , label="angle")
                ax2.set_ylabel("angle (rad)")
                ax2.set_xlabel("Frame")
                ax2.set_title(moco_dir + "moco_slice_{}_angle".format(slice))
                fig2.savefig(save_file_angle, bbox_inches="tight", dpi=300)

            return
    printlog("Could not find xml file for scan dimensions. Skipping plot.")


def print_progress_table_moco(total_vol, complete_vol, printlog, start_time, width):
    """
    There's a very similarly named function in utils!
    :param total_vol:
    :param complete_vol:
    :param printlog:
    :param start_time:
    :param width:
    :return:
    """
    fraction_complete = complete_vol / total_vol

    ### Get elapsed time ###
    elapsed = time.time() - start_time
    elapsed_hms = sec_to_hms(elapsed)

    ### Get estimate of remaining time ###
    try:
        remaining = elapsed / fraction_complete - elapsed
    except ZeroDivisionError:
        remaining = 0
    remaining_hms = sec_to_hms(remaining)

    ### Get progress bar ###
    complete_vol_str = f"{complete_vol:04d}"
    total_vol_str = f"{total_vol:04d}"
    length = (
        len(elapsed_hms)
        + len(remaining_hms)
        + len(complete_vol_str)
        + len(total_vol_str)
    )
    bar_string = utils.progress_bar(complete_vol, total_vol, width - length - 10)

    full_line = (
        "| "
        + elapsed_hms
        + "/"
        + remaining_hms
        + " | "
        + complete_vol_str
        + "/"
        + total_vol_str
        + " |"
        + bar_string
        + "|"
    )
    printlog(full_line)


def sec_to_hms(t):
    secs = f"{np.floor(t%60):02.0f}"
    mins = f"{np.floor((t/60)%60):02.0f}"
    hrs = f"{np.floor((t/3600)%60):02.0f}"
    return ":".join([hrs, mins, secs])
