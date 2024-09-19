import sys
import smtplib
import re
import os
import h5py
import math
import json
import datetime
import pyfiglet
from time import time
from time import strftime
from time import sleep
from functools import wraps
from email.mime.text import MIMEText
import numpy as np
import nibabel as nib
from xml.etree import ElementTree as ET
import subprocess
import traceback
import natsort
import pathlib
# For logfile
WIDTH = 120 # can go into a config file as well.

try:
    # only imports on linux
    import fcntl
    FCNTL_EXISTS = True
except ImportError:
    # import windows analog to being able to lock a file!
    FCNTL_EXISTS = False
    import threading

def convert_list_of_string_to_posix_path(list_of_strings):
    """
    Unfortunately when input is passed by snakemake, the PosixPaths are converted to strings.
    Convert back with this function
    :param list_of_strings:
    :return:
    """
    list_with_posix_paths = []
    for current_path in list_of_strings:
        print(current_path)
        try:
            list_with_posix_paths.append(pathlib.Path(current_path))
        except TypeError:
            # This happens because we sometimes pass empty lists for the unused channels.
            pass
    return list_with_posix_paths


def get_new_fly_number(
    target_path, first_fly_with_genotype_this_run, already_created_folders
):
    """
    Function to identify new fly number
    :param target_path: a string pointing to a path, i.e. /oak/stanford/groups/trc/data/David/Bruker/preprocessed
    :return: three digit number
    """
    # Check if target path (the gentype folder) already exists. If not, by definition
    # we need to create 'fly_001' (see below, else statement)
    if target_path.is_dir():
        # loop through target path and collect all files and folders that contain 'fly'
        fly_folders = [
            s for s in (target_path).iterdir() if "fly" in s.name and s.is_dir()
        ]
        # fly folders should then be sorted like this: ['fly_999', 'fly_998',.., 'fly_001']
        sorted_fly_folder = natsort.natsorted(fly_folders, reverse=True)
        oldest_fly = None
        # If we have not yet created any folders with this genotype in this run...
        if first_fly_with_genotype_this_run:
            # ...check if the 'incomplete' flag is present in the most recent
            # 'fly_XXX' folder
            for current_fly_folder in sorted_fly_folder:
                print("current_fly_folder " + repr(current_fly_folder))
                if pathlib.Path(current_fly_folder, "incomplete").exists():
                    # If incomplete exists, continue 'down' the folders. i.e. if
                    # the incomplete file exists in 'fly_999', go to 'fly_998' and
                    # check there etc. until a folder is found where the file does not
                    # exists!

                    # TODOCheck if other files are in the folder and warn user
                    pass
                else:
                    # if no 'incomplete' file exists, we have found the oldest fly folder
                    oldest_fly = current_fly_folder.name.split("_")[-1]
                    break
            # if even 'fly_001' has a 'incomplete' file, the 'new_fly_number' is '001'.
            if oldest_fly is None:
                new_fly_number = str(1).zfill(3)  # fly number 001
                return new_fly_number  # stop function here.
        # if we have created other paths in the same genotype file during this run before
        # use this information to higher the next higher number
        else:
            # First sort all already_created_folder and find the ones with the correct
            # genotype
            relevant_paths = []
            for current_path in already_created_folders:
                if target_path.name in current_path.parts:
                    relevant_paths.append(current_path)
            # Next, sort the relevant paths
            sorted_relevant_paths = natsort.natsorted(relevant_paths, reverse=True)
            # and pick out the oldest fly
            oldest_fly = sorted_relevant_paths[0].name.split("_")[-1]
        new_fly_number = str(int(oldest_fly) + 1).zfill(3)
    else:
        new_fly_number = str(1).zfill(3)  # fly number 001

    return new_fly_number


def write_error(logfile, error_stack, width=WIDTH):
    with open(logfile, "a+") as file:
        file.write(f"\n{'     ERROR     ':!^{width}}")
        file.write("\nTraceback (most recent call last): " + str(error_stack) + "\n\n")
        file.write("Full traceback below: \n\n")
        file.write(traceback.format_exc())
        file.write(f"\n{'     ERROR     ':!^{width}}")


def parse_true_false(true_false_string):
    if true_false_string in ["True", "true"]:
        return True
    elif true_false_string in ["False", "false"]:
        return False
    else:
        return False


def load_user_settings(user):
    current_path = pathlib.Path(__file__).parents[1].resolve()
    user_file = pathlib.Path(current_path, "users", user + ".json")
    # print(user_file)
    # user_file = os.path.join(os.path.dirname(scripts_path), 'users', user + '.json')
    with open(user_file) as file:
        settings = json.load(file)
    return settings


def get_json_data(file_path):
    with open(file_path) as f:
        data = json.load(f)
    return data


def create_logfile(fly_folder, function_name):
    pathlib.Path(fly_folder, "logs").mkdir(exist_ok=True, parents=True)
    # Have one log file per rule placed into the destination folder!
    # This will make everything super traceable!
    logfile = (
        str(fly_folder) + "/logs/" + function_name + "_" + fly_folder.name + ".txt"
    )
    # Not sure what this does exactly, from Bella's code
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    # Pipe all errors to the logfile
    sys.stderr = LoggerRedirect(logfile)
    # Pipe all print statements (and other console output) to the logfile
    sys.stdout = LoggerRedirect(logfile)
    # Problem: Snakemake runs twice. Seems to be a bug: https://github.com/snakemake/snakemake/issues/2350
    # Only print title and fly if logfile doesn't yet exist
    #width = 120  # can go into a config file as well.
    print_function_start(logfile, function_name)
    return logfile


class LoggerRedirect(object):
    """
    for redirecting stderr to a central log file.
    note, locking did not work fir fcntl, but seems to work fine without it
    keep in mind it could be possible to get errors from not locking this
    but I haven't seen anything, and apparently linux is "atomic" or some shit...
    # Renamed by David to 'LoggerRedirect' from 'Logger_stderr_sherlock' and
    # is now used for both error messages (stdeer) and for print outputs (stout).
    """

    def __init__(self, logfile):
        self.logfile = logfile

    def write(self, message):
        with open(self.logfile, "a+") as f:
            f.write(message)

    def flush(self):
        pass


'''class LoggerStdout(object):
    """
    for redirecting stdout to a central log file.
    """
    def __init__(self, logfile):
        self.logfile = logfile

    def write(self, message):
        with open(self.logfile, 'a+') as f:
            f.write(message)

    def flush(self):
        pass'''


class Printlog:
    """
    for printing all processes into same log file on sherlock
    """

    def __init__(self, logfile):
        self.logfile = logfile

        if not FCNTL_EXISTS:
            self.lock = threading.Lock()

    def write_to_win_file(self, message):
        """
        Small re-write of the commented parts below for readability (and
        less def function and class calls!
        To be tested!
        """
        with self.lock:
            with open(self.logfile, "a+") as file:
                file.write(message)
                file.write('\n')

    def print_to_log(self, message):
        # Linux/Mac
        if FCNTL_EXISTS:

            with open(self.logfile, "a+") as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                f.write(message)
                f.write("\n")
                fcntl.flock(f, fcntl.LOCK_UN)
        # Windows
        else:
            # The commented stuff works!
            #lock = threading.Lock()
            #def write_to_file():
            #    with lock:
            #        with open(self.logfile, "a+") as file:
            #            file.write(message)
            #            file.write("\n")

            #write_to_file()
            #print(message)
            #print('\n')
            self.write_to_win_file(message)

def sbatch(
    jobname,
    script,
    modules,
    args,
    logfile,
    time=1,
    mem=1,
    dep="",
    nice=False,
    silence_print=False,
    nodes=2,
    begin="now",
    global_resources=False,
):
    if dep != "":
        dep = "--dependency=afterok:{} --kill-on-invalid-dep=yes ".format(dep)

    command = f"ml {modules}; python3 {script} {json.dumps(json.dumps(args))}"

    if nice:  # For lowering the priority of the job
        nice = 1000000

    if nodes == 1:
        node_cmd = "-w sh02-07n34 "
    else:
        node_cmd = ""

    #width = 120
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    script_name = os.path.basename(os.path.normpath(script))
    print_big_header(logfile, script_name, WIDTH)

    if global_resources:
        sbatch_command = "sbatch -J {} -o ./com/%j.out -e {} -t {}:00:00 --nice={} {}--open-mode=append --cpus-per-task={} --begin={} --wrap='{}' {}".format(
            jobname, logfile, time, nice, node_cmd, mem, begin, command, dep
        )
    else:
        sbatch_command = "sbatch -J {} -o ./com/%j.out -e {} -t {}:00:00 --nice={} --partition=trc {}--open-mode=append --cpus-per-task={} --begin={} --wrap='{}' {}".format(
            jobname, logfile, time, nice, node_cmd, mem, begin, command, dep
        )
    sbatch_response = subprocess.getoutput(sbatch_command)

    if not silence_print:
        printlog(f"{sbatch_response}{jobname:.>{WIDTH-28}}")
    job_id = sbatch_response.split(" ")[-1].strip()
    return job_id


def get_job_status(job_id, logfile, should_print=False):
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    temp = subprocess.getoutput(
        "sacct -n -P -j {} --noconvert --format=State,Elapsed,MaxRSS,NCPUS,JobName".format(
            job_id
        )
    )
    if temp == "":
        return None  # is empty if the job is too new
    status = temp.split("\n")[0].split("|")[0]

    if should_print:
        if status != "PENDING":
            try:
                duration = temp.split("\n")[0].split("|")[1]
                jobname = temp.split("\n")[0].split("|")[4]
                num_cores = int(temp.split("\n")[0].split("|")[3])
                memory_used = float(temp.split("\n")[1].split("|")[2])  # in bytes
            except (IndexError, ValueError) as e:
                printlog(str(e))
                printlog(f"Failed to parse sacct subprocess: {temp}")
                return status
            core_memory = 7.77 * 1024 * 1024 * 1024  # GB to MB to KB to bytes

            if memory_used > 1024**3:
                memory_to_print = f"{memory_used/1024 ** 3 :0.1f}" + "GB"
            elif memory_used > 1024**2:
                memory_to_print = f"{memory_used/1024 ** 2 :0.1f}" + "MB"
            elif memory_used > 1024**1:
                memory_to_print = f"{memory_used/1024 ** 1 :0.1f}" + "KB"
            else:
                memory_to_print = f"{memory_used :0.1f}" + "B"

            percent_mem = memory_used / (core_memory * num_cores) * 100
            percent_mem = f"{percent_mem:0.1f}"

            #width = 120
            pretty = "+" + "-" * (WIDTH - 2) + "+"
            sep = " | "
            printlog(
                f"{pretty}\n"
                f"{'| SLURM | '+jobname+sep+job_id+sep+status+sep+duration+sep+str(num_cores)+' cores'+sep+memory_to_print+' (' + percent_mem + '%)':{WIDTH-1}}|\n"
                f"{pretty}"
            )
        else:
            printlog("Job {} Status: {}".format(job_id, status))

    return status


def wait_for_job(job_id, logfile, com_path):
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    # printlog(f'Waiting for job {job_id}')
    while True:
        status = get_job_status(job_id, logfile)
        if status in ["COMPLETED", "CANCELLED", "TIMEOUT", "FAILED", "OUT_OF_MEMORY"]:
            status = get_job_status(job_id, logfile, should_print=True)
            com_file = os.path.join(com_path, job_id + ".out")
            try:
                with open(com_file, "r") as f:
                    output = f.read()
            except:
                output = None
            return output
        else:
            sleep(5)


def print_progress_table(
    progress, logfile, start_time, print_header=False, print_footer=False
):
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")

    fly_print, expt_print, total_vol, complete_vol = [], [], [], []
    for funcanat in progress:
        fly_print.append(funcanat.split("/")[-2])
        expt_print.append(funcanat.split("/")[-1])
        total_vol.append(progress[funcanat]["total_vol"])
        complete_vol.append(progress[funcanat]["complete_vol"])
        # printlog("{}, {}".format(progress[funcanat]['total_vol'], progress[funcanat]['complete_vol']))

    total_vol_sum = np.sum([int(x) for x in total_vol])
    complete_vol_sum = np.sum([int(x) for x in complete_vol])
    # printlog("{}, {}".format(total_vol_sum, complete_vol_sum))
    fraction_complete = complete_vol_sum / total_vol_sum
    num_columns = len(fly_print)
    column_width = int((WIDTH - 20) / num_columns)
    if column_width < 9:
        column_width = 9

    if print_header:
        printlog(
            (" " * 9) + "+" + "+".join([f"{'':-^{column_width}}"] * num_columns) + "+"
        )
        printlog(
            (" " * 9)
            + "|"
            + "|".join([f"{fly:^{column_width}}" for fly in fly_print])
            + "|"
        )
        printlog(
            (" " * 9)
            + "|"
            + "|".join([f"{expt:^{column_width}}" for expt in expt_print])
            + "|"
        )
        printlog(
            (" " * 9)
            + "|"
            + "|".join([f"{str(vol)+' vols':^{column_width}}" for vol in total_vol])
            + "|"
        )
        printlog(
            "|ELAPSED "
            + "+"
            + "+".join([f"{'':-^{column_width}}"] * num_columns)
            + "+"
            + "REMAININ|"
        )

    # This exists further below, no need to re-defind
    # def sec_to_hms(t):
    #    secs=F"{np.floor(t%60):02.0f}"
    #    mins=F"{np.floor((t/60)%60):02.0f}"
    #    hrs=F"{np.floor((t/3600)%60):02.0f}"
    #    return ':'.join([hrs, mins, secs])

    elapsed = time() - start_time
    elapsed_hms = sec_to_hms(elapsed)
    try:
        remaining = elapsed / fraction_complete - elapsed
    except ZeroDivisionError:
        remaining = 0
    remaining_hms = sec_to_hms(remaining)

    single_bars = []
    for funcanat in progress:
        bar_string = progress_bar(
            progress[funcanat]["complete_vol"],
            progress[funcanat]["total_vol"],
            column_width,
        )
        single_bars.append(bar_string)
    fly_line = "|" + "|".join(single_bars) + "|"
    # fly_line = '|' + '|'.join([F"{bar_string:^{column_width}}"]*num_columns) + '|'
    fly_line = "|" + elapsed_hms + fly_line + remaining_hms + "|"
    printlog(fly_line)

    if print_footer:
        printlog(
            "|--------+"
            + "+".join([f"{'':-^{column_width}}"] * num_columns)
            + "+--------|"
        )
        # for funcanat in progress:
        #    printlog("{} {} {}".format(funcanat, progress[funcanat]['complete_vol'], progress[funcanat]['total_vol']))


def progress_bar(iteration, total, length, fill="█"):
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    fraction = f"{str(iteration):^4}" + "/" + f"{str(total):^4}"
    bar_string = f"{bar}"
    return bar_string


def moco_progress(progress_tracker, logfile, com_path):
    ##############################################################################
    ### Printing a progress bar every min until all moco_partial jobs complete ###
    ##############################################################################
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    print_header = True
    start_time = time()
    while True:
        stati = []
        ###############################
        ### Get progress and status ###
        ###############################
        ### For each expt_dir, for each moco_partial job_id, get progress from slurm.out files, and status ###
        for funcanat in progress_tracker:
            complete_vol = 0
            for job_id in progress_tracker[funcanat]["job_ids"]:
                # Read com file
                com_file = os.path.join(com_path, job_id + ".out")
                try:
                    with open(com_file, "r") as f:
                        output = f.read()
                        # complete_vol_partial = int(max(re.findall(r'\d+', output)))
                        complete_vol_partial = max(
                            [int(x) for x in re.findall(r"\d+", output)]
                        )
                except:
                    complete_vol_partial = 0
                complete_vol += complete_vol_partial
            progress_tracker[funcanat]["complete_vol"] = complete_vol
            stati.append(get_job_status(job_id, logfile))  # Track status

        ############################
        ### Print progress table ###
        ############################
        print_progress_table(progress_tracker, logfile, start_time, print_header)
        print_header = False

        ###############################################
        ### Return if all jobs are done, else sleep ###
        ###############################################
        finished = ["COMPLETED", "CANCELLED", "TIMEOUT", "FAILED", "OUT_OF_MEMORY"]
        if all([status in finished for status in stati]):
            print_progress_table(progress_tracker, logfile, start_time)
            print_progress_table(
                progress_tracker, logfile, start_time, print_footer=True
            )  # print final 100% complete line
            return
        else:
            sleep(int(60 * 5))


def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    return [tryint(c) for c in re.split("([0-9]+)", s)]


def sort_nicely(x):
    x.sort(key=alphanum_key)


def get_resolution(xml_file):
    """Gets the x,y,z resolution of a Bruker recording.

    Units of microns.

    Parameters
    ----------
    xml_file: string to full path of Bruker xml file.

    Returns
    -------
    x: float of x resolution
    y: float of y resolution
    z: float of z resolution"""

    tree = ET.parse(xml_file)
    root = tree.getroot()
    statevalues = root.findall("PVStateShard")[0].findall("PVStateValue")
    for statevalue in statevalues:
        key = statevalue.get("key")
        if key == "micronsPerPixel":
            indices = statevalue.findall("IndexedValue")
            for index in indices:
                axis = index.get("index")
                if axis == "XAxis":
                    x = float(index.get("value"))
                elif axis == "YAxis":
                    y = float(index.get("value"))
                elif axis == "ZAxis":
                    z = float(index.get("value"))
                else:
                    print("Error")
    return x, y, z


# def load_timestamps(directory, file='recording_metadata'): # file='functional.xml'):
def load_timestamps(path_to_metadata):
    """
    NEW:
    Parses a Bruker xml file to get the times of each frame
    :param path_to_metadata: a full filepath to a recording_metadata.xml file
    return: timestamps: [t,z] numpy array of times (in ms) of Bruker imaging frames.

    ORIGINAL:

    Parses a Bruker xml file to get the times of each frame, or loads h5py file if it exists.

    First tries to load from 'timestamps.h5' (h5py file). If this file doesn't exist
    it will load and parse the Bruker xml file, and save the h5py file for quick loading in the future.

    Parameters
    ----------
    directory: full directory that contains xml file (str).
    file: Defaults to 'functional.xml'

    Returns
    -------
    timestamps: [t,z] numpy array of times (in ms) of Bruker imaging frames.

    """
    # Not sure if h5py really needs to be made...We only need the timestamps once per call
    # and it's different for each recording...
    # if '.h5' in path_to_metadata:
    # try:
    #    print('Trying to load timestamp data from hdf5 file.')
    #    #with h5py.File(os.path.join(directory, 'timestamps.h5'), 'r') as hf:
    #    with h5py.File(path_to_metadata, 'r') as hf:
    #        timestamps = hf['timestamps'][:]

    # except:
    # else:
    # print('Failed. Extracting frame timestamps from bruker xml file.')
    # xml_file = os.path.join(directory, file)
    # xml_file = pathlib.Path(directory, file)

    tree = ET.parse(path_to_metadata)
    root = tree.getroot()
    timestamps = []
    # In case of an aborted scan, 'frames' would return the
    # wrong number of z-slices. Keep track of initial z-slices
    # with this variable
    initial_frames = None

    sequences = root.findall("Sequence")
    for sequence in sequences:
        if initial_frames is None:
            initial_frames = sequence.findall("Frame")
        frames = sequence.findall("Frame")
        for frame in frames:
            # filename = frame.findall('File')[0].get('filename')
            time = float(frame.get("relativeTime"))
            timestamps.append(time)
    timestamps = np.multiply(timestamps, 1000)

    if len(sequences) > 1:
        try:
            timestamps = np.reshape(timestamps, (len(sequences), len(frames)))
        except ValueError:
            print('This might be an aborted scan because the timestamps can be assigned')
            print('Attempting to create partial timestamps (slower)')
            temp = np.zeros((len(sequences), len(initial_frames)))
            temp.fill(np.nan)
            counter = 0
            for seq in range(len(sequences)):
                for fr in range(len(initial_frames)):
                    try:
                        temp[seq,fr] = timestamps[counter]
                        counter +=1
                    except IndexError:
                        break
            timestamps = temp
    else:
        timestamps = np.reshape(timestamps, (len(frames), len(sequences)))

    ### Save h5py file ###
    # with h5py.File(os.path.join(directory, 'timestamps.h5'), 'w') as hf:
    #    hf.create_dataset("timestamps", data=timestamps)

    print("Success loading timestamps.")
    return timestamps


def print_big_header(logfile, message, width=WIDTH):
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    message_and_space = "   " + message.upper() + "   "
    printlog("\n")
    printlog("=" * width)
    printlog(f"{message_and_space:=^{width}}")
    printlog("=" * width)
    print_datetime(logfile, width)


def print_title(logfile, width=WIDTH, fly_id=False):
    """
    Currently not being used
    """
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    title = pyfiglet.figlet_format("snake-brainsss", font="doom")
    title_shifted = ("\n").join([" " * 35 + line for line in title.split("\n")][:-2])
    printlog("\n")
    printlog(title_shifted)
    # print_datetime(logfile, width)


def print_datetime(logfile, width=WIDTH):
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    day_now = datetime.datetime.now().strftime("%B %d, %Y")
    time_now = datetime.datetime.now().strftime("%I:%M:%S %p")
    printlog(f"{day_now+' | '+time_now:^{width}}")


def print_footer(logfile, width=WIDTH):
    """
    Currently not being used
    """
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    sleep(3)  # to allow any final printing
    day_now = datetime.datetime.now().strftime("%B %d, %Y")
    time_now = datetime.datetime.now().strftime("%I:%M:%S %p")
    printlog("=" * width)
    printlog(f"{day_now+' | '+time_now:^{width}}")


def print_function_done(logfile, function_name, width=WIDTH, ):
    """
    Currently not in use
    """
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    day_now = datetime.datetime.now().strftime("%B %d, %Y")
    time_now = datetime.datetime.now().strftime("%I:%M:%S %p")
    day_time = str(day_now) + " | " + str(time_now)
    printlog(
        f"\n{'   ' + str(function_name) + ' finished at:  ' + str(day_time) + '   ':*^{width}}"
    )


def print_function_start(logfile, function_name):
    printlog = getattr(Printlog(logfile=logfile), "print_to_log")
    day_now = datetime.datetime.now().strftime("%B %d, %Y")
    time_now = datetime.datetime.now().strftime("%I:%M:%S %p")
    day_time = str(day_now) + " | " + str(time_now)
    printlog(
        f"\n{'   ' + str(function_name) + ' called at:  ' + str(day_time) + '   ':=^{WIDTH}}"
    )


def append_json(path, key, value):
    """
    Function that reads a json file as a dict, adds a single key and value and OVERWRITES
    the original json file.
    :param path:
    :param key:
    :param value:
    :return:
    """
    with open(path, "r") as openfile:
        dict = json.load(openfile)
    dict[key] = value
    with open(path, "w") as openfile:
        json.dump(dict, openfile, indent=4)


def sec_to_hms(t):
    secs = f"{np.floor(t%60):02.0f}"
    mins = f"{np.floor((t/60)%60):02.0f}"
    hrs = f"{np.floor((t/3600)%60):02.0f}"
    return ":".join([hrs, mins, secs])


def check_for_nan_and_inf_func(array):
    """
    Check if there are any nan or inf in the array that is being passed.

    :param array:
    :return:
    """
    try:
        if np.isnan(array).any():
            print("!!!!! WARNING - THERE ARE NAN IN THE ARRAY !!!!!")
            print(
                "The position(s) of np.nan is/are: " + repr(np.where(np.isnan(array)))
            )
    except:
        print("Could not check for nan because:\n\n")
        print(traceback.format_exc())
        print("\n")

    try:
        if np.isinf(array).any():
            print("!!!!! WARNING - THERE ARE INF IN THE ARRAY !!!!!")
            print("The position(s) of np.inf is/are " + repr(np.where(np.isnan(array))))
    except:
        print("Could not check for inf because:\n\n")
        print(traceback.format_exc())
        print("\n")
