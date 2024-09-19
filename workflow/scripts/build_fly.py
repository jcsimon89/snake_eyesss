from xml.etree import ElementTree as ET
from lxml import etree, objectify
import json
import sys
import pathlib
import traceback
import shutil
import pandas as pd

# To import files (or 'modules') from the brainsss folder, define path to scripts!
# path of workflow i.e. /Users/dtadres/snake_brainsss/workflow
scripts_path = pathlib.Path(
    __file__
).parent.resolve()
sys.path.insert(0, pathlib.Path(scripts_path, "workflow"))

from brainsss import utils

####################
# GLOBAL VARIABLES #
####################
WIDTH = 120  # This is used in all logging files

def fly_builder(autotransferred_stimpack,
                fictrac_folder_path,
                import_dirs,
                dataset_dirs):
    """
    Move folders from imports to fly dataset - need to restructure folders.
    This is based on Bella's 'fly_builder.py' script

    # Note: I removed the discrepancy between 'anat' and 'func' folders. All files
    are now just called 'channel_1.nii' or 'channel_2.nii' as it makes handling filenames
    much, much simpler, especially when going through snakemake.

    In addition, channel # is conserved and does reflect the channel # assigned by the
    Bruker software!

    :param logfile: logfile to be used for all errors (stderr) and console outputs (stdout)
    :param user: your SUnet ID as a string
    :param dirs_to_build: a list of folders to build. e.g. dir_to_build = ['20231301'] or  dir_to_build = ['20232101', '20231202']
    :param target_folder:
    :return:
    """
    try:
        # Loop through the list of dirs
        for current_import_dir, current_dataset_dir in zip(import_dirs, dataset_dirs):
            ###
            # Logging
            ###
            logfile = utils.create_logfile(
                current_dataset_dir, function_name="fly_builder"
            )
            printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
            #utils.print_function_start(logfile, "fly_builder")
            printlog(f"Building flies from: {str(current_import_dir):.>{WIDTH - 22}}")

            printlog(
                f"Building fly directory:{str(current_dataset_dir):.>{WIDTH - 22}}"
            )
            printlog(
                f"\n{'   Building ' + current_import_dir.name + ' as ' + str(current_dataset_dir.name) + '   ':-^{WIDTH}}"
            )

            ###
            # Dict to keep track of filepaths
            ###
            fly_dirs_dict = {}
            fly_dirs_dict["fly ID"] = current_dataset_dir.name

            ###
            # Copy the fly.json file
            ###
            fly_json_import_path = pathlib.Path(current_import_dir, 'fly.json')
            fly_json_fly_path = pathlib.Path(current_dataset_dir, 'fly.json')
            shutil.copyfile(fly_json_import_path, fly_json_fly_path)

            ###
            # Add date & time to fly.json
            ###
            try:
                add_date_to_fly(current_dataset_dir)
            except Exception as e:
                printlog(str(e))
                printlog(str(e))
                printlog(traceback.format_exc())

            ###
            # Copy fly data
            ####
            fly_dirs_dict = copy_fly(
                current_import_dir, current_dataset_dir, printlog,
                autotransferred_stimpack,
                fictrac_folder_path,
                fly_dirs_dict
            )
            ###
            # Done copying fly data!
            ###

            # Save json file with all relevant paths
            with open(pathlib.Path(current_dataset_dir, current_dataset_dir.name + "_dirs.json"), "w", ) as outfile:
                json.dump(fly_dirs_dict, outfile, indent=4)

            # If we are here it should mean that everything above has been copied as expected.
            # We can therefore delete the 'incomplete' file in this folder
            try:
                pathlib.Path(current_dataset_dir, "incomplete").unlink()
                printlog(
                    f'Deleted "incomplete" file in :{str(current_dataset_dir):.>{WIDTH}}'
                )
            except FileNotFoundError:
                printlog(
                    f'"Incomplete" file not found in! :{str(current_dataset_dir):.>{WIDTH}}'
                )
            # In case part of the import folder is being copied and the program crashes, it might
            # lead to accidentally copying all flies from import to new flies in data which might lead
            # to data duplication.
            # To avoid this, we write in the import folder a file called 'data_transfered_to.txt'
            # with the path of the target
            already_copied_path = pathlib.Path(
                current_import_dir, "data_transfered_to.txt"
            )
            with open(already_copied_path, "w") as outputfile:
                outputfile.write(str(current_dataset_dir))
            printlog(
                f"Wrote data_already_transfered_to.txt file in :{str(current_import_dir):.>{WIDTH}}"
            )

    except Exception as error_stack:
        printlog('!!! ERROR !!! -> check error file in "log" folder')
        logfile = utils.create_logfile(
            current_dataset_dir, function_name="ERROR_fly_builder"
        )
        utils.write_error(logfile=logfile, error_stack=error_stack, width=WIDTH)


def add_date_to_fly(fly_folder):
    """get date from xml file and add to fly.json"""

    # Check if there are func folders:
    candidate_folders = [
        pathlib.Path(fly_folder, x) for x in fly_folder.iterdir() if "func" in x.name
    ]
    # if not...
    if len(candidate_folders) == 0:
        # ...check if there are anat folders:
        candidate_folders = [
            pathlib.Path(fly_folder, x)
            for x in fly_folder.iterdir()
            if "anat" in x.name
        ]

    if len(candidate_folders) > 0:
        candidate_folder = candidate_folders[0]
        xml_file = pathlib.Path(candidate_folder, "imaging", "recording_metadata.xml")
        print("xml_file path" + repr(xml_file))
        # Extract datetime
        datetime_str, _, _ = get_datetime_from_xml(xml_file)
        # Get just date
        date = datetime_str.split("-")[0]
        time = datetime_str.split("-")[1]

        ### Add to fly.json
        json_file = pathlib.Path(fly_folder, "fly.json")
        # json_file = os.path.join(destination_fly, 'fly.json')
        with open(json_file, "r+") as f:
            metadata = json.load(f)
            metadata["Date"] = str(date)
            metadata["Time"] = str(time)
            f.seek(0)
            json.dump(metadata, f, indent=4)
            f.truncate()

    else:
        print(
            'Unable to find folder called "anat" or "func" to read "recording_metadata.xml'
        )


def copy_fly(import_dir,
             dataset_dir,
             printlog,
             autotransferred_stimpack,
             fictrac_folder_path,
             fly_dirs_dict,
             ):
    """
    There will be two types of folders in a fly folder.
    1) func_x folder
    2) anat_x folder
    For functional folders, need to copy fictrac and visual as well
    For anatomy folders, only copy folder. There will also be
    3) fly json data
    4) the recording metadata was called 'anatomy.xml' and 'functional.xml' in the past.
       this is changed to always be called 'recording_metadata.xml'
    """
    # look at every item in source fly folder
    # This is a folder like: '/oak/stanford/groups/trc/data/David/Bruker/imports/20231207/fly1'
    for current_import_file_or_folder in import_dir.iterdir():
        # This should be e.g. directory such as anat1 or func0 or fly.json
        print(
            "Currently looking at source_item: {}".format(
                current_import_file_or_folder.name
            )
        )

        # Handle folders
        if current_import_file_or_folder.is_dir():
            # Call this folder source expt folder
            current_import_imaging_folder = current_import_file_or_folder
            # Make the same folder in destination fly folder
            current_dataset_folder = pathlib.Path(
                dataset_dir, current_import_file_or_folder.name
            )
            # current_dataset_folder.mkdir(parents=True, exist_ok=True)

            # Is this folder an anatomy or functional folder?
            if "anat" in current_import_imaging_folder.name:
                # If anatomy folder, just copy everything
                # Make imaging folder and copy
                # imaging_destination = os.path.join(expt_folder, 'imaging')
                imaging_destination = pathlib.Path(current_dataset_folder, "imaging")
                # os.mkdir(imaging_destination)
                imaging_destination.mkdir(parents=True, exist_ok=True)
                copy_bruker_data(
                    current_import_imaging_folder, imaging_destination, "anat", printlog
                )
                current_fly_dir_dict = str(imaging_destination).split(
                    imaging_destination.parents[1].name
                )[-1]
                # json_key = current_imaging_folder.name + ' Imaging'
                # utils.append_json(path=fly_dirs_dict_path, key=json_key, value=current_fly_dir_dict)
                fly_dirs_dict[
                    current_import_imaging_folder.name + " Imaging"
                ] = pathlib.Path(current_fly_dir_dict).as_posix()

                ###
                # write anat info to csv file
                ###
                # Define date and especially time of this particular recording
                current_recording_metadata_path = pathlib.Path(imaging_destination, 'recording_metadata.xml')
                # Extract datetime
                datetime_str, _, _ = get_datetime_from_xml(current_recording_metadata_path)
                # Get just date
                current_date = datetime_str.split("-")[0]
                current_time = datetime_str.split("-")[1]
                # Finally write anat info to csv file
                add_fly_to_csv(import_dir, dataset_dir, current_import_imaging_folder, current_date, current_time, printlog)
                ######################################################################
                print(
                    f"anat:{current_dataset_folder}"
                )  # IMPORTANT - FOR COMMUNICATING WITH MAIN
                ######################################################################
            elif "func" in current_import_imaging_folder.name:
                # Make imaging folder and copy
                # imaging_destination = os.path.join(expt_folder, 'imaging')
                imaging_destination = pathlib.Path(current_dataset_folder, "imaging")
                # os.mkdir(imaging_destination)
                imaging_destination.mkdir(parents=True, exist_ok=True)
                copy_bruker_data(
                    current_import_imaging_folder, imaging_destination, "func", printlog
                )
                # Update fly_dirs_dict
                current_fly_dir_dict = str(imaging_destination).split(
                    imaging_destination.parents[1].name
                )[-1]
                # json_key = current_imaging_folder.name + ' Imaging'
                # utils.append_json(path=fly_dirs_dict_path, key=json_key, value=current_fly_dir_dict)
                fly_dirs_dict[
                    current_import_imaging_folder.name + " Imaging"
                ] = pathlib.Path(current_fly_dir_dict).as_posix()
                # Copy fictrac data
                # Automatic fictrac assignment (done on davidtadres/brukerbridge) where
                # each func folder would have a folder with stimpack/#/loco/*.dat
                print('autotransferred_stimpack: ' + repr(autotransferred_stimpack))
                if autotransferred_stimpack:
                    try:
                        # while data is PROBABLY TEST automatically transferred, the name
                        # of the *.dat file needs to be changed for downstream analysis!
                        automatic_copy_stimpack(import_folder=current_import_imaging_folder,
                                                target_folder=current_dataset_folder,
                                                fly_dirs_dict=fly_dirs_dict,
                                               )
                    except Exception as error:
                        printlog('***** ERROR *****')
                        printlog(str(error))
                        printlog('\n')

                # Manual fictrac assignment
                elif fictrac_folder_path is not None:
                    try:
                        fly_dirs_dict = manual_copy_fictrac(
                            current_dataset_folder,
                            printlog,
                            fictrac_folder_path,
                            current_import_imaging_folder,
                            fly_dirs_dict,
                        )
                        # printlog('Fictrac data copied')
                    except Exception as e:
                        printlog("Could not copy fictrac data because of error:")
                        printlog(str(e))
                        printlog(traceback.format_exc())
                # Copy visual data based on timestamps, and create visual.json
                try:
                    copy_visual(current_dataset_folder, printlog)
                except Exception as e:
                    printlog("Could not copy visual data because of error:")
                    printlog(str(e))
                ###
                # write func info to csv file
                ###
                # Define date and especially time of this particular recording
                current_recording_metadata_path = pathlib.Path(imaging_destination, 'recording_metadata.xml')
                # Extract datetime
                datetime_str, _, _ = get_datetime_from_xml(current_recording_metadata_path)
                # Get just date
                current_date = datetime_str.split("-")[0]
                current_time = datetime_str.split("-")[1]
                # Finally write anat info to csv file
                add_fly_to_csv(import_dir, dataset_dir, current_import_imaging_folder, current_date, current_time,
                               printlog)
                ######################################################################
                # print(f"func:{expt_folder}")  # IMPORTANT - FOR COMMUNICATING WITH MAIN
                ######################################################################
                # REMOVED TRIGGERING

            else:
                printlog(
                    "Invalid directory in fly folder (skipping): {}".format(
                        current_import_imaging_folder.name
                    )
                )

        # Copy fly.json file
        #else:
        #    current_import_file = current_import_file_or_folder
        #    if current_import_file_or_folder.name == "fly.json":
        #        target_path = pathlib.Path(dataset_dir, current_import_file.name)
        #        shutil.copyfile(current_import_file, target_path)
        else:
            # There shouldn't be any files in the top folder except fly.json.
            # If there are other files that need to be copied, make an elif above.
            printlog(
                "Invalid file in fly folder (skipping): {}".format(
                    current_import_file_or_folder.name
                )
            )

    return fly_dirs_dict


def copy_bruker_data(source, destination, folder_type, printlog, fly_dirs_dict=None):
    # Do not update destination - download all files into that destination
    # for item in os.listdir(source):
    for source_path in source.iterdir():
        # Check if item is a directory
        if source_path.is_dir():
            # Do not update destination - download all files into that destination
            copy_bruker_data(source_path, destination, folder_type, printlog)
            # In my case this leads to /oak/stanford/groups/trc/data/David/Bruker/imports/20231201/fly2/func1
            # The code then calls itself with this path which then goes one folder deeper

        # If the item is a file
        else:
            target_path = None
            # Don't copy these files
            if "SingleImage" in source_path.name:
                continue  # skip rest of the 'else' term
            # elif '.nii' in source_path.name and folder_type == 'func': # Shouldn't be necessary with '_s' check!
            #    continue  # do not copy!! Else we'll copy all the split nii files as well.
            #    # This is an artifact of splitting the nii file on Brukerbridge and might not be
            #    # relevant in the future/for other users!
            # each source path file can only be a single file - why if..if instead of if..elif?
            ### Change file names and filter various files
            # This is for split files from brukerbridge
            # Sometimes we have huge recordings. In that case, the nii files on Brukerbridge are split
            # To combine them, run the stitch snakefile. It'll write a file like 'ch1_concat.nii'
            elif "concat.nii" in source_path.name and folder_type == "func": # <this should be 'func'
                target_name = (
                    "channel_" + source_path.name.split("ch")[1].split("_")[0] + ".nii"
                )
                target_path = pathlib.Path(destination, target_name)
            # This is for non-split files from Brukerbridge - they MUST not have a '_s' string
            # which is indicative of a split file (see elif just above).
            elif (
                ".nii" in source_path.name
                and "_s" not in source_path.name
                and folder_type == "func" # <this should be 'func'
            ):
                target_name = (
                    "channel_" + source_path.name.split("channel")[1].split("_")[1]
                )  # this already gets us the '.nii'!
                print("target name" + repr(target_name))
                target_path = pathlib.Path(destination, target_name)
            # I don't completely understand whether anatomy files are never or only sometimes split.
            # To be safe, have an elif which first checks for concat
            elif ("concat.nii" in source_path.name
                  and folder_type == "anat"):  #
                target_name = (
                        "channel_" + source_path.name.split("ch")[1].split("_")[0] + ".nii"
                )
                target_path = pathlib.Path(destination, target_name)
            # Finally, similar to the func case above, if we are in split directory (indicated by .._s0.nii')
            # we ignore all nii files if they have the '_s' string.
            # Else, make that file the 'channel_x.nii' file.
            elif (".nii" in source_path.name
                  and "_s" not in source_path.name
                  and folder_type == "anat"):
                target_name = (
                    "channel_" + source_path.name.split("channel")[1].split("_")[1]
                )
                if ".nii" not in target_name:
                    target_name += "nii"  # Data I got from Yandan had double .nii!
                target_path = pathlib.Path(destination, target_name)
            # Special copy for photodiode since it goes in visual folder
            # To be tested once I have such data!!
            elif ".csv" in source_path.name:
                source_name = "photodiode.csv"
                visual_folder_path = pathlib.Path(destination, "visual")
                visual_folder_path.mkdir(exist_ok=True)
                target_path = pathlib.Path(visual_folder_path, source_name)
            # Special copy for visprotocol metadata since it goes in visual folder
            # To be tested once I have such data!!
            elif ".hdf5" in source_path.name:
                # Create folder 'visual'
                visual_folder_path = pathlib.Path(destination.name, "visual")
                visual_folder_path.mkdir(exist_ok=True)
                target_path = pathlib.Path(visual_folder_path, source_path.name)
            # Rename to recording_metadata.xml if appropriate
            elif ".xml" in source_path.name and "Voltage" not in source_path.name:
                target_name = "recording_metadata.xml"
                target_path = pathlib.Path(destination, target_name)
                copy_file_func(source_path, target_path, printlog)
                # Create json file
                create_imaging_json(target_path, printlog)
                #if folder_type == "func":
                continue

            elif ".xml" in source_path.name and "VoltageOutput" in source_path.name:
                target_path = pathlib.Path(destination, "voltage_output.xml")

            if target_path is not None:
                # Actually copy the file
                copy_file_func(source_path, target_path, printlog)


def copy_file_func(source, target, printlog):
    # printlog('Transfering file {}'.format(target))
    # to_print = ('/').join(target.split('/')[-4:])
    # print('source: ' + str(source))
    # print('target: ' + str(target))
    to_print = str(source.name + " to " + target.name)
    # width = 120
    printlog(f"Transfering file{to_print:.>{WIDTH - 16}}")
    ##sys.stdout.flush()
    shutil.copyfile(source, target)


def copy_visual(destination_region, printlog):
    printlog("copy_visual NOT IMPLEMENTED YET")
    """width = 120
    printlog(F"Copying visual stimulus data{'':.^{width - 28}}")
    visual_folder = '/oak/stanford/groups/trc/data/Brezovec/2P_Imaging/imports/visual'
    visual_destination = os.path.join(destination_region, 'visual')

    # Find time of experiment based on functional.xml
    true_ymd, true_total_seconds = get_expt_time(os.path.join(destination_region, 'imaging'))

    # Find visual folder that has the closest datetime
    # First find all folders with correct date, and about the correct time
    folders = []
    for folder in os.listdir(visual_folder):
        test_ymd = folder.split('-')[1]
        test_time = folder.split('-')[2]
        test_hour = test_time[0:2]
        test_minute = test_time[2:4]
        test_second = test_time[4:6]
        test_total_seconds = int(test_hour) * 60 * 60 + \
                             int(test_minute) * 60 + \
                             int(test_second)

        if test_ymd == true_ymd:
            time_difference = np.abs(true_total_seconds - test_total_seconds)
            if time_difference < 3 * 60:
                folders.append([folder, test_total_seconds])
                printlog('Found reasonable visual folder: {}'.format(folder))

    # if more than 1 folder, use the oldest folder
    if len(folders) == 1:
        correct_folder = folders[0]
    # if no matching folder,
    elif len(folders) == 0:
        printlog(F"{'No matching visual folders found; continuing without visual data':.<{width}}")
        return
    else:
        printlog('Found more than 1 visual stimulus folder within 3min of expt. Picking oldest.')
        correct_folder = folders[0]  # set default to first folder
        for folder in folders:
            # look at test_total_seconds entry. If larger, call this the correct folder.
            if folder[-1] > correct_folder[-1]:
                correct_folder = folder

    # now that we have the correct folder, copy it's contents
    printlog('Found correct visual stimulus folder: {}'.format(correct_folder[0]))
    try:
        os.mkdir(visual_destination)
    except:
        pass
        ##print('{} already exists'.format(visual_destination))
    source_folder = os.path.join(visual_folder, correct_folder[0])
    printlog('Copying from: {}'.format(source_folder))
    for file in os.listdir(source_folder):
        target_path = os.path.join(visual_destination, file)
        source_path = os.path.join(source_folder, file)
        ##print('Transfering from {} to {}'.format(source_path, target_path))
        ##sys.stdout.flush()
        shutil.copyfile(source_path, target_path)

    # Create visual.json metadata
    # Try block to prevent quiting if visual stimuli timing is wonky (likely went too long)
    try:
        unique_stimuli = brainsss.get_stimuli(visual_destination)
    except:
        unique_stimuli = 'brainsss.get_stimuli failed'
    with open(os.path.join(visual_destination, 'visual.json'), 'w') as f:
        json.dump(unique_stimuli, f, indent=4)"""

def automatic_copy_stimpack(import_folder, target_folder, fly_dirs_dict):
    """
    When using stimpack to trigger fictrac, davidtadres/brukerbridge has the option
    to automatically assign a given fictrac dataset to the corresponding imaging series.

    This function makes sure data is copied as expected

    :param import_folder: current_dataset_folder, pathlib object pointing to i.e. \\oak-smb-trc.stanford.edu\groups\trc\data\David\Bruker\imports\20240619\fly1\func0
    :param target_folder: pathlib object pointing to i.e. \\oak-smb-trc.stanford.edu\groups\trc\data\David\Bruker\preprocessed\FS144_x_FS61\fly_001\func0
    :param fly_dirs_dict: dict with relative path used by snakemake
    """

    # Check if fictrac data exists:
    fictrac_import_path = pathlib.Path(import_folder, 'stimpack/loco')
    fictrac_target_path = pathlib.Path(target_folder, 'stimpack/loco')
    # Go through import folder
    for current_file in fictrac_import_path.iterdir():
        # Rename this file because it needs to be consistent for snake_brainsss to work
        if '.dat' in current_file.name:
            current_target_path = pathlib.Path(fictrac_target_path, 'fictrac_behavior_data.dat')
        else:
            # Copy the rest of the content as well
            current_target_path = pathlib.Path(fictrac_target_path, current_file.name)
        # Make folder structure if not existing yet
        current_target_path.parent.mkdir(exist_ok=True, parents=True)
        # Copy file
        shutil.copyfile(current_file, current_target_path)

    # To keep track of where files are, ass to fly_dirs_dict
    relative_path = pathlib.Path('/' + target_folder.name, 'stimpack/loco/fictrac_behavior_data.dat')
    fly_dirs_dict[import_folder.name + " Fictrac "] = relative_path.as_posix() # Check if this correct!

    # Add more stuff if needed such as renaming visual data!

def manual_copy_fictrac(destination_region, printlog, fictrac_folder, source_fly, fly_dirs_dict):
    """
    Attempt to get correct fictrac data into the corresponding 'func' folder.
    This is called 'manual' because for each experiment (i.e. each TSeries) the
    user has to prepare a folder containing the .dat file for fictrac based on
    'fictrac_path' in user.json (i.e. david.json)
    For example for fly 20231201\fly2\func1 imaging data, fictrac data must
    be in the folder 20231201_fly2_func1. There must only be a single dat file in that folder!

    Consider using the stimpack autotransfer option.
    """
    # The target file will be called 'fictrac_behavior_data.dat' because it makes
    # handling files much easier in the snakefile.
    # Make fictrac folder
    fictrac_destination = pathlib.Path(destination_region, "fictrac")
    fictrac_destination.mkdir(exist_ok=True)
    # Fictrac source folder is defined in user.json
    # when doing post-hoc fictrac, Bella's code where one compare the recording
    # timestamps of imaging and fictrac doesn't work anymore.
    # I instead use a deterministic file structure:
    # for example for fly 20231201\fly2\func1 imaging data, fictrac data must
    # be in the folder 20231201_fly2_func1. There must only be a single dat file in that folder.
    source_path = pathlib.Path(
        fictrac_folder,
        source_fly.parts[-3]
        + "_"
        + source_fly.parts[-2]
        + "_"
        + source_fly.parts[-1],
    )
    for current_file in source_path.iterdir():
        if "dat" in current_file.name:
            # width = 120
            # source_path = os.path.join(source_path, file)
            dat_path = current_file
            # target_path = pathlib.Path(fictrac_destination, current_file.name)
            target_path = pathlib.Path(
                fictrac_destination, "fictrac_behavior_data.dat"
            )
            to_print = str(target_path)
            printlog(f"Transfering file{to_print:.>{WIDTH - 16}}")

            # put fictrac file path in into fly_dirs_dict
            current_fly_dir_dict = str(target_path).split(
                fictrac_destination.parents[1].name
            )[-1]
            fly_dirs_dict[destination_region.name + " Fictrac "] = pathlib.Path(current_fly_dir_dict).as_posix()
            shutil.copyfile(dat_path, target_path)
        else:
            # Copy rest of files such as videos and logs
            shutil.copyfile(current_file, pathlib.Path(fictrac_destination, current_file.name))


    """# OLD
    # Different users have different rule on what to do with the data
    if user == "dtadres":
        # TODO!!!
        fictrac_folder = pathlib.Path(
            "/oak/stanford/groups/trc/data/David/Bruker/Fictrac"
        )
        # when doing post-hoc fictrac, Bella's code where one compare the recording
        # timestamps of imaging and fictrac doesn't work anymore.
        # I instead use a deterministic file structure:
        # for example for fly 20231201\fly2\func1 imaging data, fictrac data must
        # be in the folder 20231201_fly2_func1. There must only be a single dat file in that folder.
        source_path = pathlib.Path(
            fictrac_folder,
            source_fly.parts[-3]
            + "_"
            + source_fly.parts[-2]
            + "_"
            + source_fly.parts[-1],
        )
        for current_file in source_path.iterdir():
            if "dat" in current_file.name:
                # width = 120
                # source_path = os.path.join(source_path, file)
                dat_path = current_file
                # target_path = pathlib.Path(fictrac_destination, current_file.name)
                target_path = pathlib.Path(
                    fictrac_destination, "fictrac_behavior_data.dat"
                )
                to_print = str(target_path)
                printlog(f"Transfering file{to_print:.>{WIDTH - 16}}")

                # put fictrac file path in into fly_dirs_dict
                current_fly_dir_dict = str(target_path).split(
                    fictrac_destination.parents[1].name
                )[-1]
                fly_dirs_dict[
                    destination_region.name + " Fictrac "
                ] = current_fly_dir_dict

    if user == 'jcsimon':
        fictrac_folder = pathlib.Path(
            "/oak/stanford/groups/trc/data/Jacob/Fictrac"
        )
        # when doing post-hoc fictrac, Bella's code where one compare the recording
        # timestamps of imaging and fictrac doesn't work anymore.
        # I instead use a deterministic file structure:
        # for example for fly 20231201\fly2\func1 imaging data, fictrac data must
        # be in the folder 20231201_fly2_func1. There must only be a single dat file in that folder.
        source_path = pathlib.Path(
            fictrac_folder,
            source_fly.parts[-3]
            + "_"
            + source_fly.parts[-2]
            + "_"
            + source_fly.parts[-1],
        )
        for current_file in source_path.iterdir():
            if "dat" in current_file.name:
                # width = 120
                # source_path = os.path.join(source_path, file)
                dat_path = current_file
                # target_path = pathlib.Path(fictrac_destination, current_file.name)
                target_path = pathlib.Path(
                    fictrac_destination, "fictrac_behavior_data.dat"
                )
                to_print = str(target_path)
                printlog(f"Transfering file{to_print:.>{WIDTH - 16}}")

                # put fictrac file path in into fly_dirs_dict
                current_fly_dir_dict = str(target_path).split(
                    fictrac_destination.parents[1].name
                )[-1]
                fly_dirs_dict[
                    destination_region.name + " Fictrac "
                ] = current_fly_dir_dict
    else:
        # fictrac_folder = os.path.join("/oak/stanford/groups/trc/data/fictrac", user)
        fictrac_folder = pathlib.Path("/oak/stanford/groups/trc/data/fictrac", user)

        # Find time of experiment based on functional.xml
        # true_ymd, true_total_seconds = get_expt_time(os.path.join(destination_region, 'imaging'))
        true_ymd, true_total_seconds = get_expt_time(
            pathlib.Path(destination_region, "imaging")
        )

        # printlog(f'true_ymd: {true_ymd}; true_total_seconds: {true_total_seconds}')

        # Find .dat file of 1) correct-ish time, 2) correct-ish size
        correct_date_and_size = []
        time_differences = []
        # for file in os.listdir(fictrac_folder):
        for file in fictrac_folder.iterdir():
            file = str(file)  # To be changed in the future
            # but I'm currently to lazy to change everything to
            # pathlib object below.

            # must be .dat file
            if ".dat" not in file:
                continue

            # Get datetime from file name
            datetime = file.split("-")[1][:-4]
            test_ymd = datetime.split("_")[0]
            test_time = datetime.split("_")[1]
            test_hour = test_time[0:2]
            test_minute = test_time[2:4]
            test_second = test_time[4:6]
            test_total_seconds = (
                int(test_hour) * 60 * 60 + int(test_minute) * 60 + int(test_second)
            )

            # Year/month/day must be exact
            if true_ymd != test_ymd:
                continue
            # printlog('Found file from same day: {}'.format(file))

            # Must be correct size
            # fp = os.path.join(fictrac_folder, file)
            fp = pathlib.Path(fictrac_folder, file)
            file_size = fp.stat().st_size
            # file_size = os.path.getsize(fp)
            if (
                file_size < 1000000
            ):  # changed to 1MB to accomidate 1 min long recordings. #30000000: #30MB
                # width = 120
                # printlog(F"Found correct .dat file{file:.>{width-23}}")
                # datetime_correct = datetime
                # break
                continue

            # get time difference from expt
            time_difference = np.abs(true_total_seconds - test_total_seconds)
            # Time must be within 10min
            if time_difference > 10 * 60:
                continue

            # if correct date and size append to list of potential file
            correct_date_and_size.append(file)
            time_differences.append(time_difference)

        # now that we have all potential files, pick the one with closest time
        # except clause will happen if empty list
        try:
            datetime_correct = correct_date_and_size[np.argmin(time_differences)]
        except:
            # width = 120
            printlog(
                f"{'   No fictrac data found --- continuing without fictrac data   ':*^{WIDTH}}"
            )
            return

        # Collect all fictrac files with correct datetime
        correct_time_files = [
            file for file in fictrac_folder.iterdir() if datetime_correct in file.name
        ]
        # correct_time_files = [file for file in os.listdir(fictrac_folder) if datetime_correct in file]

        # correct_time_files = []
        # for file in os.listdir(fictrac_folder):
        #     if datetime_correct in file:
        #         correct_time_files.append(file)

        # printlog('Found these files with correct times: {}'.format(correct_time_files))
        ##sys.stdout.flush()

        # Now transfer these 4 files to the fly
        fictrac_folder.mkdir()
        # os.mkdir(fictrac_destination)
        for file in correct_time_files:
            # width = 120
            target_path = pathlib.Path(fictrac_folder, file)
            source_path = pathlib.Path(fictrac_folder, file)
            # target_path = os.path.join(fictrac_destination, file)
            # source_path = os.path.join(fictrac_folder, file)
            # to_print = ('/').join(target_path.split('/')[-4:])
            to_print = str(target_path)
            printlog(f"Transfering file{to_print:.>{WIDTH - 16}}")
            # printlog('Transfering {}'.format(target_path))
            ##sys.stdout.flush()
    """


    ### Create empty xml file.
    # Update this later
    root = etree.Element("root")
    fictrac = objectify.Element("fictrac")
    root.append(fictrac)
    objectify.deannotate(root)
    etree.cleanup_namespaces(root)
    tree = etree.ElementTree(fictrac)
    # with open(os.path.join(fictrac_destination, 'fictrac.xml'), 'wb') as file:
    with open(pathlib.Path(fictrac_destination, "fictrac.xml"), "wb") as file:
        tree.write(file, pretty_print=True)

    return fly_dirs_dict


def create_imaging_json(xml_source_file, printlog):
    # Make empty dict
    source_data = {}

    # Get datetime
    try:
        datetime_str, _, _ = get_datetime_from_xml(xml_source_file)
    except:
        printlog("No xml or cannot read.")
        ##sys.stdout.flush()
        return
    date = datetime_str.split("-")[0]
    time = datetime_str.split("-")[1]
    source_data["date"] = str(date)
    source_data["time"] = str(time)

    # Get rest of data
    tree = objectify.parse(xml_source_file)
    source = tree.getroot()
    statevalues = source.findall("PVStateShard")[0].findall("PVStateValue")
    for statevalue in statevalues:
        key = statevalue.get("key")
        if key == "micronsPerPixel":
            indices = statevalue.findall("IndexedValue")
            for index in indices:
                axis = index.get("index")
                if axis == "XAxis":
                    source_data["x_voxel_size"] = float(index.get("value"))
                elif axis == "YAxis":
                    source_data["y_voxel_size"] = float(index.get("value"))
                elif axis == "ZAxis":
                    source_data["z_voxel_size"] = float(index.get("value"))
        if key == "laserPower":
            # This is not great - this is just the first pockels value
            indices = statevalue.findall("IndexedValue")
            laser_power_overall = int(float(indices[0].get("value")))
            source_data["laser_power"] = laser_power_overall
        if key == "laserWavelength":
            index = statevalue.findall("IndexedValue")
            laser_wavelength = int(float(index[0].get("value")))
            source_data["laser_wavelength"] = laser_wavelength
        if key == "pmtGain":
            indices = statevalue.findall("IndexedValue")
            for index in indices:
                index_num = index.get("index")
                # I changed this from 'red' and 'green' to the actual description used by the
                # microscope itself! Since we now have 2 Brukers, this seems safer!
                if index_num == "0":
                    source_data[index.get("description")] = int(float(index.get("value")))
                if index_num == "1":
                    source_data[index.get("description")] = int(float(index.get("value")))
                if index_num == "2":
                    source_data[index.get("description")] = int(float(index.get("value")))
        if key == "pixelsPerLine":
            source_data["x_dim"] = int(float(statevalue.get("value")))
        if key == "linesPerFrame":
            source_data["y_dim"] = int(float(statevalue.get("value")))
    sequence = source.findall("Sequence")[0]
    last_frame = sequence.findall("Frame")[-1]
    source_data["z_dim"] = int(last_frame.get("index"))

    # Save data
    # with open(os.path.join(os.path.split(xml_source_file)[0], 'scan.json'), 'w') as f:
    with open(pathlib.Path(xml_source_file.parent, "scan.json"), "w") as f:
        json.dump(source_data, f, indent=4)

def get_datetime_from_xml(xml_file):
    ##print('Getting datetime from {}'.format(xml_file))
    ##sys.stdout.flush()
    tree = ET.parse(xml_file)
    root = tree.getroot()
    datetime = root.get("date")
    # will look like "4/2/2019 4:16:03 PM" to start

    # Get dates
    date = datetime.split(" ")[0]
    month = date.split("/")[0]
    day = date.split("/")[1]
    year = date.split("/")[2]

    # Get times
    time = datetime.split(" ")[1]
    hour = time.split(":")[0]
    minute = time.split(":")[1]
    second = time.split(":")[2]

    # Convert from 12 to 24 hour time
    am_pm = datetime.split(" ")[-1]
    if am_pm == "AM" and hour == "12":
        hour = str(00)
    elif am_pm == "AM":
        pass
    elif am_pm == "PM" and hour == "12":
        pass
    else:
        hour = str(int(hour) + 12)

    # Add zeros if needed
    if len(month) == 1:
        month = "0" + month
    if len(day) == 1:
        day = "0" + day
    if len(hour) == 1:
        hour = "0" + hour

    # Combine
    datetime_str = year + month + day + "-" + hour + minute + second
    datetime_int = int(year + month + day + hour + minute + second)
    datetime_dict = {
        "year": year,
        "month": month,
        "day": day,
        "hour": hour,
        "minute": minute,
        "second": second,
    }

    return datetime_str, datetime_int, datetime_dict


def load_json(file):
    with open(file, "r") as f:
        data = json.load(f)
    return data


def load_xml(file):
    tree = objectify.parse(file)
    root = tree.getroot()
    return root

def add_fly_to_csv(import_folder, fly_folder, current_import_imaging_folder,
                   current_date, current_time, printlog):
    """
    brainsss originally had a xlsx file. However, it seemed to be a bit sensitive to slight
    changes in the keywords.
    I hope to address this (so that every built fly does have an entry, even if incomplete).
    """
    printlog('Adding fly to master_2P csv log')

    try:
        csv_path = pathlib.Path(fly_folder.parent, 'master_2P.csv')
        # Read csv, explicity state that first column is index
        csv_file = pd.read_csv(csv_path, index_col=0)
        printlog('Successfully opened master_2P log')
    except FileNotFoundError:
        try:
            # This should work if I don't move it out of that folder
            csv_path = pathlib.Path(fly_folder.parent, 'master_2P.csv')
            # Todo, could also have it pulled from github
            empty_csv_path = pathlib.Path('/oak/stanford/groups/trc/data/David/shared_files/master_2P.csv')
            shutil.copyfile(empty_csv_path, csv_path)
            # Read csv, explicity state that first column is index
            csv_file = pd.read_csv(csv_path, index_col=0)
            printlog('Successfully opened master_2P log')
        except Exception as e:
            # Better to create it from scratch!
            printlog('Error while trying to move or open the csv file:')
            print(e)

    # Load scan.json with the scanning parameters
    scan_json_path = pathlib.Path(fly_folder, current_import_imaging_folder.name,  "imaging/scan.json")
    scan_json = load_json(scan_json_path)
    printlog("Successfully loaded scan.json file")

    # Load fly.json.
    # At the moment this file is essential else snakebrainsss just doens't work.
    # So no need to wrap in try...except
    fly_json_path = pathlib.Path(fly_folder, "fly.json")
    fly_data = load_json(fly_json_path)
    printlog("Successfully loaded fly.json file")

    # Prepare dict for csv
    dict_for_csv = {}

    for column in csv_file.columns:
        if column == 'Fly ID':
            dict_for_csv[column] = fly_folder.name
        elif column == 'Expt ID':
            dict_for_csv[column] = current_import_imaging_folder.name
        elif column == 'Import folder':
            dict_for_csv[column] = import_folder.as_posix()
        elif column == 'Dataset folder':
            dict_for_csv[column] = fly_folder.as_posix()
        elif column == 'Date':
            dict_for_csv[column] = current_date
        elif column == 'Time':
            dict_for_csv[column] = current_time
        else:
            # Some data is in 'fly_data'
            if column in fly_data:
                dict_for_csv[column] = fly_data.get(column)
            # And other data is in 'scan_json'
            elif column in scan_json:
                dict_for_csv[column] = scan_json.get(column)
            else:
                printlog("master.csv file contains column: " + column + " which couldn't be assigned.")

    # It would be grand to do this also the other way around: If there are
    # fields in the json that do not have a column, create a column!
    for json_key in fly_data:
        if json_key not in csv_file.columns:
            dict_for_csv[json_key] = fly_data[json_key]

    for json_key in scan_json:
        # we already have date and time but because of capitalized/not capitalized would
        # make a new column. Hence explicitely exclude!
        if not json_key == 'date' and not json_key == 'time':
            if json_key not in csv_file.columns:
                dict_for_csv[json_key] = scan_json[json_key]
    # This of course has the slight downside that users have to be careful to
    # not slightly change their json file (i.e. from Genotype to genotype).
    # Upside is that the json file becomes very flexible in combination with this
    # csv: User x can just always add a key: visual stimulus and call it 'rectangle'
    # or 'circle' or whatever for easy sorting in the csv later on

    csv_file = pd.concat([csv_file, pd.DataFrame([dict_for_csv])], ignore_index=True)
    # Include an index as the first column!
    csv_file.to_csv(csv_path)
"""
def add_fly_to_xlsx(fly_folder, printlog):
    printlog("Adding fly to master_2P excel log")


    ### TRY TO LOAD ELSX ###
    try:
        xlsx_path = pathlib.Path(fly_folder.parent, "master_2P.xlsx")
        wb = load_workbook(filename=xlsx_path, read_only=False)
        ws = wb.active
        printlog("Successfully opened master_2P log")
    except Exception as e:
        printlog(
            "FYI you have no excel metadata sheet found, so unable to append metadata for this fly."
        )
        printlog(traceback.format_exc())
        return

    ### TRY TO LOAD FLY METADATA ###
    try:
        fly_file = pathlib.Path(fly_folder, "fly.json")
        # fly_file = os.path.join(fly_folder, 'fly.json')
        fly_data = load_json(fly_file)
        printlog("Successfully loaded fly.json")
    except:
        printlog(
            "FYI no *fly.json* found; this will not be logged in your excel sheet."
        )
        fly_data = {}
        fly_data["circadian_on"] = None
        fly_data["circadian_off"] = None
        fly_data["gender"] = None
        fly_data["age"] = None
        fly_data["temp"] = None
        fly_data["notes"] = None
        fly_data["date"] = None
        fly_data["genotype"] = None

    # Write in master xlsx only if we have a func folder. Ignore anatomical data!
    expt_folders = [
        pathlib.Path(fly_folder, x) for x in fly_folder.iterdir() if "func" in x.name
    ]
    # brainsss.sort_nicely(expt_folders)
    expt_folders = natsort.natsorted(expt_folders)
    for expt_folder in expt_folders:
        print("expt_folder" + repr(expt_folder))

        ### TRY TO LOAD EXPT METADATA ###
        try:
            # expt_file = os.path.join(expt_folder, 'expt.json')
            expt_file = pathlib.Path(expt_folders, "expt.json")
            expt_data = load_json(expt_file)
            printlog("Sucessfully loaded expt.json")
        except:
            printlog(
                "FYI no *expt.json* found; this will not be logged in your excel sheet."
            )
            expt_data = {}
            expt_data["brain_area"] = None
            expt_data["notes"] = None
            expt_data["time"] = None

        ### TRY TO LOAD SCAN DATA ###
        try:
            scan_file = pathlib.Path(expt_folder, "imaging", "scan.json")
            # scan_file = os.path.join(expt_folder, 'imaging', 'scan.json')
            scan_data = load_json(scan_file)
            scan_data["x_voxel_size"] = "{:.1f}".format(scan_data["x_voxel_size"])
            scan_data["y_voxel_size"] = "{:.1f}".format(scan_data["y_voxel_size"])
            scan_data["z_voxel_size"] = "{:.1f}".format(scan_data["z_voxel_size"])
            printlog("Sucessfully loaded scan.json")
        except:
            printlog(
                "FYI no *scan.json* found; this will not be logged in your excel sheet."
            )
            scan_data = {}
            scan_data["laser_power"] = None
            scan_data["PMT_green"] = None
            scan_data["PMT_red"] = None
            scan_data["x_dim"] = None
            scan_data["y_dim"] = None
            scan_data["z_dim"] = None
            scan_data["x_voxel_size"] = None
            scan_data["y_voxel_size"] = None
            scan_data["z_voxel_size"] = None

        visual_file = pathlib.Path(expt_folder, "visual", "visual.json")
        # visual_file = os.path.join(expt_folder, 'visual', 'visual.json')
        try:
            visual_data = load_json(visual_file)
            visual_input = visual_data[0]["name"] + " ({})".format(len(visual_data))
        except:
            visual_input = None

        # Get fly_id
        # fly_folder = expt_folder.parent
        # fly_folder = os.path.split(os.path.split(expt_folder)[0])[-1]
        fly_id = fly_folder.name.split("_")[-1]
        printlog(f"Got fly ID as {fly_id}")

        # Get expt_id
        expt_id = expt_folder.name  # probably 'func1' etc.
        printlog(f"Got expt ID as {expt_id}")
        # expt_id = 'NA' # Not sure what this is, NA for now

        # Append the new row
        # new_row = []
        new_row = [
            int(fly_id),
            str(expt_id),
            fly_data["date"],
            expt_data["brain_area"],
            fly_data["genotype"],
            visual_input,
            None,
            fly_data["notes"],
            expt_data["notes"],
            expt_data["time"],
            fly_data["circadian_on"],
            fly_data["circadian_off"],
            fly_data["gender"],
            fly_data["age"],
            fly_data["temp"],
            scan_data["laser_power"],
            scan_data["PMT_green"],
            scan_data["PMT_red"],
            scan_data["x_dim"],
            scan_data["y_dim"],
            scan_data["z_dim"],
            scan_data["x_voxel_size"],
            scan_data["y_voxel_size"],
            scan_data["z_voxel_size"],
        ]

        ws.append(new_row)
        printlog(f"Appended {fly_id} {expt_id}")

    # Save the file
    wb.save(xlsx_path)
    printlog("master_2P successfully updated")"""
