import natsort
import pathlib
import json
import datetime
scripts_path = workflow.basedir # Exposes path to this file
from brainsss import utils
from scripts import preprocess
from scripts import snake_utils
import os
import sys

print(os.getcwd())
current_user = config['user']
settings = utils.load_user_settings(current_user)
dataset_path = pathlib.Path(settings['dataset_path'])

########################################################

#>>>>
#fictrac_fps = 100 # AUTOMATE THIS!!!! ELSE FOR SURE A MISTAKE WILL HAPPEN IN THE FUTURE!!!!
# TODO!!!! Instead of just believing a framerate, use the voltage signal recorded during imaging
# that defines the position of a given frame!
#<<<<

# First n frames to average over when computing mean/fixed brain | Default None
# (average over all frames).
meanbrain_n_frames =  None

##########################################################

# On sherlock this is usually python3 but on a personal computer can be python
shell_python_command = str(settings.get('shell_python_command', "python3"))
print("shell_python_command " + shell_python_command)
moco_temp_folder = str(settings.get('moco_temp_folder', "/scratch/groups/trc"))

# Define path to imports to find fly.json!
#fly_folder_to_process_oak = pathlib.Path(dataset_path,fly_folder_to_process)
fly_folder_to_process_oak = pathlib.Path(os.getcwd())
print('Analyze data in ' + repr(fly_folder_to_process_oak.as_posix()))

# Read channel information from fly.json file
# If fails here, means the folder specified doesn't exist. Check name.
# Note: Good place to let user know to check folder and exit!
with open(pathlib.Path(fly_folder_to_process_oak, 'fly.json'), 'r') as file:
    fly_json = json.load(file)

# fictrac settings from fly.json (if not present, default values are used)
fictrac_fps = fly_json.get('fictrac_fps', 100)
# reg/moco settings from fly.json (if not present, default values are used)
moco_transform_type = fly_json.get('moco_transform_type', 'StackReg.RIGID_BODY')
moco_smooth = fly_json.get('moco_smooth', False)
moco_avg_wid = fly_json.get('moco_avg_wid', 5) # only relevant if smooth is True
moco_mean_frames = fly_json.get('moco_mean_frames', 40)
cores = fly_json.get('cores', 40)
reg_transform_type = fly_json.get('reg_transform_type', 'StackReg.RIGID_BODY')

# anatomy specific moco settings from fly.json (if not present, default values are used)
moco_smooth_anat = fly_json.get('moco_smooth_anat', False)
moco_avg_wid_anat = fly_json.get('moco_avg_wid_anat', 1) # only relevant if smooth is True
moco_mean_frames_anat = fly_json.get('moco_mean_frames_anat', 40)

# This needs to come from some sort of json file the experimenter
# creates while running the experiment. Same as genotype.
FUNCTIONAL_CHANNELS = fly_json['functional_channel']
# It is probably necessary to forcibly define STRUCTURAL_CHANNEL if not defined
# Would be better to have an error to be explicit!

# Throw an error if missing! User must provide this!
STRUCTURAL_CHANNEL = fly_json['structural_channel']
if STRUCTURAL_CHANNEL != 'channel_1' and \
    STRUCTURAL_CHANNEL != 'channel_2' and \
        STRUCTURAL_CHANNEL != 'channel_3':
    print('!!! ERROR !!!')
    print('You must provide "structural_channel" in the "fly.json" file for snake_brainsss to run!')
    sys.exit()
    # This would be a implicit fix. Not great as it'll
    # hide potential bugs. Better explicit
    #STRUCTURAL_CHANNEL = FUNCTIONAL_CHANNELS[0]


# Bool for which channel exists in this particular recording.
# We currently assume that for folders with 'func' in the folder name
# all channels are recorded from for a given fly! Else this will break.
CH1_EXISTS_FUNC = snake_utils.ch_exists("1", FUNCTIONAL_CHANNELS)
CH2_EXISTS_FUNC = snake_utils.ch_exists("2", FUNCTIONAL_CHANNELS)
CH3_EXISTS_FUNC = snake_utils.ch_exists("3", FUNCTIONAL_CHANNELS)

# IMPORTANT: For the folder with 'anat' in the name, we only consider
#
CH1_EXISTS_STRUCT = snake_utils.ch_exists("1", STRUCTURAL_CHANNEL)
CH2_EXISTS_STRUCT = snake_utils.ch_exists("2", STRUCTURAL_CHANNEL)
CH3_EXISTS_STRUCT = snake_utils.ch_exists("3", STRUCTURAL_CHANNEL)

####
# Load fly_dir.json
####
fly_dirs_dict_path = pathlib.Path(fly_folder_to_process_oak, fly_folder_to_process_oak.name + '_dirs.json')
with open(pathlib.Path(fly_dirs_dict_path),'r') as file:
    fly_dirs_dict = json.load(file)

#####
# Prepare filepaths to be used
#####
# Snakemake runs acyclical, meaning it checks which input depends on the output of which rule
# in order to parallelize a given snakefile.
# I'll therefore keep variables with lists of paths that will be feed into a given rule.
# These lists of paths are created here.

# Imaging data paths
imaging_file_paths = []
fictrac_file_paths = []
for key in fly_dirs_dict:
    if 'Imaging' in key:
        imaging_file_paths.append(fly_dirs_dict[key][1::])
        # this yields for example 'func2/imaging'
    elif 'Fictrac' in key:
        fictrac_file_paths.append(fly_dirs_dict[key][1::])
        # This yields for example 'func1/fictrac/fictrac_behavior_data.dat'
        # With automatic stimpack transfer it'll return "/func0/stimpack/loco/fictrac_behavior_data.dat"

if len(fictrac_file_paths)>0:
    fictrac_exists = True
else:
    fictrac_exists = False

#######
# Data path on OAK
#######

FICTRAC_PATHS = []
if fictrac_exists:
    for current_path in fictrac_file_paths:
        FICTRAC_PATHS.append(current_path.split('/fictrac_behavior_data.dat')[0])
    # Fictrac data can be in different folders! For correlation, need to know the
    # relative path following 'funcX'.
    # IN ONE EXPERIMENT I ASSUME THAT THE FICTRAC STRUCTURE IS CONSISTENT!
    fictrac_rel_path_correlation = None
    current_fictrac_rel_path = FICTRAC_PATHS[0]
    # Remove the first folder which is going to be likely 'func0'
    rel_path_parts = pathlib.Path(current_fictrac_rel_path).parts[1::]
    # Then put the parts back together
    fictrac_rel_path_correlation = pathlib.Path(*rel_path_parts)

# For wildcards we need lists of elements of the path for each folder.
list_of_paths = []
for current_path in imaging_file_paths:
    list_of_paths.append(current_path.split('/imaging')[0])
# This is a list of all imaging paths so something like this 
# ['anat0', 'func0', 'func1']
print('list_of_paths ' +repr(list_of_paths) )

list_of_paths_func = []
for current_path in imaging_file_paths:
    if 'func' in current_path:
        list_of_paths_func.append(current_path.split('/imaging')[0])

print("list_of_paths_func " + repr(list_of_paths_func))

# list of func paths for registration to first series (same as list_of_paths_func but without first series)
list_of_paths_func_reg = natsort.natsorted(list_of_paths_func) [1:]
fixed_path_func_reg = natsort.natsorted(list_of_paths_func)[0]
print("fixed_path_func_reg " + repr(fixed_path_func_reg))
print("list_of_paths_func_reg " + repr(list_of_paths_func_reg))

list_of_paths_struct = []
for current_path in imaging_file_paths:
    if 'anat' in current_path:
        list_of_paths_struct.append(current_path.split('/imaging')[0])
        print("list_of_paths_struct" + repr(list_of_paths_struct))
if len(list_of_paths_struct) > 1:
    print('!!!WARNING!!!')
    print('The following folders have the "anat" keyword:')
    print(list_of_paths_struct)
    print('The folder ' + repr(natsort.natsorted(list_of_paths_struct)[0]) + ' will be treated as the "main" anat folder.')
    print('The other folder(s) will be ignored for this analysis. To get moco for the other folders, use the "misc_imaging" keyword!')
    list_of_paths_struct = natsort.natsorted(list_of_paths_struct)[0]
print('list_of_paths_struct' + repr(list_of_paths_struct))

# Folders with this keyword are run through the pipeline up to moco_mean
list_of_paths_misc_imaging = []
for current_path in imaging_file_paths:
    if 'misc' in current_path:
        list_of_paths_misc_imaging.append(current_path.split('/imaging')[0])
print('list_of_paths_misc_imaging ' + repr(list_of_paths_misc_imaging))

list_of_channels_func = []
if CH1_EXISTS_FUNC:
    list_of_channels_func.append("1")
if CH2_EXISTS_FUNC:
    list_of_channels_func.append("2")
if CH3_EXISTS_FUNC:
    list_of_channels_func.append("3")
print("list_of_channels_func" + repr(list_of_channels_func))

# For func moco we want to use the structural channel if it exists!
CH1_EXISTS_FUNC_MOCO = False
CH2_EXISTS_FUNC_MOCO = False
CH3_EXISTS_FUNC_MOCO = False
list_of_channels_for_func_moco = []
if CH1_EXISTS_FUNC or CH1_EXISTS_STRUCT:
    CH1_EXISTS_FUNC_MOCO = True
    list_of_channels_for_func_moco.append("1")
if CH2_EXISTS_FUNC or CH2_EXISTS_STRUCT:
    CH2_EXISTS_FUNC_MOCO = True
    list_of_channels_for_func_moco.append("2")
if CH3_EXISTS_FUNC or CH3_EXISTS_STRUCT:
    CH3_EXISTS_FUNC_MOCO = True
    list_of_channels_for_func_moco.append("2")

list_of_channels_struct = []
if CH1_EXISTS_STRUCT:
    list_of_channels_struct.append("1")
if CH2_EXISTS_STRUCT:
    list_of_channels_struct.append("2")
if CH3_EXISTS_STRUCT:
    list_of_channels_struct.append("3")
print("list_of_channels_struct" + repr(list_of_channels_struct))

list_of_channels_misc = []
# Since we are only looking at a single fly, go into the misc folder
# to collect the channels that exist.
list_of_misc_channels = []
for counter, current_misc_folder in enumerate(list_of_paths_misc_imaging):
    temp = []
    for current_file in pathlib.Path(fly_folder_to_process_oak, current_misc_folder, 'imaging').iterdir():
        if 'channel' in current_file.name:
            temp.append(current_file.name.split('.nii')[0].split('channel_')[-1])
    if counter == 0:
        list_of_misc_channels = temp
    else:
        if list_of_misc_channels != temp:
            import sys
            print('!!!!!!!!!!!!ERROR!!!!!!!!!!!')
            print('Your misc_imaging channels have different channel settings:')
            print('The current misc folder: ' + repr(current_misc_folder) + ' has these channels: ' + repr(temp))
            print('The previous misco folder has: ' + repr(list_of_misc_channels))
            print('snakemake can not currently handle mixture of channels for misc.')
            print('Exiting')
            sys.exit(1)

if '1' in list_of_misc_channels:
    CH1_EXISTS_MISC = True
else:
    CH1_EXISTS_MISC = False
if '2' in list_of_misc_channels:
    CH2_EXISTS_MISC = True
else:
    CH2_EXISTS_MISC = False
if '3' in list_of_misc_channels:
    CH3_EXISTS_MISC = True
else:
    CH3_EXISTS_MISC = False

print("list_of_misc_channels" + repr(list_of_misc_channels))


rule all:
    """
    See: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
        By default snakemake executes the first rule in the snakefile. This gives rise to pseudo-rules at the beginning 
        of the file that can be used to define build-targets similar to GNU Make
    Or in other words: Here we define which file we want at the end. Snakemake checks which one is there and which 
    one is missing. It then uses the other rules to see how it can produce the missing files.
    """
    threads: 1 # should be sufficent
    resources: mem_mb=1000 # should be sufficient
    input:
        ###
        # Bleaching QC
        # Since func and struct can have different channels, seperate the two
        ###
        expand(str(fly_folder_to_process_oak)
               + "/{bleaching_imaging_paths}/imaging/bleaching_func.png",
            bleaching_imaging_paths=list_of_paths_func),
        ##
        expand(str(fly_folder_to_process_oak)
               + "/{bleaching_imaging_paths}/imaging/bleaching_struct.png",
            bleaching_imaging_paths=list_of_paths_struct),

        expand(str(fly_folder_to_process_oak)
               + "/{bleaching_imaging_paths}/imaging/bleaching_struct.png",
               bleaching_imaging_paths=list_of_paths_misc_imaging),

        ###
        # Fictrac QC
        ###
        
        expand(str(fly_folder_to_process_oak)
               + "/{fictrac_paths}/fictrac_2d_hist_fixed.png",# not sure if necessary: if fictrac_exists else [],
            fictrac_paths=FICTRAC_PATHS),
        # data in fly_dirs.json!

        ###
        # Meanbrain
        ###
        expand(str(fly_folder_to_process_oak)
               + "/{meanbr_imaging_paths_func}/imaging/channel_{meanbr_ch_func}_mean_func.nii",
            meanbr_imaging_paths_func=list_of_paths_func,
            meanbr_ch_func=list_of_channels_for_func_moco),
        ##
        expand(str(fly_folder_to_process_oak)
               + "/{meanbr_imaging_paths_struct}/imaging/channel_{meanbr_ch_struct}_mean_struct.nii",
            meanbr_imaging_paths_struct=list_of_paths_struct,
            meanbr_ch_struct=list_of_channels_struct),
        expand(str(fly_folder_to_process_oak)
               + "/{meanbr_imaging_paths_misc}/imaging/channel_{meanbr_ch_misc}_mean_misc.nii",
               meanbr_imaging_paths_misc=list_of_paths_misc_imaging,
               meanbr_ch_misc=list_of_misc_channels),

        ###
        # Meanbrain BG
        ###
        expand(str(fly_folder_to_process_oak)
                + "/{meanbr_imaging_paths_func}/imaging/bg/channel_2_bg_mean_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
                meanbr_imaging_paths_func=list_of_paths_func),

        ###
        # Motion correction output FUNC
        # The idea is to use the structural channel for moco!
        ###
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func}/moco/tmats_func.npy",
               moco_imaging_paths_func=list_of_paths_func),

        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func}/moco/channel_1_moco_func.nii" if CH1_EXISTS_FUNC_MOCO else [],
            moco_imaging_paths_func=list_of_paths_func),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func}/moco/channel_2_moco_bg_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
            moco_imaging_paths_func=list_of_paths_func),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func}/moco/channel_3_moco_func.nii" if CH3_EXISTS_FUNC_MOCO else [],
            moco_imaging_paths_func=list_of_paths_func),

        ###
        # Motion correction output STRUCT
        ###
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_struct}/moco/tmats_struct.npy",
               moco_imaging_paths_struct=list_of_paths_struct),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_struct}/moco/channel_1_moco_struct.nii" if CH1_EXISTS_STRUCT else [],
               moco_imaging_paths_struct=list_of_paths_struct),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_struct}/moco/channel_2_moco_struct.nii" if CH2_EXISTS_STRUCT else [],
               moco_imaging_paths_struct=list_of_paths_struct),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_struct}/moco/channel_3_moco_struct.nii" if CH3_EXISTS_STRUCT else [],
               moco_imaging_paths_struct=list_of_paths_struct),

        ###
        # Motion correction output MISC
        ###
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_misc}/moco/tmats_misc.npy",
               moco_imaging_paths_misc=list_of_paths_misc_imaging),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_misc}/moco/channel_1_moco_misc.nii" if CH1_EXISTS_MISC else [],
               moco_imaging_paths_misc=list_of_paths_misc_imaging),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_misc}/moco/channel_2_moco_misc.nii" if CH2_EXISTS_MISC else [],
               moco_imaging_paths_misc=list_of_paths_misc_imaging),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_misc}/moco/channel_3_moco_misc.nii" if CH3_EXISTS_MISC else [],
               moco_imaging_paths_misc=list_of_paths_misc_imaging),

        ###
        # Meanbrain of moco brain
        ###
        expand(str(fly_folder_to_process_oak)
               + "/{moco_meanbr_imaging_paths_func}/moco/channel_1_moco_mean_func.nii" if CH1_EXISTS_FUNC_MOCO else [],
            moco_meanbr_imaging_paths_func=list_of_paths_func),

        expand(str(fly_folder_to_process_oak)
               + "/{moco_meanbr_imaging_paths_func}/moco/channel_2_moco_bg_mean_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
            moco_meanbr_imaging_paths_func=list_of_paths_func),

        expand(str(fly_folder_to_process_oak)
               + "/{moco_meanbr_imaging_paths_func}/moco/channel_3_moco_mean_func.nii" if CH3_EXISTS_FUNC_MOCO else [],
            moco_meanbr_imaging_paths_func=list_of_paths_func),

        #
        expand(str(fly_folder_to_process_oak)
               + "/{moco_meanbr_imaging_paths_struct}/moco/channel_{meanbr_moco_ch_struct}_moco_mean_struct.nii",
            moco_meanbr_imaging_paths_struct=list_of_paths_struct,
            meanbr_moco_ch_struct=list_of_channels_struct),
        #
        expand(str(fly_folder_to_process_oak)
               + "/{moco_meanbr_imaging_paths_misc}/moco/channel_{meanbr_moco_ch_misc}_moco_mean_struct.nii",
               moco_meanbr_imaging_paths_misc=list_of_paths_misc_imaging,
               meanbr_moco_ch_misc=list_of_channels_misc),
        ####

        ###
        # Background Subtraction (line by line)
        ###
        # jcs run background_subtract on ch2 only?
        expand(str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_bg_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
                moco_imaging_paths_func=list_of_paths_func),
        expand(str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_before_bg_removal.png" if CH2_EXISTS_FUNC_MOCO else [],
                moco_imaging_paths_func=list_of_paths_func),
        expand(str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_after_bg_removal.png" if CH2_EXISTS_FUNC_MOCO else [],
                moco_imaging_paths_func=list_of_paths_func),
        expand(str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_bg_selection.tif" if CH2_EXISTS_FUNC_MOCO else [],
                moco_imaging_paths_func=list_of_paths_func),
    
        ###
        # register_series
        ###
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func_reg}/moco/tmats_func_reg.npy",
               moco_imaging_paths_func_reg=list_of_paths_func_reg),

        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func_reg}/moco/channel_1_moco_func_reg.nii" if CH1_EXISTS_FUNC_MOCO else [],
            moco_imaging_paths_func_reg=list_of_paths_func_reg),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func_reg}/moco/channel_2_moco_bg_func_reg.nii" if CH2_EXISTS_FUNC_MOCO else [],
            moco_imaging_paths_func_reg=list_of_paths_func_reg),
        expand(str(fly_folder_to_process_oak)
               + "/{moco_imaging_paths_func_reg}/moco/channel_3_moco_func_reg.nii" if CH3_EXISTS_FUNC_MOCO else [],
            moco_imaging_paths_func_reg=list_of_paths_func_reg),

rule fly_builder_rule:
    threads:
        1
    resources:
        mem_mb=snake_utils.mem_mb_times_threads,
        runtime='10m' # should generally be sufficient

    run:
        build_fly.fly_builder(
            autotransferred_stimpack=autotransferred_stimpack,
            fictrac_folder_path=fictrac_folder_path,
            import_dirs= all_imports_paths,
            dataset_dirs = all_fly_dataset_paths
                                                 )

rule fictrac_qc_rule:
    """
    """
    threads: snake_utils.threads_per_memory
    resources:
        mem_mb=snake_utils.mem_mb_times_threads,
        runtime='10m'
    input:
        str(fly_folder_to_process_oak) + "/{fictrac_paths}/fictrac_behavior_data.dat"
    output:
        str(fly_folder_to_process_oak) + "/{fictrac_paths}/fictrac_2d_hist_fixed.png"
    run:
        try:
            preprocess.fictrac_qc(fly_folder_to_process_oak,
                                    fictrac_file_path= input,
                                    fictrac_fps=fictrac_fps # AUTOMATE THIS!!!! ELSE BUG PRONE!!!!
                                    )
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_fictrac_qc_rule')
            utils.write_error(logfile=logfile,
                                 error_stack=error_stack)

rule bleaching_qc_rule_func:
    """
    """
    threads: snake_utils.threads_per_memory_less
    resources:
        mem_mb=snake_utils.mem_mb_less_times_input, # This is probably overkill todo decrease!
        runtime='10m' # In my test cases it was never more than 5 minutes!
    input:
        brains_paths_ch1=str(fly_folder_to_process_oak) + "/{bleaching_imaging_paths}/imaging/channel_1.nii" if CH1_EXISTS_FUNC else [],
        brains_paths_ch2=str(fly_folder_to_process_oak) + "/{bleaching_imaging_paths}/imaging/channel_2.nii" if CH2_EXISTS_FUNC else [],
        brains_paths_ch3=str(fly_folder_to_process_oak) + "/{bleaching_imaging_paths}/imaging/channel_3.nii" if CH3_EXISTS_FUNC else [],
    output:
        str(fly_folder_to_process_oak) +"/{bleaching_imaging_paths}/imaging/bleaching_func.png"
    run:
        try:
            preprocess.bleaching_qc(fly_directory=fly_folder_to_process_oak,
                                        path_to_read=[input.brains_paths_ch1, input.brains_paths_ch2, input.brains_paths_ch3], #imaging_paths_by_folder_scratch, # {input} didn't work, I think because it destroyed the list of list we expect to see here #imaging_paths_by_folder_scratch,
                                        path_to_save=output, # can't use output, messes things up here! #imaging_paths_by_folder_oak
                                        #print_output = output
            )
            print('Done with bleaching_qc_func')
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_bleaching_qc_rule')
            utils.write_error(logfile=logfile,
                error_stack=error_stack)
            print('Error with bleaching_qc' )

rule bleaching_qc_rule_struct:
    """
    """
    threads: snake_utils.threads_per_memory_less
    resources:
        mem_mb=snake_utils.mem_mb_less_times_input, # This is probably overkill todo decrease!
        runtime='10m' # In my test cases it was never more than 5 minutes!
    input:
        brains_paths_ch1=str(fly_folder_to_process_oak) + "/{bleaching_imaging_paths}/imaging/channel_1.nii" if CH1_EXISTS_STRUCT else [],
        brains_paths_ch2=str(fly_folder_to_process_oak) + "/{bleaching_imaging_paths}/imaging/channel_2.nii" if CH2_EXISTS_STRUCT else [],
        brains_paths_ch3=str(fly_folder_to_process_oak) + "/{bleaching_imaging_paths}/imaging/channel_3.nii" if CH3_EXISTS_STRUCT else [],
    output:
        str(fly_folder_to_process_oak) +"/{bleaching_imaging_paths}/imaging/bleaching_struct.png"
    run:
        try:
            preprocess.bleaching_qc(fly_directory=fly_folder_to_process_oak,
                                        path_to_read=[input.brains_paths_ch1, input.brains_paths_ch2, input.brains_paths_ch3], #imaging_paths_by_folder_scratch, # {input} didn't work, I think because it destroyed the list of list we expect to see here #imaging_paths_by_folder_scratch,
                                        path_to_save=output, # can't use output, messes things up here! #imaging_paths_by_folder_oak
                                        #print_output = output
            )
            print('Done with bleaching_qc')
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_bleaching_qc_rule')
            utils.write_error(logfile=logfile,
                error_stack=error_stack)
            print('Error with bleaching_qc' )


rule make_mean_brain_rule_func:
    """
    """
    threads: snake_utils.threads_per_memory_less
    resources:
        mem_mb=snake_utils.mem_mb_less_times_input,  #snake_utils.mem_mb_times_input #mem_mb=snake_utils.mem_mb_more_times_input
        runtime='10m' # should be enough
    input:
            str(fly_folder_to_process_oak) + "/{meanbr_imaging_paths_func}/imaging/channel_{meanbr_ch_func}.nii"
    output:
            str(fly_folder_to_process_oak) + "/{meanbr_imaging_paths_func}/imaging/channel_{meanbr_ch_func}_mean_func.nii"
    run:
        try:
            preprocess.make_mean_brain(fly_directory=fly_folder_to_process_oak,
                meanbrain_n_frames=meanbrain_n_frames,
                path_to_read=input,
                path_to_save=output,
                rule_name='make_mean_brain_rule')
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_make_mean_brain')
            utils.write_error(logfile=logfile,
                error_stack=error_stack)

rule make_mean_brain_rule_bg_func:
    """
    """
    threads: snake_utils.threads_per_memory_less
    resources:
        mem_mb=snake_utils.mem_mb_less_times_input,  #snake_utils.mem_mb_times_input #mem_mb=snake_utils.mem_mb_more_times_input
        runtime='10m' # should be enough
    input:
            str(fly_folder_to_process_oak) + "/{meanbr_imaging_paths_func}/imaging/bg/channel_2_bg_func.nii"
    output:
            str(fly_folder_to_process_oak) + "/{meanbr_imaging_paths_func}/imaging/bg/channel_2_bg_mean_func.nii"
    run:
        try:
            preprocess.make_mean_brain(fly_directory=fly_folder_to_process_oak,
                meanbrain_n_frames=meanbrain_n_frames,
                path_to_read=input,
                path_to_save=output,
                rule_name='make_mean_brain_rule_bg')
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_make_mean_brain_bg')
            utils.write_error(logfile=logfile,
                error_stack=error_stack)

rule make_mean_brain_rule_struct:
    """
    """
    threads: snake_utils.threads_per_memory_less
    resources:
        mem_mb=snake_utils.mem_mb_less_times_input,  #snake_utils.mem_mb_times_input #mem_mb=snake_utils.mem_mb_more_times_input
        runtime='10m' # should be enough
    input:
            str(fly_folder_to_process_oak) + "/{meanbr_imaging_paths_struct}/imaging/channel_{meanbr_ch_struct}.nii"
    output:
            str(fly_folder_to_process_oak) + "/{meanbr_imaging_paths_struct}/imaging/channel_{meanbr_ch_struct}_mean_struct.nii"
    run:
        try:
            preprocess.make_mean_brain(fly_directory=fly_folder_to_process_oak,
                meanbrain_n_frames=meanbrain_n_frames,
                path_to_read=input,
                path_to_save=output,
                rule_name='make_mean_brain_rule')
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_make_mean_brain')
            utils.write_error(logfile=logfile,
                error_stack=error_stack)

rule motion_correction_parallel_slice_func:
    # separate slices, 2d moco, restitch
    threads:
        40
    resources:
        mem_mb=snake_utils.mb_for_moco_input_slice, #.mem_mb_much_more_times_input,
        runtime=snake_utils.time_for_moco_input # runtime takes input as seconds!
    input:
        # Only use the Channels that exists - this organizes the anatomy and functional paths inside the motion correction
        # module.
        brain_paths_ch1=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/channel_1.nii" if CH1_EXISTS_FUNC_MOCO else [],
        brain_paths_ch2=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_bg_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
        brain_paths_ch3=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/channel_3.nii" if CH3_EXISTS_FUNC_MOCO else [],

    output:
        moco_path_ch1=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/moco/channel_1_moco_func.nii" if CH1_EXISTS_FUNC_MOCO else[],
        moco_path_ch2=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/moco/channel_2_moco_bg_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
        moco_path_ch3=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/moco/channel_3_moco_func.nii" if CH3_EXISTS_FUNC_MOCO else [],
        par_output=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/moco/tmats_func.npy",

    shell: shell_python_command + " " + scripts_path + "/scripts/moco_parallel_slice_turboreg_dev.py "
        "--fly_directory {fly_folder_to_process_oak} "
        "--dataset_path {dataset_path} "
        "--brain_paths_ch1 {input.brain_paths_ch1} "
        "--brain_paths_ch2 {input.brain_paths_ch2} "
        "--brain_paths_ch3 {input.brain_paths_ch3} "
        "--STRUCTURAL_CHANNEL {STRUCTURAL_CHANNEL} "
        "--FUNCTIONAL_CHANNELS {FUNCTIONAL_CHANNELS} "
        "--moco_path_ch1 {output.moco_path_ch1} "
        "--moco_path_ch2 {output.moco_path_ch2} "
        "--moco_path_ch3 {output.moco_path_ch3} "
        "--par_output {output.par_output} "
        "--moco_temp_folder {moco_temp_folder} "
        "--moco_transform_type {moco_transform_type} "
        "--moco_smooth {moco_smooth} "
        "--moco_avg_wid {moco_avg_wid} "
        "--moco_mean_frames {moco_mean_frames} "
        "--cores {cores} "

rule motion_correction_parallel_slice_struct:
    # separate slices, 2d moco, restitch
    threads:
        40
    resources:
        mem_mb=snake_utils.mb_for_moco_input_slice, #.mem_mb_much_more_times_input,
        runtime=snake_utils.time_for_moco_input # runtime takes input as seconds!
    input:
        # Only use the Channels that exists - this organizes the anatomy and functional paths inside the motion correction
        # module.
        brain_paths_ch1=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_struct}/imaging/channel_1.nii" if CH1_EXISTS_STRUCT else [],
        brain_paths_ch2=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_struct}/imaging/channel_2.nii" if CH2_EXISTS_STRUCT else [],
        brain_paths_ch3=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_struct}/imaging/channel_3.nii" if CH3_EXISTS_STRUCT else [],

    output:
        moco_path_ch1=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_struct}/moco/channel_1_moco_struct.nii" if CH1_EXISTS_STRUCT else [],
        moco_path_ch2=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_struct}/moco/channel_2_moco_struct.nii" if CH2_EXISTS_STRUCT else [],
        moco_path_ch3=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_struct}/moco/channel_3_moco_struct.nii" if CH3_EXISTS_STRUCT else [],
        par_output=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_struct}/moco/tmats_struct.npy",

    shell: shell_python_command + " " + scripts_path + "/scripts/moco_parallel_slice_turboreg_dev.py "
                                       "--fly_directory {fly_folder_to_process_oak} "
                                       "--dataset_path {dataset_path} "
                                       "--brain_paths_ch1 {input.brain_paths_ch1} "
                                       "--brain_paths_ch2 {input.brain_paths_ch2} "
                                       "--brain_paths_ch3 {input.brain_paths_ch3} "
                                       "--STRUCTURAL_CHANNEL {STRUCTURAL_CHANNEL} "
                                       "--FUNCTIONAL_CHANNELS {FUNCTIONAL_CHANNELS} "
                                       "--moco_path_ch1 {output.moco_path_ch1} "
                                       "--moco_path_ch2 {output.moco_path_ch2} "
                                       "--moco_path_ch3 {output.moco_path_ch3} "
                                       "--par_output {output.par_output} "
                                       "--moco_temp_folder {moco_temp_folder} "
                                       "--moco_transform_type {moco_transform_type} "
                                       "--moco_smooth {moco_smooth_anat} " # use anat settings
                                       "--moco_avg_wid {moco_avg_wid_anat} " # use anat settings
                                       "--moco_mean_frames {moco_mean_frames_anat} " # use anat settings
                                       "--cores {cores} "

rule moco_mean_brain_rule_func:
    """
    """
    threads: snake_utils.threads_per_memory
    resources:
        mem_mb=snake_utils.mem_mb_times_input,
        runtime='10m'# should be enough
    input:
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_func}/moco/channel_1_moco_func.nii" if CH1_EXISTS_FUNC_MOCO else [],
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_func}/moco/channel_2_moco_bg_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_func}/moco/channel_3_moco_func.nii" if CH3_EXISTS_FUNC_MOCO else [],
    output:
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_func}/moco/channel_1_moco_mean_func.nii" if CH1_EXISTS_FUNC_MOCO else [],
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_func}/moco/channel_2_moco_bg_mean_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_func}/moco/channel_3_moco_mean_func.nii" if CH3_EXISTS_FUNC_MOCO else [],
    run:
        try:
            preprocess.make_mean_brain(fly_directory=fly_folder_to_process_oak,
                                            meanbrain_n_frames=meanbrain_n_frames,
                                            path_to_read=input,
                                            path_to_save=output,
                                            rule_name='moco_mean_brain_rule',)
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_make_moco_mean_brain')
            utils.write_error(logfile=logfile,
                error_stack=error_stack)

rule moco_mean_brain_rule_struct:
    """
    """
    threads: snake_utils.threads_per_memory
    resources:
        mem_mb=snake_utils.mem_mb_times_input,
        runtime='10m'# should be enough
    input:
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_struct}/moco/channel_{meanbr_moco_ch_struct}_moco_struct.nii"
    output:
        str(fly_folder_to_process_oak) + "/{moco_meanbr_imaging_paths_struct}/moco/channel_{meanbr_moco_ch_struct}_moco_mean_struct.nii"
    run:
        try:
            preprocess.make_mean_brain(fly_directory=fly_folder_to_process_oak,
                                            meanbrain_n_frames=meanbrain_n_frames,
                                            path_to_read=input,
                                            path_to_save=output,
                                            rule_name='moco_mean_brain_rule',)
        except Exception as error_stack:
            logfile = utils.create_logfile(fly_folder_to_process_oak,function_name='ERROR_make_moco_mean_brain')
            utils.write_error(logfile=logfile,
                error_stack=error_stack)

rule background_subtract_func:
    # run alex's line-by-line background subtraction (on channel 2 only?)
    threads: snake_utils.threads_per_memory
    resources: mem_mb=snake_utils.mem_mb_times_input
    input:
        # run background_subtract on ch2 only?
        brain_paths_ch2=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/channel_2.nii" if CH2_EXISTS_FUNC_MOCO else [],
    output:
        bg_path_ch2 = str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_bg_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
        bg_img_before = str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_before_bg_removal.png" if CH2_EXISTS_FUNC_MOCO else [],
        bg_img_after = str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_after_bg_removal.png" if CH2_EXISTS_FUNC_MOCO else [],
        bg_img_selection = str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func}/imaging/bg/channel_2_bg_selection.tif" if CH2_EXISTS_FUNC_MOCO else [],
    shell: shell_python_command + " " + scripts_path + "/scripts/background_subtract.py "
        "--fly_directory {fly_folder_to_process_oak} "
        "--dataset_path {dataset_path} "
        "--brain_paths_ch2 {input.brain_paths_ch2} "
        "--FUNCTIONAL_CHANNELS {FUNCTIONAL_CHANNELS} "
        "--bg_path_ch2 {output.bg_path_ch2} "

rule register_series:
    # register each functional series to the first series so that rois are aligned and apply to all series (correct for drift)
    threads: snake_utils.threads_per_memory
    resources: mem_mb=snake_utils.mem_mb_times_input
    input:
        # fixed mean brain to register to (first series), structural channel
        fixed_moco_path=expand(str(fly_folder_to_process_oak) + "/{fixed_moco_imaging_path_func}/moco/channel_{moco_ch_struct}_moco_bg_func.nii" if STRUCTURAL_CHANNEL=='channel_2'
            else str(fly_folder_to_process_oak) + "/{fixed_moco_imaging_path_func}/moco/channel_{moco_ch_struct}_moco_func.nii",
            fixed_moco_imaging_path_func=fixed_path_func_reg,
            moco_ch_struct=list_of_channels_struct),
        # moving brains to apply the registration to (all series but first)
        moco_path_ch1=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func_reg}/moco/channel_1_moco_func.nii" if CH1_EXISTS_FUNC_MOCO else[],
        moco_path_ch2=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func_reg}/moco/channel_2_moco_bg_func.nii" if CH2_EXISTS_FUNC_MOCO else [],
        moco_path_ch3=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func_reg}/moco/channel_3_moco_func.nii" if CH3_EXISTS_FUNC_MOCO else [],
        
    output:
        reg_path_ch1=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func_reg}/moco/channel_1_moco_func_reg.nii" if CH1_EXISTS_FUNC_MOCO else [],
        reg_path_ch2=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func_reg}/moco/channel_2_moco_bg_func_reg.nii" if CH2_EXISTS_FUNC_MOCO else [],
        reg_path_ch3=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func_reg}/moco/channel_3_moco_func_reg.nii" if CH3_EXISTS_FUNC_MOCO else [],
        reg_par_output=str(fly_folder_to_process_oak) + "/{moco_imaging_paths_func_reg}/moco/tmats_func_reg.npy",
    
    shell: shell_python_command + " " + scripts_path + "/scripts/register_series.py "
        "--fly_directory {fly_folder_to_process_oak} "
        "--dataset_path {dataset_path} "
        "--STRUCTURAL_CHANNEL {STRUCTURAL_CHANNEL} "
        "--FUNCTIONAL_CHANNELS {FUNCTIONAL_CHANNELS} "
        "--fixed_moco_path {input.fixed_moco_path} "
        "--moco_path_ch1 {input.moco_path_ch1} "
        "--moco_path_ch2 {input.moco_path_ch2} "
        "--moco_path_ch3 {input.moco_path_ch3} "
        "--reg_path_ch1 {output.reg_path_ch1} "
        "--reg_path_ch2 {output.reg_path_ch2} "
        "--reg_path_ch3 {output.reg_path_ch3} "
        "--reg_par_output {output.reg_par_output} "        
        "--reg_transform_type {reg_transform_type} "