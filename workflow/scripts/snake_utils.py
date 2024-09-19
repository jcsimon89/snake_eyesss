# Functions to be called from snakefiles
import numpy as np
import nibabel as nib
import pathlib



def ch_exists(channel, CHANNELS):
    """
    Check if a given channel exists in global variables STRUCTURAL_CHANNEL and FUNCTIONAL_CHANNELS
    :param channel:
    :return:
    """
    if 'channel_' + str(channel) in CHANNELS:
        ch_exists = True
    else:
        ch_exists = False
    return(ch_exists)


def mem_mb_times_threads(wildcards, threads):
    """
    Returns memory in mb as 7500Mb/thread (I think we have ~8Gb/thread? to be confirmed)
    Note: wildcards is essential here!
    :param threads:
    :return:
    """
    return threads * 7500


def mem_mb_less_times_input(wildcards, input):
    """
    Returns memory in mb as 1.5*input memory size or 1.5Gb, whichever is larger
    :param wildcards:
    :param input:
    :return:
    """
    return max(input.size_mb * 1.5, 1500)
def threads_per_memory_less(wildcards, input):
    """
    It seems I now have to define the threads based on the memory requirements.
    We get 8Gb per core so we need e.g. 2 threads if we want to use 15Gb of RAM
    :param wildcards:
    :param input:
    :return:
    """
    calculated_memory = max(input.size_mb * 1.5, 10000)
    cores = int(np.ceil(calculated_memory/8000))
    return(cores)

def mem_mb_times_input(wildcards, input):
    """
    Returns memory in mb as 2.5*input memory size or 2Gb, whichever is larger
    :param wildcards:
    :param input:
    :return:
    """
    return max(input.size_mb * 2.5, 2000)

def threads_per_memory(wildcards, input):
    """
    It seems I now have to define the threads based on the memory requirements.
    We get 8Gb per core so we need e.g. 2 threads if we want to use 15Gb of RAM
    :param wildcards:
    :param input:
    :return:
    """
    calculated_memory = max(input.size_mb * 2.5, 10000)
    cores = int(np.ceil(calculated_memory/8000))
    return(cores)

def mem_mb_more_times_input(wildcards, input):
    """
    Returns memory in mb as 3.5*input memory size or 4Gb, whichever is larger,
    but not more than 256Gb which is the submission limit for sherlock.
    :param wildcards:
    :param input:
    :return:
    """
    calculated_memory = max(input.size_mb * 3.5, 4000)

    if calculated_memory < 256000:
        return(calculated_memory)
    else:
        return(256000)

def threads_per_memory_more(wildcards, input):
    """
    It seems I now have to define the threads based on the memory requirements.
    We get 8Gb per core so we need e.g. 2 threads if we want to use 15Gb of RAM
    :param wildcards:
    :param input:
    :return:
    """
    calculated_memory = max(input.size_mb * 3.5, 10000)
    cores = int(np.ceil(calculated_memory/8000))
    return(cores)

def mem_mb_much_more_times_input(wildcards, input):
    """
    Returns memory in mb as 5.5*input memory size or 1Gb, whichever is larger,
    but not more than 256Gb which is the submission limit for sherlock.
    :param wildcards:
    :param input:
    :return:
    """
    calculated_memory = max(input.size_mb * 5.5, 10000)

    if calculated_memory < 256000:
        return(calculated_memory)
    else:
        return(256000)


def threads_per_memory_much_more(wildcards, input):
    """
    It seems I now have to define the threads based on the memory requirements.
    We get 8Gb per core so we need e.g. 2 threads if we want to use 15Gb of RAM
    :param wildcards:
    :param input:
    :return:
    """
    calculated_memory = max(input.size_mb * 5.5, 10000)
    cores = int(np.ceil(calculated_memory/8000))
    return(cores)

def threads_8x_more_input(wildcards, input):
    """
    It seems I now have to define the threads based on the memory requirements.
    We get 8Gb per core so we need e.g. 2 threads if we want to use 15Gb of RAM
    :param wildcards:
    :param input:
    :return:
    """
    calculated_memory = max(input.size_mb * 8, 10000)
    cores = int(np.ceil(calculated_memory/8000))
    return(cores)

def disk_mb_times_input(wildcards, input):
    """
    Returns memory in mb as 2.5*input memory size or 1Gb, whichever is larger
    :param wildcards:
    :param input:
    :return:
    """
    return max(input.size_mb * 2.5, 1000)

def mb_for_moco_input(wildcards, input):
    """
    Note on memory requirements:
    https://sourceforge.net/p/advants/discussion/840261/thread/9a2a668a/?limit=25#aeff
    ANTs converts all input images to floating point. Also, the SyN algorithm requires the use
    of four displacement fields. Thus, in addition to the two input images, you also have an
    additional set of 4 displacement fields * 3 component (i.e., x, y, z) images / displacement
    field = 12 scalar images stored in memory. That's probably where most of your memory is going.
    The number of threads shouldn't matter.
    Also check here:
    https://sourceforge.net/p/advants/discussion/840261/thread/907e5819/
    :param wildcards:
    :param input:
    :return:
    """
    # There's a huge discrepancy of memory requirement between anat (1024x512) and func (256x128) files.
    # Can we read the file dimensions here?
    if input.brain_paths_ch1 != []:
        # Make a pointer to the nii file!
        print(input.brain_paths_ch1)
        brain_proxy = nib.load(input.brain_paths_ch1)
        # Without reading the image into memory ,just look at the header
        # to get image dimensions.
        brain_dims = brain_proxy.header.get_data_shape()
    # If channel_1 doesn't exist, try 2. All channels should have the same
    # dimension, so checking one is sufficient.
    elif input.brain_paths_ch2 !=[]:
        brain_proxy = nib.load(input.brain_paths_ch2)
        brain_dims = brain_proxy.header.get_data_shape()
    elif input.brain_paths_ch3 !=[]:
        brain_proxy = nib.load(input.brain_paths_ch2)
        brain_dims = brain_proxy.header.get_data_shape()
    # This is the size of the brain volume at one 'timepoint', meaning
    # for x,y and z but NOT t(ime)
    volume_size = brain_dims[0] * brain_dims[1] * brain_dims[2]
    # When I use a 'standard' recording with (256,128,49) I need ~ 2.5 times the memory
    # compared to the input. 256*128*49=1,605,632
    # In contrast, I was just barely able to avoid OOM error for (1024,512,100) with
    # Yandan's data when using 6.5x the memory. 1024*512*100=52,428,800
    # And I get OOM error when doing (1024,512,241). 1024*512*241=126,353,408
    ###
    # There's clearly a non-linear relationship between memory requirement and size of the
    # first three dimensions...
    if volume_size < 5e6:
        # Probably a functional dataset with something like (256,128,49)
        multiplier = 2.5
        # That's a pretty normal anat recording like (1024,512,100)
    elif volume_size < 6e7:
        multiplier = 7
        # that's a larger anat recording like (1024,512,275)
    elif volume_size < 1.5e8:
        multiplier = 16

    mem_mb = input.size_mb*multiplier

    # We need pretty much exactly 4 times the input file size OF ONE CHANNEL.
    # If we have two channels we only need ~2 times.
    if (
            input.brain_paths_ch1 != []
            and input.brain_paths_ch2 != []
            and input.brain_paths_ch3 != []
    ):
        # if all three channels are used
        #mem_mb = int((input.size_mb*6.5)/3)
        mem_mb /= 3
    elif (
        (input.brain_paths_ch1 != [] and input.brain_paths_ch2 != [])
        or (input.brain_paths_ch1 != [] and input.brain_paths_ch3 != [])
        or (input.brain_paths_ch2 != [] and input.brain_paths_ch3 != [])
    ):
        # if only two channels are in use
        #mem_mb = int((input.size_mb * 6.5)/2)
        mem_mb /= 2
    #else:
    #    # only one channel is provided:
    #    mem_mb = (input.size_mb * 6.5)
    # Sherlock only allows us to request up to 256 Gb per job
    if mem_mb > 256000:
        mem_mb = 256000

    return(max(mem_mb, 5000))

# if all three channels are used

def time_for_moco_input(wildcards, input):
    """
    Note that all benchmarks were done with 'threads=32' and 'cores=8'.
    Returns time in minutes based on input. Use this with the multiprocessing motion
    correction code.
    We needs about 5 minutes for 1 Gb of for two channels. Double it to be on the safe
    side for now

    :param wildcards: Snakemake requirement
    :param input: intput from snakefile. Needed to access filesize
    :return:
    """
    if input == "<TBD>":  # This should ONLY happen during a -np call of snakemake.
        time_in_minutes = 120
    elif (
        input.brain_paths_ch1 != []
        and input.brain_paths_ch2 != []
        and input.brain_paths_ch3 != []
    ):
        # if all three channels are used
        time_in_minutes = (input.size_mb / 1000) * 2.5  # /1000 to get Gb, then *minutes
    elif (
        (input.brain_paths_ch1 != [] and input.brain_paths_ch2 != [])
        or (input.brain_paths_ch1 != [] and input.brain_paths_ch3 != [])
        or (input.brain_paths_ch2 != [] and input.brain_paths_ch3 != [])
    ):
        # if only two channels are in use
        time_in_minutes = (input.size_mb / 1000) * 5 # /1000 to get Gb, then *minutes
    else:
        # only one channel is provided:
        time_in_minutes = (input.size_mb / 1000) * 10  # /1000 to get Gb, then *minutes

    # hours = int(np.floor(time_in_minutes / 60))
    # minutes = int(np.ceil(time_in_minutes % 60))
    # string_to_return = str(hours) + ':' + str(minutes) + ':00'

    # Define a minimum time of 10 minutes
    time_in_minutes = max(time_in_minutes, 10)

    # https: // snakemake.readthedocs.io / en / stable / snakefiles / rules.html
    # If we want minutes we just add a 'm' after the number
    string_to_return = str(time_in_minutes) + "m"
    return(string_to_return)
def OLDtime_for_moco_input(wildcards, input):
    """
    ### This was used for the chunk based NOT-multiprocessed code ####
    Returns time in minutes based on input.
    I know that at 5 minute recording should take at least 30 minutes with input size ~2Gb
    I also know that a 30 minute recording should take ~7 hours with ~11Gb input size.
    And an anatomical scan also ~25Gb should take ~16 hours
    The slowest seems to be the functional 7 hours scan but even that is taken
    We'll do 45min per anatomy channel
    Note that we can have 1, 2 or 3 input files...Assume that only one of the files is the
    anatomy channel!
    :param wildcards:
    :param input:
    :return:
    """
    if input == "<TBD>":  # This should ONLY happen during a -np call of snakemake.
        string_to_return = str(2) + "h"
    elif (
        input.brain_paths_ch1 != []
        and input.brain_paths_ch2 != []
        and input.brain_paths_ch3 != []
    ):
        # if all three channels are used
        time_in_minutes = (input.size_mb / 1000) * 15  # /1000 to get Gb, then *minutes
    elif (
        (input.brain_paths_ch1 != [] and input.brain_paths_ch2 != [])
        or (input.brain_paths_ch1 != [] and input.brain_paths_ch3 != [])
        or (input.brain_paths_ch2 != [] and input.brain_paths_ch3 != [])
    ):
        # if only two channels are in use
        time_in_minutes = (input.size_mb / 1000) * 30  # /1000 to get Gb, then *minutes
    else:
        # only one channel is provided:
        time_in_minutes = (input.size_mb / 1000) * 45  # /1000 to get Gb, then *minutes

    # hours = int(np.floor(time_in_minutes / 60))
    # minutes = int(np.ceil(time_in_minutes % 60))
    # string_to_return = str(hours) + ':' + str(minutes) + ':00'

    # https: // snakemake.readthedocs.io / en / stable / snakefiles / rules.html
    # If we want minutes we just add a 'm' after the number - TEST!!!mem_mb_more_times_input
    string_to_return = str(time_in_minutes) + "m"
    return string_to_return


'''def time_for_correlation(wildcards, input):
    """
    returns time in based on input - for a 5 min test case we need just under 5 minutes.
    That's about 4 Gb of h5 file for a single functional channel. Lets go with 10 minutes for 4Gb so 2.5 minutes
    for each Gb

    # >DOESNT WORKThat's about 10 Mb of fictrac data. So each Mb of fictrac data is ~1 minute of compute time
    :param wildcards:
    :param input
    :return: time in seconds
    """
    # return(input.fictrac_path.size_mb*60) DOESNT WORK!
    time_in_minutes = (input.size_mb/1000)*2.5
    return(time_in_minutes*60) # This turned out to give minutes...'''

