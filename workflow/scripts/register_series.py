"""
register each functional series to the first series
"""
import os
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
#from para_stack_reg import ParaReg


parent_path = str(pathlib.Path(pathlib.Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, parent_path)
from brainsss import moco_utils
from brainsss import utils
###
# Global variable
###
WIDTH = 120 # Todo: make a global parameter class somewhere to keep track of this variable!

def reg_slice(
    moving_path,
    fixed_path,
    functional_channel_paths,
    par_output,
    reg_settings
):
    """
    Loop doing the registration for each slice.
    This is the function that is doing the heavy lifting of the multiprocessing.
    Saves the transformed files as nii (float32).
    Saves transformation matrices as tmats_func_reg.npy (float32)
    """
    # Unpack functional paths
    if functional_channel_paths is None:
        functional_path_one = None
        functional_path_two = None
    elif len(functional_channel_paths) == 1:
        functional_path_one = functional_channel_paths[0]
        functional_one_proxy = nib.load(functional_path_one)
        functional_data_one_final = np.empty(functional_one_proxy.shape)
        functional_path_two = None
    elif len(functional_channel_paths) == 2:
        functional_path_one = functional_channel_paths[0]
        functional_one_proxy = nib.load(functional_path_one)
        functional_data_one_final = np.empty(functional_one_proxy.shape)
        functional_path_two = functional_channel_paths[1]
        functional_two_proxy = nib.load(functional_path_two)
        functional_data_two_final = np.empty(functional_two_proxy.shape)
    else:
        # Fix this, should be identical to if!
        functional_path_one = None
        functional_path_two = None

    # Load moving proxy
    moving_proxy = nib.load(moving_path) #xy(z)t

    # Determine number of slices
    if len(moving_proxy.shape)==4: # xyzt
        n_slices = moving_proxy.shape[2]
    elif len(moving_proxy.shape)==3: # xyt
        n_slices = 1
    else:
        print('Error: Could not determine number of slices in moving data.')

    # Load fixed proxy
    fixed_proxy = nib.load(fixed_path) #xy(z)t

    moving_data_final = np.empty(moving_proxy.shape)

    n_timepoints = moving_proxy.shape[-1] # source xy(z)t

    tmats_final = np.empty([3,3,n_slices]) # for rigid, 3x3 tmat for each slice and timepoint
    # for a rigid transformation, tmat is a 3x3 matrix [[cos(r) −sin(r) tx],[sin(r) cos(r) ty],[0 0 1]] where r is angle, t is translation
    
    # register moving mean data to fixed mean data

    for slice in range(n_slices):
        # Keeping track of time
        t_function_start = time.time()
        if n_slices==1: # data is single plane
            moving_data_mean = np.mean(np.asarray(moving_proxy.dataobj[:,:,:],dtype='float32'),axis=-1) #xy, source data xyt
            fixed_data_mean = np.mean(np.asarray(fixed_proxy.dataobj[:,:,:],dtype='float32'),axis=-1) #xy, source data xyt
        else: # data is volume
            moving_data_mean = np.mean(np.squeeze(np.asarray(moving_proxy.dataobj[:,:,slice,:],dtype='float32')),axis=-1) #xy, source data xyzt
            fixed_data_mean = np.mean(np.squeeze(np.asarray(fixed_proxy.dataobj[:,:,slice,:],dtype='float32')),axis=-1) #xy, source data xyzt
        #moving_data_mean = np.expand_dims(moving_data_mean, axis=0) # add dummy time axis=0, txy
        #fixed_data_mean = np.expand_dims(fixed_data_mean, axis=0) # add dummy time axis=0, txy
        print('moving_data_mean shape: ' + repr(moving_data_mean.shape))
        print('fixed_data_mean shape: ' + repr(fixed_data_mean.shape))
        
        # pr = ParaReg(reg_mode=reg_settings['reg_mode'],
        #              smooth=reg_settings['smooth'],
        #              avg_wid=reg_settings['avg_wid'],
        #              n_proc=reg_settings['n_proc'],
        #              )
        sr=StackReg(reg_settings['reg_mode'])
        tmat=sr.register(ref=fixed_data_mean, mov=moving_data_mean)
        #tmats=sr.register_stack(reference=fixed_data_mean, img=moving_data_mean)
        # apply transforms
        
        if n_slices==1: # data is single plane
            moving_data = np.asarray(moving_proxy.dataobj[:,:,:],dtype='float32') #xyt
            moving_data = np.moveaxis(moving_data, -1, 0) #rearrange moving axes to t,x,y
            moving_data = sr.transform_stack(img=moving_data,tmats=np.repeat(np.expand_dims(tmat,axis=0), n_timepoints, axis=0))
            moving_data = np.moveaxis(moving_data, 0, -1)#rearrange moving axes back to x,y,t
            moving_data_final[:,:,:] = moving_data

            if functional_path_one is not None:
                functional_data_one = np.asarray(functional_one_proxy.dataobj[:,:,:],dtype='float32') #xyt
                functional_data_one = np.moveaxis(functional_data_one, -1, 0) #rearrange moving axes to t,x,y
                functional_data_one = sr.transform_stack(img=functional_data_one,tmats=np.repeat(np.expand_dims(tmat,axis=0), n_timepoints, axis=0))
                functional_data_one = np.moveaxis(functional_data_one, 0, -1) #rearrange moving axes back to x,y,t
                functional_data_one_final[:,:,:] = functional_data_one

                if functional_path_two is not None:
                    functional_data_two = np.asarray(functional_two_proxy.dataobj[:,:,:],dtype='float32') #xyt
                    functional_data_two = np.moveaxis(functional_data_two, -1, 0) #rearrange moving axes to t,x,y
                    functional_data_two = sr.transform_stack(img=functional_data_two,tmats=np.repeat(np.expand_dims(tmat,axis=0), n_timepoints, axis=0))
                    functional_data_two = np.moveaxis(functional_data_two, 0, -1) #rearrange moving axes back to x,y,t
                    functional_data_two_final[:,:,:] = functional_data_two

        else: # data is volume
            moving_data = np.squeeze(np.asarray(moving_proxy.dataobj[:,:,slice,:],dtype='float32')) #xyt
            moving_data = np.moveaxis(moving_data, -1, 0) #rearrange moving axes to t,x,y
            #moving_data = sr.transform(mov=moving_data,tmat=tmat)
            moving_data = sr.transform_stack(img=moving_data,tmats=np.repeat(np.expand_dims(tmat,axis=0), n_timepoints, axis=0))
            moving_data = np.moveaxis(moving_data, 0, -1)#rearrange moving axes back to x,y,t
            moving_data_final[:,:,slice,:] = moving_data
 
            if functional_path_one is not None:
                functional_data_one = np.squeeze(np.asarray(functional_one_proxy.dataobj[:,:,slice,:],dtype='float32')) #xyt
                functional_data_one = np.moveaxis(functional_data_one, -1, 0) #rearrange moving axes to t,x,y
                functional_data_one = sr.transform_stack(img=functional_data_one,tmats=np.repeat(np.expand_dims(tmat,axis=0), n_timepoints, axis=0))
                functional_data_one = np.moveaxis(functional_data_one, 0, -1) #rearrange moving axes back to x,y,t
                functional_data_one_final[:,:,slice,:] = functional_data_one

                if functional_path_two is not None:
                    functional_data_two = np.squeeze(np.asarray(functional_two_proxy.dataobj[:,:,slice,:],dtype='float32')) #xyt
                    functional_data_two = np.moveaxis(functional_data_two, -1, 0) #rearrange moving axes to t,x,y
                    functional_data_two = sr.transform_stack(img=functional_data_two,tmats=np.repeat(np.expand_dims(tmat,axis=0), n_timepoints, axis=0))
                    functional_data_two = np.moveaxis(functional_data_two, 0, -1) #rearrange moving axes back to x,y,t
                    functional_data_two_final[:,:,slice,:] = functional_data_two
        
        
        # save transform params for slice to tmats_final
        tmats_final[:,:,slice] = np.asarray(tmat,dtype='float32')
            # for a rigid transformation, tmat is a 3x3 matrix [[cos(r) −sin(r) tx],[sin(r) cos(r) ty],[0 0 1]] where r is angle, t is translation

        print('series registration for ' + moving_path.as_posix()
            + 'at slice ' + str(slice) + ' took : '
            + repr(round(time.time() - t_function_start, 1))
            + 's\n')
    
    # save tmats_final to tmats.npy file
    np.save(pathlib.Path(par_output), np.squeeze(tmats_final))

    # save nifti files
    if n_slices==1:
        aff = np.eye(3)
    else:
        aff = np.eye(4)

    nib.Nifti1Image(np.squeeze(moving_data_final), aff).to_filename(moving_output_path)
    
    if functional_path_one is not None:
        nib.Nifti1Image(np.squeeze(functional_data_one_final), aff).to_filename(functional_channel_output_paths[0])
        
        if functional_path_two is not None:
            nib.Nifti1Image(np.squeeze(functional_data_two_final), aff).to_filename(functional_channel_output_paths[1])
    

if __name__ == '__main__':
    ############################
    ### Organize shell input ###
    ############################
    parser = argparse.ArgumentParser()
    parser.add_argument("--fly_directory", help="Folder of fly to save log")
    parser.add_argument("--dataset_path", nargs="?", help="Folder pointing 'preprocessed'")

    parser.add_argument("--STRUCTURAL_CHANNEL", nargs="?", help="variable with string containing the structural channel")
    parser.add_argument("--FUNCTIONAL_CHANNELS", nargs="?", help="list with strings containing the functional channel")

    parser.add_argument("--fixed_moco_path", nargs="?", help="Path to fixed moco mean file for series registration")

    parser.add_argument("--moco_path_ch1", nargs="?", help="Path to ch1 moco corrected file, if Ch1 exists")
    parser.add_argument("--moco_path_ch2", nargs="?", help="Path to ch2 moco corrected file, if Ch2 exists")
    parser.add_argument("--moco_path_ch3", nargs="?", help="Path to ch3 moco corrected file, if Ch3 exists")

    parser.add_argument("--reg_path_ch1", nargs="?", help="Path to registered ch1 moco corrected file, if Ch1 exists")
    parser.add_argument("--reg_path_ch2", nargs="?", help="Path to registered ch2 moco corrected file, if Ch2 exists")
    parser.add_argument("--reg_path_ch3", nargs="?", help="Path to registered ch3 moco corrected file, if Ch3 exists")

    parser.add_argument("--reg_par_output", nargs="?", help="Path to registration parameter output")

    # reg settings
    parser.add_argument("--reg_transform_type", nargs="?", help="Type of transformation to use for registration") # default is rigid

    args = parser.parse_args()

    print('args: ' + repr(args))

    reg_par_output = args.reg_par_output
    fixed_moco_path = args.fixed_moco_path

    #####################
    ### SETUP LOGGING ###
    #####################
    fly_directory = pathlib.Path(args.fly_directory)
    logfile = utils.create_logfile(fly_directory, function_name="register_series")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    ####################################
    ### Identify the structural channel ###
    ####################################
    # Normally, we would have one structural channel
    if args.STRUCTURAL_CHANNEL is not None:
        if 'channel_1' == args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.moco_path_ch1)
            moving_output_path = pathlib.Path(args.reg_path_ch1)
        elif 'channel_2' ==  args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.moco_path_ch2)
            moving_output_path = pathlib.Path(args.reg_path_ch2)
        elif 'channel_3' ==  args.STRUCTURAL_CHANNEL:
            moving_path = pathlib.Path(args.moco_path_ch3)
            moving_output_path = pathlib.Path(args.reg_path_ch3)
        # Convert the string represenation of a list to a list - it's either
        # ['channel_1',] or ['channel_1','channel_2'] or similar
        functional_channel_paths = []
        functional_channel_output_paths = []
        if 'channel_1' in args.FUNCTIONAL_CHANNELS and 'channel_1' not in args.STRUCTURAL_CHANNEL:
            # It possible to record i.e. anat scan with ONLY the structural
            # marker. Hence, check here if the functional channel exists.
            if args.moco_path_ch1 is not None:
                functional_channel_paths.append(pathlib.Path(args.moco_path_ch1))
                functional_channel_output_paths.append(pathlib.Path(args.reg_path_ch1))
            else:
                print('Information: Channel 1 in func channel, but this channel does not exist')
        if 'channel_2' in args.FUNCTIONAL_CHANNELS and 'channel_2' not in args.STRUCTURAL_CHANNEL:
            if args.moco_path_ch2 is not None:
                functional_channel_paths.append(pathlib.Path(args.moco_path_ch2))
                functional_channel_output_paths.append(pathlib.Path(args.reg_path_ch2))
            else:
                print('Information: Channel 2 in func channel, but this channel does not exist')
        if 'channel_3' in args.FUNCTIONAL_CHANNELS and 'channel_3' not in args.STRUCTURAL_CHANNEL:
            if args.moco_path_ch3 is not None:
                functional_channel_paths.append(pathlib.Path(args.moco_path_ch3))
                functional_channel_output_paths.append(pathlib.Path(args.reg_path_ch3))
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
    # the channel to be used as the structural channel (previously 'anatomy_channel').
    # Even if you don't have a dedicated structural marker/channel, you must define
    # one of the functional channels as the structural channel.  Moco will be done on the structural
    # channel and the resulting transforms will be applied to any other functional channels
    """

    print("moving_path" + repr(moving_path))
    print("moving_output_path" + repr(moving_output_path))


    # Put moving anatomy image into a proxy for nibabel
    moving_proxy = nib.load(moving_path)
    # Read the header to get dimensions
    brain_shape = moving_proxy.header.get_data_shape()
    n_timepoints = brain_shape[-1]
    if len(brain_shape)==4: # xyzt
        n_slices = brain_shape[2]
    elif len(brain_shape)==3: # xyt
        n_slices = 1


    # reg settings
        # unpack transform type from string
    if args.reg_transform_type == 'StackReg.RIGID_BODY':
        reg_transform_type = StackReg.RIGID_BODY
    reg_settings = {}
    reg_settings['reg_mode'] = reg_transform_type
    reg_settings['smooth'] = False
    reg_settings['avg_wid'] = 1
    reg_settings['n_proc'] = 1

    print('reg_settings: ' + repr(reg_settings))

    print('Will perform registration on a total of ' + repr(n_timepoints) + ' timepoints and ' + repr(n_slices) + ' slice(s).')

    print('Starting series registration')
    time_start = time.time()

    # DO registration
    reg_slice(
             moving_path=moving_path,
             fixed_path=fixed_moco_path,
             functional_channel_paths=functional_channel_paths,
             par_output = reg_par_output,
             reg_settings = reg_settings,
             )
    
    print('Took: ' + repr(round(time.time() - time_start,1)) + 's\n to register series')
    print('Series registration done.')