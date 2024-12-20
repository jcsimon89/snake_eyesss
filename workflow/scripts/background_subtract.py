import os
import glob
import numpy as np
import nibabel as nib
from skimage import io
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy import signal
import argparse
import sys
import pathlib
import traceback

# To import files (or 'modules') from the brainsss folder, define path to scripts!
# path of workflow i.e. /Users/dtadres/snake_brainsss/workflow
scripts_path = pathlib.Path(
    __file__
).parent.resolve()
workflow_path = pathlib.Path(scripts_path).parent
sys.path.insert(0, os.path.join(workflow_path))

from brainsss import utils

def BgRemover3D(args, half_wid=20):
    
    # LOGGING
    ####
    logfile = utils.create_logfile(pathlib.Path(args.fly_directory), function_name="background_subtract_func")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, rule_name)

    path = args.brain_paths_ch2
    
    #img shoud have dimension x, y, z, t here, x is along the line scan direction
    # Doesn't load anything, just points to a given location
    img_proxy = nib.load(path)
    # Load data, it's float32 at this point
    img = np.asarray(img_proxy.dataobj, dtype='float32')
    dir = os.path.dirname(path)
    file_head = path.split('.')[0].split('/')[-1]
    
    # save before fig
    half_wid = 5
    half_y = 15 
    fs = 180
    kernel2d = np.ones((half_wid*2, half_y*2))/(4*half_wid*half_y)
    conv_template = signal.convolve2d(img.mean(-1).mean(-1), kernel2d, boundary='symm', mode='valid')
    test_x, test_y = np.unravel_index(np.argmin(conv_template), conv_template.shape)
    test_x += half_wid
    test_y += half_y
    
    test_patch = img[test_x-half_wid:test_x+half_wid, test_y-half_y:test_y+half_y, :, :]
    test = test_patch.mean(axis=(0,1))
    test = test.flatten(order='F') 
    test = (test-test.mean())/test.std()
    f, Pxx_den = signal.periodogram(test, fs)
    plt.semilogy(f, Pxx_den)
    plt.ylim([1e-7, 1000])
    plt.savefig(os.path.join(dir, file_head +'_before_removal.png'))
    plt.close()

    ### draw bg
    wid = 2*half_wid
    kernel = np.ones(wid)/wid
    template = np.mean(img, axis=-1)
    template = np.moveaxis(template, (0, 2), (2, 0))
    bg_ind = []
    for patch in template:
        bg_ind_tmp = []
        for line in patch:
            tmp = np.convolve(line, kernel, 'valid')
            bg_center = np.argmin(tmp) + half_wid
            bg_ind_tmp.append([bg_center-half_wid, bg_center+half_wid])
        bg_ind.append(bg_ind_tmp)
    
    ### show bg
    show_bg = np.mean(img, axis=-1)
    mv = np.round(np.max(show_bg))
    show_bg = np.moveaxis(show_bg, (0, 2), (2, 0))
    for i in range(show_bg.shape[0]):
        for j in range(show_bg.shape[1]):
            show_bg[i, j, bg_ind[i][j][0]:bg_ind[i][j][1]] = mv
    
    selection_save_name = os.path.join(dir, file_head + '_bg_selection.tif')
    io.imsave(selection_save_name, np.round(show_bg).astype('int16'))
    
    ### remove bg
    img_temp = np.moveaxis(img, (0,1,2,3), (3,1,2,0))
    out = np.zeros_like(img_temp)
    for ind_y in range(img_temp.shape[1]):
        for ind_z in range(img_temp.shape[2]):
            patch = img_temp[:, ind_y, ind_z, :]
            bg_patch = img_temp[:, ind_y, ind_z, bg_ind[ind_z][ind_y][0]:bg_ind[ind_z][ind_y][1]]
            bg = bg_patch.mean(axis=-1)
            patch = patch-bg[None].T
            out[:, ind_y, ind_z, :] = patch
    out = np.moveaxis(out, (0,1,2,3), (3,1,2,0))

    ### save after
    test_patch = out[test_x-half_wid:test_x+half_wid, test_y-half_y:test_y+half_y, :, :]
    test = test_patch.mean(axis=(0,1))
    test = test.flatten(order='F') 
    test = (test-test.mean())/test.std()
    f, Pxx_den = signal.periodogram(test, fs)
    plt.semilogy(f, Pxx_den)
    plt.ylim([1e-7, 1000])
    plt.savefig(os.path.join(dir, file_head +'_after_removal.png'))
    plt.close()

    ### save out
    try:
        assert img.shape == out.shape
    except Exception as e:
                printlog(str(e))
                printlog(str(e))
                printlog(traceback.format_exc())
    
    save_name = args.moco_path_ch2
    nib.Nifti1Image(out.astype('float32'), np.eye(4)).to_filename(save_name)

if __name__ == '__main__':
    # parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--fly_directory", help="Folder of fly to save log")
    parser.add_argument("--dataset_path", nargs="?", help="Folder pointing 'preprocessed'")
    parser.add_argument("--brain_paths_ch2", nargs="?", help="Path to ch2 file, if it exists")
    parser.add_argument("--FUNCTIONAL_CHANNELS", nargs="?", help="list with strings containing the functional channel")
    parser.add_argument("--moco_path_ch2", nargs="?", help="Path to ch2 moco corrected file, if Ch2 exists")
    args = parser.parse_args()

    BgRemover3D(args, half_wid=20) #original setting: 20

