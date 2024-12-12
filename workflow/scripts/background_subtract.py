import os
import glob
import numpy as np
import nibabel as nib
from skimage import io
import matplotlib.pyplot as plt
from scipy import signal
import argparse
import sys
import pathlib

# To import files (or 'modules') from the brainsss folder, define path to scripts!
# path of workflow i.e. /Users/dtadres/snake_brainsss/workflow
scripts_path = pathlib.Path(
    __file__
).parent.resolve()
print("scripts path: " + repr(scripts_path))
workflow_path = pathlib.Path(scripts_path).parent
sys.path.insert(0, os.path.join(workflow_path))
print("sys path: " + repr(sys.path))

from brainsss import utils

def BgRemover3D(args, half_wid=20):
    
    # LOGGING
    ####
    logfile = utils.create_logfile(args.fly_directory, function_name="background_subtract_func")
    printlog = getattr(utils.Printlog(logfile=logfile), "print_to_log")
    #utils.print_function_start(logfile, rule_name)

    path = args.brain_paths_ch2
    #img shoud have dimension x, y, z, t here, x is along the line scan direction
    img = np.asarray(nib.load(path).get_data().squeeze(), dtype='float32')
    dir = os.path.dirname(path)
    file_head = path.split('.')[0].split('/')[-1]
    
    ### save before
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
    img = np.moveaxis(img, (0,1,2,3), (3,1,2,0))
    out = np.zeros_like(img)
    for ind_y in range(img.shape[1]):
        for ind_z in range(img.shape[2]):
            patch = img[:, ind_y, ind_z, :]
            bg_patch = img[:, ind_y, ind_z, bg_ind[ind_z][ind_y][0]:bg_ind[ind_z][ind_y][1]]
            bg = bg_patch.mean(axis=-1)
            patch = patch-bg[None].T
            out[:, ind_y, ind_z, :] = patch
    out = np.moveaxis(out, (0,1,2,3), (3,1,2,0))
    
    ### show spectrum ???
    # half_wid_test = 5
    # half_y_test = 15 
    # kernel2d = np.ones((half_wid_test*2, half_y_test*2))/(4*half_wid_test*half_y_test)
    # conv_template = signal.convolve2d(img.mean(-1).mean(-1), kernel2d, boundary='symm', mode='valid')
    # test_x, test_y = np.unravel_index(np.argmin(conv_template), conv_template.shape)
    # test_x += half_wid_test
    # test_y += half_y_test
    # test_patch = self.img[test_x-half_wid:test_x+half_wid, test_y-half_y:test_y+half_y, :, :]
    # test = test_patch.mean(axis=(0,1))
    # test = test.flatten(order='F') 
    # test = (test-test.mean())/test.std()
    # f, Pxx_den = signal.periodogram(test, fs)
    # plt.semilogy(f, Pxx_den)
    # plt.ylim([1e-7, 1000])

    #def show_spectrum(self, fs=180):


    ### save out

    assert img.shape == out.shape
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

    br.draw_bg()
    br.show_bg()
    br.remove_bg()
    br.show_spectrum(fs=180)
    br.save_out()
