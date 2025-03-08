from functools import partial
from multiprocessing.pool import Pool
from time import time

from tqdm import tqdm
import numpy as np
import tifffile as tif
from pystackreg import StackReg 


def single_slice_reg(sr, ref, mov_slice):
    return sr.register(ref, mov_slice)

def single_slice_trans(sr, mov_slice, tmat):
    return sr.transform(mov_slice, tmat)

def smooth(x, axis=0, wid=5):
    # this is way faster than convolve
    if wid < 2:
        return x
    cumsum_vec = np.cumsum(np.insert(x, 0, 0, axis=axis), axis=axis)
    ma_vec = (cumsum_vec[wid:] - cumsum_vec[:-wid]) / wid
    y = x.copy()
    start_ind = int(np.floor((wid-1)/2))
    end_ind = wid-1-start_ind
    y[start_ind:-end_ind] = ma_vec
    return y

class ParaReg(object):
    def __init__(self, reg_mode, smooth, avg_wid, n_proc, mean_frames=None):
        super().__init__()
        self.sr = StackReg(reg_mode)
        self.smooth = smooth
        self.avg_wid = avg_wid
        self.n_proc=n_proc
        self._tmats=[]
        self.mean_frames = mean_frames
    
    def register(self, img, ref=None):
        if self.smooth:
            img = smooth(img, wid=self.avg_wid)
        if ref is None:
            ref = img[:self.mean_frames, :, :].mean(axis=0)
        register_worker = partial(single_slice_reg, self.sr, ref)
        with Pool(processes=self.n_proc) as p:
            res = p.imap(register_worker, img, 128)
            tmats = []
            for r in tqdm(res):
                tmats.append(r)
        self._tmats = tmats
    
    def transform(self, img):
        transform_worker = partial(single_slice_trans, self.sr)
        p = Pool(processes=self.n_proc)
        assert len(self._tmats) == img.shape[0]
        with Pool(processes=self.n_proc) as p:
            iters = [(b, t) for b, t in zip(img, self._tmats)]
            res = p.starmap(transform_worker, iters, 128)
            registered = []
            for r in res:
                registered.append(r)
        return np.array(registered, dtype='float32')
