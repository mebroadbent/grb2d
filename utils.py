"""
Utility functions for image modeling and simulation.

Created 2013-05-11 by Tom Loredo (drawing from other modules)
"""

import os, pickle
import os.path as path

import matplotlib as mpl
import matplotlib.pyplot as pl
from numpy import random, linspace, diag, sqrt


__all__ = ['save_rng', 'restore_rng', 'pixel_image']


# Make images pixelized with origin at lower left.
mpl.rc('image', origin='lower', interpolation='nearest', aspect='equal', cmap='jet')


MT_id = 'MT19937'  # NumPy's RNG as of 1.0.3


def save_rng(fname='.numpy-rng-state'):
    """
    Save the state of NumPy's RNG to a file in the CWD ('.numpy-rng-state' by
    default).  Backup the previous two saved states if present:  If the RNG
    state file exists (from a previous save), rename the previous one with a
    '.1' suffix.  If a '.1' file exists, rename it with a '.2' suffix.
    """
    state = random.get_state()
    if state[0] == MT_id:
        id, state = state[0], (state[1], state[2]) # (ID, (key, pos)) for MT
    else:
        raise RuntimeError, 'numpy.random using unrecognized RNG type!'
    if path.exists(fname):
        fname1 = fname + '.1'
        if path.exists(fname1):
            fname2 = fname + '.2'
            os.rename(fname1, fname2)
        os.rename(fname, fname1)
    ofile = open(fname, 'w')
    pickle.dump((id, state), ofile)
    ofile.close()
    
def restore_rng(fname='.numpy-rng-state', notify=True):
    """
    Restore the state of NumPy's RNG from the contents of a file in the CWD
    if the file exists; otherwise use (and save) the default initialization.
    The default file name is '.numpy-rng-state'.
    """
    if os.access(fname, os.R_OK | os.W_OK):
        rng_file = open(fname, 'r')
        id, state = pickle.load(rng_file)
        rng_file.close()
        if id == MT_id:
            # Note key is numpy,uint32 -> need repr() to see value.
            if notify:
                print 'Recovered RNG state:  %s [%s %s ...] %i' %\
                    (id, repr(state[0][0]), repr(state[0][1]), state[1])
            random.set_state((id, state[0], state[1]))
        else:
            raise ValueError, 'Invalid ID for RNG in %s!' % fname
    else:
        print 'No accessible RNG status file; using (and saving) default initialization.'
        save_rng(fname)

def pixel_image(data, slices=False, xscale=None, units='arc sec',
               log_rng=None, norm=False, title=None, axes=None):
    """
    Plot pixelized image data, with a formatter showing row, col, and
    intensity.

    If `log_rng` is specified, the base-10 log of the image data is
    used for the plot, with the color map spanning from the max value
    to the max - `log_rng`.

    If `norm` is true, the data are normalized to their peak value
    before plotting.  This only affects the colorbar scale.
    """
    if slices:
        if axes:
            raise ValueError('Cannot use existing axes for slices!')
        fig = pl.figure(figsize=(15,6))
        ax = fig.add_subplot(121)
    else:
        if axes:
            ax = axes
        else:
            fig = pl.figure(figsize=(10,8))
            ax = fig.add_subplot(111)

    if norm:
        ndata = data/data.max()
    else:
        ndata = data

    if log_rng:
        dmin = ndata.min()
        dmax = ndata.max()
        if dmin <= 0.:
            ndata = ndata + abs(dmin) + 1.e-20*dmax
        ndata = log10(ndata)
        dmax = ndata.max()
        im = ax.imshow(ndata, vmax=dmax, vmin=dmax-log_rng,
                       interpolation='nearest')
    else:
        im = ax.imshow(ndata, interpolation='nearest')
    ax.get_figure().colorbar(im)
    if title:
        pl.title(title)

    # This is based on mpl's example image_zcoord.py.
    numrows, numcols = ndata.shape
    def format_coord(x, y, data=ndata):
        col = int(x+0.5)
        row = int(y+0.5)
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = data[row,col]
            return 'x=%i, y=%i, z=%1.4f'%(col, row, z)
        else:
            return 'x=%1.4f, y=%1.4f'%(x, y)

    ax.format_coord = format_coord

    if slices:
        ax = fig.add_subplot(122)
        # *** Note the handling of the image scale assumes the same
        # units in x and y.
        if xscale is None:
            xscale = 1.  # just px units
            xlbl = 'offset (px)'
        else:
            xlbl = 'offset (%s)' % units
        n, m = ndata.shape
        if n != m:
            raise ValueError('Slice plot only for square images!')
        center = n/2
        # Cross-sections:
        offsets = xscale*(linspace(0, n-1, n) - center)
        # Diagonal:
        ax.plot(sqrt(2)*offsets, diag(ndata), 'b-', lw=2, label='diagonal')
        # x and y:
        ax.plot(offsets, ndata[center][:], 'g--', lw=2, label='x cut')
        ax.plot(offsets, ndata[:][center], 'r.', lw=2, label='y cut', alpha=.5)
        ax.legend()
        ax.set_xlabel(xlbl)
        y_l, y_u = ax.get_ylim()
        if log_rng is None:
            ax.set_ylim(0, y_u)
        else:
            ax.set_ylim(y_u-log_rng-1, y_u)
