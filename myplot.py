from matplotlib.pyplot import *
from matplotlib import rc

__all__ = ['tex_on', 'tex_off', 'Args', 'Fig1x2', 'Fig2x1', 'crosshair',
    'thick', 'labels', 'axlabels', 'markers', 'rabove', 'cabove', 'cbelow',
    'xhair', 'ebar_k']

rc('figure.subplot', bottom=.125, top=.95, right=.95)  # left=0.125
#rc('lines', linewidth=2.0) # doesn't affect frame lines
rc('font', size=14)  # default for labels (not axis labels)
rc('font', family='serif')  # default for labels (not axis labels)
rc('axes', labelsize=18)
rc('xtick.major', pad=8)
rc('xtick', labelsize=14)
rc('ytick.major', pad=8)
rc('ytick', labelsize=14)
rc('savefig', dpi=150)
rc('axes.formatter', limits=(-4,4))

# Turn TeX processing of labels on/off.
def tex_on():
    rc('text', usetex=True)
    rc('font',**{'family':'serif','serif':['Computer Modern Roman']})

def tex_off():
    rc('text', usetex=False)
    rc('font',**{'family':'serif','serif':['Times']})


# Dict-like object that allows easy appending of collections of plot
# keyargs.

class Args(dict):
    """
    A dict-like object for keyword args that supports "addition" of Args
    instances.
    """

    def __init__(self, *args, **kwds):
        """
        Create a dict-like object from the arguments and keywords.
        
        Each non-keyword argument must be a dict instance or an object
        with a __dict__ attribute.  The dictionaries from these objects
        will all be merged into a single mapping (with later args taking
        precedence when keys overlap).
        
        Keyword arguments will also be merged into the mapping.
        """
        for arg in args:
            if isinstance(arg, dict):
                self.update(arg)
            elif hasattr(arg,'__dict__'):
                self.update(arg.__dict__)
            else:
                raise TypeError("Mapping or object with a __dict__ required!")
        if kwds:
            self.update(kwds)

    def __add__(self, arg):
        result = self.copy()
        result.update(arg)
        return result

#    __radd__ = __add__

    def __radd__(self, arg):
        # Don't simply duplicate __add__, or key precedence is violated.
        result = arg.copy()
        result.update(self)
        return result


# Black error bars:
ebar_k = Args({
    'markersize' : 5,
    'markeredgewidth' : 1,
    'markerfacecolor' : 'k',  # mfc
    'markeredgecolor' : 'k',  # mec
    'ecolor' : 'k'
    })

# Thick lines:
thick = Args(linewidth=2.0)

labels = Args(fontsize=14)
axlabels = Args(fontsize=18)
markers = Args(markersize=5, markeredgewidth=1)  # markeredgecolor='k'

# Text to the right and above a location.
rabove = Args(fontsize=12, ha='left', va='bottom')
# use transform=ax.transAxes for (0,1) coords

# Text centered above or below a location.
cabove = Args(ha='center', va='bottom', fontsize=12)
cbelow = Args(ha='center', va='top', fontsize=12)

# Thin-lined crosshair:
xhair = Args({ 'color' : '0.5', 'linestyle' : ':' , 'linewidth' : '1.5'})

def crosshair(ax, x, y):
    """
    Plot a crosshair.
    
    Plot a thin-lined crosshair (vertical and horizontal lines spanning the
    axes) on the specified axes, crossing at (x,y).
    """
    ax.axvline(x, **xhair)
    ax.axhline(y, **xhair)

def axis_text(x, y, s, **kwds):
    """
    Add text to the current figure using axis coords (in [0,1] for x & y).
    """
    ax = gcf().gca()
    ax.text(x, y, s, transform=ax.transAxes, **kwds)

class Fig1x2:
    """
    A figure with 1 row of 2 axes, side-by-side.
    """

    def __init__(self, figsize=(12,4)):
        self.fig = figure(figsize=figsize)
        # Left and right axes:
        self.leftax = self.fig.add_axes([.09, .15, .385, .8])  # l, b, w, h
        self.rightax = self.fig.add_axes([.59, .15, .385, .8])
        # Use thicker frame lines.
        self.leftax.get_frame().set_lw(1.25)  # thicker frame lines
        self.rightax.get_frame().set_lw(1.25)  # thicker frame lines
        # Leave with the left axes as current.
        self.fig.sca(self.leftax)

    def left(self):
        self.fig.sca(self.leftax)
        return self.leftax

    def right(self):
        self.fig.sca(self.rightax)
        return self.rightax


# Unfortunately I originally named it col X row instead of row X col.
# For compatibility with earlier use:
Fig2x1 = Fig1x2


class FigLRAxes:
    """
    A figure with two ordinate axes (left and right) sharing a common
    abscissa axis.
    
    In matplotlib lingo, this is a two-scale plot using twinx().
    """

    def __init__(self, figsize=(8,6), l=0.175, r=0.825):
        self.fig = figure(figsize=figsize)
        # Left and right axes:
        self.leftax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=l, right=r)
        self.rightax = self.leftax.twinx()
        # Use thicker frame lines.
        self.leftax.get_frame().set_lw(1.25)  # thicker frame lines
        self.rightax.get_frame().set_lw(1.25)  # thicker frame lines
        # Leave with the left axes as current.
        self.fig.sca(self.leftax)

    def left(self):
        self.fig.sca(self.leftax)
        return self.leftax

    def right(self):
        self.fig.sca(self.rightax)
        return self.rightax


def shelves(hts, lo=0, hi=1, fmt='b-', ends=False, dy=False):
    """
    Plot piecewise-constant "shelves" with ordinates in hts, evenly spaced
    from lo to hi.  Plot the lines with the specified format.
    
    If ends==True, plot end segments to 0.
    
    If dy != False, shift all heights by dy.
    """
    n = len(hts)
    if dy:
        hts = hts + dy
    step = 1./n
    left = lo
    xvals = []
    yvals = []
    if ends:
        xvals.append(left)
        yvals.append(0.)
    for i in range(n):
        xvals.append(left)
        yvals.append(hts[i])
        xvals.append(left+step)
        yvals.append(hts[i])
        left += step
    if ends:
        xvals.append(left)
        yvals.append(0.)
    plot(xvals, yvals)

# *** Just use axhline?
def hline(y, fmt='k-', **kwds):
    lo, hi = xlim()
    plot([lo, hi], [y, y], fmt, **kwds)

# *** Just use axvline?
def vline(x, fmt='k-', **kwds):
    lo, hi = ylim()
    plot([x, x], [lo, hi], fmt, **kwds)
    ylim(lo, hi)  # sometimes the line changes the limits(?)


if __name__ == '__main__':
    class C:
        def __init__(self, x, y):
            self.x = x
            self.y = y
    c = C(1,2)
    d = dict(a='a', b='b')
    args1 = Args(c, d, k1=5, k2='hello', b='B')
    args2 = Args(key1=1, key2='two')
    e = dict(e1=100, e2=200, key1='one')
    
    def kwds(**kwds):
        for key, val in kwds.iteritems():
            print '%s = %s' % (key, val)
    
        print 'Defined objects for testing Args class.'
