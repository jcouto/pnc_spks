import pylab as plt

def plot_raster(spks,offset=0.0,height=1.0,colors='black',ax = None):
    ''' Plot a raster from sets of spiketrains.
            - "spks" is a list of spiketrains
            - "colors" can be an array of colors
        Joao Couto - January 2016
    '''
    if ax is None:
        ax = plt.gca()
    nspks = len(spks)
    if type(colors) is str:
        colors = [colors]*nspks
    for i,(sp,cc) in enumerate(zip(spks,colors)):
        ax.vlines(sp,offset+(i*height),
                  offset+((i+1)*height),
                  colors=cc)
