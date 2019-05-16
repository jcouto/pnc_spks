import numpy as np
from .utils import unpackbits

def unpack_npix3a_sync(syncdat,srate=1,onsets = True):
    '''Unpacks neuropixels phase 3a external input data
events = unpack_npix3a_sync(trigger_data_channel)    
    Inputs:
        syncdat       : trigger data channel to unpack (pass the last channel of the memory mapped file)
        srate (1)     : sampling rate of the data; to convert to time - meta['imSampRate']
        onsets (True) : if true extracts onsets, if false extracts offsets 
    Outputs
        events        : dictionary of events. the keys are the channel number, the items the sample times of the events.

    Usage:
Load and get trigger times in seconds:
    dat,meta = load_spikeglx('test3a.imec.lf.bin')
    srate = meta['imSampRate']
    events = unpack_npix3a_sync(dat[-1,:],srate);
Plot events:
    plt.figure(figsize = [10,4])
    for ichan,times in events.items():
        plt.vlines(times,ichan,ichan+.8,linewidth = 0.5)
    plt.ylabel('Sync channel number'); plt.xlabel('time (s)')
    '''
    dd = unpackbits(syncdat.flatten(),16)
    mult = 1
    if not onsets:
        mult = -1
    sync_idx = np.where(mult*np.diff(dd,axis = 0)>0)
    events = {}
    for ichan in np.unique(sync_idx[1]):
        events[ichan] = sync_idx[0][sync_idx[1] == ichan]/srate
    return events
