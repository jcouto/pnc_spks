import pandas as pd
import os
import numpy as np
from os.path import join as pjoin
from glob import glob
from .io import *

def read_phy_data(sortfolder,srate = 30000,bin_file=None):
    ''' 
    Read the spike times saved as phy format.
    Does not use template data.
    Computes mean waveforms from the binary file for cluster depth.
    TODO: Add waveform stats and unit depths.
    '''
    keys = ['spike_times',
            'spike_clusters']
    tmp = dict() 
    def read_npy(key):
        res = None
        if os.path.isfile(pjoin(sortfolder,key+'.npy')):
            res = np.load(pjoin(sortfolder,key+'.npy'))
        return res
    for k in keys:
        tmp[k] = read_npy(k)
    sp = tmp['spike_times']
    clu = tmp['spike_clusters']
    uclu = np.unique(clu)
    cgroupsfile = pjoin(sortfolder,'cluster_groups.csv')
    if not os.path.isfile(cgroupsfile):
        print('Labels are from KS.')
        cgroupsfile = pjoin(sortfolder,'cluster_KSLabel.tsv')
    cgroups = pd.read_csv(cgroupsfile,sep='\t',
                          names = ['cluster_id','label'],
                         header = 0)
    spks = [(u,
             cgroups[cgroups.cluster_id == u].label.iloc[0],
             sp[clu == u]/srate) for u in uclu]
    res = pd.DataFrame(spks,columns = ['cluster_id',
                                        'cluster_label',
                                        'ts'])
    if not bin_file is None:
        print('Reading mean waveforms from the binary file.')
        assert os.path.isfile(bin_file), 'File {0} not found.'.format(bin_file)
        mwaves = par_get_mean_waveforms(bin_file,[s[2] for s in spks])
        chmap = read_phy_channelmap(sortfolder)
        # discard unused channels
        mwaves = mwaves[:,:,np.array(chmap.ichan)]
        res['mean_waveforms'] = [m for m in mwaves]
        # peak to baseline per channel 
        ptb = [mwave[20:50,:].max(axis=0) - mwave[20:50,:].min(axis=0)
               for mwave in mwaves]
        # "active" channels
        activeidx = [np.where(np.abs((p-np.mean(p))/np.std(p))>1)[0]
                     for p in ptb]
        peakchan =  [chmap.ichan.iloc[np.argmax(np.abs(p))] for p in ptb]
        res['peak_channel'] = peakchan
        res['active_channels'] = activeidx
    return res


def read_phy_channelmap(sortfolder):
    ''' 
    
    chmap = read_phy_channelmap(sortfolder)
    
    Reads channelmap from phy (template-gui).

    Needed files: channel_map.npy, channel_positions.npy
    '''
    fname = pjoin(sortfolder,'channel_map.npy')
    assert os.path.isfile(fname),'Could not find {0}'.format(fname)
    chidx = np.load(pjoin(sortfolder,fname))
    fname = pjoin(sortfolder,'channel_positions.npy')
    assert os.path.isfile(fname),'Could not find {0}'.format(fname)
    channel_positions = np.load(pjoin(sortfolder,fname))
    chx = [x[0] for x in channel_positions]
    chy = [x[1] for x in channel_positions]
    return pd.DataFrame(zip(chidx.flatten(),chx,chy),
                        columns = ['ichan','x','y'])

# helper functions to get mean waveforms
vars = {}
def _init_spikeglx_bin(bin_file):
    vars['dat'],vars['meta'] = load_spikeglx_binary(bin_file)
    vars['srate'] = vars['meta']['imSampRate']
def _work_mwaves(ts,chmap = None):
    if chmap is None:
        chmap = np.arange(vars['dat'].shape[1])
    res = get_random_waveforms(
        vars['dat'],
        chmap,
        (ts*vars['srate']).astype(int),
        nwaves=100,
        npre=30,
        npost=30,
        dofilter=True)
    return np.median(res,axis=0)

def par_get_mean_waveforms(bin_file,ts):
    '''
    mwaves = par_get_mean_waveforms(bin_file,ts)
    Gets the mean waveforms from a binary file given the timestamps.
    Usage:
    
    '''
    from multiprocessing import Pool,cpu_count
    with Pool(processes=cpu_count(),
              initializer=_init_spikeglx_bin,
              initargs=([bin_file])) as pool:
        res = pool.map(_work_mwaves,ts)
    return np.stack(res)



# Old/deprecated... don't look at this.

def get_cluster_groups(sortfolder,verbose = False,create=False):
    df = None
    sortfolder = os.path.abspath(sortfolder)
    clustergroups = glob(pjoin(sortfolder,'cluster_groups.csv'))
    KSlabel = glob(pjoin(sortfolder,'cluster_KSLabel.tsv'))
    if len(clustergroups):
        df = pd.read_csv(clustergroups[0],
                         header=0,
                         index_col=None,
                         sep='\t')
    elif len(KSlabel):
        df = pd.read_csv(KSlabel[0],
                         header=0,
                         index_col=None,
                         sep='\t')
        if verbose:
            print('Loaded (KS auto) '+ KSlabel[0])
    
    if df is None and create:
        clus = None
        for root, dirs, filenames in os.walk(sortfolder):
            for filename in filenames:
                if not '/.' in root:
                    if 'spike_templates.npy' in filename:
                        folder = root
                        clus = np.load(pjoin(folder,filename))
                    if 'spike_clusters.npy' in filename:
                        folder = root
                        print(folder)
                        clus = np.load(pjoin(folder,filename))
        assert not clus is None, 'Could not find clu file.'
        unit_id = np.unique(clus)
        unitscladic = {'cluster_id':unit_id,
              'group':np.array(['unsorted' for u in unit_id])}
        df = pd.DataFrame.from_dict(unitscladic)
        df.to_csv(folder+'/cluster_groups.csv',
                                 index=False,
                                 sep='\t')
        print('cluster_groups.csv created in ' + folder)
    return df

def save_cluster_groups(df,sortfolder):
    '''
    Saves the cluster groups, backs-up previous.
    '''
    olddf = get_cluster_groups(sortfolder,verbose = False)
    folder = None
    sortfolder = os.path.abspath(sortfolder)
    for root, dirs, filenames in os.walk(sortfolder):
        for filename in filenames:
            if not '/.' in root:
                if 'cluster_groups.csv' in filename:
                    folder = root

                    df.to_csv(folder+'/cluster_groups.csv',
                                    index=False,
                                    sep='\t')
                    olddf.to_csv(folder+'/.cluster_groups.bak',
                                index=False,
                                sep='\t')
                    print('Saved cluster_groups in:' + folder)
                    break
    return folder

def get_units_ts(sortfolder,selection='all'):
    '''
    Get timestamps from a sort folder.

    spks,unit_ids = get_units_ts(sortfolder,selection='good')

    Joao Couto - March 2016
    '''
    for root, dirs, filenames in os.walk(sortfolder):
        if not '/.' in root:
            for filename in filenames:
                if 'spike_times.npy' in filename:
                    ts = np.load(pjoin(root,'spike_times.npy'))
                if 'spike_templates.npy' in filename:
                    clus = np.load(pjoin(root,'spike_templates.npy'))
    # Use spike clusters if existing.
    for root, dirs, filenames in os.walk(sortfolder):
        if not '/.' in root:
            for filename in filenames:
                if 'spike_clusters.npy' in filename:
                    clus = np.load(pjoin(root,'spike_clusters.npy'))
    if selection == 'all':
        unit_ids = np.unique(clus)
    else:
        cluinfo = get_cluster_groups(sortfolder)
        unit_ids = np.array(cluinfo[cluinfo.group == selection].cluster_id)
    spks = [ts[clus==i] for i in unit_ids]
    return spks,unit_ids


def get_random_waveforms(data, datchannels, timestamps, nwaves = 100,
                         srate = 30000, npre = 15, npost=25,dofilter = True):
    '''Gets waveforms sampled randomly from a set of timestamps.'''
    spks2extract = np.random.choice(timestamps,
                          np.clip(nwaves,
                                  1,len(timestamps)),
                          replace=False)
    indexes = np.arange(-npre,npost,dtype=np.int32)
    waveforms = np.zeros((len(spks2extract),npre+npost,len(datchannels)),
                         dtype=np.float32)
    for i,s in enumerate(np.sort(spks2extract)):
        waveforms[i,:,:] = np.take(data[indexes+s,:].astype(np.float32),datchannels,axis=1)
    if dofilter:
        from scipy import signal
        b,a = signal.butter(3,(500 / (srate / 2.), 5000 / (srate / 2.)),'pass')
        waveforms = signal.filtfilt(b,a,waveforms,axis = 1)
    return waveforms


    waves = []
    from tqdm import tqdm
    for i,sp in tqdm(enumerate(spks)):
        waves.append(np.mean(get_random_waveforms(data, datchannels,
                                          sp.flatten(),
                                          nwaves = nwaves,
                                          npre = npre,
                                          npost=npost),axis=0))
    return np.array(waves)

