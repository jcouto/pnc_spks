import pandas as pd
import os
import numpy as np

def get_cluster_groups(sortfolder,verbose = False):
    df = None
    sortfolder = os.path.abspath(sortfolder)
    for root, dirs, filenames in os.walk(sortfolder):
        for filename in filenames:
            if not '/.' in root:
                if 'cluster_groups.csv' in filename:
                    folder = root
                    df = pd.DataFrame.from_csv(pjoin(folder,'cluster_groups.csv'),
                                           header=0,
                                           index_col=None,
                                           sep='\t')
                    if verbose:
                        print('Loaded '+ pjoin(folder,'cluster_groups.csv'))

    if df is None:
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
                         npre = 15,npost=25,dofilter = True):
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
        b,a = signal.butter(3,(500 / (25000. / 2.), 5000 / (25000 / 2.)),'pass')
        waveforms = signal.filtfilt(b,a,waveforms,axis = 1)
    return waveforms

def get_mean_waveforms(data,datchannels,spks,nwaves = 100,npre = 15,npost=25):
    waves = []
    from tqdm import tqdm
    for i,sp in tqdm(enumerate(spks)):
        waves.append(np.mean(get_random_waveforms(data, datchannels,
                                          sp.flatten(),
                                          nwaves = nwaves,
                                          npre = npre,
                                          npost=npost),axis=0))
    return np.array(waves)

