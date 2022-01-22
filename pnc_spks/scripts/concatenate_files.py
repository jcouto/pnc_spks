from ..io import *
from natsort import natsorted
from glob import glob

import numpy as np


from glob import glob
from os.path import join as pjoin
from natsort import natsorted
import numpy as np
import shutil
from pnc_spks import load_spikeglx_binary


usagemsg = '''

concat-bin [-h] SESSIONS
         --output OUTPUT_FILENAME
         [--type .ap.bin] [--no-metadata FALSE]
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
'''


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Concatenate spikeGLX binary files.')
    parser.add_argument('sessions',metavar='sessions',
                        type=str,
                        help='Sessions or search string ("207988-1_Fmr1KO_M_4mo*")',
                        nargs='+')
    parser.add_argument('--no-metadata',
                        action='store_true',
                        help='do not write the metadata',
                        default=False,required=False)
    parser.add_argument('-o','--output',
                        action='store',
                        help='output filename',
                        type=str,required=True)
    parser.add_argument('-t','--type',
                        action='store',
                        help='file type to concatenate',
                        type=str,default = '.ap.bin',required=False)

    opts = parser.parse_args()

    print(opts)
    # careful with the star in the end
    folders = opts.sessions
    file_type = opts.type
    output_file = os.path.abspath(opts.output)
    fix_metadata = not opts.no_metadata 
    
    
    files = []
    for folder in folders:
        files.extend(glob(pjoin(folder,'*','*'+file_type)))
    if len(files):
        files = nat_sorted(files)
    else:
        sys.exit()

    if not os.path.exists(os.path.dirname(output_file)):
        print('Creating {0}'.format(os.path.dirname(output_file)))
        os.makedirs(os.path.dirname(output_file))
        
    concatenate_binary_files(files, output_file, fix_metadata = fix_metadata)    
