'''
This code creates a text file that contains code for MadGraph.
MadGraph is then run and the data is extracted and stored.
The following is implemented :
    - Combinatorix
    - Generation Level Cuts
    - Matching
'''

from distutils.util import strtobool
import os
import subprocess
import sys

# Constants
HOME_DIR = '/home/blizzard/Research/'
CODE_DIR = 'codes/bin/'
DELPHES_FILE = '/Events/run_01/tag_1_delphes_events.root'
MADGRAPH_DIR = '/home/blizzard/Research/tools/madGraph/'
DATASET_DIR = 'datasets/csvdata/'
TXT_DIR = '/home/blizzard/Research/scratch/'

# Specific Parameters
EVENT_NAME = 'ttbar'
SIGNAL = False

# Hyper Parameters
NUM_RUNS = 1
EVENT_NUM = int(1e4)

def proc_to_gen(proc, signal = 'n2n2'):
    ret_val = ''
    if proc == 'n2n2' and signal == 'n2n2':

        # Original Decay
        ret_val += 'generate pb pb > zp > n2 n2, '
        ret_val += '(n2 > z vl~,z > j j), (n2 > w la, w > j j)\n'

        # Combinatorics on the Decay
        ret_val += 'add process pb pb > zp > n2 n2, '
        ret_val += '(n2 > w la, w > j j), (n2 > z vl~,z > j j)\n'

    elif proc == 'ttbar' and signal == 'n2n2':

        # Original Decay
        ret_val += 'generate pb pb > t t~, '
        ret_val += '(t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > j j)\n'

        # Alternative Decay
        ret_val += 'add process pb pb > t t~, '
        ret_val += '(t > w+ b, w+ > j j), (t~ > w- b~, w- > l- vl~)\n'

        # Original Decay with One Jet
        ret_val += 'add process pb pb > t t~ j, '
        ret_val += '(t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > j j)\n'

        # Alternative Decay with One Jet
        ret_val += 'add process pb pb > t t~ j, '
        ret_val += '(t > w+ b, w+ > j j), (t~ > w- b~, w- > l- vl~)\n'

    elif proc == 'wmp' and signal == 'n2n2':
        
        # Decay
        ret_val += 'generate pb pb > w, w > la vla\n'

        # Decay with One Jet
        ret_val += 'add process pb pb > w j, w > la vla\n'

        # Decay with Two Jets
        ret_val += 'add process pb pb > w j j, w > la vla\n'

    elif proc == 'wpwm' and signal == 'n2n2':
        
        # Decay
        ret_val += 'generate pb pb > w+ w-, w+ > l+ vl, w- > j j\n'

        ret_val += 'add process pb pb > w+ w-, w+ > j j, w- > l- vl~\n' 

        # Decay with One Jet
        ret_val += 'add process pb pb > w+ w- j, w+ > l+ vl, w- > j j\n'
        
        ret_val += 'add process pb pb > w+ w- j, w+ > j j, w- > l- vl~\n' 

        # Decay with Two Jets
        ret_val += 'add process pb pb > w+ w- j j, w+ > l+ vl, w- > j j\n'

        ret_val += 'add process pb pb > w+ w- j j, w+ > j j, w- > l- vl~\n' 

    elif proc == 'zwpm' and signal == 'n2n2':

        # Decay
        ret_val += 'generate pb pb > z w,(w > la vla), (z > j j)\n'

        # Decay with One Jet
        ret_val += 'add process pb pb > z w j,(w > la vla), (z > j j)\n'

        # Decay with Two Jets
        ret_val += 'add process pb pb > z w j j,(w > la vla), (z > j j)\n'

    return ret_val

def gen_cuts():
    ret_val = ''

    ret_val += 'set ptl      120 \n'
    ret_val += 'set ptj1min  120 \n'
    ret_val += 'set etaj     2   \n'

    return ret_val

def jet_matching(proc):
    # Reference : https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/IntroMatching
    ret_val = ''

    # Matching for ttbar at 1 Jet
    if proc == 'ttbar':
        ret_val += 'set xqcut 20\n'
        ret_val += 'set JetMatching:qCut 30\n'
        ret_val += 'set JetMatching:nJetMax 1\n'

    # Matching for W&Z at 2 Jets
    elif proc in ['wmp','wpwm','zwpm']:
        ret_val += 'set xqcut 10\n'
        ret_val += 'set JetMatching:qCut 15\n'
        ret_val += 'set JetMatching:nJetMax 2\n'

    return ret_val

def main(proc_name,sig_flag,gen_proc = True):
    # The loop starts at 1 as default seed (0) takes a random values of seed
    for i in range(5,NUM_RUNS+5):

        # Making File for MadGraph
        f = open(TXT_DIR + proc_name + '.txt','w')

        # Default Import and Variable Definitions
        f.write('import model ' + MADGRAPH_DIR + 'models/zprhn_leptophobic_UFO\n')
        f.write('define pb = p b b~\n')
        f.write('define w = w+ w-\n')
        f.write('define la = l+ l-\n')
        f.write('define vla = vl vl~\n')

        # Generation for a Particular Channel
        f.write(proc_to_gen(proc_name))

        # General Output (Same for all the Channels)
        f.write('output ' + HOME_DIR + 'results/' + proc_name + '\n')
        f.write('launch\n1\n2\n3\ndone\n')
        f.write('set nevents ' + str(EVENT_NUM) + '\n')
        f.write('set ebeam1 7000.0\nset ebeam2 7000.0\n')
        f.write(jet_matching(proc_name))
        f.write(gen_cuts())
        f.write('set iseed ' + str(i) + '\n')

        # True only for Background
        if not sig_flag:
            f.write('set cut_decays True\n')
        
        # Closing the file
        f.close()

        # Flag checking if process is to be generated
        if gen_proc:
            # madGraph Operations
            p = subprocess.Popen([MADGRAPH_DIR + 'bin/mg5_aMC',
                                 TXT_DIR + proc_name + '.txt'])
            p.wait()

            p = subprocess.Popen([HOME_DIR + CODE_DIR + 'dtset', 
                                  HOME_DIR + 'results/' + proc_name + DELPHES_FILE,
                                  HOME_DIR + DATASET_DIR + proc_name + str(i) + '.csv'])
            
            p.wait()
            ## Deleting Garbage
            #os.system('rm -rf ' + HOME_DIR + 'results/' + proc_name)
            os.system('rm ' + TXT_DIR + proc_name + '.txt')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        proc = EVENT_NAME
        sig = SIGNAL
    
    elif len(sys.argv) == 3:
        proc = sys.argv[1]
        sig = strtobool(sys.argv[2])

    else:
        raise ValueError('Incorrect arguments. Format : python program.py proc_name is_signal')

    main(gen_proc=True,proc_name = proc,sig_flag = sig)
