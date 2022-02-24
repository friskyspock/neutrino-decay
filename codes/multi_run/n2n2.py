import os
from data_class import GenDataset

# hyperparameters
num_runs = 6
for i in range(2,num_runs):
    
    # File Operations
    f = open('n2n2.txt','w')
    f.write('import model ~/tools/madGraph/models/zprhn_leptophobic_UFO\n')
    f.write('generate p p > zp > n2 n2, (n2 > z vl~,z > j j), (n2 > w+ l-, w+ > j j)\n')
    f.write('add process p p > zp > n2 n2, (n2 > z vl~, z > j j), (n2 > w- l+, w- > j j)\n')
    f.write('output ~/neutrinoDecay/results/ppTon2n2/n2n2_1M\n')
    f.write('launch\n1\n2\n3\ndone\n')
    f.write('set nevents 100000\nset ebeam1 7000.0\nset ebeam2 7000.0\n')
    f.write('set iseed ' + str(i) + '\n')
    f.close()

    # madGraph Operations
    os.system('python3.7 ~/tools/madGraph/bin/mg5_aMC n2n2.txt')

    # Getting Dataset
    s = GenDataset()
    s.load_from_root_file('~/neutrinoDecay/results/ppTon2n2/n2n2_1M/Events/run_01/tag_1_delphes_events.root')
    print("INFO : The amount of events generated are : " + str(len(s.dataset_dict)))
    s.sanitize_data()
    s.create_dataset('/home/student03/neutrinoDecay/datasets/sanitize/n2n2_100k_' + str(i+1) + '.h5')
    #print("INFO : The amount of events that survived are " + str(len(s.dataset_dict)))
    del(s)

    os.system('~/neutrinoDecay/codes/getnsub ~/neutrinoDecay/results/ppTon2n2/n2n2_1M/Events/run_01/tag_1_delphes_events.root ~/neutrinoDecay/datasets/jetdata/n2n2_100k_' + str(i+1) + '.csv')

    # Deleting Garbage
    os.system('rm -rf ~/neutrinoDecay/results/ppTon2n2/n2n2_1M')
    os.system('rm n2n2.txt')