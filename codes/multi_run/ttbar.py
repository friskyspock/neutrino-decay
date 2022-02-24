import os
from data_class import GenDataset

# hyperparameters
num_runs = 6
for i in range(2,num_runs):
    
    # File Operations
    f = open('ttbar.txt','w')
    f.write('import model ~/tools/madGraph/models/zprhn_leptophobic_UFO\n')
    f.write('generate p p > t t~, (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > j j)\n')
    f.write('add process p p > t t~, (t > w+ b, w+ > j j), (t~ > w- b~, w- > l- vl~)\n')
    f.write('output ~/neutrinoDecay/results/ppTon2n2/ttbar_1M\n')
    f.write('launch\n1\n2\n3\ndone\n')
    f.write('set nevents 100000\nset ebeam1 7000.0\nset ebeam2 7000.0\n')
    f.write('set iseed ' + str(i) + '\n')
    f.close()

    # madGraph Operations
    os.system('python3.7 ~/tools/madGraph/bin/mg5_aMC ttbar.txt')

    # Getting Dataset
    s = GenDataset()
    s.load_from_root_file('~/neutrinoDecay/results/ppTon2n2/ttbar_1M/Events/run_01/tag_1_delphes_events.root')
    print("INFO : The amount of events generated are : " + str(len(s.dataset_dict)))
    s.sanitize_data()
    s.create_dataset('/home/student03/neutrinoDecay/datasets/ttbar_100k_' + str(i+1) + '.h5')
    #print("INFO : The amount of events that survived are " + str(len(s.dataset_dict)))
    del(s)

    os.system('~/neutrinoDecay/codes/getnsub ~/neutrinoDecay/results/ppTon2n2/ttbar_1M/Events/run_01/tag_1_delphes_events.root ~/neutrinoDecay/datasets/jetdata/ttbar_100k_' + str(i+1) + '.csv')

    # Deleting Garbage
    os.system('rm -rf ~/neutrinoDecay/results/ppTon2n2/ttbar_1M')
    os.system('rm ttbar.txt')