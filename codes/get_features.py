import ROOT
from ROOT import TLorentzVector
from math import sqrt
import h5py
import numpy as np
import awkward as ak
import pandas as pd
import json

class GetFeatures():
    def __init__(self, procname = 'hh', file_location = '~/Tests/bbWW/', num_events = '100k'):
        self.procname = procname
        self.filepath = file_location + procname + '_' + num_events + '.h5'
        self.part = []
        self.eta = []
        self.phi = []
        self.angle = []

    def get_data(self):
        '''
        Accessor function that gives back lists from hdf5 file. Due to the various conversions from 
        List --> Awkward --> HDF5 file, this function is made as to make the conversion back to list 
        simple.

        Parameters :
        ------------
        None

        Returns :
        ---------
        None
        '''

        print("INFO : Started Getting Data From File")

        hf = h5py.File(self.filepath,'r')
        partArray = hf.get("ParticleArray")
        azimuthalArray = hf.get("AzimuthalAngle")
        etaArray = hf.get("EtaAngle")
        phiArray = hf.get("PhiAngle")

        reconstitutedPartArray = ak.from_buffers(
            ak.forms.Form.fromjson(partArray.attrs["form"]),
            json.loads(partArray.attrs["length"]),
            {k: np.asarray(v) for k, v in partArray.items()},
        )

        reconstitutedAzAngle = ak.from_buffers(
            ak.forms.Form.fromjson(azimuthalArray.attrs["form"]),
            json.loads(azimuthalArray.attrs["length"]),
            {k: np.asarray(v) for k, v in azimuthalArray.items()},
        )

        reconstitutedEtaAngle = ak.from_buffers(
            ak.forms.Form.fromjson(etaArray.attrs["form"]),
            json.loads(etaArray.attrs["length"]),
            {k: np.asarray(v) for k, v in etaArray.items()},
        )

        reconstitutedPhiAngle = ak.from_buffers(
            ak.forms.Form.fromjson(phiArray.attrs["form"]),
            json.loads(phiArray.attrs["length"]),
            {k: np.asarray(v) for k, v in phiArray.items()},
        )

        self.part = ak.to_list(reconstitutedPartArray)
        self.angle = ak.to_list(reconstitutedAzAngle)
        self.eta = ak.to_list(reconstitutedEtaAngle)
        self.phi = ak.to_list(reconstitutedPhiAngle)

        print("INFO : Done Getting Data from File")

    def get_mbb(self):
        mbbTemp = []

        for i in range(len(self.part)):
            ptJetsb = []        # List containing b-tagged Jets
            flag = True         # Flag that True -> Survive the Cut,
                                # False -> Does not Survive the Cut

            for momenta in self.part[i]:
                # Getting a List of Jets to  
                # compute cuts on them
                
                # B-Tagged Jet
                if momenta[2] == 1:
                    ptJetsb.append(momenta[4])
            
            # Record invalid if the two leading momenta are not of b-tagged Jets
            ptJetsb.sort(reverse=True)
            
            # Getting the Index of Leading Jets
            jetIndex = []

            for momenta in self.part[i]:
                for j in range(2):
                    if momenta[2] == 1 and momenta[4] == ptJetsb[j]:
                        jetIndex.append(self.part[i].index(momenta))

            Jet1 = TLorentzVector()
            Jet2 = TLorentzVector()
            
            Jet1.SetPtEtaPhiM(self.part[i][jetIndex[0]][4],
                              self.eta[i][jetIndex[0]],
                              self.phi[i][jetIndex[0]],
                              self.part[i][jetIndex[0]][-1])
            
            Jet2.SetPtEtaPhiM(self.part[i][jetIndex[1]][4],
                              self.eta[i][jetIndex[1]],
                              self.phi[i][jetIndex[1]],
                              self.part[i][jetIndex[1]][-1])
            
            mbbTemp.append((Jet1 + Jet2).M())
        return mbbTemp

    def get_mll(self):
        mllTemp = []
    
        for i in range(len(self.part)):
            leptonIndex = []
            for momenta in self.part[i]:
                if momenta[1] != 0:
                        leptonIndex.append(self.part[i].index(momenta))
                
            for j in range(len(leptonIndex)):
                for k in range(j+1,len(leptonIndex)):
                        
                    Lepton1 = TLorentzVector()
                    Lepton2 = TLorentzVector()
    
                    Lepton1.SetPtEtaPhiM(self.part[i][leptonIndex[j]][4],
                                      self.eta[i][leptonIndex[j]],
                                      self.phi[i][leptonIndex[j]],
                                      self.part[i][leptonIndex[j]][-1])
                                      
                    Lepton2.SetPtEtaPhiM(self.part[i][leptonIndex[k]][4],
                                      self.eta[i][leptonIndex[k]],
                                      self.phi[i][leptonIndex[k]],
                                      self.part[i][leptonIndex[k]][-1])
                    
                    mllTemp.append((Lepton1 + Lepton2).M())
        
        return mllTemp

    def get_rbb(self):
        rbbTemp = []

        for i in range(len(self.part)):
            ptJetsb = []        # List containing b-tagged Jets
            flag = True         # Flag that True -> Survive the Cut,
                                # False -> Does not Survive the Cut

            for momenta in self.part[i]:
                # Getting a List of Jets to  
                # compute cuts on them
                
                # B-Tagged Jet
                if momenta[2] == 1:
                    ptJetsb.append(momenta[4])
            
            # Record invalid if the two leading momenta are not of b-tagged Jets
            ptJetsb.sort(reverse=True)

            # Getting the Index of Leading Jets
            jetIndex = []

            for momenta in self.part[i]:
                for j in range(2):
                    if momenta[2] == 1 and momenta[4] == ptJetsb[j]:
                        jetIndex.append(self.part[i].index(momenta))

            rbbTemp.append(self.angle[i][jetIndex[0]][jetIndex[1]])

        return rbbTemp

    def get_rll(self):
        rllTemp = []
    
        for i in range(len(self.part)):
            leptonIndex = []
            for momenta in self.part[i]:
                if momenta[1] != 0:
                        leptonIndex.append(self.part[i].index(momenta))
                
            for j in range(len(leptonIndex)):
                for k in range(j+1,len(leptonIndex)):
                    rllTemp.append(self.angle[i][leptonIndex[j]][leptonIndex[k]])
                          
        return rllTemp

    def get_ptj(self):
        ''' Gets the transverse momenta of the leading jet'''
        ptjTemp = []

        for i in range(len(self.part)):
            ptAllJets = []      # List of all the Jets

            for momenta in self.part[i]:
                if momenta[2] != 0:
                    ptAllJets.append(momenta[4])
            
            ptAllJets.sort(reverse=True)
            ptjTemp.append(ptAllJets[0])

        return ptjTemp

    def get_ptl(self):
        ptlTemp = []
    
        for i in range(len(self.part)):
            leptonpt = []
            for momenta in self.part[i]:
                if momenta[1] != 0:
                        leptonpt.append(momenta[4])
        
            leptonpt.sort(reverse=True)    
            ptlTemp.append(leptonpt[0])
                          
        return ptlTemp

    def get_El(self):
        ElTemp = []
    
        for i in range(len(self.part)):
            leptonpt = []
            for momenta in self.part[i]:
                if momenta[1] != 0:
                        leptonpt.append(momenta[5])
        
            leptonpt.sort(reverse=True)    
            ElTemp.append(leptonpt[0])
                          
        return ElTemp

    def get_MET(self):
        METTemp = []
    
        for i in range(len(self.part)):
            leptonpt = []
            for momenta in self.part[i]:
                if momenta[3] != 0:
                        leptonpt.append(momenta[4])
        
            METTemp.append(leptonpt[0])
                          
        return METTemp

    def get_jetE(self):
        EjTemp = []
    
        for i in range(len(self.part)):
            jetE = []
            for momenta in self.part[i]:
                if momenta[2] != 0:
                        jetE.append(momenta[5])
        
            jetE.sort(reverse=True)    
            EjTemp.append(jetE[0])
                          
        return EjTemp

    def get_jetM(self):
        jetMTemp = []
    
        for i in range(len(self.part)):
            jetE = []
            for momenta in self.part[i]:
                if momenta[2] != 0:
                        jetE.append(momenta[6])
        
            jetE.sort(reverse=True)    
            jetMTemp.append(jetE[0])
                          
        return jetMTemp

    def save_csv(self,output_file,i):
        ptj = self.get_ptj()
        ptl = self.get_ptl()
        El = self.get_El()
        jetE = self.get_jetE()
        met = self.get_MET()
        jetM = self.get_jetM()

        if self.procname != 'n2n2':
            tag = np.zeros_like(ptl)    
        else:
            tag = np.ones_like(ptl)

        #Nsubjettiness variables
        pb = pd.read_csv('/home/student03/neutrinoDecay/datasets/jetdata/' + self.procname + '_100k_' + str(i) + '.csv')
        
        features = np.vstack((ptj,ptl,jetE,jetM,met,El,list(pb['Number N21']),list(pb['Number N32']),tag))

        df = pd.DataFrame(features.T,columns=['ptj','ptl','jetE','jetM','MET','El','n21','n32','tag'])
        df.to_csv(output_file)

if __name__ == "__main__":
    process = ['n2n2','ttbar','wmp','wpwm','zzjj','zwpm']
    for proc in process:
        for i in range(2,7):
            f1 = GetFeatures(proc,file_location='/home/student03/neutrinoDecay/datasets/sanitize/',num_events='100k_' + str(i))
            f1.get_data()
            f1.save_csv('~/neutrinoDecay/datasets/csv/' + proc + '_100k_' + str(i) + '.csv',i)
