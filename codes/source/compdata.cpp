//#ifdef __CLING__
//R__LOAD_LIBRARY (libDelphes);
//#endif
//R_ADD_INCLUDE_PATH(/home/student02/hep/Delphes-3.4.1)

/* 
This code is to get nsubjettiness of the jets of the root 
files and keep them in an ordered csv or a list to be 
imported in a python code. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <exception>
#include <cstdlib>
#include <cstdio>
#include <limits>


#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLorentzVector.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/internal/BasicRandom.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"

#include "classes/DelphesClasses.h"


using namespace std;
using namespace H5;
using namespace fastjet;
using namespace fastjet::contrib;

int main(int argc,char* argv[]){
    gSystem->Load("libDelphes");
	TChain chain("Delphes");

    double beta = 2.0; 

    if(argc < 3)
	{
		if(argc < 1) cout<<"INFO : No input root file was mentioned."<<endl;
		cout << "INFO : Exiting" << endl;
        return 0;
	}

    char *fileName = argv[1];
    chain.Add(fileName);
    cout <<"INFO : File is being processed" << endl;

    fstream myfile;
    myfile.open(argv[2],ios::out);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);			// Create class ExRootTreeReader Object
	Long64_t allEntries = treeReader->GetEntries();

    // Initialization of the arrays
    TClonesArray *branchJet = treeReader->UseBranch ("Jet");			// Get pointers to branches used in this analysis
	TClonesArray *branchElectron = treeReader->UseBranch ("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
	TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
	TClonesArray *branchTower = treeReader->UseBranch("Tower");
	TClonesArray *branchTrack = treeReader->UseBranch("Track");

    // Jet Variables
	Double_t fjet_R = 0.8;
	RecombinationScheme fjet_recomb_scheme=E_scheme;				// RecombinationScheme recomb_scheme= WTA_pt_scheme;
	Strategy fjet_strategy = Best;

	JetDefinition fjet_def(antikt_algorithm, fjet_R, fjet_recomb_scheme, fjet_strategy);

    Double_t fjet_mass_l = 70.0;
	Double_t fjet_mass_u = 90.0;
	Double_t ptmin_fjet  = 30.0;
	Double_t nsub32_cut = 0.81; // cut for checking 3-prong structure
	Double_t nsub21_cut = 0.35; // cut for checking 2-prong structure

    // Storage Vector
    vector<vector<vector<double_t>>> comp_dtset;

    // Looping through all entries
    for(Long64_t entry = 0; entry < allEntries; entry++){
        // Vector for the event level data
        vector<vector<double_t>> event_dtset;
        treeReader->ReadEntry(entry);

        vector<double_t> photon_dtset;
        for(Int_t i = 0; i < branchPhoton->GetEntriesFast(); ++i){
            Photon *photon_entry = (Photon*) branchPhoton->At(i);
            photon_dtset = {0.0,photon_entry->PT,photon_entry->P4().E(),photon_entry->P4().M(),photon_entry->Eta,photon_entry->Phi};
        }
        event_dtset.push_back(photon_dtset);

        vector<double_t> electron_dtset;
        for(Int_t i = 0; i < branchElectron->GetEntriesFast(); ++i){
            Electron *electron_entry = (Electron*) branchElectron->At(i);
            photon_dtset = {(double_t) electron_entry->Charge,electron_entry->PT,electron_entry->P4().E(),electron_entry->P4().M(),electron_entry->Eta,electron_entry->Phi};
        }
        event_dtset.push_back(electron_dtset);

        vector<double_t> muon_dtset;
        for(Int_t i = 0; i < branchMuon->GetEntriesFast(); ++i){
            Muon *muon_entry = (Muon*) branchMuon->At(i);
            photon_dtset = {(double_t) muon_entry->Charge,muon_entry->PT,muon_entry->P4().E(),muon_entry->P4().M(),muon_entry->Eta,muon_entry->Phi};
        }
        event_dtset.push_back(muon_dtset);

        vector<double_t> jet_dtset;
        for(Int_t i = 0; i < branchJet->GetEntriesFast(); ++i){
            Jet *jet_entry = (Jet*) branchJet->At(i);
            jet_dtset = {2.0,jet_entry->PT,jet_entry->P4().E(),jet_entry->P4().M(),jet_entry->Eta,jet_entry->Phi};
        }
        event_dtset.push_back(jet_dtset);
        
        // For getting the N-Subjettiness data
        vector<PseudoJet> towerjets;
        for(Int_t i = 0; i < branchTower->GetEntriesFast(); ++i){
				Tower *delphes_tower = (Tower*) branchTower->At(i);	
				TLorentzVector lorentz_jet = delphes_tower->P4();
				PseudoJet towerjet(lorentz_jet.Px(), lorentz_jet.Py(), lorentz_jet.Pz(),lorentz_jet.E());
				towerjets.push_back(towerjet);            
			}

		towerjets = sorted_by_pt(towerjets);
		
        vector<PseudoJet> input_particles = towerjets;        	        // form fjets and define filters
        ClusterSequence clust_seq(input_particles, fjet_def);
		// This makes a vector "fjets" with all fatjets with pt greater than or eq. to "ptmin_fjet"
		vector<PseudoJet> fjets = sorted_by_pt(clust_seq.inclusive_jets(ptmin_fjet));  
		
        // Defining the Subjettiness variables
		NsubjettinessRatio nSub21_beta(2, 1, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
		NsubjettinessRatio nSub32_beta(3, 2, OnePass_KT_Axes(), UnnormalizedMeasure(beta));

        vector<Double_t> lnsub32, lnsub21; 
        Double_t fjet_m, fjet_pt, nsub32, nsub21,num_n32=0,num_n21=0;

        for(unsigned ijet=0; ijet<fjets.size(); ijet++){
			fjet_m = fjets[ijet].m();
			fjet_pt = fjets[ijet].pt();
			
			nsub32 = nSub32_beta(fjets[ijet]);
			nsub21 = nSub21_beta(fjets[ijet]);
            
            lnsub32.push_back(nsub32);
            lnsub21.push_back(nsub21);
        }

        comp_dtset.push_back(event_dtset);
    }
    myfile.close();
    cout << "INFO : The dataset has been generated" << endl;
}