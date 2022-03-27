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
    cout <<"INFO : File is being processed for Nsubjettiness" << endl;

    char *outputFile = argv[2];
    std::ofstream myfile;
    myfile.open(outputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);			// Create class ExRootTreeReader Object
	Long64_t allEntries = treeReader->GetEntries();

	TClonesArray *branchTower = treeReader->UseBranch("Tower");

	Double_t fjet_R = 0.9;
	RecombinationScheme fjet_recomb_scheme=E_scheme;				// RecombinationScheme recomb_scheme= WTA_pt_scheme;
	Strategy fjet_strategy = Best;

	JetDefinition fjet_def(antikt_algorithm, fjet_R, fjet_recomb_scheme, fjet_strategy);

    Double_t fjet_mass_l = 70.0;
	Double_t fjet_mass_u = 90.0;
	Double_t ptmin_fjet  = 30.0;
	Double_t nsub32_cut = 0.81; // cut for checking 3-prong structure
	Double_t nsub21_cut = 0.35; // cut for checking 2-prong structure

    for(Long64_t entry = 0; entry < allEntries; entry++){
        treeReader->ReadEntry(entry);
        
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

        if(lnsub21.size() > 0){
            for(unsigned i=0;i<lnsub21.size() - 1;i++)
                    myfile << lnsub21[i] << ",";
            myfile << lnsub21[lnsub21.size()-1] << ";";
        }
        else
            myfile << ";";
        
        if(lnsub32.size() > 0){
            for(unsigned i=0;i<lnsub32.size() - 1;i++)
                    myfile << lnsub32[i] << ",";
            myfile << lnsub21[lnsub21.size()-1] << "\n";
        }
        else
            myfile << "\n";

        //cout << endl;
        //cout << endl;
    }
    myfile.close();
    cout << "INFO : The subjets have been added to the file" << endl;
}