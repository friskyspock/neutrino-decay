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

double_t delR(double_t eta1, double_t eta2,double_t phi1, double_t phi2){
    /*
    Calculation of the Delta R, with the condition that if the difference 
    is more than pi, then it must be subtracted by pi as the range of 
    eta and phi is [0,pi]
    */
    double_t del_eta,del_phi;
    del_eta = fabs(eta1 - eta2);
    del_phi = fabs(phi1 - phi2);
    if(del_eta > M_PI)
        del_eta -= M_PI;
    if(del_phi > M_PI)
        del_phi -= M_PI;    

    return sqrt(pow(del_eta,2.0) + pow(del_phi,2.0));
}

tuple <double_t,double_t> invFL(vector<TLorentzVector> leptons, 
                                vector<TLorentzVector> fatjets,double_t sourceMass){
        /*
        Returns the min(m_wl - m_vr) over all w and l.
        This gives the lepton and W pair which are most likely 
        to come from the source particle (In this case, the right 
        handed neutrino)
        */
    
    double_t max_invFL = fabs((leptons[0] + fatjets[0]).M() /*- sourceMass*/);
    double_t drfl = delR(leptons[0].Eta(),fatjets[0].Eta(),leptons[0].Phi(),fatjets[0].Phi());
    
    for(Int_t i=0;i<leptons.size();i++){
        for(Int_t j=0;j<fatjets.size();j++){
            
            if(max_invFL < fabs((leptons[i] + fatjets[j]).M()/*- sourceMass*/));{
                max_invFL = fabs((leptons[i] + fatjets[j]).M()/* - sourceMass*/);
                drfl = delR(leptons[i].Eta(),fatjets[j].Eta(),leptons[i].Phi(),fatjets[j].Phi());
            }
        }
    }
    return make_tuple(max_invFL,drfl);
}

vector<Double_t> removeNan(vector<Double_t> orig_vec){
    /*
    Removes the NaN from the nsubjetiness pairs as 
    I don't know what a NaN is supposed to remember and 
    no statistical method would work on it.
    */
    vector<Double_t> ret_val;
    for(Int_t i=0;i<orig_vec.size();i++){
        if(!isnan(orig_vec[i]))
            ret_val.push_back(orig_vec[i]);
    }
    return ret_val;
}

int main(int argc,char* argv[]){
    gSystem->Load("libDelphes");
	TChain chain("Delphes");

    double beta = 2.0; 

    if(argc < 3){
		if(argc < 2) cout<<"INFO : The format is ./dtset input_file output_file"<<endl;
		cout << "INFO : Exiting" << endl;
        return 0;
	}

    /* Input the ROOT file that comes from Delphes for analysis */
    char *fileName = argv[1];
    chain.Add(fileName);
    cout <<"INFO : File is being processed for event level data" << endl;

    /* The name of the output CSV */
    char *outputFile = argv[2];
    std::ofstream myfile;
    myfile.open(outputFile);

    // Create class ExRootTreeReader Object
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t allEntries = treeReader->GetEntries();

    // Get pointers to branches used in this analysis
    TClonesArray *branchJet = treeReader->UseBranch ("Jet");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMET = treeReader->UseBranch("MissingET");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");

    // RecombinationScheme recomb_scheme= WTA_pt_scheme;
	Double_t fjet_R = 0.8;
	RecombinationScheme fjet_recomb_scheme=E_scheme;
	Strategy fjet_strategy = Best;

    /* Different radii and algorithms were tested but did not provide any 
    difference in the results and as the peaks remained the same (for 
    NSubJettiness), R = 0.8 and Anti-KT algorithm is chosen */
	JetDefinition fjet_def(antikt_algorithm, fjet_R, fjet_recomb_scheme, fjet_strategy);

    Double_t fjet_mass_l = 70.0;
	Double_t fjet_mass_u = 90.0;
	Double_t ptmin_fjet  = 30.0;
	
    /*
    The Features that are calculted are :
    - ptl : PT of the leading Lepton
    - etal : Eta of the leading Lepton
    - energyl : Energy of the leading Lepton
    - ptj : PT of the leading Jet
    - etaj : Eta of the leading Jet
    - energyj : Energy of the leading Jet
    - massj : Mass of the leading Jet
    - mjj : Invariant Mass of the two leading Jets
    - rjj : Delta R of the two leading Jets
    - rjl : Delta R of the leading Jet and leading Lepton
    - met : Missing Transverse Energy
    - n21_1 : The least valued nsub21
    - n21_2 : The second least valued nsub21
    - n32_1 : The least valued nsub32
    - n32_2 : The second least valued nsub32
    - infl : min(m_wl - m_vr)
    - drfl : Delta R between the Lepton and the FatJet involved in infl
    */

    myfile << "ptl,etal,energyl,ptj,etaj,energyj,massj,";
    myfile << "mjj,rjj,rjl,met,n21_1,n21_2,n32_1,n32_2,infl,drfl\n";

    for(Long64_t entry = 0; entry < allEntries; entry++){
        treeReader->ReadEntry(entry);
        
        int jet_no=0,lep_no=0;

        /* The vector which takes all the 
           leptons for checking the infl */ 
        vector<TLorentzVector> lepton_set;
        // The leading Lepton
        TLorentzVector lepton;
        Muon *leading_muon;
        double_t top_ptm = 0.0;
        for(Int_t i = 0; i < branchMuon->GetEntriesFast(); ++i){
            lep_no++;
			Muon *muonset = (Muon*) branchMuon->At(i);
            lepton_set.push_back(muonset->P4());	
			if(muonset->PT > top_ptm){
                top_ptm = muonset->PT;
                lepton = muonset->P4();
                leading_muon = muonset;
            }         
		}

        /* Getting the PT for Electron as the leading can be 
           any of Muon or Electron and thus checkign it */
        Electron *leading_electron;
        double_t top_pte = 0.0;
        for(Int_t i = 0; i < branchElectron->GetEntriesFast(); ++i){
			lep_no++;
            Electron *electronset = (Electron*) branchElectron->At(i);
            lepton_set.push_back(electronset->P4());	
			if(electronset->PT > top_pte){
                top_pte = electronset->PT;
                lepton = electronset->P4();
                leading_electron = electronset;
            }           
		}

        /* Getting the Leading the Second to Leading 
           Jets for invariant mass */
        Jet *leading_jet, *sec_jet;
        TLorentzVector lead_jet;
        double_t top_ptj=0.0,sec_ptj=0.0;
        for(Int_t i = 0; i < branchJet->GetEntriesFast(); ++i){
            jet_no++;
            Jet *jet_entry = (Jet*) branchJet->At(i);
            if(jet_entry->PT > top_ptj){
                top_ptj = jet_entry->PT;
                leading_jet = jet_entry;
                lead_jet = jet_entry->P4();
            }
            else if(jet_entry->PT > sec_ptj && jet_entry->PT < top_ptj){
                sec_ptj = jet_entry->PT;
                sec_jet = jet_entry;
            }
        }

        /* Getting the FatJets for getting the NSubJettiness */
        vector<PseudoJet> towerjets;
        for(Int_t i = 0; i < branchTower->GetEntriesFast(); ++i){
                Tower *delphes_tower = (Tower*) branchTower->At(i);	
                TLorentzVector lorentz_jet = delphes_tower->P4();
                PseudoJet towerjet(lorentz_jet.Px(), lorentz_jet.Py(), 
                                   lorentz_jet.Pz(),lorentz_jet.E());
                towerjets.push_back(towerjet);            
            }
        towerjets = sorted_by_pt(towerjets);
		
        // form fjets and define filters
        vector<PseudoJet> input_particles = towerjets;
        ClusterSequence clust_seq(input_particles, fjet_def);

		// This makes a vector "fjets" with all 
        // fatjets with pt greater than or eq. to "ptmin_fjet"
		vector<PseudoJet> fjets = sorted_by_pt(clust_seq.inclusive_jets(ptmin_fjet));

        // Defining the Subjettiness variables
		NsubjettinessRatio nSub21_beta(2, 1, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
		NsubjettinessRatio nSub32_beta(3, 2, OnePass_KT_Axes(), UnnormalizedMeasure(beta));

        /* Going through the clustered Jets and storing them as 
           TLorentzVector and also N21 and N32 */
        vector<TLorentzVector> tower_set;
        vector<Double_t> lnsub32, lnsub21; 
        Double_t fjet_m, fjet_pt, nsub32, nsub21,num_n32=0,num_n21=0;

        for(unsigned ijet=0; ijet<fjets.size(); ijet++){
				fjet_m = fjets[ijet].m();
				fjet_pt = fjets[ijet].pt();

                TLorentzVector ljet;
                /* Making a different TLorentzVector as the P4() 
                   provides something different => Need to investigate */
                ljet.SetPtEtaPhiM(fjets[ijet].pt(),fjets[ijet].eta(),
                                  fjets[ijet].phi(),fjets[ijet].m());
                tower_set.push_back(ljet);

				nsub32 = nSub32_beta(fjets[ijet]);
				nsub21 = nSub21_beta(fjets[ijet]);

                lnsub32.push_back(nsub32);
                lnsub21.push_back(nsub21);
            }

        // Removing NaNs from n21 and n32
        lnsub21 = removeNan(lnsub21);
        lnsub32 = removeNan(lnsub32);

        if(lnsub21.size() < 2 || lnsub32.size() < 2)
            continue;

        if(lep_no < 1 || jet_no < 2)
            continue;

        double_t lead_lep_eta,lead_lep_phi;

        /* Storing the ptl,eta_l and also getting eta_l and phi_l
           for the calculation of delR_jl */
        
        if(top_ptm >= top_pte){
            myfile << leading_muon->PT << "," << leading_muon->Eta << ",";
            lead_lep_eta = leading_muon->Eta;
            lead_lep_phi = leading_muon->Phi;
        }
        else{
            myfile << leading_electron->PT << "," << leading_electron->Eta << ",";
            lead_lep_eta = leading_electron->Eta;
            lead_lep_phi = leading_electron->Phi;
        }

        // Adding Lepton Energy
        myfile << lepton.E() << ",";
        myfile << leading_jet->PT << "," << leading_jet->Eta << ",";

        // Adding Jet Energy and Mass
        myfile << lead_jet.E() << "," << lead_jet.M() << ",";

        // The invariant mass of two jets
        myfile << (leading_jet->P4() + sec_jet->P4()).M() << ",";

        // The del_R
        myfile << delR(leading_jet->Eta,sec_jet->Eta,leading_jet->Phi,sec_jet->Phi) << ",";
        myfile << delR(leading_jet->Eta,lead_lep_eta,leading_jet->Phi,lead_lep_phi) << ",";
		
        TLorentzVector MET;
        for(Int_t i = 0; i < branchMET->GetEntriesFast(); ++i){
			MissingET *metset = (MissingET*) branchMET->At(i);	
			MET = metset->P4();
		}
        // Adding MET
        myfile << MET.E() << ",";


        sort(lnsub21.begin(),lnsub21.end());
        sort(lnsub32.begin(),lnsub32.end());

        for(Int_t i=0;i<2;i++)
            myfile << lnsub21[i] << ",";

        for(Int_t i=0;i<2;i++)
            myfile << lnsub32[i] << ",";

        double_t inFL,drFL;
        tie(inFL,drFL) = invFL(lepton_set,tower_set,(double_t) 1000.0);
        myfile << inFL << "," << drFL << "\n";

    }
    myfile.close();
    cout << "INFO : The subjets have been added to the file" << endl;
    return 0;
}