//#ifdef __CLING__
//R__LOAD_LIBRARY (libDelphes);
//#endif
//R_ADD_INCLUDE_PATH(/home/student02/hep/Delphes-3.4.1)
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
//#include <fastjet/tools/JHTopTagger.hh>

//#include <fastjet/tools/MassDropTagger.hh>  	// needed for HEPTopTagger
//#include <fastjet/tools/Filter.hh>   		// needed for HEPTopTagger
//#include <fastjet/tools/Pruner.hh>   		// needed for HEPTopTagger
//#include <fastjet/tools/HEPTopTagger.hh>

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"

#include "classes/DelphesClasses.h"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//utility function to sort 2d vectors
bool sortcol(const TLorentzVector &v1, const TLorentzVector &v2)
{
	return v1.Perp() > v2.Perp();			// perp() : returns the scalar transverse momentum 
}

int main(int argc, char *argv[])
{

	const char *inputFileName;	
	int check1(0);
	double beta = 2.0;

	//reading Delphes tree from input file

	gSystem->Load("libDelphes");
	//  gSystem->Load ("mt2_bisect_cpp");  // root -l; .L mt2_bisect.cpp++ ==> ++ at the end, creates .so file to link for further use
	TChain chain("Delphes"); 
	cout<<endl<<"#----------------------------------------------------"<<endl;
	if(argc < 2)
	{
		if(argc < 1) cout<<"# [ERROR] No input root file was mentioned."<<endl;
		cout<<"# Stopping."<<endl<<"#"<<endl<<"# Run with:"<< endl<<"# $"<<argv[0]<<" <input file 1> [<input file 2> ...]"<<endl;
		cout<<"#----------------------------------------------------"<<endl;
		return 0;
	}
	else
	{	

		cout<<"# HEPTopTagger will be used."<<endl<<"#"<<endl;     
		for(int i = 1; i < argc; ++i)
		{
			inputFileName = argv[i];
			chain.Add (inputFileName);
			cout<<"# Input File ("<<i-1<<"): \""<<inputFileName<<"\""<<endl;
		}

	}
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);			// Create class ExRootTreeReader Object
	Long64_t allEntries = treeReader->GetEntries();

	cout<<"#"<<endl<<"# Processing "<<allEntries<<" events..."<<endl;
	cout<<"#----------------------------------------------------"<<endl<<endl;

	TClonesArray *branchJet = treeReader->UseBranch ("Jet");			// Get pointers to branches used in this analysis
	TClonesArray *branchElectron = treeReader->UseBranch ("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
	TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
	TClonesArray *branchTower = treeReader->UseBranch("Tower");
	TClonesArray *branchTrack = treeReader->UseBranch("Track");

	//----------------------------------------------------------------------
	// edit cuts here
	//----------------------------------------------------------------------

	/* -----------------Fat jet variables----------------------------------*/

	Double_t fjet_R = 1.5;
	RecombinationScheme fjet_recomb_scheme=E_scheme;				// RecombinationScheme recomb_scheme= WTA_pt_scheme;
	Strategy fjet_strategy = Best;

	JetDefinition fjet_def(cambridge_algorithm, fjet_R, fjet_recomb_scheme, fjet_strategy);

	//double ptmin_fjet(180.0);
	Double_t fjet_mass_l = 135.0;
	Double_t fjet_mass_u = 215.0;
	Double_t ptmin_fjet  = 135.0;//180.0;

	/* -----------------HTT variables-------------------------------------*/

	Double_t SubjetMass_max 	= 50.0;
	Double_t MassDropThreshold 	= 0.8;
	bool 	 CheckMassratioFirst	= true;
	Double_t topmass_lower 		= 135.0;	/* top mass range */
	Double_t topmass_upper 		= 215.0;	/* top mass range */
	Double_t PrunRcut		= 0.75;
	Double_t PrunZcut		= 0.1;

	/* -------------- lepton nuber-------------------------------------- */
	Int_t	reqlepton_no = 2; //total number of e+mu+tau-h

	/* -----------lepton invariant mass cut | Z-veto-------------------- */

	Double_t mmll_cut = 120.0;

	/* ------------------------pt cut---------------------------------- */

	//Double_t min_pt_lep = 40.0;
	Double_t min_pt_el = 30.0;
	Double_t min_pt_mu = 30.0;
	Double_t min_pt_tau_h = 100.0;
	//Double_t min_pt_jt = 30.0;

	Double_t min_pt_top = 300.0;

	/* ------------------------eta cut--------------------------------- */

	Double_t max_eta_el = 2.5;
	Double_t max_eta_mu = 2.5;
	Double_t max_eta_jt = 5.0;
	Double_t max_eta_tau_h = 5.0;
	Double_t max_eta_fj = 5.0;


	/* --------------------barrel-endcap cut--------------------------- */  

	Double_t eta_BE_ll = 1.37;
	Double_t eta_BE_ul = 1.52;

	/* ---------------------missing ET cut---------------------------- */

	Double_t min_met = 0.0;
	
	/* ---------------------Scalar HT cut---------------------------- */

	Double_t min_ht = 1.2*1000.0;

	/*----------------------LQ cuts-------------------------------------*/
	Double_t min_lq_m = 0.8*1000.0;

	//---------------------------------------------------------------------- 
	// end cuts
	//----------------------------------------------------------------------



	gROOT->Reset();

	// TCanvas* c1 = new TCanvas("c1", "pt", 200, 10, 600, 400);

	// TH1D* h1 = new TH1D("LQ from Fatjet", "", 100, 0.0, 2500.0);
	// TH1D* h3 = new TH1D("Fatjet Mass", "", 100, 0.0, 2500.0);

	// TH1D* h2 = new TH1D("LQ from Tagged Tops", "", 100, 0.0, 2500.0);
	// TH1D* h4 = new TH1D("Tagged-top Mass", "", 100, 0.0, 2500.0);

	// TH1D* h5 = new TH1D("M_{ll}", "", 100, 0.0, 500.0);

	Electron* el[50]; 
	Muon* 	  mu[50];
	Jet*      jt[50]; 
	ScalarHT  *ht;

	TLorentzVector t_top, f_top;
	int total_tops=0, lq_cand=0, singl_top=0, dble_top=0, n_top=0, lep_top_events=0, leptop_lq_cand=0, hadtop_lq_cand=0, less_lep_events=0;
	int lq_cand_type[5];
	memset(lq_cand_type, 0, sizeof(lq_cand_type));		// Filling the whole array lq_cand_type with 0.
	int event_type = -1;

	// ofstream myfile;
    // 	myfile.open ("Lqtat_output_Sig.dat",ios::app); 


	for(Long64_t entry = 0; entry < allEntries; entry++)				// loop over all events
	{
		int tops_per_event = 0;
		Bool_t had_top = 0, lep_top = 0;                                // if there's an hadronic/leptonic top

		treeReader->ReadEntry(entry);  					// load selected branches with data from specified event

		vector<TLorentzVector> lepton, tau_h;			// defining a vector named lepton and tau_h with class TLorentzVector
		TLorentzVector v;								// TLorentzVector class with classname v

		// checking electron
		for(int i = 0; i < branchElectron->GetEntriesFast(); i++)
		{
			el[i] = (Electron *) branchElectron->At(i);		// Putting all the electron in the branch electron into array el[]
			if(el[i]->PT >= min_pt_el && (fabs(el[i]->Eta) <= eta_BE_ul || fabs(el[i]->Eta) >= eta_BE_ll) && fabs(el[i]->Eta) <= max_eta_el) 
			{
				v.SetPxPyPzE(el[i]->P4().Px(), el[i]->P4().Py(), el[i]->P4().Pz(), el[i]->P4().E());  
				lepton.push_back(v);		// It push the elements into a vector from the back.
			}
		}

		// checking muon
		for(int i = 0; i < branchMuon->GetEntriesFast(); i++)
		{
			mu[i] = (Muon *) branchMuon->At(i);
			if(mu[i]->PT >= min_pt_mu && fabs(mu[i]->Eta) <= max_eta_mu )   
			{
				v.SetPxPyPzE(mu[i]->P4().Px(), mu[i]->P4().Py(), mu[i]->P4().Pz(), mu[i]->P4().E());  		
				lepton.push_back(v);
			}
		}

		sort(lepton.begin(), lepton.end(), sortcol);

		// checking jets for hadronic tau
		for(int i = 0; i < branchJet->GetEntriesFast(); i++)
		{
			jt[i] = (Jet*) branchJet->At(i);
			if(jt[i]->TauTag==1 && jt[i]->PT >= min_pt_tau_h && fabs(jt[i]->Eta) <= max_eta_tau_h)
			{
				v.SetPxPyPzE(jt[i]->P4().Px(), jt[i]->P4().Py(), jt[i]->P4().Pz(), jt[i]->P4().E());  //P4 is a class which returns components of 4 momentum
				tau_h.push_back(v);	
			}
			if(tau_h.size() >= reqlepton_no) break;			// The vector tau_h cannot have more than 2 leptons.
		}

		// analyse scalar HT
		if(branchScalarHT->GetEntriesFast() > 0)
		{
			ht = (ScalarHT*) branchScalarHT->At(0);
		}

		// clustering using fastjet
		vector<PseudoJet> towerjets;					// jets from towers

		for(Int_t i = 0; i < branchTower->GetEntriesFast(); ++i)
		{
			Tower *delphes_tower = (Tower*) branchTower->At(i);			//Tower is the class.
			TLorentzVector lorentz_jet = delphes_tower->P4();
			PseudoJet towerjet(lorentz_jet.Px(), lorentz_jet.Py(), lorentz_jet.Pz(),lorentz_jet.E());
			towerjets.push_back(towerjet);            
		}

		towerjets = sorted_by_pt(towerjets);

		vector<PseudoJet> input_particles = towerjets;        	        // form fjets and define filters
		ClusterSequence clust_seq(input_particles, fjet_def);
		// This makes a vector "fjets" with all fatjets with pt greater than or eq. to "ptmin_fjet"
		vector<PseudoJet> fjets = sorted_by_pt(clust_seq.inclusive_jets(ptmin_fjet));  
		NsubjettinessRatio nSub21_beta(2, 1, OnePass_KT_Axes(), UnnormalizedMeasure(beta));


		// Top tagging using HEP toptagger
	
		// cout << "#----- HEPTopTagger definition -----#" << endl;
		fastjet::HEPTopTagger heptop1(SubjetMass_max); //Maximum subjet mass
		heptop1.set_mass_drop_threshold(MassDropThreshold);
		heptop1.check_massratio_first(CheckMassratioFirst);
		heptop1.set_top_range(topmass_lower,topmass_upper);
		//heptop1.set_mass_ratio_cut(0.35,0.2,1.3);
		//heptop1.set_max_subjet_mass(50.);
		heptop1.set_prun_rcutfac(PrunRcut);
		heptop1.set_prun_zcut(PrunZcut);
		heptop1.set_filter_jetalgorithm(kt_algorithm); //kt performs better than c/a
		heptop1.set_recluster_jetalgorithm(kt_algorithm);
		heptop1.set_prune_jetalgorithm(kt_algorithm); //added 
		heptop1.return_only_tagged(true);



		// loop over all jets
		for(unsigned ijet=0; ijet<fjets.size(); ijet++)
		{
			PseudoJet taggedtop=heptop1(fjets[ijet]);

			if(taggedtop != PseudoJet()) // heptop1.is_maybe_top();
			{
				bool istop_bool = taggedtop.structure_of<HEPTopTagger>().is_tagged();
				// test with "min_pt_top" is not really needed if "min_pt_top" 
				// is set to be less than or equal to "ptmin_fjet"
				if(istop_bool && taggedtop.perp()>=min_pt_top) 					// perp returns scalar transverse momentum.
				{
					had_top = 1;

					t_top.SetPxPyPzE(taggedtop.px(), taggedtop.py(), taggedtop.pz(), taggedtop.e());
					tops_per_event++;

					total_tops++;
				}

			}
		}
		switch(tops_per_event)
		{
			case 0: break; 
			case 1: singl_top++; break;
			case 2: dble_top++; break;
			default: n_top++; break;
		}


		double total_m = 0;
		event_type = -1;

		double lq_pt = 0.0;
		vector<TLorentzVector> lq(2);

		/*------Hadronic tops------*/

		// 2 hadronic taus and a hadronic top (hhH - 2 had-tau + top)
		if(tau_h.size() == reqlepton_no && had_top>0 && lepton.size() == 0)
		{
			less_lep_events++;
			event_type = 0;
			for(int i = 0; i < 2; i++)
			{
				lq[i] = t_top + tau_h[i];
				if(lq[i].M() > total_m) 
				{
					total_m = lq[i].M();
					lq_pt = lq[i].Pt();
				}
			}
			if( (tau_h[0] + tau_h[1]).M() > mmll_cut && ht->HT >= min_ht && total_m >= min_lq_m)
			{
				lq_cand++;
				lq_cand_type[0]++;
			}
		}

	}

	cout << endl << "#----------------------------------------------------" << endl;

	cout << "# Processed " << allEntries << " events." << endl; 
	cout << "# Total single hadronic top events: " << singl_top << " (" << (100.0*singl_top/allEntries) << " %)" << endl;
	cout << "# Total double hadronic top events: " << dble_top << " (" << (100.0*dble_top/allEntries) << " %)"<< endl;  
	//cout << "# Total n (>2) hadronic top events: " << n_top << " (" << (100.0*n_top/allEntries) << " %)" << endl;        
	cout << "# Total number of hhH events with a LQ candidate: " << lq_cand_type[0] << " (" << 100.0*lq_cand_type[0]/allEntries << "%)" << endl;
	// cout << "# Total number of hlH events with a LQ candidate: " << lq_cand_type[1] << " (" << 100.0*lq_cand_type[1]/allEntries << "%)" << endl;
	cout << "# Total number of events with a LQ candidate: " << lq_cand << " (" << (100.0*lq_cand/allEntries) <<" %)"<<endl;
	cout << "# Total number of events with at least "<< reqlepton_no<<" leptons (including had-taus): " << less_lep_events << " (" << 100.0*(less_lep_events)/allEntries << " %)"<< endl;
	cout << "#----------------------------------------------------" << endl << endl;
	//myfile << "For value min_pt_lep :" << min_pt_lep << endl;
	//myfile << "For value min_pt_jt :" << min_pt_jt << endl;
	// myfile << "For value max_eta_jt :" <<max_eta_jt << endl;
	// myfile << "For value max_eta_tau_h :"<< max_eta_tau_h << endl;
	// myfile << "For value min_ht:" << min_ht << endl;
	// myfile << "For value min_pt_tau_h :"<< min_pt_tau_h << endl;
	// myfile << "For value min_pt_el:" << min_pt_el << endl;
	// myfile << "For value min_pt_mu:" << min_pt_mu << endl;
	// myfile << "For value min_pt_top:" << min_pt_top << endl;
	// myfile << "For value min_lq_m:" << min_lq_m << endl;
  	// myfile << "Total entries: " << allEntries << "| hhH events: " << (Double_t)lq_cand_type[0]/allEntries << "| hlH events: " << (Double_t)lq_cand_type[1]/allEntries << "| llH events: " << (Double_t)lq_cand_type[2]/allEntries << "| Total LQ candidates: " << (Double_t)lq_cand/allEntries <<endl;
  	// myfile << "# of signal events = " << 49395* ((Double_t)lq_cand/allEntries) <<endl;
  	// myfile << "------------------------------------------------------------------------------------------------------------------------------" <<endl;
  	// myfile << endl;
  	// myfile.close();

	return 0;
}


