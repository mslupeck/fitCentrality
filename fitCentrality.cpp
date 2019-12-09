#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/stat.h> // check for existence of file

#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>

#include "TreeStorage.h"

inline bool fileExists (const std::string& name) {
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

// Reads root files after ESD distiller and merges them into one file
int ReadAndMerge(UInt_t nFiles){
	std::string inputBaseFolder = "inputs/B500_FullImpact/";
	std::vector<std::string> vInputPath;
	for(UInt_t ifile=0; ifile<nFiles; ifile++){
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(3) << ifile+1 << "/";
		vInputPath.push_back(inputBaseFolder + ss.str() + "AnalysisResults.root");
	}

	std::vector<TFile*> vInFile;
	std::cout << "Open all files..." << std::endl;
	for(UInt_t ifile=0; ifile<vInputPath.size(); ifile++){
		if(fileExists(vInputPath[ifile].c_str())){
			vInFile.push_back(new TFile(vInputPath[ifile].c_str(),"READ"));
			if(vInFile.at(vInFile.size()-1)->IsOpen()){
				std::cout << "File opened: " <<  vInputPath[ifile] << std::endl;
			}
		}
	}

	TFile *fout = new TFile("outputs/merged.root","RECREATE");
	TTree *mtr = new TTree("mtr","mtr");
	TreeStorage *mts = new TreeStorage();
	mtr->Branch("MTreeStorage",mts);

	TreeStorage *ts = NULL;
	TTree *tr = NULL;
	for(UInt_t ifile=0; ifile<vInFile.size(); ifile++){
		TFile *f = vInFile.at(ifile);
		if(vInFile.at(ifile)->IsOpen()){
			f->cd();
			TList* l = (TList*)(f->Get("mssFIT"));
			tr = (TTree*)(l->At(0));
			if(tr==NULL){
				std::cout << "<E> Problems opening tree." << std::endl;
				return -1;
			}
			tr->SetBranchAddress("TreeStorage", &ts);
			if(ts==NULL){
				std::cout << "<E> Problems opening TreeStorage." << std::endl;
				return -2;
			}
			// Copy the contents of the tree read from file to the in-memory merged one
			UInt_t nEntries = tr->GetEntries();
			std::cout << nEntries << std::endl;
			for(UInt_t ientry=0; ientry<nEntries; ientry++){
				tr->GetEntry(ientry);
				(*mts) = (*ts);
				mtr->Fill();
			}
			f->Close();
		}
	}

	fout->cd();
	mtr->Write();
	fout->Close();

	return 0;
}

void InitGraphNames(std::vector<std::string> &vGrName){
	vGrName.push_back("FIT_T0A");
	vGrName.push_back("FIT_T0C");
	vGrName.push_back("FIT_V0A");
	vGrName.push_back("etaFIT_T0A");
	vGrName.push_back("etaFIT_T0C");
	vGrName.push_back("etaFIT_V0A");
	vGrName.push_back("FIT_T0A+T0C");
	vGrName.push_back("FIT_T0A+T0C+V0A");
	vGrName.push_back("etaFIT");
	vGrName.push_back("etaV0A_Run2");
	vGrName.push_back("etaV0C_Run2");
	vGrName.push_back("etaV0AC_Run2");
}

void InitGraphs(std::vector<TGraph*> &vGr, UInt_t nEntries){
	std::vector<std::string> vGrNames;
	InitGraphNames(vGrNames);
	for(UInt_t igr=0; igr<vGrNames.size(); igr++){
		TGraph* gr = new TGraph(nEntries);
		gr->SetName(vGrNames.at(igr).c_str());
		vGr.push_back(gr);
	}
}

void RemoveLastPoints(std::vector<TGraphErrors*> &vge){
	for(UInt_t igr=0; igr<vge.size(); igr++){
		vge.at(igr)->Set(8);
	}
}

int Draw(int argc, char* argv[]){
	TApplication theApp("App", &argc, argv);

	std::cout << "Drawing..." << std::endl;
	TFile *f = new TFile("outputs/centrRes.root","READ");
	if(!f->IsOpen()){
		std::cout << "<E> Root file not opened. Exiting" << std::endl;
		return -1;
	}
	std::vector<std::string> vGrNames;
	InitGraphNames(vGrNames);
	std::vector<TGraphErrors*> vge;
	for(UInt_t igr=0; igr<vGrNames.size(); igr++){
		std::stringstream ss;
		ss << vGrNames.at(igr) << "-uncert-" << igr;
		vge.push_back((TGraphErrors*)f->Get(ss.str().c_str()));
		if(vge.at(igr) == NULL){
			std::cout << "<E> Problem reading graph. Exiting" << std::endl;
			return -2;
		}
	}


	std::vector<Color_t> vCol;
	vCol.push_back(kRed);
	vCol.push_back(kBlue);
	vCol.push_back(kGreen+2);
	vCol.push_back(kRed-9);
	vCol.push_back(kBlue-9);
	vCol.push_back(kGreen-9);
	vCol.push_back(kMagenta);
	vCol.push_back(kBlack);
	vCol.push_back(kGray);
	vCol.push_back(kOrange);
	vCol.push_back(kOrange+2);
	vCol.push_back(kOrange+4);

	RemoveLastPoints(vge);

	TCanvas *c = new TCanvas("c","c",1700,920);
	gPad->SetGridx();
	gPad->SetGridy();
	TLegend *leg = new TLegend(0.11,0.88, 0.35, 0.55);

	for(UInt_t igr=0; igr<vge.size(); igr++){
		std::string sDraw="";
		if(igr>0){
			sDraw = "same";
		}
		vge.at(igr)->Draw(sDraw.c_str());
		vge.at(igr)->SetTitle("");
		vge.at(igr)->SetLineWidth(2);
		vge.at(igr)->SetLineColor(vCol.at(igr));
		vge.at(igr)->SetMarkerSize(1.4);
		leg->AddEntry(vge.at(igr), vge.at(igr)->GetName(), "lep");
	}
	leg->Draw();

	TCanvas *cFinal = new TCanvas("cFinal","cFinal",900,900);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.04);
	gPad->SetTopMargin(0.02);
	gPad->SetBottomMargin(0.09);

	TLegend *legFinal = new TLegend(0.42, 0.13, 0.92, 0.33);
	vge.at(0)->GetXaxis()->SetRangeUser(0,80);
	vge.at(0)->Draw("");
	vge.at(0)->GetXaxis()->SetLabelSize(0.05);
	vge.at(0)->GetXaxis()->SetTitleOffset(1.0);
	vge.at(0)->GetYaxis()->SetLabelSize(0.04);
	vge.at(0)->GetYaxis()->SetTitleOffset(1.4);

	vge.at(1)->Draw("same LP");
	vge.at(2)->Draw("same LP");
//	vge.at(6)->Draw("same LP");
	vge.at(7)->Draw("same LP");

	vge.at(0)->SetMarkerStyle(22);
	vge.at(1)->SetMarkerStyle(23);
	vge.at(2)->SetMarkerStyle(21);
	vge.at(7)->SetMarkerStyle(20);
	vge.at(0)->SetMarkerColor(kBlack);
	vge.at(1)->SetMarkerColor(kBlack);
	vge.at(2)->SetMarkerColor(kBlack);
	vge.at(7)->SetMarkerColor(kBlack);

	legFinal->AddEntry(vge.at(0), "T0A+", "lep");
	legFinal->AddEntry(vge.at(1), "T0C+", "lep");
	legFinal->AddEntry(vge.at(2), "V0A+", "lep");
//	legFinal->AddEntry(vge.at(6), "T0A+ & T0C+", "le");
	legFinal->AddEntry(vge.at(7), "FIT (T0A+ & T0C+ & V0A+)", "lep");
	legFinal->Draw();

	TLatex *t = new TLatex(0.14, 0.9, "#bf{#splitline{ALICE simulation}{Pb-Pb, #sqrt{s_{NN}} = 5.5TeV}}");
	t->SetNDC();
	t->SetTextSize(0.04);
	t->SetTextAttributes();
	t->Draw();

	std::cout << "Done." << std::endl;
	theApp.SetIdleTimer(10000,"exit()"); // exit after 30 min of being idle
	theApp.SetReturnFromRun(true);
	theApp.Run();
	return 0;
}

int main(int argc, char* argv[]){
	// Simple command line handler
	if(argc > 1){
		if(argv[1][0] == 'd'){
			return Draw(argc, argv);
		}
		if(argv[1][0] == 'm'){
			// Merge ESD distiller files
			int err = ReadAndMerge(36);
			if(err < 0){
				std::cout << "<E> Error while merging files" << std::endl;
				return err;
			}
		}
	}

	// A-side
	const float etaMinV0A_R2 = 2.8;
	const float etaMaxV0A_R2 = 5.1;
	const float etaMinV0A = 2.15;
	const float etaMaxV0A = 5.03;
	const float etaMinT0A = 3.75;
	const float etaMaxT0A = 5.04;

	// C-side
	const float etaMinV0C_R2 = -3.7;
	const float etaMaxV0C_R2 = -1.7;
	const float etaMinT0C = -3.35;
	const float etaMaxT0C = -2.26;

	const UInt_t chMaxT0A = 96;
	const UInt_t chMaxT0C = 208;
	const UInt_t chMaxV0A = 288;

	TApplication theApp("App", &argc, argv);
	gStyle->SetOptStat("oueni");

	// Open the merged file and setup tree reading
	TFile *f = new TFile("outputs/merged.root","READ");
	if(!f->IsOpen()){
		std::cout << "<E> Merged file failed to open." << std::endl;
		return -10;
	}
	TTree *mtr = (TTree*)(f->Get("mtr"));
	if(mtr == NULL){
		std::cout << "<E> Merged tree failed to open." << std::endl;
		return -11;
	}
	TreeStorage *mts = NULL;
	mtr->SetBranchAddress("MTreeStorage",&mts);
	if(mts == NULL){
		std::cout << "<E> Error while assigning address of the merged branch." << std::endl;
		return -12;
	}

	// Run through tree entries, copy selected data to graphs (to make it easier to sort)
	UInt_t nEntries = mtr->GetEntries();
//	nEntries = 2000;
	std::cout << "  <I> nEntries = " << nEntries << std::endl;
	std::vector<TGraph*> vgr;
	InitGraphs(vgr, nEntries);
	for(UInt_t ientry=0; ientry<nEntries; ientry++){
		mtr->GetEntry(ientry);
		TH1S* hEtaPrimCh = mts->hEtaPrimCharged;
		TH1S* hEtaAllCh = mts->hEtaAllCharged;

		float bTrue = mts->trueImpactParam;
		int nPrimCharged = hEtaPrimCh->Integral();
		int nAllCharged = hEtaAllCh->Integral();
		int mV0A_R2 = 0;
		int mV0C_R2 = 0;
		int mEtaT0A = 0;
		int mEtaT0C = 0;
		int mEtaV0A = 0;
		int mT0A = 0;
		int mT0C = 0;
		int mV0A = 0;

		// Extract old V0 multiplicities in the ideal case
		TAxis* ax = hEtaAllCh->GetXaxis();
		mV0A_R2 = hEtaAllCh->Integral(ax->FindBin(etaMinV0A_R2), ax->FindBin(etaMaxV0A_R2));
		mV0C_R2 = hEtaAllCh->Integral(ax->FindBin(etaMinV0C_R2), ax->FindBin(etaMaxV0C_R2));
		mEtaT0A = hEtaAllCh->Integral(ax->FindBin(etaMinT0A), ax->FindBin(etaMaxT0A));
		mEtaT0C = hEtaAllCh->Integral(ax->FindBin(etaMinT0C), ax->FindBin(etaMaxT0C));
		mEtaV0A = hEtaAllCh->Integral(ax->FindBin(etaMinV0A), ax->FindBin(etaMaxV0A));

		// Go through all channels and count number of MIPs for each detector array
		for(UInt_t ich=0; ich<mts->nCh; ich++){
			if(ich<chMaxT0A){ // T0A
				mT0A += round(mts->fitAmp[ich]);
			}
			else if(ich<chMaxT0C){ // T0C
				mT0C += round(mts->fitAmp[ich]);
			}
			else if(ich<chMaxV0A){ // V0A
				mV0A += round(mts->fitAmp[ich]);
			}
			else{
				std::cout << " <W> Wrong channel number !!! Check the ranges within the code." << std::endl;
			}
		}

		vgr.at(0)->SetPoint(ientry, mT0A, nPrimCharged);
		vgr.at(1)->SetPoint(ientry, mT0C, nPrimCharged);
		vgr.at(2)->SetPoint(ientry, mV0A, nPrimCharged);
		vgr.at(3)->SetPoint(ientry, mEtaT0A, nPrimCharged);
		vgr.at(4)->SetPoint(ientry, mEtaT0C, nPrimCharged);
		vgr.at(5)->SetPoint(ientry, mEtaV0A, nPrimCharged);
		vgr.at(6)->SetPoint(ientry, mT0A+mT0C, nPrimCharged);
		vgr.at(7)->SetPoint(ientry, mT0A+mT0C+mV0A, nPrimCharged);
		vgr.at(8)->SetPoint(ientry, mEtaT0A+mEtaT0C+mEtaV0A, nPrimCharged);
		vgr.at(9)->SetPoint(ientry, mV0A_R2, nPrimCharged);
		vgr.at(10)->SetPoint(ientry, mV0C_R2, nPrimCharged);
		vgr.at(11)->SetPoint(ientry, mV0A_R2+mV0C_R2, nPrimCharged);
	}

	// Make copies of the graphs to be sorted
	std::vector<TGraph*> vGrSortX, vGrSortY;
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		vGrSortX.push_back(new TGraph(*(vgr.at(igr))));
		vGrSortY.push_back(new TGraph(*(vgr.at(igr))));
	}

	// Sort all the graphs
	std::cout << "  <I> Sort: ";
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		std::cout << igr << " ";
		std::cout.flush();
		vGrSortX.at(igr)->Sort(TGraph::CompareX,kFALSE);
		vGrSortY.at(igr)->Sort(TGraph::CompareY,kFALSE);
	}
	std::cout << std::endl;

	// Make graphs with centrality fraction assigned to y-axis
	const int minCentr = 0;		//  0%
	const int maxCentr = 100;	// 100%
	const int widthCentr = maxCentr - minCentr;
	std::vector<TGraph*> vGrX2Cent, vGrY2Cent;
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		vGrX2Cent.push_back(new TGraph(vgr.at(igr)->GetN()));
		vGrY2Cent.push_back(new TGraph(vgr.at(igr)->GetN()));
		for(int ipoint=0; ipoint<vGrSortX.at(igr)->GetN(); ipoint++){
			double x,y;
			vGrSortX.at(igr)->GetPoint(ipoint,x,y);
			vGrX2Cent.at(igr)->SetPoint(ipoint, x, (1.0+ipoint)/vGrSortX.at(igr)->GetN()*widthCentr + minCentr);
		}
		for(int ipoint=0; ipoint<vGrSortY.at(igr)->GetN(); ipoint++){
			double x,y;
			vGrSortY.at(igr)->GetPoint(ipoint,x,y);
			vGrY2Cent.at(igr)->SetPoint(ipoint, y, (1.0+ipoint)/vGrSortY.at(igr)->GetN()*widthCentr + minCentr);
		}
	}

	// Graphs with centrality correlation and difference between the detector and reference (all nChargedPart)
	std::vector<TGraph*> vGrCentCorr;
	std::vector<TGraph*> vGrDiff;
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		TGraph* g = vgr.at(igr);
		if(true){
			std::cout << "  <I> ";
			std::cout << std::setfill(' ') << std::setw(2);
			std::cout << igr << "/" << vgr.size() << " (";
			std::cout << std::setw(2) << (int)(100.0*igr/vgr.size()) << "%) - ";
		}
		vGrCentCorr.push_back(new TGraph(g->GetN()));
		vGrDiff.push_back(new TGraph(g->GetN()));
		for(int ipoint=0; ipoint<g->GetN(); ipoint++){
			if(ipoint%5000==0){
				std::cout << std::setprecision(2);
				std::cout << 100.0*ipoint/g->GetN() << "% ";
				std::cout.flush();
			}
			double x,y;
			g->GetPoint(ipoint,x,y);

			double xshits=-1,yshits=-1;
			for(int ipointx=0; ipointx<g->GetN(); ipointx++){
				vGrX2Cent.at(igr)->GetPoint(ipointx,xshits,yshits);
				if(fabs(x-xshits)<1e-5) break;
			}

			double xref=-1,yref=-1;
			for(int ipointy=0; ipointy<g->GetN(); ipointy++){
				vGrY2Cent.at(igr)->GetPoint(ipointy,xref,yref);
				if(fabs(y-xref)<1e-5) break;
			}
			vGrCentCorr.at(igr)->SetPoint(ipoint, yshits, yref);
			vGrDiff.at(igr)->SetPoint(ipoint, yshits, yshits-yref);
		}
		std::cout << std::endl;
	}

	// Final histos (binned into 10 centrality classes)
	const int nDivisions = 10;
	std::vector<TH2D*> vh2diff;
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		std::stringstream ss;
		ss << vgr[igr]->GetName() << "-diff";
		vh2diff.push_back(new TH2D(ss.str().c_str(), ss.str().c_str(), nDivisions,minCentr,maxCentr, 200*nDivisions,-100,100));
		for(int ipoint=0; ipoint<vGrDiff.at(igr)->GetN(); ipoint++){
			double x,y;
			vGrDiff.at(igr)->GetPoint(ipoint,x,y);
			vh2diff.at(igr)->Fill(x,y);
		}
	}
	std::vector<TH2D*> vh2corr;
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		std::stringstream ss;
		ss << vgr[igr]->GetName() << "-corr";
		vh2corr.push_back(new TH2D(ss.str().c_str(), ss.str().c_str(), nDivisions,minCentr,maxCentr, nDivisions,minCentr,maxCentr));
		for(int ipoint=0; ipoint<vGrCentCorr[igr]->GetN(); ipoint++){
			double x,y;
			vGrCentCorr.at(igr)->GetPoint(ipoint,x,y);
			vh2corr.at(igr)->Fill(x,y);
		}
	}
	std::vector<TH1D*> vh1spectr;
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		std::stringstream ss;
		ss << vgr[igr]->GetName() << "-spectr";
		vh1spectr.push_back(new TH1D(ss.str().c_str(), ss.str().c_str(), 1000,0,1000));
		for(int ipoint=0; ipoint<vgr.at(igr)->GetN(); ipoint++){
			double x,y;
			vgr.at(igr)->GetPoint(ipoint,x,y);
			vh1spectr.at(igr)->Fill(x);
		}
	}

	// Calculate centrality resolution
	std::vector<TGraphErrors*> vGrUncertainties;
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		//std::cout << " " << igr << std::endl;
		vGrUncertainties.push_back(new TGraphErrors(nDivisions));
		std::stringstream ss;
		ss << vgr.at(igr)->GetName() << "-uncert-" << igr;
		vGrUncertainties.at(igr)->SetName(ss.str().c_str());
		vGrUncertainties.at(igr)->SetTitle(ss.str().c_str());
		for(int idiv=0; idiv<nDivisions; idiv++){
			//std::cout << "   " << idiv << std::endl;
			std::stringstream ss;
			ss << vh2diff.at(igr)->GetName() << "-diff-" << idiv;

			TH1D* hProjY = vh2diff.at(igr)->ProjectionY((ss.str()+"-pY").c_str(), idiv+1, idiv+2);
			TF1* fFit = new TF1((ss.str()+"-f-gaus").c_str(),"gaus");
			hProjY->Fit(fFit,"QN","");

			float shiftedX = (float)minCentr + (0.5+idiv)*((float)widthCentr/nDivisions);
			vGrUncertainties.at(igr)->SetPoint(idiv, shiftedX, fFit->GetParameter(2));
			vGrUncertainties.at(igr)->SetPointError(idiv, 0, fFit->GetParError(2));
			delete fFit;
		}
	}


	// Draw
	TCanvas *cSort = new TCanvas("cSort","cSort",1800,950);
	cSort->Divide(3,4);
	for(UInt_t igr=0; igr<vGrX2Cent.size(); igr++){
		cSort->cd(igr+1);
		vGrX2Cent.at(igr)->SetTitle(vgr.at(igr)->GetName());
		vGrX2Cent.at(igr)->GetXaxis()->SetTitle("nChargedPart");
		vGrX2Cent.at(igr)->GetYaxis()->SetTitle("Reconstructed centrality [%]");
		vGrX2Cent.at(igr)->GetYaxis()->SetRangeUser(0,1.1*maxCentr);
		vGrX2Cent.at(igr)->Draw("AP");
	}

	TCanvas *cCorr = new TCanvas("cCorr","cCorr",1800,950);
	cCorr->Divide(3,4);
	for(UInt_t igr=0; igr<vh2corr.size(); igr++){
		cCorr->cd(igr+1);
		gPad->SetLogz();
		vh2corr.at(igr)->SetXTitle("Measured centr. (from detector) [%]");
		vh2corr.at(igr)->SetYTitle("Reference centr. (from nPrimCharged) [%]");
		vh2corr.at(igr)->GetZaxis()->SetRangeUser(1,1e4);
		vh2corr.at(igr)->Draw("colz");
	}

	TCanvas *cDiff = new TCanvas("cDiff","cDiff",1800,950);
	cDiff->Divide(3,4);
	for(UInt_t igr=0; igr<vh2diff.size(); igr++){
		cDiff->cd(igr+1);
		gPad->SetLogz();
		vh2diff.at(igr)->SetXTitle("Measured centr. (from detector) [%]");
		vh2diff.at(igr)->SetYTitle("Difference between measured and reference [%]");
		vh2diff.at(igr)->GetYaxis()->SetRangeUser(-10, 10);
		vh2diff.at(igr)->GetZaxis()->SetRangeUser(1, 1e3);
		vh2diff.at(igr)->Draw("colz");
	}

	TCanvas *cUncert = new TCanvas("cUncert","cUncert",1800,950);
	cUncert->Divide(3,4);
	for(UInt_t igr=0; igr<vGrUncertainties.size(); igr++){
		cUncert->cd(igr+1);
		gPad->SetGridx();
		gPad->SetGridy();
		vGrUncertainties.at(igr)->SetMarkerColor(kBlue);
		vGrUncertainties.at(igr)->SetMarkerStyle(8);
		vGrUncertainties.at(igr)->SetMarkerSize(1);
		vGrUncertainties.at(igr)->GetXaxis()->SetTitle("Centrality class");
		vGrUncertainties.at(igr)->GetXaxis()->Set(round(widthCentr/nDivisions), minCentr, maxCentr);
		vGrUncertainties.at(igr)->GetXaxis()->SetNdivisions(nDivisions, kFALSE);
		vGrUncertainties.at(igr)->GetYaxis()->SetTitle("Absolute centrality resolution [%]");
		vGrUncertainties.at(igr)->GetYaxis()->SetRangeUser(0, 3.5);
		for(int ibin=1; ibin<=vGrUncertainties.at(igr)->GetXaxis()->GetNbins(); ibin++){
			std::stringstream ss;
			ss << minCentr + ((float)ibin)/nDivisions*widthCentr << "%";
			vGrUncertainties.at(igr)->GetXaxis()->SetBinLabel(ibin, ss.str().c_str());
		}
		vGrUncertainties.at(igr)->Draw("AP");
	}

	TFile *fout = new TFile("outputs/centrRes.root","RECREATE");
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		vGrX2Cent.at(igr)->Write();
	}
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		vh2corr.at(igr)->Write();
	}
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		vh2diff.at(igr)->Write();
	}
	for(UInt_t igr=0; igr<vgr.size(); igr++){
		vGrUncertainties.at(igr)->Write();
	}
	fout->Close();

	std::cout << "Done." << std::endl;
	theApp.SetIdleTimer(10000,"exit()"); // exit after 30 min of being idle
	theApp.SetReturnFromRun(true);
	theApp.Run();
	return 0;
}

