# In this .py file I will try to write my own version of the Zpeak analysis provided in RDataFrames_Zpeak_alt

import ROOT
import os, argparse
import numpy as np
import sys
import re
import time
import yaml

#########################################################################################
# C++ Helper Functions
#########################################################################################

cpp = R"""

using FourVector = ROOT::Math::PxPyPzMVector;
using RVecFloat = ROOT::RVec<Float_t>;

using FourVectorPtEtaPhiM = ROOT::Math::PtEtaPhiMVector;
using FourVectorPtEtaPhiE = ROOT::Math::PtEtaPhiEVector;
using FourVectorPxPyPzM = ROOT::Math::PxPyPzMVector;
using namespace ROOT::VecOps;

float getLeading(RVecFloat vec){
    auto idxmax = ROOT::VecOps::ArgMax(vec);
    return vec[idxmax];

}

float getTrailing(RVecFloat vec){
    auto idxmin = ROOT::VecOps::ArgMin(vec);
    return vec[idxmin];
}
auto myDiLep(RVecFloat& pt, RVecFloat& eta, RVecFloat& phi, float mass){
    Float_t px1 = pt[0]*cos(phi[0]);Float_t py1 = pt[0]*sin(phi[0]);Float_t pz1 = pt[0]*cos(eta[0]);
    Float_t px2 = pt[1]*cos(phi[1]);Float_t py2 = pt[1]*sin(phi[1]);Float_t pz2 = pt[1]*cos(eta[1]);
    FourVector P1 {px1, py1, pz1, mass};
    FourVector P2 {px2, py2, pz2, mass};
    FourVector dilep = P1 + P2;
    return dilep;
}

auto myLep(RVecFloat& pt, RVecFloat& eta, RVecFloat& phi, float mass, int tag){
    Float_t px1 = pt[0]*cos(phi[0]);Float_t py1 = pt[0]*sin(phi[0]);Float_t pz1 = pt[0]*cos(eta[0]);
    Float_t px2 = pt[1]*cos(phi[1]);Float_t py2 = pt[1]*sin(phi[1]);Float_t pz2 = pt[1]*cos(eta[1]);
    FourVector P1 {px1, py1, pz1, mass};
    FourVector P2 {px2, py2, pz2, mass};
    if (tag == 0)
        return P1;
    else
        return P2;
}

RVecFloat getDeltaR(RVecFloat& eta, RVecFloat& phi){
/** Get all DeltaR from all possible pairs in the data
 * Input should be the relevant eta and phi columns respectively
 * Returns Vec with length equal to amount of combinations where each entry is the DeltaR for a pair
 */
    if (eta.size() != phi.size()) {
        std::cout << "Eta and phi are not the same size!" << std::endl;
        return RVecFloat();
    } else {
        auto idx = ROOT::VecOps::Combinations(eta, 2);
        RVecFloat deltaR(idx[0].size()); // Initialize deltaR with the correct size
                               
        auto phi1 = Take(phi, idx[0]); auto phi2 = Take(phi, idx[1]);
        auto eta1 = Take(eta, idx[0]); auto eta2 = Take(eta, idx[1]);

        deltaR = ROOT::VecOps::DeltaR(eta1, eta2, phi1, phi2);
        return deltaR;
    }

}

Int_t MatchGenbis(FourVectorPtEtaPhiM &gen_lv, RVec<double>& packedCandsPx,RVec<double>& packedCandsPy,RVec<double>& packedCandsPz,Double_t mass,RVec<int>& packedCandsCharge,Int_t genPCharge) {
    double minDR = 100;
    Int_t index = 0;
    double dr_max = 0.03;
    double delR=minDR;
    //std::cout<<"Stuck at index:"<<gen_ind<<std::endl;
    for (int i=0; i<packedCandsPx.size(); i++) {
        FourVectorPxPyPzM recoP {packedCandsPx[i],packedCandsPy[i],packedCandsPz[i],mass}; //packedCandsE[i]
        delR = ROOT::Math::VectorUtil::DeltaR(gen_lv,recoP);
        //if (recoP.Pt() < 5) continue;
        //std::cout<<"GenIndex, Pt, Charge, ID: "<<gen_ind<<", "<<event.genParPt[gen_ind]<<", "<<event.genParCharge[gen_ind]<<", "<<event.genParId[gen_ind]<<";"<<std::endl;
        //std::cout<<"PCdIndex, Pt, Charge, ID: "<<j<<", "<<packedCand.Pt()<<",  "<<event.packedCandsCharge[j]<<", "<<event.packedCandsPdgId[j]<<";"<<std::endl;
        //std::cout<<"Check the dR on this: "<<delR<<std::endl;
        //std::cout<<"GenIndex: "<<gen_ind<<", PackedCandIndex: "<<j<<std::endl;
        if (genPCharge!=packedCandsCharge[i]) continue;
        // if (abs(packedCandsPdgId[i])==13) continue;
        if (minDR < delR) continue;
        minDR = delR;
        index = i;
    }
    //std::cout<<"Index, PDGID, MPDGID: "<<index<<", "<<genParId[index]<<", "<<genParMotherId[index]<<";"<<std::endl;

    if (minDR < dr_max) 
        return index; //abs(genParId[index])
    else
        return -1;
}

Int_t MatchGen(FourVectorPtEtaPhiM &gen_lv, RVec<Int_t>& packedCandsPdgId, RVec<Float_t>& packedCandsPx,RVec<Float_t>& packedCandsPy,RVec<Float_t>& packedCandsPz,Double_t mass,RVec<Int_t>& packedCandsCharge,Int_t genPCharge) {
    double minDR = 100;
    Int_t index = 0;
    double dr_max = 0.03;
    double delR=minDR;
    //std::cout<<"Stuck at index:"<<gen_ind<<std::endl;
    for (int i=0; i<packedCandsPx.size(); i++) {
        if (abs(packedCandsPdgId[i]) != 211) continue;
        FourVectorPxPyPzM recoP {packedCandsPx[i],packedCandsPy[i],packedCandsPz[i],mass}; //packedCandsE[i]
        delR = ROOT::Math::VectorUtil::DeltaR(gen_lv,recoP);
        //if (recoP.Pt() < 5) continue;
        //std::cout<<"GenIndex, Pt, Charge, ID: "<<gen_ind<<", "<<event.genParPt[gen_ind]<<", "<<event.genParCharge[gen_ind]<<", "<<event.genParId[gen_ind]<<";"<<std::endl;
        //std::cout<<"PCdIndex, Pt, Charge, ID: "<<j<<", "<<packedCand.Pt()<<",  "<<event.packedCandsCharge[j]<<", "<<event.packedCandsPdgId[j]<<";"<<std::endl;
        //std::cout<<"Check the dR on this: "<<delR<<std::endl;
        //std::cout<<"GenIndex: "<<gen_ind<<", PackedCandIndex: "<<j<<std::endl;
        if (genPCharge!=packedCandsCharge[i]) continue;
        // if (abs(packedCandsPdgId[i])==13) continue;
        if (minDR < delR) continue;
        minDR = delR;
        index = i;
    }
    //std::cout<<"Index, PDGID, MPDGID: "<<index<<", "<<genParId[index]<<", "<<genParMotherId[index]<<";"<<std::endl;

    if (minDR < dr_max) 
        return index; //abs(genParId[index])
    else
        return -1;
}

Int_t matchReco(FourVectorPxPyPzM &reco_lv,RVec<Float_t>& genParPt,RVec<Float_t>& genParEta,RVec<Float_t>& genParPhi,RVec<Float_t>& genParE,RVec<Int_t>& genParId,RVec<Int_t>& genParMotherId) {
    double minDR = 100;
    unsigned index = 0;
    double dr_max = 0.03;
    double delR=minDR;
    //std::cout<<"Stuck at index:"<<gen_ind<<std::endl;
    for (int i=0; i<genParPt.size(); i++) {
        FourVectorPtEtaPhiE genP {genParPt[i],genParEta[i],genParPhi[i],genParE[i]};
        delR = ROOT::Math::VectorUtil::DeltaR(reco_lv,genP);
        //std::cout<<"GenIndex, Pt, Charge, ID: "<<gen_ind<<", "<<event.genParPt[gen_ind]<<", "<<event.genParCharge[gen_ind]<<", "<<event.genParId[gen_ind]<<";"<<std::endl;
        //std::cout<<"PCdIndex, Pt, Charge, ID: "<<j<<", "<<packedCand.Pt()<<",  "<<event.packedCandsCharge[j]<<", "<<event.packedCandsPdgId[j]<<";"<<std::endl;
        //std::cout<<"Check the dR on this: "<<delR<<std::endl;
        //std::cout<<"GenIndex: "<<gen_ind<<", PackedCandIndex: "<<j<<std::endl;
        if (minDR < delR) continue;
        minDR = delR;
        index = i;
        }
    //std::cout<<"Index, PDGID, MPDGID: "<<index<<", "<<genParId[index]<<", "<<genParMotherId[index]<<";"<<std::endl;
    if (minDR < dr_max) 
        return index; //abs(genParId[index])
    else
        return -1;
}

RVec<Float_t> getKinematics(const RVec<FourVector> &tracks, TString var = "pt"){
    auto pt = [](const FourVector &v) { return v.Pt(); };
    auto px = [](const FourVector &v) { return v.Px(); };
    auto py = [](const FourVector &v) { return v.Py(); };
    auto pz = [](const FourVector &v) { return v.Pz(); };
    auto eta = [](const FourVector &v) { return v.Eta(); };
    auto phi = [](const FourVector &v) { return v.Phi(); };
    if (var.Contains("eta")) return Map(tracks, eta);
    else if (var.Contains("phi")) return Map(tracks, phi);
    else if (var.Contains("px")) return Map(tracks, px);
    else if (var.Contains("py")) return Map(tracks, py);
    else if (var.Contains("pz")) return Map(tracks, pz);
    else return Map(tracks, pt);
}

ROOT::VecOps::RVec<ROOT::Math::PxPyPzMVector> makeLVs(ROOT::VecOps::RVec<Int_t>& packedCandsPdgId,ROOT::VecOps::RVec<Int_t>& packedCandsCharge,ROOT::VecOps::RVec<Int_t>& packedCandsHasTrackDetails,ROOT::VecOps::RVec<Float_t>& packedCandsPx,ROOT::VecOps::RVec<Float_t>& packedCandsPy,ROOT::VecOps::RVec<Float_t>& packedCandsPz, double chsMass_, bool flag = true) {
    ROOT::VecOps::RVec<ROOT::Math::PxPyPzMVector> lvs;
    // std::cout<<"check sizes:"<<packedCandsPdgId.size()<<","<<packedCandsCharge.size()<<","<<packedCandsHasTrackDetails.size()<<","<<packedCandsPx.size()<<","<<packedCandsPy.size()<<","<<packedCandsPz.size()<<std::endl;
    for (Int_t k = 0; k < packedCandsPx.size(); k++) {
        ROOT::Math::PxPyPzMVector lVec {packedCandsPx[k], packedCandsPy[k], packedCandsPz[k], chsMass_};
        lvs.push_back(lVec);
    }
    // std::cout<<"Gets out"<<std::endl;
    return lvs;
}

RVec<unsigned long> getIndices(const int num) {
    RVec<unsigned long> idx;
    for (int i = 0;i<num;i++) idx.emplace_back(i);
    return idx;
}


int checkTrkPair(Int_t leadingidx, Int_t subleadingidx, Float_t num, const RVec<Float_t>& mu1_idx, const RVec<Float_t>& mu2_idx) { 
    int TkPairIdx = -1;
    for (int i{0}; i < num; i++) {
    if ((mu1_idx[i] == leadingidx && mu2_idx[i] == subleadingidx) || (mu1_idx[i] == subleadingidx && mu2_idx[i] == leadingidx)) {
        TkPairIdx = i;
        break;
        }
    }
    int idx = TkPairIdx;
    if (idx < 0){
        return -1;
    } else {
        return TkPairIdx;
    }
}
int getMuonTrackPairIndex(int leadingidx, int subleadingidx, Float_t num, const RVec<Float_t>& mu1_idx, const RVec<Float_t>& mu2_idx,const RVec<Float_t>& mu1_px, const RVec<Float_t>& mu1_py, const RVec<Float_t>& mu1_pz, const RVec<Float_t>& mu2_px, const RVec<Float_t>& mu2_py, const RVec<Float_t>& mu2_pz,  double muonMass_, std::vector<ROOT::Math::PxPyPzMVector> &refit_trks) { 
    int TkPairIdx = -1;
    for (int i{0}; i < num; i++) {
    if ((mu1_idx[i] == leadingidx && mu2_idx[i] == subleadingidx) || (mu1_idx[i] == subleadingidx && mu2_idx[i] == leadingidx)) {
        TkPairIdx = i;
        break;
        }
    }
    int idx = TkPairIdx;
    if (idx < 0) return -1;

    if ( std::isnan(mu1_px[idx])  || std::isnan(mu2_px[idx]) ) return -1;
    if ( std::isnan(mu1_py[idx])  || std::isnan(mu2_py[idx]) ) return -1;
    if ( std::isnan(mu1_pz[idx]) || std::isnan(mu2_pz[idx]) ) return -1;

    ROOT::Math::PxPyPzMVector muTrk1{mu1_px[idx], mu1_py[idx], mu1_pz[idx], muonMass_};
    ROOT::Math::PxPyPzMVector muTrk2{mu2_px[idx], mu2_py[idx], mu2_pz[idx], muonMass_};
    if (mu1_idx[idx] == leadingidx && mu2_idx[idx] == subleadingidx) {
        refit_trks[0] = muTrk1;
        refit_trks[1] = muTrk2;
    }
    else {
        refit_trks[0]  = muTrk2;
        refit_trks[1]  = muTrk1;
    }
    return idx;
}


std::vector<Int_t> getDileptonCand2(const RVec<unsigned long>& mu_pt_sorted_idx, const RVec<Int_t>& mu_ch, const RVec<Float_t>& mu_px, const RVec<Float_t>& mu_py,const RVec<Float_t>& mu_pz, const Float_t mass, Float_t num, const RVec<Float_t>& mu1_idx, const RVec<Float_t>& mu2_idx, const RVec<Float_t>& mu1_px, const RVec<Float_t>& mu1_py, const RVec<Float_t>& mu1_pz, const RVec<Float_t>& mu2_px, const RVec<Float_t>& mu2_py, const RVec<Float_t>& mu2_pz, TString ptype, const RVec<Float_t>& packedCandsPx, const RVec<Float_t>& packedCandsPy, const RVec<Float_t>& packedCandsPz, const RVec<Float_t>& packedCandsE, const RVec<Int_t>& packedCandsCharge, const RVec<Int_t>& packedCandsPdgId, const RVec<Int_t>& packedCandsFromPV, Int_t numPackedCands, ROOT::Math::PxPyPzMVector& leadinglv, ROOT::Math::PxPyPzMVector& subleadinglv, Float_t diMuonPt_=0, double dr_max = 0.4) {

    std::vector<Int_t> objidx(3);
    objidx[0]=-1;
    objidx[1]=-1;
    objidx[2]=-1;
    float maxDileptonDeltaR_ = 0.4;
    for (int i=0;i<mu_pt_sorted_idx.size();i++) {
        for (int j=i+1;j<mu_pt_sorted_idx.size();j++) { //it is Pt sorted array
            if (mu_ch[mu_pt_sorted_idx[i]] * mu_ch[mu_pt_sorted_idx[j]] >= 0) continue;
            FourVector lepton1{mu_px[mu_pt_sorted_idx[i]], mu_py[mu_pt_sorted_idx[i]], mu_pz[mu_pt_sorted_idx[i]], mass};
            FourVector lepton2{mu_px[mu_pt_sorted_idx[j]], mu_py[mu_pt_sorted_idx[j]], mu_pz[mu_pt_sorted_idx[j]], mass};
            double delR = ROOT::Math::VectorUtil::DeltaR(lepton1,lepton2);
            if ( delR > maxDileptonDeltaR_ ) continue;
            std::vector<FourVector> refittedTrks(2);
            int TkPairIdx = getMuonTrackPairIndex(mu_pt_sorted_idx[i], mu_pt_sorted_idx[j], num, mu1_idx, mu2_idx,mu1_px, mu1_py,mu1_pz, mu2_px, mu2_py, mu2_pz,  mass, refittedTrks);
            if (TkPairIdx < 0) continue; //checking if refitted tracks are present -> necessary for SV construction
            

            double pT { (lepton1+lepton2).Pt() };
            if (pT < diMuonPt_) continue;
            
            objidx[0] = mu_pt_sorted_idx[i];
            objidx[1] = mu_pt_sorted_idx[j];
            objidx[2] = TkPairIdx;
            leadinglv = lepton1;
            subleadinglv = lepton2;
            return objidx;
        }
    }
    return objidx;
}

std::vector<Int_t> getDileptonCand(const RVec<unsigned long>& mu_pt_sorted_idx, const RVec<Int_t>& mu_ch, const RVec<Float_t>& mu_px, const RVec<Float_t>& mu_py,const RVec<Float_t>& mu_pz, const Float_t mass, Float_t num, const RVec<Float_t>& mu1_idx, const RVec<Float_t>& mu2_idx, const RVec<Float_t>& mu1_px, const RVec<Float_t>& mu1_py, const RVec<Float_t>& mu1_pz, const RVec<Float_t>& mu2_px, const RVec<Float_t>& mu2_py, const RVec<Float_t>& mu2_pz, TString ptype, const RVec<Float_t>& packedCandsPx, const RVec<Float_t>& packedCandsPy, const RVec<Float_t>& packedCandsPz, const RVec<Float_t>& packedCandsE, const RVec<Int_t>& packedCandsCharge, const RVec<Int_t>& packedCandsPdgId, const RVec<Int_t>& packedCandsFromPV, Int_t numPackedCands, ROOT::Math::PxPyPzMVector& leadinglv, ROOT::Math::PxPyPzMVector& subleadinglv, Float_t diMuonPt_=0, double dr_max = 0.4) {

    std::vector<Int_t> objidx(3);
    objidx[0]=-1;
    objidx[1]=-1;
    objidx[2]=-1;
    float maxDileptonDeltaR_ = 0.4;
    for (int i=0;i<mu_pt_sorted_idx.size();i++) {
        for (int j=i+1;j<mu_pt_sorted_idx.size();j++) { //it is Pt sorted array
            if (mu_ch[mu_pt_sorted_idx[i]] * mu_ch[mu_pt_sorted_idx[j]] >= 0) continue;
            FourVector lepton1{mu_px[mu_pt_sorted_idx[i]], mu_py[mu_pt_sorted_idx[i]], mu_pz[mu_pt_sorted_idx[i]], mass};
            FourVector lepton2{mu_px[mu_pt_sorted_idx[j]], mu_py[mu_pt_sorted_idx[j]], mu_pz[mu_pt_sorted_idx[j]], mass};
            double delR = ROOT::Math::VectorUtil::DeltaR(lepton1,lepton2);
            if ( delR > maxDileptonDeltaR_ ) continue;
            std::vector<FourVector> refittedTrks(2);
            int TkPairIdx = getMuonTrackPairIndex(mu_pt_sorted_idx[i], mu_pt_sorted_idx[j], num, mu1_idx, mu2_idx,mu1_px, mu1_py,mu1_pz, mu2_px, mu2_py, mu2_pz,  mass, refittedTrks);
            if (TkPairIdx < 0) continue; //checking if refitted tracks are present -> necessary for SV construction

            int leadingidx = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? mu_pt_sorted_idx[i] : mu_pt_sorted_idx[j];   
            int subleadingidx = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? mu_pt_sorted_idx[j] : mu_pt_sorted_idx[i];
            FourVector lepton1_refit = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? refittedTrks[0] : refittedTrks[1];
            FourVector lepton2_refit = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? refittedTrks[1] : refittedTrks[0];
            double pT { (lepton1_refit+lepton2_refit).Pt() };
            if (pT < diMuonPt_) continue;
            
            objidx[0] = leadingidx;
            objidx[1] = subleadingidx;
            objidx[2] = TkPairIdx;
            leadinglv = lepton1_refit;
            subleadinglv = lepton2_refit;
            return objidx;
        }
    }
    return objidx;
}

std::vector<Int_t> getDileptonCand1(const RVec<unsigned long>& mu_pt_sorted_idx, const RVec<Int_t>& mu_ch, const RVec<Float_t>& mu_px, const RVec<Float_t>& mu_py,const RVec<Float_t>& mu_pz, const Float_t mass, Float_t num, const RVec<Float_t>& mu1_idx, const RVec<Float_t>& mu2_idx, const RVec<Float_t>& mu1_px, const RVec<Float_t>& mu1_py, const RVec<Float_t>& mu1_pz, const RVec<Float_t>& mu2_px, const RVec<Float_t>& mu2_py, const RVec<Float_t>& mu2_pz, TString ptype, const RVec<Float_t>& packedCandsPx, const RVec<Float_t>& packedCandsPy, const RVec<Float_t>& packedCandsPz, const RVec<Float_t>& packedCandsE, const RVec<Int_t>& packedCandsCharge, const RVec<Int_t>& packedCandsPdgId, const RVec<Int_t>& packedCandsFromPV, Int_t numPackedCands, ROOT::Math::PxPyPzMVector& leadinglv, ROOT::Math::PxPyPzMVector& subleadinglv, Float_t diMuonPt_=0, double dr_max = 0.4) {

    std::vector<Int_t> objidx(3);
    std::vector<std::tuple<int, int, double>> delR_values;
    objidx[0]=-1;
    objidx[1]=-1;
    objidx[2]=-1;

    //iterate over all (i,j) combinations to compute delR
    for (int i=0;i<mu_pt_sorted_idx.size();i++) {
        for (int j=i+1;j<mu_pt_sorted_idx.size();j++) { //it is Pt sorted array
            if (mu_ch[mu_pt_sorted_idx[i]] * mu_ch[mu_pt_sorted_idx[j]] >= 0) continue;

            FourVector lepton1{mu_px[mu_pt_sorted_idx[i]], mu_py[mu_pt_sorted_idx[i]], mu_pz[mu_pt_sorted_idx[i]], mass};
            FourVector lepton2{mu_px[mu_pt_sorted_idx[j]], mu_py[mu_pt_sorted_idx[j]], mu_pz[mu_pt_sorted_idx[j]], mass};
            double delR = ROOT::Math::VectorUtil::DeltaR(lepton1,lepton2);
            
            //Set up pairs alongside delR
            delR_values.emplace_back(i, j, delR);
        }
    }
    
    //sort according to ascending delR
    std::sort(delR_values.begin(), delR_values.end(), [](const auto &a, const auto &b) {
        return std::get<2>(a) < std::get<2>(b);
    });
    
    //iterate over acending delR pairs to check other conditions
    for (int k=0; k < delR_values.size(); k++){
        int i = std::get<0>(delR_values[k]);
        int j = std::get<1>(delR_values[k]);
    
        std::vector<FourVector> refittedTrks(2);
        int TkPairIdx = getMuonTrackPairIndex(mu_pt_sorted_idx[i], mu_pt_sorted_idx[j], num, mu1_idx, mu2_idx,mu1_px, mu1_py,mu1_pz, mu2_px, mu2_py, mu2_pz,  mass, refittedTrks);
        if (TkPairIdx < 0) continue; //checking if refitted tracks are present -> necessary for SV construction

        int leadingidx = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? mu_pt_sorted_idx[i] : mu_pt_sorted_idx[j];   
        int subleadingidx = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? mu_pt_sorted_idx[j] : mu_pt_sorted_idx[i];
        FourVector lepton1_refit = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? refittedTrks[0] : refittedTrks[1];
        FourVector lepton2_refit = refittedTrks[0].Pt() > refittedTrks[1].Pt() ? refittedTrks[1] : refittedTrks[0];
        double pT { (lepton1_refit+lepton2_refit).Pt() };
        if (pT < diMuonPt_) continue;
     
        //If all conditions passed return out of function, hence lowest delR now favoured
        objidx[0] = leadingidx;
        objidx[1] = subleadingidx;
        objidx[2] = TkPairIdx;
        leadinglv = lepton1_refit;
        subleadinglv = lepton2_refit;
        return objidx;

    }
    //if no condition is passed in previous for loop objidx will be unaltered [-1,-1,-1]
    return objidx;
}
"""


#########################################################################################
# Setting up the environment for loading in the data from configs
#########################################################################################

ROOT.gInterpreter.Declare(cpp)
start_time = time.time()
cpu_count = 16 # give the same for request_cpus on condor script
#ROOT.ROOT.EnableImplicitMT(cpu_count)

#Setting up the loading in of the appropriate dataset along with the right cuts, weights, etc...
parser = argparse.ArgumentParser(description='Analysis of Zpeak using Z to muons')
parser.add_argument("-c", "--config", dest="config",   help="Enter config file to process", type=str)
parser.add_argument("-o","--output", dest="output", help="Destination directory", type=str)
parser.add_argument("-y", "--year", help="Enter the year of the dataset you want to analyse", type=str)
# For now the cuts will be defined in this .py file 
# parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)

args = parser.parse_args()
#Open the config file and set a variable to be dictionary corresponding to content of the config file
fconfig = open(args.config, 'r')
conf_pars = yaml.safe_load(fconfig)

#Get the cuts
fcuts = open(conf_pars['cuts'])
cuts_pars = yaml.safe_load(fcuts)

#Get data location from the config file (in my config file location is directory path where multiple .root are stored)
location = conf_pars['locations']
directory = location if (location[-1] == '/') else location + '/'
listOfFiles = [os.path.join(directory, file) for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file))]

#This is the  luminosity for the total 2017UL run (see file name). 41474
#These values can be found in the config file and are idealy taken from here in an automated way dependent on which sampleset is called in the commandline to analyze.
#givenLuminosity = conf_pars['luminosity'][args.year]

#Get the weights, cross section, and luminosity from MonteCarlo. Set to 1 if Data is not MC (Non_MC will always contain 'Run' in name?)
sumWeights = 1 if 'Run' in args.config else conf_pars['sum_weights']
crossSection = 1 if 'Run' in args.config  else conf_pars['cross_section']
luminosity = 1 if 'Run' in args.config else 1

#Lets perform a check
if not 'Run' in args.config:
    f = ROOT.TFile(listOfFiles[0])
    weightPlot = f.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
    weightPlot.SetDirectory(0)
    f.Close()
    for fnum, fstr in enumerate(listOfFiles):
        print(fstr)
        if fnum == 0:
            continue
        else:
            f = ROOT.TFile(fstr)
            tmpPlot = f.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
            weightPlot.Add(tmpPlot)
            f.Close()
    sumWeightsPlot = weightPlot.GetBinContent(2) - weightPlot.GetBinContent(3) # bins filled from 1, but bins available from 0
    sumWeightsPlotbis = weightPlot.GetBinContent(1)
    print(sumWeightsPlot, sumWeightsPlotbis)

    if sumWeights != sumWeightsPlot:
        sumWeights = sumWeightsPlot
        print("\nSum of weights in config file did not match actual sum of weights, change the value in the config file accordingly")
        print("\ncorrect sum of weights is: {} ".format(sumWeights))
    else:
        print("\nSum of weights in config matches sum of weights calculated from simulation")
        print("\ncorrect sum of weights is: {} ".format(sumWeights))


#########################################################################################
# Helper Functions
#########################################################################################

ROOT.gStyle.SetOptStat(111111)


def genTurnOn(df, dataTag, triggerName, mmin, mmax, steps):
    dataStr = 'Data_' + args.year if dataTag else 'MC'
    #Filter out the muons from the dataset before applying trigger to get a better gauge on trigger efficiency
    #If this is not done, max efficiency plateau at 50% due to presence of electrons
    df_muons = df.Filter('_nEles == 0')
    df_muontrigger = df_muons.Filter(triggerName)
    
    #Get the leading muon pT since this is most likely the muon the event trigger triggered on
    df_muons = df_muons.Define('leadingPt', 'getLeading(_lPt)')\
            .Define('subleadingPt', 'getTrailing(_lPt)')
    df_muontrigger = df_muontrigger.Define('leadingPt', 'getLeading(_lPt)')\
            .Define('subleadingPt', 'getTrailing(_lPt)')
    
    pTBins = np.arange(mmin,mmax,steps)

    eff_hist = df_muontrigger.Histo1D(('hist_triggerleading', 'Triggered muons', len(pTBins)-1, pTBins), 'leadingPt')
    eff_hist_total = df_muons.Histo1D(('hist_totalleading', 'Total  muons', len(pTBins)-1, pTBins), 'leadingPt')

    turnon = ROOT.TEfficiency(eff_hist.GetValue(), eff_hist_total.GetValue())
    turnon.SetName("turnon_" + triggerName + '_' + dataStr)
    turnon.SetTitle("Turn on curve " + triggerName + dataStr + ";Muon pT (GeV);Efficiency")
    
    #Print the maximum reached efficiency
    max_efficiency = 0
    for i in range(1, turnon.GetTotalHistogram().GetNbinsX() + 1):
        if turnon.GetEfficiency(i) > max_efficiency:
                max_efficiency = turnon.GetEfficiency(i)
    print(triggerName + " Maximum Efficiency:" + str(max_efficiency))
    
    turnon.Write()

def genPtPlot(df, datatag, mmin, mmax, bins):
    dataStr = 'Data_' + args.year if dataTag else 'MC_'

    df_leading_subleading = df.Define('leadingPt', 'getLeading(mu_pt)')\
        .Define('subleadingPt', 'getTrailing(mu_pt)')
    # Histograms
    if dataTag:

        mu_leadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_leading_' + dataStr, 'Leading Muon PTs ' + dataStr, bins, mmin, mmax), 'leadingPt')
        mu_subleadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_subleading_' + dataStr, 'Subleading Muon PT ' + dataStr, bins, mmin, mmax), 'subleadingPt')

        mu_leadingPt_hist.GetXaxis().SetTitle("Muon pT (GeV)")  # X-axis label
        mu_leadingPt_hist.GetYaxis().SetTitle("Events")         # Y-axis label
        mu_leadingPt_hist.SetLineColor(ROOT.kBlue)  # Set color for the leading pt histogram
        mu_subleadingPt_hist.SetLineColor(ROOT.kRed)  # Set color for the subleading pt histogram

        mu_leadingPt_hist.Write()
        mu_subleadingPt_hist.Write()
    else:
        for year in years:
            mu_leadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_leading_' + dataStr + year, 'Leading Muon PTs ' + dataStr, bins, mmin, mmax), 'leadingPt', 'evt_weight_' + year)
            mu_subleadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_subleading_' + dataStr + year, 'Subleading Muon PT ' + dataStr, bins, mmin, mmax), 'subleadingPt', 'evt_weight_' + year)

            mu_leadingPt_hist.GetXaxis().SetTitle("Muon pT (GeV)")  # X-axis label
            mu_leadingPt_hist.GetYaxis().SetTitle("Events")         # Y-axis label
            mu_leadingPt_hist.SetLineColor(ROOT.kBlue)  # Set color for the leading pt histogram
            mu_subleadingPt_hist.SetLineColor(ROOT.kRed)  # Set color for the subleading pt histogram

            mu_leadingPt_hist.Write()
            mu_subleadingPt_hist.Write()


def genXPlot(df, dataTag, varName, mmin, mmax, bins):
    dataStr = 'Data_'+ args.year if dataTag else 'MC_'
    if dataTag:
        hist = df.Histo1D(('hist_' + varName + '_' + dataStr, varName + dataStr , bins, mmin, mmax), varName)
        hist.GetXaxis().SetTitle(varName)
        hist.GetYaxis().SetTitle("Events")
        hist.Write()
    else:
        for year in years:
            hist = df.Histo1D(('hist_' + varName + '_' + dataStr + year, varName + dataStr , bins, mmin, mmax), varName, 'evt_weight_' + year)
            hist.GetXaxis().SetTitle(varName)
            hist.GetYaxis().SetTitle("Events")
            hist.Write()


def genInvMassPlot(df, dataTag, mmin, mmax, bins):
    dataStr = 'Data_' + args.year if dataTag else 'MC_'
    if dataTag:
        hist = df.Histo1D(('hist_mu_invmass_' + dataStr, 'Invariant mass distribution of Z', bins, mmin, mmax), 'dimuon_mass')
        hist.GetXaxis().SetTitle('m#_{\mu\mu} (GeV)')
        hist.GetYaxis().SetTitle('Events')
        hist.Write()

    else:
        for year in years:
            hist = df.Histo1D(('hist_mu_invmass_' + dataStr + year, 'Invariant mass distribution of Z', bins, mmin, mmax), 'dimuon_mass', 'evt_weight_' + year)
            hist.GetXaxis().SetTitle('m#_{\mu\mu} (GeV)')
            hist.GetYaxis().SetTitle('Events')
            hist.Write()

def genXHist(df, varName,  mmin, mmax, bins):
    if 'invmass' in varName:
        if 'Z' in varName:
            nameStr = 'm#_{\mu\mu} (GeV)'
        else:
            nameStr = 'm_{hh} (GeV)'
    elif 'pT' in varName:
        nameStr = 'pT (GeV)'
    else:
        nameStr = varName
    dataStr = 'MC'
    hist = df.Histo1D(('hist_' + varName + '_' + dataStr, varName + dataStr , bins, mmin, mmax), varName)
    hist.GetXaxis().SetTitle(nameStr)
    hist.GetYaxis().SetTitle("Events")
    hist.Write()

def gen2DHist(df, varName1, varName2, mmin, mmax, bins, nameStr1, nameStr2):
    dataStr = 'MC'
    hist = df.Histo2D(('hist2D_' + varName1+'|'+varName2 + '_' + dataStr, '2D Hist'+ dataStr + ';' +varName1 + ';' + varName2, bins, mmin, mmax, bins, mmin, mmax), varName1, varName2)
    hist.GetXaxis().SetTitle(nameStr1)
    hist.GetYaxis().SetTitle(nameStr2)
    hist.Write()

#########################################################################################
# Loading in the Tree and performing selection on Dataframe
#########################################################################################

#Load in the tree from the data in the .root file
treeName = "makeTopologyNtupleMiniAOD/tree"

df = ROOT.RDataFrame(treeName, listOfFiles)
print('\nTree loaded in succesfully')

totalEntries = df.Count().GetValue()
print(f"\nTotal Entries : {totalEntries}")

df.Describe().Print()
#displayList = ['genParMotherId', 'genParDaughterId1', 'genParDaughterId2', 'muonPF2PATPt', 'numMuonPF2PAT']
#displayList1 = ['genParDaughterId1', 'genParDaughterId2', 'genParDaughter1Index', 'genParDaughter2Index']
#df.Display(displayList1, 10).Print()


##################
#######Cuts#######
##################

mu_cuts = cuts_pars['muons']
ch_cuts = cuts_pars['hadrons']
dihadron_cuts = cuts_pars['dihadrons']

diChPt_ = str(dihadron_cuts['pt'])
diChdR_ = str(dihadron_cuts['dR'])

#############################
#######Generated Muons#######
#############################

#Define muon pt of Z muons and all muons
isMu = 'abs(genParId) == 13 && abs(genParMotherId) != 13'
isMuNotZ = 'abs(genParId) == 13 && genParMotherId != 23 && abs(genParMotherId) != 13'
isMuZ = 'genParId == 13 && genParMotherId == 23'
isAntiMuZ = 'genParId == -13 && genParMotherId == 23'
isMuZComb = 'abs(genParId) == 13 && genParMotherId == 23'

df_genMu = df.Define('genMu_pT', f'genParPt[{isMu}]')\
        .Define('genMu_size', 'genMu_pT.size()')\
        .Define('genMuNotZ_pT', f'genParPt[{isMuNotZ}]')\
        .Define('genMuNotZ_size', 'genMuNotZ_pT.size()')

df_genMuZ = df.Define('genNegMuZ_pT', f'genParPt[{isMuZ}]')\
        .Define('genPosMuZ_pT', f'genParPt[{isAntiMuZ}]')
#Set up kinematics for muonZ sorted according to descending pT
df_genMuZ = df_genMuZ.Define('genMuZ_pT_unsorted',f'genParPt[{isMuZComb}]')\
        .Define('sort_indices', 'ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(genMuZ_pT_unsorted))')\
        .Define('genMuZ_pT', 'Take(genMuZ_pT_unsorted, sort_indices)')\
        .Define('genMuZ_sizeZ', 'genMuZ_pT.size()')\
        .Define('genMuZ_eta', f'Take(genParEta[{isMuZComb}], sort_indices)')\
        .Define('genMuZ_phi', f'Take(genParPhi[{isMuZComb}], sort_indices)')

#Define leading and subleading muon Z pt
df_genMuZ = df_genMuZ.Define('genMuZ_pT_leading', 'genNegMuZ_pT[0] > genPosMuZ_pT[0] ? genNegMuZ_pT[0]:genPosMuZ_pT[0]')\
        .Define('genMuZ_pT_subleading', 'genNegMuZ_pT[0] > genPosMuZ_pT[0] ? genPosMuZ_pT[0]:genNegMuZ_pT[0]')



#Filter out the events with 2 or more genMuons from Z
print('Total events', df.Count().GetValue())

print('Events with minimum 1 muons', df_genMu.Filter('genMu_size >= 1').Count().GetValue())
print('Events with minimum 1 muons and muons are from Z', df_genMuZ.Filter('genMuZ_sizeZ >= 1').Count().GetValue())

print('Events with minimum 2 muons', df_genMu.Filter('genMu_size >= 2').Count().GetValue())
print('Events with minimum 2 muons and muons are from Z', df_genMuZ.Filter('genMuZ_sizeZ >= 2').Count().GetValue())
print('Events with minimum 2 muons not from a Z', df_genMu.Filter('genMuNotZ_size >= 2').Count().GetValue())

print('Events with exactly 2 muons from Z', df_genMuZ.Filter('genMuZ_sizeZ == 2').Count().GetValue())
print('\n genparMotherId of all muons')
df_genMu.Filter('genMu_size >= 2')\
        .Define('genMuMotherId', f'genParMotherId[{isMu}]').Define('genMuId', f'genParId[{isMu}]')\
        .Display(['genMuId', 'genMuMotherId'], 50).Print()


df_genMuZ = df_genMuZ.Filter('genMuZ_sizeZ >= 2')


#Set up the invariant Z mass
muon_massVal = 0.105 #GeV
df_genMuZ = df_genMuZ.Define('mu_mass', str(muon_massVal))\
        .Define('genZ_invmass','myDiLep(genMuZ_pT, genMuZ_eta, genMuZ_phi, mu_mass).M()')\
        .Define('genZ_pT', 'myDiLep(genMuZ_pT, genMuZ_eta, genMuZ_phi, mu_mass).Pt()')\
        .Define('genZ_deltaR', 'getDeltaR(genMuZ_eta, genMuZ_phi)')
df_genMuZ_invcut = df_genMuZ.Filter('genZ_invmass < 110 && genZ_invmass > 70')

genDisplayList1 = ['genMu_size', 'genMu_pT']
genDisplayList2 = ['genMuZ_sizeZ', 'genMuZ_pT', 'genMuZ_pT_leading', 'genMuZ_pT_subleading']
genDisplayList3 = ['genZ_invmass', 'genZ_pT']
#df_genMu.Display(genDisplayList1).Print()
#df_genMuZ.Display(genDisplayList2).Print()
#df_genMuZ.Display(genDisplayList3).Print()

####################################
#######Generated Scalar to hh#######
####################################

isScalar = 'abs(genParId) == 9000006'
isK = 'abs(genParId) == 321'
isKfromS1 = 'abs(genParId) == 321 && genParMotherId == 9000006'
isKfromS2 = 'abs(genParId) == 321 && genParMotherId == -9000006'

#Check genId
df_genScalar = df_genMuZ.Define('genKId', f'genParId[{isK}]').Define('genK_size', 'genKId.size()')\
        .Define('genKS1Id', f'genParId[{isKfromS1}]')\
        .Define('genS1', f'genParMotherId[{isKfromS1}]')\
        .Define('genKS2Id', f'genParId[{isKfromS2}]')\
        .Define('genS2', f'genParMotherId[{isKfromS2}]')\
        .Define('genKS1_size', 'genKS1Id.size()').Define('genKS2_size', 'genKS2Id.size()')

genSDisplayList = ['genKS1Id', 'genS1', 'genKS2Id', 'genS2']
#df_genScalar.Display(genSDisplayList, 20).Print()

#Setup gen Kinematics
df_genScalar = df_genScalar.Define('genKS1_pT_unsorted', f'genParPt[{isKfromS1}]').Define('genKS2_pT_unsorted', f'genParPt[{isKfromS2}]')\
        .Define('indicesKS1', 'ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(genKS1_pT_unsorted))')\
        .Define('indicesKS2', 'ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(genKS2_pT_unsorted))')\
        .Define('genKS1_pT', 'Take(genKS1_pT_unsorted, indicesKS1)').Define('genKS2_pT', 'Take(genKS2_pT_unsorted, indicesKS2)')\
        .Define('genKS1_eta', f'Take(genParEta[{isKfromS1}], indicesKS1)').Define('genKS2_eta', f'Take(genParEta[{isKfromS2}], indicesKS2)')\
        .Define('genKS1_phi', f'Take(genParPhi[{isKfromS1}], indicesKS1)').Define('genKS2_phi', f'Take(genParPhi[{isKfromS2}], indicesKS2)')\
        .Define('genKS1_pT_leading', 'genKS1_pT[0]').Define('genKS2_pT_leading', 'genKS2_pT[0]')\
        .Define('genKS1_pT_subleading', 'genKS1_pT[1]').Define('genKS2_pT_subleading', 'genKS2_pT[1]')

genSDisplayList1 = ['genKS1_pT', 'genKS2_pT', 'genKS1_eta', 'genKS2_eta', 'genKS1_phi', 'genKS2_phi']
#df_genScalar.Display(genSDisplayList1, 20).Print()
#df_genScalar.Display(['genKS1_pT_leading', 'genKS1_pT_subleading', 'genKS2_pT_leading', 'genKS2_pT_subleading'],20).Print()

#Filter out the events which have at least 2 hadrons from both scalars
df_genScalar.Filter('genKS1_size >= 2 && genKS2_size >= 2')
print('scalar to hh events and Z to mumu in event: ', df_genScalar.Count().GetValue())

Kmass_val = 0.493677 #GeV
Kmass_val = float(Kmass_val)

df_genScalar = df_genScalar.Define('genKS1_dPhi', 'ROOT::VecOps::DeltaPhi(genKS1_phi[0], genKS1_phi[1])')\
        .Define('genKS2_dPhi', 'ROOT::VecOps::DeltaPhi(genKS2_phi[0], genKS2_phi[1])')\
        .Define('genKS1_dR', 'ROOT::VecOps::DeltaR(genKS1_eta[0], genKS1_eta[1], genKS1_phi[0], genKS1_phi[1])')\
        .Define('genKS2_dR', 'ROOT::VecOps::DeltaR(genKS2_eta[0], genKS2_eta[1], genKS2_phi[0], genKS2_phi[1])')


#Setup invariant mass
df_genScalar = df_genScalar.Define('K_mass', str(Kmass_val))\
        .Define('genK1_lv', f'ROOT::Math::PtEtaPhiMVector(genKS1_pT[0], genKS1_eta[0], genKS1_phi[0], {Kmass_val})')\
        .Define('genK2_lv', f'ROOT::Math::PtEtaPhiMVector(genKS1_pT[1], genKS1_eta[1], genKS1_phi[1], {Kmass_val})')\
        .Define('genK3_lv', f'ROOT::Math::PtEtaPhiMVector(genKS2_pT[0], genKS2_eta[0], genKS2_phi[0], {Kmass_val})')\
        .Define('genK4_lv', f'ROOT::Math::PtEtaPhiMVector(genKS2_pT[1], genKS2_eta[1], genKS2_phi[1], {Kmass_val})')\
        .Define('genS1_lv', 'genK1_lv + genK2_lv')\
        .Define('genS2_lv', 'genK3_lv + genK4_lv')\
        .Define('genS1_invmass1', 'genS1_lv.M()')\
        .Define('genS2_invmass1', 'genS2_lv.M()')\
        .Define('genS1_invmass', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, K_mass).M()')\
        .Define('genS1_pT', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, K_mass).Pt()')\
        .Define('genS1_eta', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, K_mass).Eta()')\
        .Define('genS1_phi', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, K_mass).Phi()')\
        .Define('genS2_invmass', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, K_mass).M()')\
        .Define('genS2_pT', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, K_mass).Pt()')\
        .Define('genS2_eta', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, K_mass).Eta()')\
        .Define('genS2_phi', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, K_mass).Phi()')\
        .Define('genS12_dPhi', 'ROOT::VecOps::DeltaPhi(genS1_phi, genS2_phi)')\
        .Define('genS12_dR', 'ROOT::VecOps::DeltaR(genS1_eta, genS2_eta, genS1_phi, genS2_phi)')

genSDisplayList2 = ['genS1_invmass', 'genS1_pT', 'genS2_invmass', 'genS2_pT']
#df_genScalar.Display(genSDisplayList2, 20).Print()


#############################
#####Reconstructed Muons#####
#############################

df_recoMu = df_genMuZ.Filter('numMuonPF2PAT >= 2', '2 or more recoMu')
recoMu_qualitycut = 'abs(muonPF2PATEta) < ' + str(mu_cuts['eta']) + ' && muonPF2PATLooseCutId &&  muonPF2PATPt > ' + str(mu_cuts['pt'])
recoMu_qualitycutEta = 'abs(muonPF2PATEta) < ' + str(mu_cuts['eta'])
recoMu_qualitycutId = 'muonPF2PATLooseCutId'
recoMu_qualitycutPt = 'muonPF2PATPt > ' + str(mu_cuts['pt'])
#Take first two entries which correspond to highest pT muons, and some quality cuts on the muons
df_recoMu = df_recoMu.Define("indices", "ROOT::RVec<int>{0, 1}")\
        .Define('recoMu_pT', f'Take(muonPF2PATPt[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_eta', f'Take(muonPF2PATEta[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_phi', f'Take(muonPF2PATPhi[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_ch', f'Take(muonPF2PATCharge[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_RelIsoLoose', f'Take(muonPF2PATPfIsoLoose[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_RelIsoTight', f'Take(muonPF2PATPfIsoTight[{recoMu_qualitycut}], indices)')
print('PfRelIsoLoose pass: ', df_recoMu.Filter('recoMu_RelIsoLoose[0] == 1 && recoMu_RelIsoLoose[1] == 1').Count().GetValue())
print('PfRelIsoTight pass: ', df_recoMu.Filter('recoMu_RelIsoTight[0] == 1 && recoMu_RelIsoTight[1] == 1').Count().GetValue())
#print('TEST')
#df_recoMu.Display(['muonPF2PATPt', 'recoMu_pT', 'muonPF2PATLooseCutId', 'muonPF2PATTightCutId'], 20).Print()
#The definition of recoMu_pT makes it that muons not passing the conditions have their pt set to 0
df_recoMu = df_recoMu.Filter('Min(recoMu_pT) > ' + str(mu_cuts['pt']), 'Quality filter')
df_recoMu = df_recoMu.Define('recoMu_pT_leading', 'recoMu_pT[0]').Define('recoMu_pT_subleading', 'recoMu_pT[1]')
#print(df_recoMu.Count().GetValue())

#Only keep events with leading pT above cut and opposite charged muons
recoMu_cut = 'Max(recoMu_pT) > ' + str(mu_cuts['leadingPt']) + ' && recoMu_ch[0]*recoMu_ch[1] < 0'
df_recoMu_cut = df_recoMu.Filter(recoMu_cut, 'Muon cuts')


recoDisplayList = ['recoMu_pT', 'recoMu_eta', 'recoMu_phi', 'recoMu_ch']
recoDisplayList1 = ['recoMu_pT', 'recoMu_pT_leading', 'recoMu_pT_subleading']
#print('UNCUT')
#df_recoMu.Display(recoDisplayList).Print()
#df_recoMu.Display(recoDisplayList1).Print()
#print('CUT')
#df_recoMu_cut.Display(recoDisplayList).Print()
#df_recoMu_cut.Display(recoDisplayList1).Print()

#Set up the invariant Z mass
df_recoMu_cut = df_recoMu_cut.Define('recoZ_invmass', 'myDiLep(recoMu_pT, recoMu_eta, recoMu_phi, mu_mass).M()')\
        .Define('recoZ_pT', 'myDiLep(recoMu_pT, recoMu_eta, recoMu_phi, mu_mass).Pt()')\
        .Define('recoZ_deltaR', 'getDeltaR(recoMu_eta, recoMu_phi)')

recoDisplayList2 = ['recoZ_invmass', 'recoZ_pT']
#df_recoMu_cut.Display(recoDisplayList2).Print()

invMassCut = 'recoZ_invmass < 110 && recoZ_invmass > 70'
df_recoMu_cut_invcut = df_recoMu_cut.Filter(invMassCut, 'Invariant mass cuts')
print('CUTFLOWREPORT MUONS BELOW')
df_recoMu_cut_invcut.Report().Print()

#Matching the recoMu to genMu
df_recoMu_cut_matched = df_recoMu_cut.Define('leadingRecoMu_lv', 'myLep(recoMu_pT, recoMu_eta, recoMu_phi, mu_mass, 0)')\
        .Define('subleadingRecoMu_lv', 'myLep(recoMu_pT, recoMu_eta, recoMu_phi, mu_mass, 1)')\
        .Define('leadingRecoMu_genParIndex', 'matchReco(leadingRecoMu_lv, genParPt, genParEta, genParPhi,genParE,genParId,genParMotherId)')\
        .Define('subleadingRecoMu_genParIndex', 'matchReco(subleadingRecoMu_lv, genParPt, genParEta, genParPhi,genParE,genParId,genParMotherId)')\
        .Filter('leadingRecoMu_genParIndex != -1 || subleadingRecoMu_genParIndex != -1')\
        .Define('genParIdLeading', 'genParId[leadingRecoMu_genParIndex]').Define('genParIdSubleading', 'genParId[subleadingRecoMu_genParIndex]')\
        .Define('genParMotherIdLeading', 'genParMotherId[leadingRecoMu_genParIndex]').Define('genParMotherIdSubleading', 'genParMotherId[subleadingRecoMu_genParIndex]')

#df_recoMu_cut_matched.Display(['leadingRecoMu_genParIndex', 'subleadingRecoMu_genParIndex', 'genParIdLeading']).Print()
#df_recoMu_cut_matched.Describe().Print()
#print(df_recoMu_cut_matched.Count().GetValue())

isleadingMufromZCond = 'abs(genParIdLeading) == 13 && genParMotherIdLeading == 23'
issubleadingMufromZCond = 'abs(genParIdSubleading) == 13 && genParMotherIdSubleading == 23'
test = 'abs(genParIdSubleading) == 13 && genParMotherIdSubleading == 23 && abs(genParIdLeading) == 13 && genParMotherIdLeading == 23'
df_recoMu_cut_matched = df_recoMu_cut_matched.Filter(test)

#print('MATCHING')
#df_recoMu_cut_matched.Display(['recoMu_pT_leading', 'recoMu_pT_subleading']).Print()
#df_recoMu_cut_matched.Display(['genParIdLeading', 'genParMotherIdLeading']).Print()

#############################
#####Reconstructed Kaons#####
#############################

chsMass_ = Kmass_val

#df.Display(['packedCandsPdgId'], 10).Print()

#Get the Kaons from the packedCands
isRecoCh = 'abs(packedCandsPdgId) == 211 && packedCandsCharge!=0 && packedCandsHasTrackDetails==1' #reconstructed charged hadrons 

#cuts from muons
df_recoCh = df_recoMu_cut.Define('pCandChId', f'packedCandsPdgId[{isRecoCh}]')\
        .Define('recoCh_charge', f'packedCandsCharge[{isRecoCh}]')\
        .Define('recoCh_trkDetail', f'packedCandsHasTrackDetails[{isRecoCh}]')\
        .Define('K_mass', str(Kmass_val))\
        .Define('recoCh_px', f'packedCandsPx[{isRecoCh}]')\
        .Define('recoCh_py', f'packedCandsPy[{isRecoCh}]')\
        .Define('recoCh_pz', f'packedCandsPz[{isRecoCh}]')
        
#Check how many of the 3000 events have well reconstructed ch, keeping in mind all 3000 events have s->KK at gen level
df_recoCh_noMu = df.Define('pCandChIdNoMu', f'packedCandsPdgId[{isRecoCh}]')\
        .Define('recoChNoMu_charge', f'packedCandsCharge[{isRecoCh}]')

#Cuts on charged hadrons
ch_cuts =  'lv_trk_pt > ' + str(ch_cuts['pt']) + ' && ' + 'abs(lv_trk_eta) < ' + str(ch_cuts['eta'])

df_recoCh_cut = df_recoCh.Define('LVs', f'makeLVs(pCandChId, recoCh_charge, recoCh_trkDetail, recoCh_px, recoCh_py, recoCh_pz, {chsMass_}, true)')\
        .Define('lv_trk_pt','getKinematics(LVs, "pt")')\
        .Define('lv_trk_eta','getKinematics(LVs, "eta")')\
        .Define('ch_trk_pt',f'lv_trk_pt[{ch_cuts}]')\
        .Define('ch_trk_eta',f'lv_trk_eta[{ch_cuts}]')\
        .Define('LVs_sel',f'LVs[{ch_cuts}]')\
        .Define('ch_trk_sel','LVs_sel.size()>=4')


#Filter on the events which have 4 or more Charged Hadrons
df_recoCh_sel = df_recoCh_cut.Filter('ch_trk_sel')



#Set up the 2 hadron pairs, related to the two scalars
df_diCh = df_recoCh_sel.Define('cand_globalidx_hadronsonly',f'getIndices(numPackedCands)[{isRecoCh}]')\
        .Define('ch_globalidx1',f'cand_globalidx_hadronsonly[{ch_cuts}]')\
        .Define('ch1_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch2_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch_pair_idx1',f'''getDileptonCand1(ch_globalidx1, packedCandsCharge, packedCandsPx, packedCandsPy, packedCandsPz, {chsMass_}, numChsTrackPairs,chsTkPairIndex1,chsTkPairIndex2,chsTkPairTk1Px,chsTkPairTk1Py,chsTkPairTk1Pz,chsTkPairTk2Px,chsTkPairTk2Py,chsTkPairTk2Pz, "hadron", packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, ch1_lv, ch2_lv, {diChPt_}, {diChdR_})''')\
        .Define('ch_globalidx2', 'ch_globalidx1[ch_globalidx1 != ch_pair_idx1[0] && ch_globalidx1 != ch_pair_idx1[1]]')\
        .Define('ch3_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch4_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch_pair_idx2',f'''getDileptonCand1(ch_globalidx2, packedCandsCharge, packedCandsPx, packedCandsPy, packedCandsPz, {chsMass_}, numChsTrackPairs,chsTkPairIndex1,chsTkPairIndex2,chsTkPairTk1Px,chsTkPairTk1Py,chsTkPairTk1Pz,chsTkPairTk2Px,chsTkPairTk2Py,chsTkPairTk2Pz, "hadron", packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, ch3_lv, ch4_lv, {diChPt_}, {diChdR_})''')\
        .Define('ch_pair_check','ch_pair_idx1[0]>=0 && ch_pair_idx1[1]>=0 && ch_pair_idx2[0]>=0 && ch_pair_idx2[1]>=0')\
        .Define('ch_pair_checksingle','(ch_pair_idx1[0]>=0 && ch_pair_idx1[1]>=0) || (ch_pair_idx2[0]>=0 && ch_pair_idx2[1]>=0)')
       
df_diCh_check_single = df_diCh.Filter('ch_pair_checksingle')



df_diCh.Display(['ch_trk_pt','ch_globalidx1', 'ch_globalidx2'], 20).Print()


#Seperate the kinematic variables for the hadron pairs
df_diCh_check = df_diCh.Filter('ch_pair_check')\
        .Define('ch1_pT_leading', 'ch1_lv.Pt()').Define('ch2_pT_subleading', 'ch2_lv.Pt()')\
        .Define('ch1_eta', 'ch1_lv.Eta()').Define('ch2_eta', 'ch2_lv.Eta()')\
        .Define('ch1_phi', 'ch1_lv.Phi()').Define('ch2_phi', 'ch2_lv.Phi()')\
        .Define('ch3_pT_leading', 'ch3_lv.Pt()').Define('ch4_pT_subleading', 'ch4_lv.Pt()')\
        .Define('ch3_eta', 'ch3_lv.Eta()').Define('ch4_eta', 'ch4_lv.Eta()')\
        .Define('ch3_phi', 'ch3_lv.Phi()').Define('ch4_phi', 'ch4_lv.Phi()')\
        .Define('ch1_px', 'ch1_lv.Px()').Define('ch1_py', 'ch1_lv.Py()').Define('ch1_pz', 'ch1_lv.Pz()')\
        .Define('ch2_px', 'ch2_lv.Px()').Define('ch2_py', 'ch2_lv.Py()').Define('ch2_pz', 'ch2_lv.Pz()')\
        .Define('ch3_px', 'ch3_lv.Px()').Define('ch3_py', 'ch3_lv.Py()').Define('ch3_pz', 'ch3_lv.Pz()')\
        .Define('ch4_px', 'ch4_lv.Px()').Define('ch4_py', 'ch4_lv.Py()').Define('ch4_pz', 'ch4_lv.Pz()')\
        .Define('ch1_charge', 'packedCandsCharge[ch_pair_idx1[0]]').Define('ch2_charge', 'packedCandsCharge[ch_pair_idx1[1]]')\
        .Define('ch3_charge', 'packedCandsCharge[ch_pair_idx2[0]]').Define('ch4_charge', 'packedCandsCharge[ch_pair_idx2[1]]')
df_diCh_check.Describe().Print()
#Check if pairs can be found as gen matched pairs
isnKfromS1 = 'genParId == -321 && genParMotherId == 9000006'
ispKfromS1 = 'genParId == 321 && genParMotherId == 9000006'
isnKfromS2 = 'genParId == -321 && genParMotherId == -9000006'
ispKfromS2 = 'genParId == 321 && genParMotherId == -9000006'

df_diCh_genmatch = df_diCh_check.Define('genParId_nKS1', f'genParId[{isnKfromS1}]').Define('genParId_pKS1', f'genParId[{ispKfromS1}]')\
        .Define('genParCharge_nKS1', f'genParCharge[{isnKfromS1}]').Define('genParCharge_pKS1', f'genParCharge[{ispKfromS1}]')\
        .Define('genParId_nKS2', f'genParId[{isnKfromS2}]').Define('genParId_pKS2', f'genParId[{ispKfromS2}]')\
        .Define('genParCharge_nKS2', f'genParCharge[{isnKfromS2}]').Define('genParCharge_pKS2', f'genParCharge[{ispKfromS2}]')
print('Filtered Kaon gen')
#df_genKMatch.Display(['genParId_nKS1', 'genParCharge_nKS1', 'genParId_pKS1', 'genParCharge_pKS1']).Print()

#Get the neg and pos kaon kinematics from the scalars
df_diCh_genmatch = df_diCh_genmatch.Define('nKS1_pT', f'genParPt[{isnKfromS1}]').Define('nKS1_eta', f'genParEta[{isnKfromS1}]').Define('nKS1_phi', f'genParPhi[{isnKfromS1}]')\
        .Define('pKS1_pT', f'genParPt[{ispKfromS1}]').Define('pKS1_eta', f'genParEta[{ispKfromS1}]').Define('pKS1_phi', f'genParPhi[{ispKfromS1}]')\
        .Define('nKS2_pT', f'genParPt[{isnKfromS2}]').Define('nKS2_eta', f'genParEta[{isnKfromS2}]').Define('nKS2_phi', f'genParPhi[{isnKfromS2}]')\
        .Define('pKS2_pT', f'genParPt[{ispKfromS2}]').Define('pKS2_eta', f'genParEta[{ispKfromS2}]').Define('pKS2_phi', f'genParPhi[{ispKfromS2}]')

#Make gen LVs
df_diCh_genmatch = df_diCh_genmatch.Define('nKS1_lv', 'ROOT::Math::PtEtaPhiMVector(nKS1_pT[0], nKS1_eta[0], nKS1_phi[0], 0.493677)')\
        .Define('pKS1_lv', 'ROOT::Math::PtEtaPhiMVector(pKS1_pT[0], pKS1_eta[0], pKS1_phi[0], 0.493677)')\
        .Define('nKS2_lv', 'ROOT::Math::PtEtaPhiMVector(nKS2_pT[0], nKS2_eta[0], nKS2_phi[0], 0.493677)')\
        .Define('pKS2_lv', 'ROOT::Math::PtEtaPhiMVector(pKS2_pT[0], pKS2_eta[0], pKS2_phi[0], 0.493677)')

df_diCh_genmatch = df_diCh_genmatch.Define('chsPx', '''ROOT::VecOps::RVec<double> v = {ch1_px, ch2_px, ch3_px, ch4_px}; return v;''')\
        .Define('chsPy', '''ROOT::VecOps::RVec<double> v = {ch1_py, ch2_py, ch3_py, ch4_py}; return v;''')\
        .Define('chsPz', '''ROOT::VecOps::RVec<double> v = {ch1_pz, ch2_pz, ch3_pz, ch4_pz}; return v;''')\
        .Define('chsCharge', '''ROOT::VecOps::RVec<int> v = {ch1_charge, ch2_charge, ch3_charge, ch4_charge}; return v;''')\
        .Define('chsIdx', '''ROOT::VecOps::RVec<Int_t> v = {1, 2, 3, 4}; return v;''')\

#Matching
df_diCh_genmatch = df_diCh_genmatch.Define('nKS1_matchedIdx', 'MatchGenbis(nKS1_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, -1)')\
        .Define('pKS1_matchedIdx', 'MatchGenbis(pKS1_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, 1)')\
        .Define('nKS2_matchedIdx', 'MatchGenbis(nKS2_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, -1)')\
        .Define('pKS2_matchedIdx', 'MatchGenbis(pKS2_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, 1)')\

#Filter out unmatched events
df_diCh_genmatch_filtered = df_diCh_genmatch.Filter('nKS1_matchedIdx != -1 && pKS1_matchedIdx != -1 && nKS2_matchedIdx != -1 && pKS2_matchedIdx != -1 ')

print('Amount of events passing muon cuts (no invmass cut)', df_recoMu_cut.Count().GetValue())
#print('Amount of events with ch reco larger than 4 per events ', df_recoCh_noMu.Filter('recoChNoMu_charge.size() >= 4').Count().GetValue())
print('Amount of events with >= 4 ch, events passed muon selection ', df_recoCh.Filter('recoCh_charge.size() >= 4').Count().GetValue())
print('Amount of events with >= 4 ch after hadron cuts (pT & Eta)', df_recoCh_sel.Count().GetValue())
print('Amount of events where single pair found', df_diCh_check_single.Count().GetValue())
print('Amount of events where two pairs found', df_diCh_check.Count().GetValue())
print('Events before matching conditions on pair selected events', df_diCh_genmatch.Count().GetValue())
print('Amount of events able to be genmatched', df_diCh_genmatch_filtered.Count().GetValue())

print('\nnum of Ch events after quality cut  with 1 or more svVertexChi2', df_diCh.Filter('svVertexChi2.size() >= 1').Count().GetValue())
print('num of Ch events after quality cut with 2 or more svVertexChi2', df_diCh.Filter('svVertexChi2.size() >= 2').Count().GetValue())

df_diCh_genmatch_filtered = df_diCh_genmatch_filtered.Define('genMatchPair1', '''ROOT::VecOps::RVec<Int_t> v = {chsIdx[nKS1_matchedIdx], chsIdx[pKS1_matchedIdx]}; return v;''')\
        .Define('genMatchPair2', '''ROOT::VecOps::RVec<Int_t> v = {chsIdx[nKS2_matchedIdx], chsIdx[pKS2_matchedIdx]}; return v;''')

#df_diCh_genmatch_filtered.Display(['nKS1_matchedIdx', 'pKS1_matchedIdx', 'nKS2_matchedIdx', 'pKS2_matchedIdx']).Print()
df_diCh_genmatch_filtered.Display(['genMatchPair1', 'genMatchPair2'], 100).Print()

#Define the angular differences
df_diCh_check = df_diCh_check.Define('ch12_dEta', 'ch1_eta - ch2_eta').Define('ch34_dEta', 'ch3_eta - ch4_eta')\
        .Define('ch12_dPhi', 'ROOT::VecOps::DeltaPhi(ch1_phi, ch2_phi)').Define('ch34_dPhi', 'ROOT::VecOps::DeltaPhi(ch3_phi, ch4_phi)')\
        .Define('ch12_dR', 'ROOT::VecOps::DeltaR(ch1_eta, ch2_eta, ch1_phi, ch2_phi)').Define('ch34_dR', 'ROOT::VecOps::DeltaR(ch3_eta, ch4_eta, ch3_phi, ch4_phi)')

#df_diCh_check.Display(['ch1_pT_leading', 'ch2_pT_subleading', 'ch_pair_idx1', 'ch3_pT_leading', 'ch4_pT_subleading', 'ch_pair_idx2'], 20).Print()


df_diCh_check = df_diCh_check.Define('s1_lv', 'ch1_lv + ch2_lv' ).Define('s2_lv', 'ch3_lv + ch4_lv')\
        .Define('s1_invmass', 's1_lv.M()').Define('s2_invmass', 's2_lv.M()')\
        .Define('s1_pT', 's1_lv.Pt()').Define('s2_pT', 's2_lv.Pt()')\
        .Define('s12_invmass', '''double s = (s1_invmass + s2_invmass)/2.0; return s;''')\
        .Define('s12_dPhi', 'ROOT::VecOps::DeltaPhi(s1_lv.Phi(), s2_lv.Phi())')\
        .Define('s12_dR', 'ROOT::VecOps::DeltaR(s1_lv.Eta(), s2_lv.Eta(), s1_lv.Phi(), s2_lv.Phi())')

#df_diCh_check.Display(['s1_invmass', 's1_pT', 's2_invmass', 's2_pT'], 20).Print()

#################################################################################################

#Alternate df with same pair selection but now kinematics are kept reco, not the track kinematics
print('Pair selection without track kinematics update')
df_diCh_noTrkKin = df_recoCh_sel.Define('cand_globalidx_hadronsonly',f'getIndices(numPackedCands)[{isRecoCh}]')\
        .Define('ch_globalidx1',f'cand_globalidx_hadronsonly[{ch_cuts}]')\
        .Define('ch1_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch2_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch_pair_idx1',f'''getDileptonCand2(ch_globalidx1, packedCandsCharge, packedCandsPx, packedCandsPy, packedCandsPz, {chsMass_}, numChsTrackPairs,chsTkPairIndex1,chsTkPairIndex2,chsTkPairTk1Px,chsTkPairTk1Py,chsTkPairTk1Pz,chsTkPairTk2Px,chsTkPairTk2Py,chsTkPairTk2Pz, "hadron", packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, ch1_lv, ch2_lv, {diChPt_}, {diChdR_})''')\
        .Define('ch_globalidx2', 'ch_globalidx1[ch_globalidx1 != ch_pair_idx1[0] && ch_globalidx1 != ch_pair_idx1[1]]')\
        .Define('ch3_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch4_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch_pair_idx2',f'''getDileptonCand2(ch_globalidx2, packedCandsCharge, packedCandsPx, packedCandsPy, packedCandsPz, {chsMass_}, numChsTrackPairs,chsTkPairIndex1,chsTkPairIndex2,chsTkPairTk1Px,chsTkPairTk1Py,chsTkPairTk1Pz,chsTkPairTk2Px,chsTkPairTk2Py,chsTkPairTk2Pz, "hadron", packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, ch3_lv, ch4_lv, {diChPt_}, {diChdR_})''')\
        .Define('ch_pair_check','ch_pair_idx1[0]>=0 && ch_pair_idx1[1]>=0 && ch_pair_idx2[0]>=0 && ch_pair_idx2[1]>=0')\
        .Define('ch_pair_checksingle','(ch_pair_idx1[0]>=0 && ch_pair_idx1[1]>=0) || (ch_pair_idx2[0]>=0 && ch_pair_idx2[1]>=0)')
       
df_diCh_noTrkKin_check_single = df_diCh_noTrkKin.Filter('ch_pair_checksingle')

#Seperate the kinematic variables for the hadron pairs
df_diCh_noTrkKin_check = df_diCh_noTrkKin.Filter('ch_pair_check')\
        .Define('ch1_pT_leading', 'ch1_lv.Pt()').Define('ch2_pT_subleading', 'ch2_lv.Pt()')\
        .Define('ch1_eta', 'ch1_lv.Eta()').Define('ch2_eta', 'ch2_lv.Eta()')\
        .Define('ch1_phi', 'ch1_lv.Phi()').Define('ch2_phi', 'ch2_lv.Phi()')\
        .Define('ch3_pT_leading', 'ch3_lv.Pt()').Define('ch4_pT_subleading', 'ch4_lv.Pt()')\
        .Define('ch3_eta', 'ch3_lv.Eta()').Define('ch4_eta', 'ch4_lv.Eta()')\
        .Define('ch3_phi', 'ch3_lv.Phi()').Define('ch4_phi', 'ch4_lv.Phi()')\

print('Amount of events passing muon cuts (no invmass cut)', df_recoMu_cut.Count().GetValue())
#print('Amount of events with ch reco larger than 4 per events ', df_recoCh_noMu.Filter('recoChNoMu_charge.size() >= 4').Count().GetValue())
print('Amount of events with >= 4 ch, events passed muon selection ', df_recoCh.Filter('recoCh_charge.size() >= 4').Count().GetValue())
print('Amount of events with >= 4 ch after hadron cuts (pT & Eta)', df_recoCh_sel.Count().GetValue())
print('Amount of events where single pair found, noTrkKin', df_diCh_noTrkKin_check_single.Count().GetValue())
print('Amount of events where two pairs found, noTrkKin', df_diCh_noTrkKin_check.Count().GetValue())

print('\nnum of Ch events after quality cut  with 1 or more svVertexChi2, no TrkKin', df_diCh_noTrkKin.Filter('svVertexChi2.size() >= 1').Count().GetValue())
print('num of Ch events after quality cut with 2 or more svVertexChi2, no TrkKin', df_diCh_noTrkKin.Filter('svVertexChi2.size() >= 2').Count().GetValue())
#Define the angular differences
df_diCh_noTrkKin_check = df_diCh_noTrkKin_check.Define('ch12_dEta', 'ch1_eta - ch2_eta').Define('ch34_dEta', 'ch3_eta - ch4_eta')\
        .Define('ch12_dPhi', 'ROOT::VecOps::DeltaPhi(ch1_phi, ch2_phi)').Define('ch34_dPhi', 'ROOT::VecOps::DeltaPhi(ch3_phi, ch4_phi)')\
        .Define('ch12_dR', 'ROOT::VecOps::DeltaR(ch1_eta, ch2_eta, ch1_phi, ch2_phi)').Define('ch34_dR', 'ROOT::VecOps::DeltaR(ch3_eta, ch4_eta, ch3_phi, ch4_phi)')

#df_diCh_check.Display(['ch1_pT_leading', 'ch2_pT_subleading', 'ch_pair_idx1', 'ch3_pT_leading', 'ch4_pT_subleading', 'ch_pair_idx2'], 20).Print()


df_diCh_noTrkKin_check = df_diCh_noTrkKin_check.Define('s1_lv', 'ch1_lv + ch2_lv' ).Define('s2_lv', 'ch3_lv + ch4_lv')\
        .Define('s1_invmass', 's1_lv.M()').Define('s2_invmass', 's2_lv.M()')\
        .Define('s1_pT', 's1_lv.Pt()').Define('s2_pT', 's2_lv.Pt()')\
        .Define('s12_dPhi', 'ROOT::VecOps::DeltaPhi(s1_lv.Phi(), s2_lv.Phi())')\
        .Define('s12_dR', 'ROOT::VecOps::DeltaR(s1_lv.Eta(), s2_lv.Eta(), s1_lv.Phi(), s2_lv.Phi())')
#############################
#####Matched Reco Kaons#####
#############################

isnKfromS1 = 'genParId == -321 && genParMotherId == 9000006'
ispKfromS1 = 'genParId == 321 && genParMotherId == 9000006'
isnKfromS2 = 'genParId == -321 && genParMotherId == -9000006'
ispKfromS2 = 'genParId == 321 && genParMotherId == -9000006'

df_genKMatch = df_recoMu_cut.Define('genParId_nKS1', f'genParId[{isnKfromS1}]').Define('genParId_pKS1', f'genParId[{ispKfromS1}]')\
        .Define('genParCharge_nKS1', f'genParCharge[{isnKfromS1}]').Define('genParCharge_pKS1', f'genParCharge[{ispKfromS1}]')\
        .Define('genParId_nKS2', f'genParId[{isnKfromS2}]').Define('genParId_pKS2', f'genParId[{ispKfromS2}]')\
        .Define('genParCharge_nKS2', f'genParCharge[{isnKfromS2}]').Define('genParCharge_pKS2', f'genParCharge[{ispKfromS2}]')
print('Filtered Kaon gen')
df_genKMatch.Display(['genParId_nKS1', 'genParCharge_nKS1', 'genParId_pKS1', 'genParCharge_pKS1']).Print()

#Get the neg and pos kaon kinematics from the scalars
df_genKMatch = df_genKMatch.Define('nKS1_pT', f'genParPt[{isnKfromS1}]').Define('nKS1_eta', f'genParEta[{isnKfromS1}]').Define('nKS1_phi', f'genParPhi[{isnKfromS1}]')\
        .Define('pKS1_pT', f'genParPt[{ispKfromS1}]').Define('pKS1_eta', f'genParEta[{ispKfromS1}]').Define('pKS1_phi', f'genParPhi[{ispKfromS1}]')\
        .Define('nKS2_pT', f'genParPt[{isnKfromS2}]').Define('nKS2_eta', f'genParEta[{isnKfromS2}]').Define('nKS2_phi', f'genParPhi[{isnKfromS2}]')\
        .Define('pKS2_pT', f'genParPt[{ispKfromS2}]').Define('pKS2_eta', f'genParEta[{ispKfromS2}]').Define('pKS2_phi', f'genParPhi[{ispKfromS2}]')

#Make gen LVs
df_genKMatch = df_genKMatch.Define('nKS1_lv', 'ROOT::Math::PtEtaPhiMVector(nKS1_pT[0], nKS1_eta[0], nKS1_phi[0], 0.493677)')\
        .Define('pKS1_lv', 'ROOT::Math::PtEtaPhiMVector(pKS1_pT[0], pKS1_eta[0], pKS1_phi[0], 0.493677)')\
        .Define('nKS2_lv', 'ROOT::Math::PtEtaPhiMVector(nKS2_pT[0], nKS2_eta[0], nKS2_phi[0], 0.493677)')\
        .Define('pKS2_lv', 'ROOT::Math::PtEtaPhiMVector(pKS2_pT[0], pKS2_eta[0], pKS2_phi[0], 0.493677)')

#Get the matched indices of the reco Kaons
df_recoCh_matched = df_genKMatch.Define('nKS1_matchedIdx', 'MatchGen(nKS1_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, -1)')\
        .Define('pKS1_matchedIdx', 'MatchGen(pKS1_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, 1)')\
        .Define('nKS2_matchedIdx', 'MatchGen(nKS2_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, -1)')\
        .Define('pKS2_matchedIdx', 'MatchGen(pKS2_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, 1)')

print('Events before matching conditions', df_recoCh_matched.Count().GetValue())

#Filter out unmatched events
df_recoCh_matched_filtered = df_recoCh_matched.Filter('nKS1_matchedIdx != -1 && pKS1_matchedIdx != -1 && nKS2_matchedIdx != -1 && pKS2_matchedIdx != -1 ')
#Create a df with the unmatched events
df_recoCh_nomatch = df_recoCh_matched.Filter('nKS1_matchedIdx == -1 || pKS1_matchedIdx == -1 || nKS2_matchedIdx == -1 || pKS2_matchedIdx == -1 ')

df_recoCh_nomatch = df_recoCh_nomatch.Define('genKS1NoMatch_absEta', 'abs(nKS1_eta[0]) > abs(pKS1_eta[0]) ? abs(nKS1_eta[0]) : abs(pKS1_eta[0])')\
        .Define('genKS2NoMatch_absEta', 'abs(nKS2_eta[0]) > abs(pKS2_eta[0]) ? abs(nKS2_eta[0]) : abs(pKS2_eta[0])')\
        .Define('genKS1NoMatch_pT_leading', 'nKS1_pT[0] > pKS1_pT[0] ? nKS1_pT[0] : pKS1_pT[0]')\
        .Define('genKS1NoMatch_pT_subleading', 'nKS1_pT[0] < pKS1_pT[0] ? nKS1_pT[0] : pKS1_pT[0]')\
        .Define('genKS2NoMatch_pT_leading', 'nKS2_pT[0] > pKS2_pT[0]? nKS2_pT[0] : pKS2_pT[0]')\
        .Define('genKS2NoMatch_pT_subleading', 'nKS2_pT[0] < pKS2_pT[0] ? nKS2_pT[0] : pKS2_pT[0]')

print('Matching idx and TkPairCheck')
df_recoCh_matched_filtered.Display(['nKS1_matchedIdx', 'pKS1_matchedIdx', 'nKS2_matchedIdx', 'pKS2_matchedIdx']).Print()

#Check if the gen pairs are track pairs, check crossover as well
df_recoCh_matched_TkCheck = df_recoCh_matched_filtered.Define('TkPair1Check', f'checkTrkPair(nKS1_matchedIdx, pKS1_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
        .Define('TkPair2Check', f'checkTrkPair(nKS2_matchedIdx, pKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
        .Define('TkPairCrossCheck_n1p2', f'checkTrkPair(nKS1_matchedIdx, pKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
        .Define('TkPairCrossCheck_n2p1', f'checkTrkPair(nKS2_matchedIdx, pKS1_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
        .Define('TkPairCrossCheck_n1n2', f'checkTrkPair(nKS1_matchedIdx, nKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
        .Define('TkPairCrossCheck_p1p2', f'checkTrkPair(pKS1_matchedIdx, pKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')
 
df_recoCh_matched_TkCheck.Display(['TkPair1Check', 'TkPair2Check'], 10).Print()
print('Number of unmatched events: ', df_recoCh_nomatch.Count().GetValue())
print('Number of matched events: ', df_recoCh_matched_filtered.Count().GetValue())
print('Number of matched events and one correct Tkpair found: ', df_recoCh_matched_TkCheck.Filter('TkPair1Check != -1 || TkPair2Check != -1').Count().GetValue())
print('Number of matched events and both correct Tkpair found: ', df_recoCh_matched_TkCheck.Filter('TkPair1Check != -1 && TkPair2Check != -1').Count().GetValue())
print('-----------CrossPairs-----------')
print('Number of matched events where n1p2 or n2p1 is a track pair: ', df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1p2 != -1 || TkPairCrossCheck_n2p1 != -1').Count().GetValue())
print('Number of matched events where n1p2 and n2p1 is a track pair: ', df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1p2 != -1 && TkPairCrossCheck_n2p1 != -1').Count().GetValue())
print('Number of matched events where n1n2 or p1p2 is a track pair: ', df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1n2 != -1 || TkPairCrossCheck_p1p2 != -1').Count().GetValue())
print('Number of matched events where n1n2 andp1p2 is a track pair: ', df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1n2 != -1 && TkPairCrossCheck_p1p2 != -1').Count().GetValue())
print('-------------------------------')
print('num of matched gen events with 1 or more svVertexChi2', df_recoCh_matched_filtered.Filter('svVertexChi2.size() >= 1').Count().GetValue())
print('num of matched gen events with 2 or more svVertexChi2', df_recoCh_matched_filtered.Filter('svVertexChi2.size() >= 2').Count().GetValue())
print('num of total events with 1 or more svVertexChi2', df.Filter('svVertexChi2.size() >= 1').Count().GetValue())
print('num of total events with 2 or more svVertexChi2', df.Filter('svVertexChi2.size() >= 2').Count().GetValue())
df_recoCh_matched_filtered.Display(['pvChi2', 'svVertexChi2'], 20).Print()
#Select the matched reco px, py, pz
df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('nKS1Matched_px', 'packedCandsPx[nKS1_matchedIdx]')\
        .Define('nKS1Matched_py', 'packedCandsPy[nKS1_matchedIdx]').Define('nKS1Matched_pz', 'packedCandsPz[nKS1_matchedIdx]')\
        .Define('pKS1Matched_px', 'packedCandsPx[pKS1_matchedIdx]').Define('pKS1Matched_py', 'packedCandsPy[pKS1_matchedIdx]').Define('pKS1Matched_pz', 'packedCandsPz[pKS1_matchedIdx]')\
        .Define('nKS2Matched_px', 'packedCandsPx[nKS2_matchedIdx]').Define('nKS2Matched_py', 'packedCandsPy[nKS2_matchedIdx]').Define('nKS2Matched_pz', 'packedCandsPz[nKS2_matchedIdx]')\
        .Define('pKS2Matched_px', 'packedCandsPx[pKS2_matchedIdx]').Define('pKS2Matched_py', 'packedCandsPy[pKS2_matchedIdx]').Define('pKS2Matched_pz', 'packedCandsPz[pKS2_matchedIdx]')

#Make matched reco LVs
df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('nKS1Matched_lv', 'ROOT::Math::PxPyPzMVector(nKS1Matched_px, nKS1Matched_py, nKS1Matched_pz, 0.493677)')\
        .Define('pKS1Matched_lv', 'ROOT::Math::PxPyPzMVector(pKS1Matched_px, pKS1Matched_py, pKS1Matched_pz, 0.493677)')\
        .Define('nKS2Matched_lv', 'ROOT::Math::PxPyPzMVector(nKS2Matched_px, nKS2Matched_py, nKS2Matched_pz, 0.493677)')\
        .Define('pKS2Matched_lv', 'ROOT::Math::PxPyPzMVector(pKS2Matched_px, pKS2Matched_py, pKS2Matched_pz, 0.493677)')

#Get pT, Eta, Phi from the mathced reco LVs
df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('nKS1Matched_pT', 'nKS1Matched_lv.Pt()')\
        .Define('nKS1Matched_eta', 'nKS1Matched_lv.Eta()').Define('nKS1Matched_phi', 'nKS1Matched_lv.Phi()')\
        .Define('pKS1Matched_pT', 'pKS1Matched_lv.Pt()').Define('pKS1Matched_eta', 'pKS1Matched_lv.Eta()').Define('pKS1Matched_phi', 'pKS1Matched_lv.Phi()')\
        .Define('nKS2Matched_pT', 'nKS2Matched_lv.Pt()').Define('nKS2Matched_eta', 'nKS2Matched_lv.Eta()').Define('nKS2Matched_phi', 'nKS2Matched_lv.Phi()')\
        .Define('pKS2Matched_pT', 'pKS2Matched_lv.Pt()').Define('pKS2Matched_eta', 'pKS2Matched_lv.Eta()').Define('pKS2Matched_phi', 'pKS2Matched_lv.Phi()')\
        .Define('S1_lv', 'nKS1Matched_lv + pKS1Matched_lv').Define('S2_lv', 'nKS2Matched_lv + pKS2Matched_lv')\
        .Define('s1Matched_invmass', 'S1_lv.M()').Define('s2Matched_invmass', 'S2_lv.M()')\
        .Define('s1Matched_pT', 'S1_lv.Pt()').Define('s2Matched_pT', 'S2_lv.Pt()')

print('unsorted matched negative positive KAON kinematics')
df_recoCh_matched_filtered.Display(['nKS1Matched_pT', 'pKS1Matched_pT', 'nKS2Matched_pT', 'pKS2Matched_pT']).Print()
df_recoCh_matched_filtered.Display(['nKS1Matched_eta', 'pKS1Matched_eta', 'nKS2Matched_eta', 'pKS2Matched_eta']).Print()
df_recoCh_matched_filtered.Display(['nKS1Matched_phi', 'pKS1Matched_phi', 'nKS2Matched_phi', 'pKS2Matched_phi']).Print()

#Seperate leading and subleading contributions
df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('ch1Matched_pT_leading', 'nKS1Matched_pT > pKS1Matched_pT ? nKS1Matched_pT:pKS1Matched_pT')\
        .Define('ch2Matched_pT_subleading', 'nKS1Matched_pT > pKS1Matched_pT ? pKS1Matched_pT:nKS1Matched_pT')\
        .Define('ch1Matched_eta', 'nKS1Matched_pT > pKS1Matched_pT ? nKS1Matched_eta:pKS1Matched_eta')\
        .Define('ch2Matched_eta', 'nKS1Matched_pT > pKS1Matched_pT ? pKS1Matched_eta:nKS1Matched_eta')\
        .Define('ch1Matched_phi', 'nKS1Matched_pT > pKS1Matched_pT ? nKS1Matched_phi:pKS1Matched_phi')\
        .Define('ch2Matched_phi', 'nKS1Matched_pT > pKS1Matched_pT ? pKS1Matched_phi:nKS1Matched_phi')\
        .Define('ch3Matched_pT_leading', 'nKS2Matched_pT > pKS2Matched_pT ? nKS2Matched_pT:pKS2Matched_pT')\
        .Define('ch4Matched_pT_subleading', 'nKS2Matched_pT > pKS2Matched_pT ? pKS2Matched_pT:nKS2Matched_pT')\
        .Define('ch3Matched_eta', 'nKS2Matched_pT > pKS2Matched_pT ? nKS2Matched_eta:pKS2Matched_eta')\
        .Define('ch4Matched_eta', 'nKS2Matched_pT > pKS2Matched_pT ? pKS2Matched_eta:nKS2Matched_eta')\
        .Define('ch3Matched_phi', 'nKS2Matched_pT > pKS2Matched_pT ? nKS2Matched_phi:pKS2Matched_phi')\
        .Define('ch4Matched_phi', 'nKS2Matched_pT > pKS2Matched_pT ? pKS2Matched_phi:nKS2Matched_phi')

print('PT sorted matched KAON kinematics')
df_recoCh_matched_filtered.Display(['ch1Matched_pT_leading', 'ch1Matched_eta', 'ch1Matched_phi']).Print()
df_recoCh_matched_filtered.Display(['ch2Matched_pT_subleading', 'ch2Matched_eta', 'ch2Matched_phi']).Print()
df_recoCh_matched_filtered.Display(['ch3Matched_pT_leading', 'ch3Matched_eta', 'ch3Matched_phi']).Print()
df_recoCh_matched_filtered.Display(['ch4Matched_pT_subleading', 'ch4Matched_eta', 'ch4Matched_phi']).Print()

############################################################

#Check quality cut on the matched reco ch
print('performing pair selection on the matched ch')
df_recoCh_matched_pairsel = df_recoCh_matched_filtered.Define('LVs','''ROOT::VecOps::RVec<ROOT::Math::PxPyPzMVector> v = {nKS1Matched_lv, pKS1Matched_lv, nKS2Matched_lv, pKS2Matched_lv}; return v;''' )\
        .Define('lv_trk_pt','getKinematics(LVs, "pt")')\
        .Define('lv_trk_eta','getKinematics(LVs, "eta")')\
        .Define('ch_trk_pt',f'lv_trk_pt[{ch_cuts}]')\
        .Define('ch_trk_eta',f'lv_trk_eta[{ch_cuts}]')\
        .Define('LVs_sel',f'LVs[{ch_cuts}]')\
        .Define('ch_trk_sel','LVs_sel.size()>=4')


#Filter on the events which have 4 or more Charged Hadrons
df_recoCh_matched_pairsel_qualcut = df_recoCh_matched_pairsel.Filter('ch_trk_sel')



#print('Amount of events with ch reco larger than 4 per events ', df_recoCh_noMu.Filter('recoChNoMu_charge.size() >= 4').Count().GetValue())
print('Amount of matched events with >= 4 ch after quality cuts (pT & Eta)', df_recoCh_matched_pairsel_qualcut.Count().GetValue())

print('\nnum of matched events after quality cut  with 1 or more svVertexChi2', df_recoCh_matched_pairsel_qualcut.Filter('svVertexChi2.size() >= 1').Count().GetValue())
print('num of matched events after quality cut with 2 or more svVertexChi2', df_recoCh_matched_pairsel_qualcut.Filter('svVertexChi2.size() >= 2').Count().GetValue())

##################################################################

#Setup the weights
#years = []
#if not 'Run' in args.config:
#    totalLumi = 0
#    luminosity_ = conf_pars.get("luminosity", {})
#    for year, lumi in luminosity_.items():
#        year = str(year)
#        totalLumi += lumi
#        df_dimuon = df_dimuon.Define("evt_weight_" + year, f'({crossSection}*{lumi}/{sumWeights})*_weight')
#        years.append(year)
#    print(totalLumi) 
#    df_dimuon = df_dimuon.Define("evt_weight_TOTAL", f'({crossSection}*{totalLumi}/{sumWeights})*_weight')
#    years.append('TOTAL')

#########################################################################################
# Extracting information and distributions
#########################################################################################

#if False:
    #Print out the efficiencies for each step:
    #totalMuEntries = df_muons.Count().GetValue()
    #eff_triggercut = df_muontrigger.Count().GetValue()/totalMuEntries
    #eff_triggercut2 = df_muontrigger2.Count().GetValue()/totalMuEntries
    #eff_triggercut3 = df_muontrigger3.Count().GetValue()/totalMuEntries
    #eff_cuts = df_muons_aftercut.Count().GetValue()/totalMuEntries
    #eff_sel = df_muons_aftercutselection.Count().GetValue()/totalMuEntries
    #eff_dimuon = df_dimuon.Count().GetValue()/totalMuEntries
    #print("\nTotal events with only muons in dataset: " + str(totalMuEntries))
    #print("\nEfficiency after HLTIsoMu27: " + str(eff_triggercut))
    #print("\nEfficiency after HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL: " + str(eff_triggercut2))
    #print("\nEfficiency after HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ: " + str(eff_triggercut3))
    #print("\nEfficiency after cuts: " + str(eff_cuts))
    #print("\nEfficiency after cuts and selection: " + str(eff_sel))
    #print("\nEfficiency after dimuon mass selection: " + str(eff_dimuon))


#Generate Histograms and put them into a Dictonary
dataTag = True if 'Run' in args.config else False

needGeneratePlots = True
if(needGeneratePlots):
    #generate the plots and save to output file
    outFile = ROOT.TFile(args.output, "RECREATE")
    #print('Generating gen Kaon pT plots')
    #genXHist(df_genScalar, 'genKS1_pT_leading', 0, 100, 30)
    #genXHist(df_genScalar, 'genKS1_pT_subleading', 0, 100, 30)
    #genXHist(df_genScalar, 'genKS2_pT_leading', 0, 100, 30)
    #genXHist(df_genScalar, 'genKS2_pT_subleading', 0, 100, 30)
    
    #print('Generating gen Kaon angular plots')
    #genXHist(df_genScalar, 'genKS1_eta', -2.5, 2.5, 30)
    #genXHist(df_genScalar, 'genKS2_eta', -2.5, 2.5, 30)
    
    #genXHist(df_genScalar, 'genKS1_phi', -4, 4, 30)
    #genXHist(df_genScalar, 'genKS2_phi', -4, 4, 30)

    #print('Generating gen Kaon/Scalar angular seperation plots')
    #genXHist(df_genScalar, 'genKS1_dPhi', -0.4, 0.4, 80)
    #genXHist(df_genScalar, 'genKS2_dPhi', -0.4, 0.4, 80)
    #genXHist(df_genScalar, 'genKS1_dR', -0.05, 0.475, 80)
    #genXHist(df_genScalar, 'genKS2_dR', -0.05, 0.475, 80)
    
    #genXHist(df_genScalar, 'genS12_dPhi', -5, 5, 30)
    #genXHist(df_genScalar, 'genS12_dR', -0.5, 5, 30)
    
    #print('Generating gen scalar invmass and pT plots')
    #genXHist(df_genScalar, 'genS1_invmass', -0.025, 2.375, 48)
    #genXHist(df_genScalar, 'genS2_invmass', -0.025, 2.375, 48)
    #genXHist(df_genScalar, 'genS1_invmass1', -0.025, 2.375, 48)
    #genXHist(df_genScalar, 'genS2_invmass1', -0.025, 2.375, 48)

    #genXHist(df_genScalar, 'genS1_pT', -0.0, 200, 50)
    #genXHist(df_genScalar, 'genS2_pT', -0.0, 200, 50)

    #print('Generating gen Muon plots')
    #print('Generating number of muons histograms')
    #genXHist(df_genMu, 'genMu_size', -0.5, 9.5, 10)
    #genXHist(df_genMuZ, 'genMuZ_sizeZ', -0.5, 9.5, 10)
    #print('Generating pT histograms')
    #genXHist(df_genMuZ, 'genMuZ_pT_leading', 0, 100, 30)
    #genXHist(df_genMuZ, 'genMuZ_pT_subleading', 0, 100, 30)
    #print('Generating angular histograms')
    #genXHist(df_genMuZ, 'genMuZ_eta', -5.5, 5.5, 30)
    #genXHist(df_genMuZ, 'genMuZ_phi', -4, 4, 30)
    #print('Generating invariant mass histogram')
    #genXHist(df_genMuZ_invcut, 'genZ_invmass', 70, 110, 50)
    
    #print('Generating reco K pT plots')
    #genXHist(df_diCh_check, 'ch1_pT_leading', 0, 100, 30)
    #genXHist(df_diCh_check, 'ch2_pT_subleading', 0, 100, 30)
    #genXHist(df_diCh_check, 'ch3_pT_leading', 0, 100, 30)
    #genXHist(df_diCh_check, 'ch4_pT_subleading', 0, 100, 30)
    
    #print('Generating reco K pT plots, no TrkKin')
    #genXHist(df_diCh_noTrkKin_check, 'ch1_pT_leading', 0, 100, 30)
    #genXHist(df_diCh_noTrkKin_check, 'ch2_pT_subleading', 0, 100, 30)
    #genXHist(df_diCh_noTrkKin_check, 'ch3_pT_leading', 0, 100, 30)
    #genXHist(df_diCh_noTrkKin_check, 'ch4_pT_subleading', 0, 100, 30)
    
    #print('Generating Unmatched gen plots')
    #genXHist(df_recoCh_nomatch, 'genKS1NoMatch_pT_leading', 0, 100, 30)
    #genXHist(df_recoCh_nomatch, 'genKS1NoMatch_pT_subleading', 0, 100, 30)
    #genXHist(df_recoCh_nomatch, 'genKS2NoMatch_pT_leading', 0, 100, 30)
    #genXHist(df_recoCh_nomatch, 'genKS2NoMatch_pT_subleading', 0, 100, 30)

    #genXHist(df_recoCh_nomatch, 'genKS1NoMatch_absEta', 0, 2.4, 40)
    #genXHist(df_recoCh_nomatch, 'genKS2NoMatch_absEta', 0, 2.4, 40)
    
    #gen2DHist(df_recoCh_nomatch, 'genKS1NoMatch_absEta', 'genKS2NoMatch_absEta', 0, 2.4, 20, 'abs($\eta$)_{KS1}', 'abs($\eta$)_{KS2}')
    #gen2DHist(df_recoCh_nomatch, 'genKS1NoMatch_pT_subleading', 'genKS2NoMatch_pT_subleading', 5, 100, 30, 'pTsubleading_{KS1}', 'pTsubleading_{KS2}')
    #gen2DHist(df_recoCh_nomatch, 'genKS1NoMatch_pT_leading', 'genKS2NoMatch_pT_leading', 5, 100, 30, 'pTleading_{KS1}', 'pTleading_{KS2}')

    #print('Generating Matched reco K pT plots')
    #genXHist(df_recoCh_matched_filtered, 'ch1Matched_pT_leading', 0, 100, 30)
    #genXHist(df_recoCh_matched_filtered, 'ch2Matched_pT_subleading', 0, 100, 30)
    #genXHist(df_recoCh_matched_filtered, 'ch3Matched_pT_leading', 0, 100, 30)
    #genXHist(df_recoCh_matched_filtered, 'ch4Matched_pT_subleading', 0, 100, 30)

    #genXHist(df_recoCh_matched_filtered, 's1Matched_invmass', 0.8, 6.6, 100)
    #genXHist(df_recoCh_matched_filtered, 's2Matched_invmass', 0.8, 6.6, 100)

    #genXHist(df_recoCh_matched_filtered, 's1Matched_pT', 0, 200, 30)
    #genXHist(df_recoCh_matched_filtered, 's2Matched_pT', 0, 200, 30)

    #print('Generating reco K angular plots')
    #genXHist(df_diCh_check, 'ch1_eta', -2.5, 2.5, 30)
    #genXHist(df_diCh_check, 'ch2_eta', -2.5, 2.5, 30)
    #genXHist(df_diCh_check, 'ch3_eta', -2.5, 2.5, 30)
    #genXHist(df_diCh_check, 'ch4_eta', -2.5, 2.5, 30)
    
    #genXHist(df_diCh_check, 'ch1_phi', -4, 4, 30)
    #genXHist(df_diCh_check, 'ch2_phi', -4, 4, 30)
    #genXHist(df_diCh_check, 'ch3_phi', -4, 4, 30)
    #genXHist(df_diCh_check, 'ch4_phi', -4, 4, 30)

    #print('Generating reco Scalar invmass, pT plots')
    #genXHist(df_diCh_check, 's1_invmass', 0.8, 1.6, 40)
    #genXHist(df_diCh_check, 's2_invmass', 0.8, 1.6, 40)
    genXHist(df_diCh_check, 's12_invmass', 0.8, 1.6, 40)
    #genXHist(df_diCh_check, 's1_pT', -0, 200, 50)
    #genXHist(df_diCh_check, 's2_pT', -0, 200, 50)
    
    print('Generating 2D scalar invmass plot')
    
    gen2DHist(df_diCh_check, 's1_invmass', 's2_invmass', 0.4, 2.0, 80, 'm_{s1} (GeV)', 'm_{s2} (GeV)')

    #print('Generating reco Scalar invmass, pT plots, no TrkKin')
    #genXHist(df_diCh_noTrkKin_check, 's1_invmass', 0.8, 1.6, 40)
    #genXHist(df_diCh_noTrkKin_check, 's2_invmass', 0.8, 1.6, 40)
    #genXHist(df_diCh_noTrkKin_check, 's1_pT', -0, 200, 50)
    #genXHist(df_diCh_noTrkKin_check, 's2_pT', -0, 200, 50)
    
    #print('Generating reco Kaon/Scalar angular difference plots')
    #genXHist(df_diCh_check, 'ch12_dPhi', -0.4, 0.4, 80)
    #genXHist(df_diCh_check, 'ch34_dPhi', -0.4, 0.4, 80)
    #genXHist(df_diCh_check, 'ch12_dR', -0.05, 0.45, 80)
    #genXHist(df_diCh_check, 'ch34_dR', -0.05, 0.45, 80)

    #genXHist(df_diCh_check, 's12_dPhi', -5, 5, 30)
    #genXHist(df_diCh_check, 's12_dR', -0.5, 5, 30)

    #print('Generating reco Muon angular plots')
    #genXHist(df_recoMu_cut, 'recoMu_pT_leading', 0, 100, 30)
    #genXHist(df_recoMu_cut, 'recoMu_pT_subleading', 0, 100, 30)
    #genXHist(df_recoMu_cut, 'recoMu_eta', -2.5, 2.5, 30)
    #genXHist(df_recoMu_cut, 'recoMu_phi', -4, 4, 30)
    genXHist(df_recoMu_cut_invcut, 'recoZ_invmass', 70, 110, 50)
    genXHist(df_recoMu_cut_invcut, 'recoZ_pT', 0, 200, 80)
    outFile.Close()

sys.stderr.write("\nTime taken: --- %s seconds ---" % (time.time() - start_time))

print("\n Program running... Press Enter to stop.")
input()  # Waits for the key to be pressed
print("Program stopped.")













