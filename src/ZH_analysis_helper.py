import ROOT

def gInterpreter_lv():
	lvCode = '''
		using FourVector = ROOT::Math::PxPyPzMVector;
		using RVecFloat = ROOT::RVec<Float_t>;
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
	'''
	ROOT.gInterpreter.Declare(lvCode)	

def gInterpreter_getIndices():
	getIndicesCode = '''
		using namespace ROOT::VecOps;
		RVec<unsigned long> getIndices(const int num) {
    		RVec<unsigned long> idx;
    		for (int i = 0;i<num;i++) idx.emplace_back(i);
    		return idx;
		}	
	'''
	ROOT.gInterpreter.Declare(getIndicesCode)

def gInterpreter_pairselection():
	TkpairCode = '''
		using FourVector = ROOT::Math::PxPyPzMVector;
		using RVecFloat = ROOT::RVec<Float_t>;
		using FourVectorPtEtaPhiM = ROOT::Math::PtEtaPhiMVector;
		using FourVectorPtEtaPhiE = ROOT::Math::PtEtaPhiEVector;
		using FourVectorPxPyPzM = ROOT::Math::PxPyPzMVector;
		using namespace ROOT::VecOps;

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
	'''

	pairselCode = '''
		using FourVector = ROOT::Math::PxPyPzMVector;
		using RVecFloat = ROOT::RVec<Float_t>;
		using FourVectorPtEtaPhiM = ROOT::Math::PtEtaPhiMVector;
		using FourVectorPtEtaPhiE = ROOT::Math::PtEtaPhiEVector;
		using FourVectorPxPyPzM = ROOT::Math::PxPyPzMVector;
		using namespace ROOT::VecOps;
		
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
			} else {
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
		'''
	ROOT.gInterpreter.Declare(TkpairCode)
	ROOT.gInterpreter.Declare(pairselCode)
def gInterpreter_matching():
	matchingCode = '''
		using FourVector = ROOT::Math::PxPyPzMVector;
		using RVecFloat = ROOT::RVec<Float_t>;
		using FourVectorPtEtaPhiM = ROOT::Math::PtEtaPhiMVector;
		using FourVectorPtEtaPhiE = ROOT::Math::PtEtaPhiEVector;
		using FourVectorPxPyPzM = ROOT::Math::PxPyPzMVector;
		using namespace ROOT::VecOps;

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
	'''
	ROOT.gInterpreter.Declare(matchingCode)

def gInterpreter_getKinematics():
	getLeadingTrailingCode = '''
		using namespace ROOT::VecOps;
		using RVecFloat = ROOT::RVec<Float_t>;
		float getLeading(RVecFloat vec){
    		auto idxmax = ROOT::VecOps::ArgMax(vec);
    		return vec[idxmax];

		}

		float getTrailing(RVecFloat vec){
    		auto idxmin = ROOT::VecOps::ArgMin(vec);
    		return vec[idxmin];
		}
	'''
	getDeltaRCode = '''
		using namespace ROOT::VecOps;
		using RVecFloat = ROOT::RVec<Float_t>;
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
	'''
	getKinematicsCode = '''
		using namespace ROOT::VecOps;
		using FourVector = ROOT::Math::PxPyPzMVector;
		using RVecFloat = ROOT::RVec<Float_t>;
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
	'''
	ROOT.gInterpreter.Declare(getLeadingTrailingCode)
	ROOT.gInterpreter.Declare(getDeltaRCode)
	ROOT.gInterpreter.Declare(getKinematicsCode)
	



