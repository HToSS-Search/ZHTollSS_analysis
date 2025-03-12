import ROOT
import uproot
import correctionlib
correctionlib.register_pyroot_binding()

def get_pileup(era,shift="nominal"): # code below taken from spark tnp
	'''
	Get the pileup distribution scalefactors to apply to simulation
	for a given era.
	'''

	# get the pileup
	dataPileup = {
		'UL2016_APV': 'pileup/UL2016/dataPileup_nominal.root',
		'UL2016': 'pileup/UL2016/dataPileup_nominal.root',
		'UL2017': 'pileup/UL2017/dataPileup_nominal.root',
		'UL2018': 'pileup/UL2018/dataPileup_nominal.root'
	}
	mcPileup = {
		'UL2016_APV': 'pileup/UL2016/mcPileup.root',
		'UL2016': 'pileup/UL2016/mcPileup.root',
		'UL2017': 'pileup/UL2017/mcPileup.root',
		'UL2018': 'pileup/UL2018/mcPileup.root'
	}
	# get absolute path
	#    baseDir = os.path.dirname(__file__)
	#    dataPileup = {k: os.path.join(baseDir, dataPileup[k]) for k in dataPileup}
	#    mcPileup = {k: os.path.join(baseDir, mcPileup[k]) for k in mcPileup}
	
	dataPileup_name=dataPileup[era].replace("nominal",shift)
	with uproot.open(dataPileup_name) as f:
		data_edges = f['pileup'].axis(0).edges()
		data_pileup = f['pileup'].values()
		data_pileup /= sum(data_pileup)
	with uproot.open(mcPileup[era]) as f:
		mc_edges = f['pileup'].axis(0).edges()
		mc_pileup = f['pileup'].values()
		mc_pileup /= sum(mc_pileup)
	pileup_edges = data_edges if len(data_edges) < len(mc_edges) else mc_edges
	pileup_ratio = [d/m if m else 1.0 for d, m in zip(
		data_pileup[:len(pileup_edges)-1], mc_pileup[:len(pileup_edges)-1])]

	return pileup_ratio, pileup_edges
def gInterpreter_SF():
    getSF_code ='''
		#include <algorithm> 
        // using FourVector = ROOT::Math::PxPyPzMVector;
        // using namespace ROOT::VecOps;
		Float_t getSF(TH1D* h2_sf, TString tag, Float_t mu1_x=-99)
        {
			// for getting entries from histogram
			//std::cout<<"Looking at: "<<tag.Data()<<std::endl;
			//std::cout<<"mu1_x,mu1_y:"<<mu1_x<<","<<mu1_y<<std::endl;

			
			int x_bin=-1;int y_bin=-1;
			int nbins=h2_sf->GetNbinsX();
			//std::cout<<"outside; y_bin,x_bin:"<<y_bin<<","<<x_bin<<std::endl;

			for (int i=1;i<=nbins;i++) {
				//std::cout<<"in x; y_bin,x_bin,lo_bin:"<<y_bin<<","<<x_bin<<","<<h2_sf->GetYaxis()->GetBinLowEdge(i+1)<<std::endl;
				if (mu1_x < h2_sf->GetXaxis()->GetBinLowEdge(i+1)) {
					x_bin = i;
					break;
				}
			}
			if (mu1_x > h2_sf->GetXaxis()->GetBinLowEdge(nbins+1)) x_bin=nbins;
			if (mu1_x < h2_sf->GetXaxis()->GetBinLowEdge(1)) x_bin=1;
			return h2_sf->GetBinContent(x_bin);
        }
		Float_t getSF(TH2D* h2_sf, TString tag, Float_t mu1_x=-99, Float_t mu1_y=-99, TString shift="nominal")
        {
			// for getting entries from histogram
			// std::cout<<"Looking at: "<<tag.Data()<<std::endl;
			// std::cout<<"mu1_x,mu1_y:"<<mu1_x<<","<<mu1_y<<std::endl;

			
			int x_bin=-1;int y_bin=-1;
			int nbins=h2_sf->GetNbinsX();
			// std::cout<<"outside; y_bin,x_bin,nbins:"<<y_bin<<","<<x_bin<<","<<nbins<<std::endl;

			for (int i=1;i<=nbins;i++) {
				// std::cout<<"in x; y_bin,x_bin,lo_bin:"<<y_bin<<","<<x_bin<<","<<h2_sf->GetYaxis()->GetBinLowEdge(i+1)<<std::endl;
				if (mu1_x < h2_sf->GetXaxis()->GetBinLowEdge(i+1)) {
					x_bin = i;
					break;
				}
			}
			if (mu1_x > h2_sf->GetXaxis()->GetBinLowEdge(nbins+1)) x_bin=nbins;
			if (mu1_x < h2_sf->GetXaxis()->GetBinLowEdge(1)) x_bin=1;
			nbins=h2_sf->GetNbinsY();
			for (int i=1;i<=nbins;i++) {
				// std::cout<<"in eta; y_bin,x_bin:"<<y_bin<<","<<x_bin<<std::endl;
				if (mu1_y < h2_sf->GetYaxis()->GetBinLowEdge(i+1)) {
					y_bin = i;
					break;
				}
			}
			if (mu1_y > h2_sf->GetYaxis()->GetBinLowEdge(nbins+1)) y_bin=nbins;
			if (mu1_y < h2_sf->GetYaxis()->GetBinLowEdge(1)) y_bin=1;
			// std::cout<<"y_bin,x_bin:"<<y_bin<<","<<x_bin<<std::endl;
			// std::cout<<"SF-"<<h2_sf->GetBinContent(x_bin,y_bin)<<std::endl;
			// if (tag.Contains("dR")) x_bin=1;
			if (h2_sf->GetBinContent(x_bin,y_bin)==0) {
				if (x_bin==1) x_bin=x_bin+1;
				else x_bin=x_bin-1;
			}
			if (shift.Contains("up")) 
				return h2_sf->GetBinContent(x_bin,y_bin)+h2_sf->GetBinError(x_bin,y_bin);
			else if (shift.Contains("down"))
				return h2_sf->GetBinContent(x_bin,y_bin)-h2_sf->GetBinError(x_bin,y_bin);
			else
				return h2_sf->GetBinContent(x_bin,y_bin);
        }
        Float_t getSF(std::shared_ptr<const correction::Correction> m_SF_map_Z_tag, std::shared_ptr<const correction::Correction> m_SF_map_JPsi_tag,  TString type, Float_t mu1_pt=-99, Float_t mu1_eta=-99, Float_t mu2_pt=-99,Float_t mu2_eta=-99)
        {
			Float_t mu1_sf=-1;
			Float_t mu2_sf=-1;
			//std::cout<<"Looking at: "<<type.Data()<<std::endl;
			//std::cout<<"mu1_pt,mu1_eta,mu2_pt,mu2_eta:"<<mu1_pt<<","<<mu1_eta<<","<<mu2_pt<<","<<mu2_eta<<std::endl;
			if (type.Contains("reco")) {
				std::cout<<"problem: function not overloading"<<std::endl;
			}
			if (type.Contains("id")) {
				mu1_pt = mu1_pt>=2 ? mu1_pt : 2;
				mu2_pt = mu2_pt>=2 ? mu2_pt : 2;
				mu1_eta = std::abs(mu1_eta)<2.4 ? mu1_eta : 2.39;
				mu2_eta = std::abs(mu2_eta)<2.4 ? mu2_eta : 2.39;
				mu1_sf = mu1_pt>=15 ? m_SF_map_Z_tag->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"}) : m_SF_map_JPsi_tag->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"});
				mu2_sf = mu2_pt>=15 ? m_SF_map_Z_tag->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"}) : m_SF_map_JPsi_tag->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"});
			}
			if (type.Contains("iso")) {
				mu1_pt=mu1_pt >= 15 ? mu1_pt : 15;
				mu2_pt=mu2_pt >= 15 ? mu2_pt : 15;
				mu1_eta = std::abs(mu1_eta)<2.4 ? mu1_eta : 2.39;
				mu2_eta = std::abs(mu2_eta)<2.4 ? mu2_eta : 2.39;
				mu1_sf = m_SF_map_Z_tag->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"});
				mu2_sf = m_SF_map_Z_tag->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"});
			}if (type.Contains("iso")) {
				mu1_pt=mu1_pt >= 15 ? mu1_pt : 15;
				mu2_pt=mu2_pt >= 15 ? mu2_pt : 15;
				mu1_pt=mu1_pt < 100 ? mu1_pt : 99.9;
				mu2_pt=mu2_pt < 100 ? mu2_pt : 99.9;
				mu1_eta = std::abs(mu1_eta)<2.4 ? mu1_eta : 2.39;
				mu2_eta = std::abs(mu2_eta)<2.4 ? mu2_eta : 2.39;
				mu1_sf = m_SF_map_Z_tag->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"});
				mu2_sf = m_SF_map_Z_tag->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"});
			}
			if (type.Contains("trg")) {
				//mu1_pt=mu1_pt > 26 ? mu1_pt : 26.1;
				if (type.Contains("2017")) {
					mu1_pt=mu1_pt > 29 ? mu1_pt : 29.1;
				}
				else {
					mu1_pt=mu1_pt > 26 ? mu1_pt : 26.1;
				}
				mu1_pt=mu1_pt < 200 ? mu1_pt : 199.9;
				mu1_eta = std::abs(mu1_eta)<2.4 ? mu1_eta : 2.39;
				mu1_sf = m_SF_map_Z_tag->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"});
			}
			if (type.Contains("trg")) return mu1_sf;
			else return mu1_sf*mu2_sf;
        }
		

		Float_t getSF(std::shared_ptr<const correction::Correction> m_SF_map_Z_tag, TH2D* h2_sf, TString type, Float_t mu1_pt=-99, Float_t mu1_eta=-99, Float_t mu2_pt=-99,Float_t mu2_eta=-99)
        {
			//std::cout<<"Looking at: "<<type.Data()<<std::endl;
			//std::cout<<"mu1_pt,mu1_eta,mu2_pt,mu2_eta:"<<mu1_pt<<","<<mu1_eta<<","<<mu2_pt<<","<<mu2_eta<<std::endl;
			Float_t mu1_sf=-1;
			Float_t mu2_sf=-1;
			if (type.Contains("reco")) {
				mu1_eta = std::abs(mu1_eta)<2.4 ? mu1_eta : 2.39;
				mu2_eta = std::abs(mu2_eta)<2.4 ? mu2_eta : 2.39;
				mu1_pt = mu1_pt<200? mu1_pt : 199.9;
				mu2_pt = mu2_pt<200? mu2_pt : 199.9;
				if (mu1_pt < 10) {
					mu1_sf = getSF(h2_sf,type,std::abs(mu1_eta),mu1_pt);
				}
				else {
					mu1_pt = mu1_pt>=40 ? mu1_pt:40; // bin is 40-infty
					mu1_sf = m_SF_map_Z_tag->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"});
				}
				if (mu2_pt<10) {
					mu2_sf = getSF(h2_sf,type,mu2_pt,mu2_eta);
				}
				else {
					mu2_pt = mu2_pt>=40 ? mu2_pt:40; // bin is 40-infty
					mu2_sf = m_SF_map_Z_tag->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"});
				}
			}
			if (type.Contains("id")) {
				std::cout<<"problem: function not overloading"<<std::endl;
			}
			if (type.Contains("iso")) {
				std::cout<<"problem: function not overloading"<<std::endl;
			}
			return mu1_sf*mu2_sf;
		}
    '''
    ROOT.gInterpreter.Declare(getSF_code)
# def get_muon_sf(era): # code below taken from spark tnp
# 	'''
# 	Get the muon scalefactors to apply to simulation
# 	for a given era.
# 	'''
# 	muonreco_files= {
# 		'UL2016_APV': {'Z':'scale_factors/UL2016_APV/muon_Z.json','JPsi':'scale_factors/UL2016_APV/muon_JPsi.json'},
# 		'UL2016': {'Z':'scale_factors/UL2016/muon_Z.json','JPsi':'scale_factors/UL2016/muon_JPsi.json'},
# 		'UL2017': {'Z':'scale_factors/UL2017/Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.json','JPsi':'scale_factors/UL2017/Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.json'},
# 		'UL2018': {'Z':'scale_factors/UL2018/muon_Z.json','JPsi':'scale_factors/UL2018/muon_JPsi.json'}
# 	}
# 	muonIDiso_files= {
# 		'UL2016_APV': {'Z':'scale_factors/UL2016_APV/muon_Z.json','JPsi':'scale_factors/UL2016_APV/muon_JPsi.json'},
# 		'UL2016': {'Z':'scale_factors/UL2016/muon_Z.json','JPsi':'scale_factors/UL2016/muon_JPsi.json'},
# 		'UL2017': {'Z':'scale_factors/UL2017/muon_Z.json','JPsi':'scale_factors/UL2017/muon_JPsi.json'},
# 		'UL2018': {'Z':'scale_factors/UL2018/muon_Z.json','JPsi':'scale_factors/UL2018/muon_JPsi.json'}
# 	}
# 	trigSF_files= {
# 		'UL2017': {'Z':'scale_factors/UL2017/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers_schemaV2.json','JPsi':'scale_factors/UL2017/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers_schemaV2.json'}
# 	}
# 	fnames=muonSF_files[era]
# 	fname2=trigSF_files[era]
# 	ROOT.gInterpreter.Declare(f'auto m_SF_map_Z = correction::CorrectionSet::from_file("{fnames['Z']}");')
# 	ROOT.gInterpreter.Declare(f'auto m_SF_map_JPsi = correction::CorrectionSet::from_file("{fnames['JPsi']}");')
# 	# ROOT.gInterpreter.Declare(f'auto m_SF_map = correction::CorrectionSet::from_file({"fname"});')
# 	ROOT.gInterpreter.Declare('auto m_SF_map_UL_ID =  m_SF_map_Z->at("NUM_LooseID_DEN_TrackerMuons");')
# 	ROOT.gInterpreter.Declare('auto m_SF_map_UL_ISO = m_SF_map_Z->at("NUM_LooseRelIso_DEN_LooseID");')
# 	ROOT.gInterpreter.Declare(f'auto trg_SF_map = correction::CorrectionSet::from_file("{fname2}");')
# 	ROOT.gInterpreter.Declare('auto trg_SF_map_UL_TRG = trg_SF_map->at("NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight");')
# 	sf_string = '''
# 		if (mu2_pt > 15) {
# 			m_SF_map_UL_ID->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"})*m_SF_map_UL_ID->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"})
# 		}
# 		else {
# 			m_SF_map_UL_ID->evaluate({std::abs(mu1_eta),mu1_pt,"nominal"})*m_SF_map_UL_ID->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"})
# 		}
# 	'''
# 	# *
# 	# 	m_SF_map_UL_ID->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"})*
# 	# 	m_SF_map_UL_ISO->evaluate({std::abs(mu2_eta),mu2_pt,"nominal"})
# 	# csetEl_2016preID->evaluate({"2016preVFP", "sf", "wp90iso", std::abs(Electron_eta[0]), Electron_pt[0]})
# 	# # '''
# 	# ('csetEl_2016preID->evaluate({"2016preVFP", "sf", "wp90iso", std::abs(Electron_eta[0]), Electron_pt[0]})')
# 	# getSF_code = '''
#     #     using FourVectorPtEtaPhiEVector = ROOT::Math::PtEtaPhiEVector;
#     #     Float_t getMuonSF(float old_lt, float new_lt, float L_scalar1, float L_scalar2, FourVectorPtEtaPhiEVector& vec_scalar1, FourVectorPtEtaPhiEVector& vec_scalar2){
            
#     #         return weight;
#     #     }
#     # '''

# 	return sf_string