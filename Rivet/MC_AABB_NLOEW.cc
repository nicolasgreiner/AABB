// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include <csignal>
#include <math.h>


namespace Rivet {

  using namespace Cuts;

  class MC_AABB_NLOEW : public Analysis {
  private:

    class NLOHisto1D : public YODA::Histo1D {
    private:

      YODA::Histo1D* _tmphist;
      int _current_event_number;

      void _syncHists() {
        for (size_t i=0; i<_tmphist->bins().size(); ++i) {
          if (_tmphist->bin(i).area()) YODA::Histo1D::fillBin(i, _tmphist->bin(i).area());
        }
        if (_tmphist->overflow().sumW())  YODA::Histo1D::overflow()+=_tmphist->overflow();
        if (_tmphist->underflow().sumW()) YODA::Histo1D::underflow()+=_tmphist->underflow();
        _tmphist->reset();
      }

    public:

      NLOHisto1D(size_t nbins, double lower, double upper, const string& path) :
        YODA::Histo1D(nbins, lower, upper, path),
        _current_event_number(-1)
      {
        _tmphist = new Histo1D(nbins, lower, upper, path+"_tmp");
      }

      NLOHisto1D(const vector<double>& binedges, const string& path) :
        YODA::Histo1D(binedges, path),
        _current_event_number(-1)
      {
        _tmphist = new Histo1D(binedges, path+"_tmp");
      }

      ~NLOHisto1D()
      {
        delete _tmphist;
      }

      void fill(double x, const Event& event)
      {
        if (_current_event_number==-1)
          _current_event_number = event.genEvent()->event_number();

        if (event.genEvent()->event_number()!=_current_event_number) {
          _syncHists();
          _current_event_number = event.genEvent()->event_number();
        }

        _tmphist->fill(x, event.weight());
      }

      void fillBin(size_t i, const Event& event, const double& fac)
      {
        if (_current_event_number==-1)
          _current_event_number = event.genEvent()->event_number();

        if (event.genEvent()->event_number()!=_current_event_number) {
          _syncHists();
          _current_event_number = event.genEvent()->event_number();
        }

        _tmphist->fillBin(i, event.weight()*fac);
      }

      void finalize()
      {
        _syncHists();
      }

    };

    typedef shared_ptr<NLOHisto1D> NLOHisto1DPtr;


    NLOHisto1DPtr bookNLOHisto1D(const string& hname,
                                      size_t nbins, double lower, double upper)
    {
      NLOHisto1DPtr hist(new NLOHisto1D(nbins, lower, upper, histoPath(hname)));
      addAnalysisObject(hist);
      return hist;
    }

    NLOHisto1DPtr bookNLOHisto1D(const string& hname,
                                 const vector<double>& binedges)
    {
      NLOHisto1DPtr hist(new NLOHisto1D(binedges, histoPath(hname)));
      addAnalysisObject(hist);
      return hist;
    }

    class NLOHisto2D : public YODA::Histo2D {

      YODA::Histo2D* _tmphist;
      int _current_event_number;

      void _syncHists() {
	for (size_t i=0; i<_tmphist->bins().size(); ++i) {
	  if (_tmphist->bin(i).volume()) YODA::Histo2D::fillBin(i, _tmphist->bin(i).volume());
	}
	// enable once YODA::Histo2D possesses overflow, underflow or outflow
	// if (_tmphist->outflow().sumW())  YODA::Histo2D::outflow()+=_tmphist->outflow();
      }

    public:

      NLOHisto2D(size_t nxbins, double xlower, double xupper,
		size_t nybins, double ylower, double yupper, const string& path) :
	YODA::Histo2D(nxbins, xlower, xupper, nybins, ylower, yupper, path),
	_current_event_number(-1)
      {
	_tmphist = new Histo2D(nxbins, xlower, xupper,
			      nybins, ylower, yupper, path+"_tmp");
      }
      
      NLOHisto2D(const vector<double>& xbinedges,
		const vector<double>& ybinedges, const string& path) :
	YODA::Histo2D(xbinedges, ybinedges, path),
	_current_event_number(-1)
      {
	_tmphist = new Histo2D(xbinedges, ybinedges, path+"_tmp");
      }

      // NLOHisto2D(const Scatter2D& refscatter, const string& path) :
      //   YODA::Histo2D(refscatter, path),
      //   _current_event_number(-1)
      // {
      //   _tmphist = new Histo2D(refscatter, path+"_tmp");
      // }

      ~NLOHisto2D()
      {
	delete _tmphist;
      }
	
      void fill(double x, double y, const Event& event, const bool& useweight=true)
      {
	if (_current_event_number==-1)
	  _current_event_number = event.genEvent()->event_number();

	if (event.genEvent()->event_number()!=_current_event_number) {
	  _syncHists();
	  _current_event_number = event.genEvent()->event_number();
	}

	_tmphist->fill(x, y, useweight?event.weight():1.);
      }

      void fillBin(size_t i, const Event& event, const double& fac)
      {
	if (_current_event_number==-1)
	  _current_event_number = event.genEvent()->event_number();

	if (event.genEvent()->event_number()!=_current_event_number) {
	  _syncHists();
	  _current_event_number = event.genEvent()->event_number();
	}

	_tmphist->fillBin(i, event.weight()*fac);
      }

      void finalize()
      {
	_syncHists();
      }
    };

    typedef shared_ptr<NLOHisto2D> NLOHisto2DPtr;


    NLOHisto2DPtr bookNLOHisto2D(const string& hname,
				size_t nxbins, double xlower, double xupper,
				size_t nybins, double ylower, double yupper)
    {
      NLOHisto2DPtr hist( new NLOHisto2D(nxbins, xlower, xupper,
					nybins, ylower, yupper, histoPath(hname)) );
      addAnalysisObject(hist);
      return hist;
    }

    NLOHisto2DPtr bookNLOHisto2D(const string& hname,
				const vector<double>& xbinedges,
				const vector<double>& ybinedges)
    {
      NLOHisto2DPtr hist( new NLOHisto2D(xbinedges, ybinedges, histoPath(hname)) );
      addAnalysisObject(hist);
      return hist;
    }

    // NLOHisto2DPtr bookNLOHisto2D(const string& hname,
    //                              const Scatter2D& refscatter) {
    //   NLOHisto2DPtr hist( new NLOHisto2D(refscatter, histoPath(hname)) );
    //   addAnalysisObject(hist);
    //   return hist;
    // }

    // NLOHisto2DPtr bookNLOHisto2D(const string& hname)
    // {
    //   const Scatter2D& refdata = refData(hname);
    //   return bookNLOHisto2D(hname, refdata);
    // }

    // NLOHisto2DPtr bookNLOHisto2D(unsigned int datasetId, unsigned int xAxisId,
    //                              unsigned int yAxisId)
    // {
    //   const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    //   return bookNLOHisto2D(axisCode);
    // }
  private:

    /// Variables
    double _dR_dress;
    double _etay_max, _pty1_min, _pty2_min, _dRyy_min, _dRyj_min, _zyj;
    double _etaj_max, _ptj_min, _dRj;
    fastjet::JetDefinition _jdef;

    /// Histograms
    std::map <string, NLOHisto1DPtr> _histos;
    std::map <string, NLOHisto2DPtr> _histos2D;

  public:
    /// Constructor
    MC_AABB_NLOEW()
      : Analysis("MC_AABB_NLOEW"),
        _dR_dress(0.1), _etay_max(2.5), _pty1_min(30.*GeV), 
        _pty2_min(30.*GeV), _dRyy_min(0.4), _dRyj_min(0.4), _zyj(0.5),
        _etaj_max(4.4), _ptj_min(20.*GeV), _dRj(0.4),
        _jdef(fastjet::antikt_algorithm, _dRj)
   {    }

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections

      //photons
      const IdentifiedFinalState photon(PID::PHOTON);
      addProjection(photon, "photon");

      //partons
      vector<PdgId> quark_ids;
      quark_ids +=  PID::DQUARK, PID::UQUARK, PID::SQUARK, PID::CQUARK, PID::BQUARK,
                   -PID::DQUARK,-PID::UQUARK,-PID::SQUARK,-PID::CQUARK,-PID::BQUARK;
      const IdentifiedFinalState quark(etaIn(-6.0, 6.0) & (pT >= 0*GeV), quark_ids);
      addProjection(quark, "quark");

      const IdentifiedFinalState bquark(etaIn(-6.0, 6.0) & (pT >= 0*GeV),  PID::BQUARK );
      addProjection(bquark, "bquark");

      const IdentifiedFinalState bbarquark(etaIn(-6.0, 6.0) & (pT >= 0*GeV),  -PID::BQUARK );
      addProjection(bbarquark, "bbarquark");
      
      vector<PdgId> b_ids;
      b_ids += PID::BQUARK, -PID::BQUARK;
      const IdentifiedFinalState b_or_bbar(etaIn(-6.0, 6.0) & (pT >= 0*GeV),  b_ids);
      addProjection(b_or_bbar, "b_or_bbar");                  

      const IdentifiedFinalState gluon(etaIn(-6.0, 6.0) & (pT >= 0*GeV),  PID::GLUON );
      addProjection(gluon, "gluon");
      
      _histos["xs_inc"]    = bookNLOHisto1D("xs_inc",1,-0.5,0.5);
      _histos["xs_excl"]   = bookNLOHisto1D("xs_excl",4,-0.5,3.5);

      _histos["pT_yy"]     = bookNLOHisto1D("pT_yy", 100, 0., 200.);
      _histos["pT_yy_log"] = bookNLOHisto1D("pT_yy_log", logspace(20, 2.e1, 2.e3));

      _histos["m_yy"]     = bookNLOHisto1D("m_yy", 50, 0., 200.);
      _histos["m_yy_log"] = bookNLOHisto1D("m_yy_log", logspace(30, 1.e1, 6.e3));

      _histos["phistar_yy_log"] = bookNLOHisto1D("phistar_yy_log", logspace(45, 1.e-3,1.e6));
      
      _histos["pT_y1"]    = bookNLOHisto1D("pT_y1", 100, 0., 500.);
      _histos["pT_y2"]    = bookNLOHisto1D("pT_y2", 100, 0., 500.);
      
      _histos["pT_j1"]    = bookNLOHisto1D("pT_j1", 100, 0., 500.);
      _histos["pT_j2"]    = bookNLOHisto1D("pT_j2", 100, 0., 500.);
      
      _histos["m_jj"]     = bookNLOHisto1D("m_jj", 50, 0., 200.);
      _histos["m_jj_log"] = bookNLOHisto1D("m_jj_log", logspace(30, 1.e1, 6.e3));      
      
      _histos["eta_y1"]    = bookNLOHisto1D("eta_y1", 20, -2.5, 2.5);
      _histos["eta_y2"]    = bookNLOHisto1D("eta_y2", 20, -2.5, 2.5);      
      
      _histos["eta_j1"]    = bookNLOHisto1D("eta_j1", 50, -5., 5.);
      _histos["eta_j2"]    = bookNLOHisto1D("eta_j2", 50, -5., 5.);      
      
      _histos["deta_y1y2"]    = bookNLOHisto1D("deta_y1y2", 50, 0., 10.);
      _histos["deta_j1j2"]    = bookNLOHisto1D("deta_j1j2", 50, 0., 10.);      
      
      _histos["dR_y1y2"]    = bookNLOHisto1D("dR_y1y2", 50, 0., 10.);
      _histos["dR_j1j2"]    = bookNLOHisto1D("dR_j1j2", 50, 0., 10.);            
      _histos["dR_y1j1"]    = bookNLOHisto1D("dR_y1j1", 50, 0., 10.);
      _histos["dR_y1j2"]    = bookNLOHisto1D("dR_y1j2", 50, 0., 10.);                  
      _histos["dR_y2j1"]    = bookNLOHisto1D("dR_y2j1", 50, 0., 10.);                  
      _histos["dR_y2j2"]    = bookNLOHisto1D("dR_y2j2", 50, 0., 10.);                  
      
      _histos["dphi_y1y2"]    = bookNLOHisto1D("dphi_y1y2", 20, 0., 3.14);
      _histos["dphi_j1j2"]    = bookNLOHisto1D("dphi_j1yj", 20, 0., 3.14);                        

      _histos2D["pT_yy_m_yy"]         = bookNLOHisto2D("pT_yy_m_yy", 100, 0., 200., 50, 0., 200.);
      _histos2D["pT_yy_log_m_yy_log"] = bookNLOHisto2D("pT_yy_log_m_yy_log", logspace(20, 2.e1, 2.e3),logspace(30, 1.e1, 6.e3));

    }


    /// Perform the per-event analysis
    void analyze(const Event& evt) {

      // get the partons
      ParticleVector cand_quarks
	= applyProjection<IdentifiedFinalState>(evt, "quark").particlesByPt();
      ParticleVector cand_gluons
	= applyProjection<IdentifiedFinalState>(evt, "gluon").particlesByPt();

      // get the photons
      ParticleVector cand_photons
	= applyProjection<IdentifiedFinalState>(evt, "photon").particlesByPt();
      ParticleVector cand_photons_cuts;
      ParticleVector cand_photons_remain;      
        

      if (cand_photons.size() < 2) vetoEvent;
      
      //determine the minimal required number of jets. Note that this is easy here as the real emission
      //contains only a photon so the minimal number of jets is give by the number of the partons
      double njets_min(cand_quarks.size()+cand_gluons.size());
      
      
      //check pT and eta cut of the first two photons 
      if ( cand_photons[0].pt() < _pty1_min || fabs(cand_photons[0].eta()) > _etay_max  \
          || cand_photons[1].pt() < _pty2_min || fabs(cand_photons[1].eta()) > _etay_max ) vetoEvent;
      else {
         cand_photons_cuts.push_back(cand_photons[0]);
         cand_photons_cuts.push_back(cand_photons[1]);
      }
      
      //it might be that there are more than two photons fulfiling the cuts. We only want to check the isolation 
      //with the ones that pass the cuts
      if (cand_photons.size() >2 ) {
        for (size_t i(2); i<cand_photons.size(); ++i){
          //std::cout<<i<<std::endl;
          if (cand_photons[i].pt() > _pty2_min && fabs(cand_photons[i].eta()) < _etay_max  )  cand_photons_cuts.push_back(cand_photons[i]);
          else cand_photons_remain.push_back(cand_photons[i]);
        }
      }

      // build isolated photons with the candidates that passed the cuts. If partons are found within a cone they
      // are removed from the list of possible jet candidates, independent whether the photon is found to be isolated or not
      ParticleVector isolatedphotons;
      ParticleVector cand_quarks_remain;
      ParticleVector cand_gluons_remain;
      double eps(0.05);
      double R_iso(0.4);
      bool isolated(true);

      foreach (const Particle & p, cand_photons_cuts) {
          isolated=true;

          for (size_t i(0); i<cand_quarks.size(); ++i) {
            if (Rivet::deltaR(p.momentum(), cand_quarks[i].momentum()) < R_iso )
             {
              double Ehad_max=eps*p.pt()*(1.0-cos(Rivet::deltaR(p.momentum(), cand_quarks[i].momentum())))/(1.0-cos(R_iso));              
              if (Ehad_max<cand_quarks[i].pt()) {
                  isolated=false;  
              //if isolation fails, both photon and parton are kept and passed to the jet algorithm 
                  cand_photons_remain.push_back(p);
                  cand_quarks_remain.push_back(cand_quarks[i]);
                  
              }
             }
            else cand_quarks_remain.push_back(cand_quarks[i]);
           }

           for (size_t i(0); i<cand_gluons.size(); ++i) {
            if (Rivet::deltaR(p.momentum(), cand_gluons[i].momentum()) < R_iso )
             {
              double Ehad_max=eps*p.pt()*(1.0-cos(Rivet::deltaR(p.momentum(), cand_gluons[i].momentum())))/(1.0-cos(R_iso));
              if (Ehad_max<cand_gluons[i].pt()) {
                isolated=false;  
              //if isolation fails, both photon and parton are kept and passed to the jet algorithm 
                  cand_photons_remain.push_back(p);
                  cand_gluons_remain.push_back(cand_gluons[i]);                
              } 
             }
            else cand_gluons_remain.push_back(cand_gluons[i]);
           } 
            
           if (isolated==true) isolatedphotons.push_back(p);

         
       }
      if (isolatedphotons.size() < 2) vetoEvent;

      if(isolatedphotons.size()==2){
       if ( Rivet::deltaR(isolatedphotons[0].momentum(),isolatedphotons[1].momentum())<0.4) vetoEvent;
      }
      if(isolatedphotons.size()==3){
        if ( Rivet::deltaR(isolatedphotons[0].momentum(),isolatedphotons[1].momentum())<0.4) vetoEvent;
        if ( Rivet::deltaR(isolatedphotons[0].momentum(),isolatedphotons[2].momentum())<0.4) vetoEvent;
        if ( Rivet::deltaR(isolatedphotons[1].momentum(),isolatedphotons[2].momentum())<0.4) vetoEvent;
      }
      
      if ( isolatedphotons[0].momentum().pt() < _pty1_min || fabs(isolatedphotons[0].momentum().eta()) > _etay_max) vetoEvent;      
      //std::cout<<cand_photons_remain.size()<<std::endl;


      // dress the quarks or gluons with the remaining photons
      // works reliably only if at most one photon to cluster
      // (otherwise first survey the closest photons for each quark first, 
      //  then dress keeping the axis for dR fixed)
      
      //if (cand_photons.size() > 1) vetoEvent;
//        foreach (const Particle & p, cand_photons_remain) {
//          size_t imin(0);
// 	double dRmin_quark(100000.);
//         double dRmin_gluon(100000.);
//         for (size_t i(0); i<cand_quarks_remain.size(); ++i) {
//           double dRi(Rivet::deltaR(p.momentum(), cand_quarks_remain[i].momentum()));
//           if (dRi < _dR_dress || dRi<dRmin_quark) {
// 	    imin=i;
// 	    dRmin_quark=dRi;
// 	  }
//         }        
//         for (size_t i(0); i<cand_gluons_remain.size(); ++i) {
//          double dRi(Rivet::deltaR(p.momentum(), cand_gluons_remain[i].momentum()));
//           if (dRi < _dR_dress || dRi<dRmin_gluon) {
// 	    imin=i;
// 	    dRmin_gluon=dRi;
// 	  }
//         }        
//         
//       //choose the closest 
//         if (dRmin_quark < dRmin_gluon){      
//           if( cand_quarks_remain.size() >0 ){
// 	    cand_quarks_remain[imin].setMomentum(cand_quarks_remain[imin].momentum()+p.momentum());
//            }
//          }
//         else {
//           if( cand_gluons_remain.size() >0 ){
// 	    cand_gluons_remain[imin].setMomentum(cand_gluons_remain[imin].momentum()+p.momentum());            
//         }   
// 
//        }

      //Jet/Photon-Finder
      std::vector<fastjet::PseudoJet> input;
      foreach (const Particle & p, cand_quarks_remain) {
	input.push_back(p.pseudojet());
	input.back().set_user_index(p.pdgId());
      }
      foreach (const Particle & p, cand_gluons_remain) {
	input.push_back(p.pseudojet());
	input.back().set_user_index(p.pdgId());
      }
      foreach (const Particle & p, cand_photons_remain) {
	input.push_back(p.pseudojet());
	input.back().set_user_index(p.pdgId());
      }

      // calculate jets
      fastjet::ClusterSequence cs(input, _jdef);
      Jets jets;
      bool contains_photon;
      bool contains_gluon;
      //bool contains_quark;


      int nr_bjets(0), nr_b(0), nr_bbar(0);
      nr_bjets=0;

      foreach (const PseudoJet& jetcand, sorted_by_pt(cs.inclusive_jets())) {
        nr_b=0;
        nr_bbar=0;
	// check eta and pt
	if (fabs(jetcand.eta()) < _etaj_max && jetcand.perp() > _ptj_min) {
	  // check constituents 
          contains_photon=false;
          contains_gluon=false;
          //contains_quark=false;
           FourMomentum py(0.,0.,0.,0.);
           foreach (const PseudoJet& con, jetcand.constituents()){
           if (con.user_index()==PID::GLUON) contains_gluon=true;
           //if (con.user_index()==quark) contains_quark=true;
           if (con.user_index()==PID::PHOTON) {
               contains_photon=true;
               py += FourMomentum(con.e(), con.px(), con.py(), con.pz());
             }
           if (con.user_index()==PID::BQUARK) nr_b+=1;
           if (con.user_index()==-PID::BQUARK) nr_bbar+=1;
           }
          
          if ( (nr_b>0 || nr_bbar>0) && !(nr_b>0 && nr_bbar>0) ) nr_bjets+=1;
          
          if (contains_photon==true && contains_gluon==true ){
 	       //if (py.E() < _zyj*jetcand.e()) jets.push_back(jetcand);
             if (py.E() > _zyj*jetcand.e()) vetoEvent;
             }           
          if ( (nr_b>0 || nr_bbar>0) && !(nr_b>0 && nr_bbar>0) ) jets.push_back(jetcand);
          }          
	}


      if (jets.size()<njets_min) vetoEvent;
      if (nr_bjets <2 ) vetoEvent;
      
      
      FourMomentum y1(isolatedphotons[0].momentum());
      FourMomentum y2(isolatedphotons[1].momentum());
      
//       if (isolatedphotons[0].momentum().pt() < isolatedphotons[1].momentum().pt() ){//|| isolatedphotons[0].pt() <isolatedphotons[2].pt() ){
//       std::cout<<isolatedphotons[0].momentum().pt()<<" "<<isolatedphotons[1].momentum().pt()<<" "<<isolatedphotons[2].momentum().pt()<<std::endl;
//       }
                

      double ptyy((y1+y2).pt());
      double myy((y1+y2).mass());
      double pty1(y1.pt());
      double pty2(y2.pt());
      double etay1(y1.eta());
      double etay2(y2.eta());
      double detay1y2(Rivet::deltaEta(y1,y2));
      double dry1y2(Rivet::deltaR(y1,y2));
      double dphiy1y2(Rivet::deltaPhi(y1,y2));
      double costhetastaryy(tanh(0.5*std::abs(etay1-etay2)));
      double sinthetastaryy((costhetastaryy>1)?0.:sqrt(1.-sqr(costhetastaryy)));
      double phistaryy(tan(0.5*(M_PI-dphiy1y2))*sinthetastaryy);
      //std::cout<<evt.weight()<<std::endl;
      
      _histos["xs_inc"]->fill(0.,evt);

      _histos["pT_yy"]->fill(ptyy,evt);
      _histos["pT_yy_log"]->fill(ptyy,evt);
      _histos["m_yy"]->fill(myy,evt);
      _histos["m_yy_log"]->fill(myy,evt);
      _histos["phistar_yy_log"]->fill(phistaryy,evt);
      _histos["pT_y1"]->fill(pty1,evt);
      _histos["pT_y2"]->fill(pty2,evt);  
      _histos["eta_y1"]->fill(etay1,evt);
      _histos["eta_y2"]->fill(etay2,evt);
      _histos["deta_y1y2"]->fill(detay1y2,evt);
      _histos["dR_y1y2"]->fill(dry1y2,evt);
      _histos["dphi_y1y2"]->fill(dphiy1y2,evt);

      _histos2D["pT_yy_m_yy"]->fill(ptyy,myy,evt);
      _histos2D["pT_yy_log_m_yy_log"]->fill(ptyy,myy,evt);
      
      
     if (jets.size()==0) _histos["xs_excl"]->fill(0.,evt); 
      
     if (jets.size()==1){
        FourMomentum j1(jets[0].momentum());
        _histos["xs_excl"]->fill(1.,evt); 
        double ptj1(j1.pt());
        double etaj1(j1.y());
        double dry1j1(Rivet::deltaR(y1,j1));
        double dry2j1(Rivet::deltaR(y2,j1));
        _histos["pT_j1"]->fill(ptj1,evt);
        _histos["eta_j1"]->fill(etaj1,evt);
        _histos["dR_y1j1"]->fill(dry1j1,evt);
        _histos["dR_y2j1"]->fill(dry2j1,evt);
          
      }
      
      if (jets.size()>1){
        if (jets.size()==2) _histos["xs_excl"]->fill(2.,evt);
        if (jets.size()==3) _histos["xs_excl"]->fill(3.,evt); 
        FourMomentum j1(jets[0].momentum());
        FourMomentum j2(jets[1].momentum());
        double ptj1(j1.pt());
        double ptj2(j2.pt());
        double etaj1(j1.y());
        double etaj2(j2.y());
        double detaj1j2(Rivet::deltaRap(j1,j2));
        double dphij1j2(Rivet::deltaPhi(j1,j2));
        double drj1j2(Rivet::deltaR(j1,j2));
        double dry1j1(Rivet::deltaR(y1,j1));
        double dry1j2(Rivet::deltaR(y1,j2));
        double dry2j1(Rivet::deltaR(y2,j1));
        double dry2j2(Rivet::deltaR(y2,j2));
        double mjj((j1+j2).mass());        
        _histos["pT_j1"]->fill(ptj1,evt);
        _histos["pT_j2"]->fill(ptj2,evt);
        _histos["eta_j1"]->fill(etaj1,evt);
        _histos["eta_j2"]->fill(etaj2,evt);
        _histos["deta_j1j2"]->fill(detaj1j2,evt);
        _histos["dphi_j1j2"]->fill(dphij1j2,evt);
        _histos["dR_j1j2"]->fill(drj1j2,evt);
        _histos["dR_y1j1"]->fill(dry1j1,evt);
        _histos["dR_y1j2"]->fill(dry1j2,evt);
        _histos["dR_y2j1"]->fill(dry2j1,evt);
        _histos["dR_y2j2"]->fill(dry2j2,evt);
        _histos["m_jj"]->fill(mjj,evt);
        _histos["m_jj_log"]->fill(mjj,evt);        
        
        
       }
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = crossSection()/picobarn/sumOfWeights();

      typedef std::map<string, NLOHisto1DPtr>::iterator it;
      typedef std::map<string, NLOHisto2DPtr>::iterator it2D;

      for(it iter = _histos.begin(); iter != _histos.end(); iter++) {
        iter->second->finalize();
        scale(iter->second,norm);
      }
       for(it2D iter = _histos2D.begin(); iter != _histos2D.end(); iter++) {
        iter->second->finalize();
        scale(iter->second,norm);
       }

    }

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_AABB_NLOEW);

}
