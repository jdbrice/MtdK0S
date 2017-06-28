#ifndef PICO_DST_SKIMMER_H
#define PICO_DST_SKIMMER_H


#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"

#include "TreeAnalyzer.h"
#include "XmlRange.h"

#include "PicoDstP16id/StPicoEvent.h"
#include "PicoDstP16id/StPicoMtdHit.h"
#include "PicoDstP16id/StPicoTrack.h"
#include "PicoDstP16id/StPicoMtdPidTraits.h"
#include "PicoDstP16id/StPicoBTofPidTraits.h"

#include "vendor/loguru.h"

#include "TLorentzVector.h"

class PicoDstSkimmer : public TreeAnalyzer
{
public:

	const float MASS_PI = 0.13957; // in GeV
	const float MASS_MU = 0.105658; // in GeV
	bool tpcOnlyPairs = false;
	bool makeLikeSign = false;

	struct TrackProxy {
		StPicoTrack * _track = nullptr;
		StPicoMtdPidTraits *_mtdPid = nullptr;
		TLorentzVector _lv;
		StThreeVectorD _p;
		StPhysicalHelixD _h;

	};


	const int DEBUG = 1;
	PicoDstSkimmer() {}
	~PicoDstSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( this->chain, "Event" );
		_rMtdHit.setup( this->chain, "MtdHit" );
		_rTrack.setup( this->chain, "Tracks" );
		_rMtdPid.setup( this->chain, "MtdPidTraits" );
		_rBTofPid.setup( this->chain, "BTofPidTraits" );

		vector<string> kids = config.childrenOf( nodePath + ".K0S", "XmlRange" );
		for ( auto p : kids ){
			XmlRange  xr( config, p );
			LOG_F( INFO, "[%s] = %s", config[ p +":name" ].c_str(), xr.toString().c_str()  );
			K0SCuts[ config[ p +":name" ] ] = xr;
		}

	}
protected:

	map< string, XmlRange > K0SCuts;

	TClonesArrayReader < StPicoEvent        > _rEvent;
	TClonesArrayReader < StPicoMtdHit       > _rMtdHit;
	TClonesArrayReader < StPicoTrack        > _rTrack;
	TClonesArrayReader < StPicoMtdPidTraits > _rMtdPid;
	TClonesArrayReader < StPicoBTofPidTraits > _rBTofPid;

	vector< shared_ptr<TrackProxy> > pip;
	vector< shared_ptr<TrackProxy> > pim;

	vector< shared_ptr<TrackProxy> > mup;
	vector< shared_ptr<TrackProxy> > mum;


	void fillMTD( shared_ptr<TrackProxy> t1, shared_ptr<TrackProxy> t2, TLorentzVector &lv ){
		bool is_K0S = fabs( lv.M() - 0.497 ) < 0.02;
		fillMTD( t1 );
		if ( is_K0S ) 
			fillMTD( t1, "k0s_" );
		fillMTD( t2 );
		if ( is_K0S ) 
			fillMTD( t2, "k0s_" );
	}
	void fillMTD( shared_ptr<TrackProxy> t, string prefix = "" ){
		if ( nullptr == t->_mtdPid )
			return;
		
		string chg = "pos";
		if ( t->_track->charge() < 0 )
			chg = "neg";
			
		book->fill( prefix + chg + "Cell", t->_track->gPt(), t->_mtdPid->cell() );
		book->fill( prefix + chg + "DeltaY", t->_track->gPt(), t->_mtdPid->deltaY() );
		book->fill( prefix + chg + "DeltaZ", t->_track->gPt(), t->_mtdPid->deltaZ() );
		book->fill( prefix + chg + "DeltaTOF", t->_track->gPt(), t->_mtdPid->deltaTimeOfFlight() );	

	}

	void fillPiPi( shared_ptr<TrackProxy> t1, shared_ptr<TrackProxy> t2, TLorentzVector &lv, string prefix = "" ){
		
		if ( tpcOnlyPairs )
			book->fill( prefix + "tpc_pt_mass", lv.M(), lv.Pt() );
		
		if ( nullptr != t1->_mtdPid && nullptr == t2->_mtdPid ){
			book->fill( prefix + "mtd_pt_mass", lv.M(), lv.Pt() );
		}
		if ( nullptr != t2->_mtdPid && nullptr == t1->_mtdPid ){
			book->fill( prefix + "mtd_pt_mass", lv.M(), lv.Pt() );
		}
		if ( nullptr != t1->_mtdPid &&  nullptr != t2->_mtdPid){
			book->fill( prefix + "mtd2_pt_mass", lv.M(), lv.Pt() );
		}
	}

	void fillMuMu( shared_ptr<TrackProxy> t1, shared_ptr<TrackProxy> t2, string prefix = "" ){
		if ( nullptr == t1->_mtdPid || nullptr == t2->_mtdPid )
			return;
		

		TLorentzVector lv1, lv2, lv;
		lv1.SetPtEtaPhiM( t1->_track->pMom().perp(), t1->_track->pMom().pseudoRapidity(), t1->_track->pMom().phi(), MASS_MU );
		lv2.SetPtEtaPhiM( t2->_track->pMom().perp(), t2->_track->pMom().pseudoRapidity(), t2->_track->pMom().phi(), MASS_MU );

		lv = lv1 + lv2;

		book->fill( "mumu_pt_mass", lv.M(), lv.Pt() );

		if ( fabs( lv.M() - 3.096 ) < 0.1 ){
			fillMTD( t1, "jpsi_" );
			fillMTD( t2, "jpsi_" );
		}
	}

	TLorentzVector analyzePair( StPicoEvent * _event, shared_ptr<TrackProxy> pos, shared_ptr<TrackProxy> neg ){
		float bField = _event->bField() * kilogauss;

		TLorentzVector lv;
		lv.SetPtEtaPhiM( 0, 0, 0, 0 );

		pos->_h = pos->_track->helix();
		neg->_h = neg->_track->helix();

		pair<double,double> pathLengths = pos->_h.pathLengths(neg->_h);

		StThreeVectorD pNegPosAtDca = pos->_h.at(pathLengths.first);
		StThreeVectorD pPosPosAtDca = neg->_h.at(pathLengths.second);

		StThreeVectorD dcaVector = pNegPosAtDca - pPosPosAtDca;
		StThreeVectorD primVtxPos = _event->primaryVertex();
		StThreeVectorD secVtxPos = pPosPosAtDca + ( dcaVector * 0.5 );
		StThreeVectorD decLenVec = secVtxPos - primVtxPos;

		if ( false == K0SCuts[ "DecayLength" ].inInclusiveRange( decLenVec.mag() ) ) return lv;
		if ( false == K0SCuts[ "MutualDCA" ].inInclusiveRange( dcaVector.mag() ) ) return lv;
		// if (decLenVec.mag() < 2.7) return lv; // cut on path length
		// if (dcaVector.mag() > 1.5) return lv; // cut on track mutual dca
		
		StThreeVectorD pNegMomAtDca = pos->_h.momentumAt(pathLengths.first, bField );	// CHECK THIS
		StThreeVectorD pPosMomAtDca = neg->_h.momentumAt(pathLengths.second, bField );

		StThreeVectorD K0sMomAtDCA = pNegMomAtDca + pPosMomAtDca;

		pos->_lv.SetPtEtaPhiM( pNegMomAtDca.perp(), pNegMomAtDca.pseudoRapidity(), pNegMomAtDca.phi(), MASS_PI );
		neg->_lv.SetPtEtaPhiM( pPosMomAtDca.perp(), pPosMomAtDca.pseudoRapidity(), pPosMomAtDca.phi(), MASS_PI );

		double K0pT = lv.Pt();
		double pointingAngle = cos(K0sMomAtDCA.angle(decLenVec) );
		if ( false == K0SCuts[ "PointingAngle" ].inInclusiveRange( pointingAngle ) ) return lv;
		
		// if (TMath::Abs(pointingAngle) > 0.106+0.056-0.1123*K0pT+0.025*K0pT*K0pT) return lv;

		double openingAngle = pNegMomAtDca.angle(pPosMomAtDca);
		
		lv = pos->_lv + neg->_lv;
		return lv;
	}

	void analyzePairs( StPicoEvent * _event ){

		float bField = _event->bField() * kilogauss;

		// ####### K0S Pion candidates #########
		for ( auto pos : pip ){
			for ( auto neg : pim ){
				if ( false == tpcOnlyPairs && nullptr == pos->_mtdPid && nullptr == neg->_mtdPid ) continue;

				TLorentzVector lv = analyzePair( _event, pos, neg );
				if ( lv.M() <= 0 ) continue;

				fillPiPi( pos, neg, lv, "" );
				
				fillMTD( pos, neg, lv );
			} // loop on pim
		} // loop on pip

		// ####### Muon candidates #########
		for ( auto pos : mup ){
			for ( auto neg : mum ){
				fillMuMu( pos, neg );
			}
		}


		// ##################### Like-sign #########################
		if ( false == makeLikeSign )
			return;

		for ( size_t i = 0; i < pip.size()-1; i++ ){
			auto p1 = pip[i];
			for ( size_t j = i+1; j < pip.size(); j++ ){
				auto p2 = pip[j];

				TLorentzVector lv = analyzePair( _event, p1, p2 );
				if ( lv.M() <= 0 ) continue;
				
				book->fill( "lsp_tpc_pt_mass", lv.M(), lv.Pt() );
				book->fill( "ls_tpc_pt_mass", lv.M(), lv.Pt() );

				if ( nullptr != p1->_mtdPid && nullptr == p2->_mtdPid ){
					book->fill( "ls_mtd_pt_mass", lv.M(), lv.Pt() );
				}
				if ( nullptr != p2->_mtdPid && nullptr == p1->_mtdPid ){
					book->fill( "ls_mtd_pt_mass", lv.M(), lv.Pt() );
				}
				if ( nullptr != p1->_mtdPid &&  nullptr != p2->_mtdPid){
					book->fill( "ls_mtd_pt_mass", lv.M(), lv.Pt() );
					book->fill( "ls_mtd2_pt_mass", lv.M(), lv.Pt() );
				}
			}
		}

		for ( size_t i = 0; i < pim.size()-1; i++ ){
			auto p1 = pim[i];
			for ( size_t j = i+1; j < pim.size(); j++ ){
				auto p2 = pim[j];

				TLorentzVector lv = analyzePair( _event, p1, p2 );
				if ( lv.M() <= 0 ) continue;
				
				book->fill( "lsn_tpc_pt_mass", lv.M(), lv.Pt() );
				book->fill( "ls_tpc_pt_mass", lv.M(), lv.Pt() );

				if ( nullptr != p1->_mtdPid && nullptr == p2->_mtdPid ){
					book->fill( "ls_mtd_pt_mass", lv.M(), lv.Pt() );
				}
				if ( nullptr != p2->_mtdPid && nullptr == p1->_mtdPid ){
					book->fill( "ls_mtd_pt_mass", lv.M(), lv.Pt() );
				}
				if ( nullptr != p1->_mtdPid &&  nullptr != p2->_mtdPid){
					book->fill( "ls_mtd_pt_mass", lv.M(), lv.Pt() );
					book->fill( "ls_mtd2_pt_mass", lv.M(), lv.Pt() );
				}
			}
		}
		// ##################### Like-sign #########################

	}


	virtual void analyzeEvent() {
		StPicoEvent *event = _rEvent.get( 0 );

		if ( nullptr == event ){
			return;
		}

		double bField = event->bField() * kilogauss;

		// if ( event->refMult() > 200 ) return;
		// if ( event->refMult() < 100 ) return;

		auto pVtx = event->primaryVertex();

		size_t nmtd = 0;
		size_t nmup = 0, nmum = 0;
		pip.clear();
		pim.clear();
		mup.clear();
		mum.clear();
		size_t ntrk = _rTrack.N();
		for ( size_t iTrack = 0; iTrack < ntrk; iTrack++ ){
			StPicoTrack *track = _rTrack.get( iTrack );
			StPicoMtdPidTraits *mtdPid = nullptr;
			if ( track->mtdPidTraitsIndex() >= 0 ){
				mtdPid = _rMtdPid.get( track->mtdPidTraitsIndex() );
			}

			float nHitsFit = track->nHitsFit() * track->charge();
			float nHitsRatio = nHitsFit / (float)track->nHitsMax() ;
			// book->fill( "nHitsFit", nHitsFit );
			// book->fill( "nHitsRatio", nHitsRatio );

			// if ( !track->isHFTTrack() ) continue; 
			if ( abs(nHitsFit ) < 10 ) continue;
			if ( abs( nHitsRatio ) < 0.52 ) continue;
			float gP = track->helix().momentum( bField ).magnitude();
			// book->fill( "gMom", gP );
			if ( gP <= 0 ) continue;
			float gDCA = track->globalDCA( bField, pVtx );
			// book->fill( "gDCA", gDCA );
			if ( abs( track->nSigmaPion() ) > 3.0 ) continue;
			


			auto tp = shared_ptr<TrackProxy>( new TrackProxy() );
			tp->_track = track;
			tp->_mtdPid = mtdPid;

			tp->_h = tp->_track->helix();
			tp->_p = tp->_h.momentum( bField );
			tp->_lv.SetPtEtaPhiM( tp->_p.perp(), tp->_p.pseudoRapidity(), tp->_p.phi(), MASS_PI );
		
			if ( gDCA < 1.5 && nullptr != mtdPid ){
				if ( track->charge() > 0 )
					mup.push_back( tp );
				else 
					mum.push_back( tp );
			}

			if ( false == K0SCuts[ "DaughterDCA" ].inInclusiveRange( gDCA ) ) continue;

			if ( track->charge() > 0 )
				pip.push_back( tp );
			else 
				pim.push_back( tp );


			if ( nullptr != mtdPid )
				nmtd++;
			if ( nullptr != mtdPid && track->charge() > 0 )
				nmup++;
			if ( nullptr != mtdPid && track->charge() < 0 )
				nmum++;
			
		}

		if ( pip.size() < 1 || pim.size() < 1 )
			return;
		if ( nmtd <= 0 )
			return;

		analyzePairs( event );

		book->fill( "n_pip", pip.size() );
		book->fill( "n_pim", pim.size() );

		book->fill( "n_mup", nmup );
		book->fill( "n_mum", nmum );
		book->fill( "n_mup_mum", nmup, nmum );
		


	}
	
};


#endif