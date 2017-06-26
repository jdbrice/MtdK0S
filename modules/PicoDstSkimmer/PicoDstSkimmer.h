#ifndef PICO_DST_SKIMMER_H
#define PICO_DST_SKIMMER_H

#include "TreeAnalyzer.h"

#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"


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

	struct TrackProxy {
		StPicoTrack * _track = nullptr;
		StPicoMtdPidTraits *_mtdPid = nullptr;
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

	}
protected:

	TClonesArrayReader < StPicoEvent        > _rEvent;
	TClonesArrayReader < StPicoMtdHit       > _rMtdHit;
	TClonesArrayReader < StPicoTrack        > _rTrack;
	TClonesArrayReader < StPicoMtdPidTraits > _rMtdPid;
	TClonesArrayReader < StPicoBTofPidTraits > _rBTofPid;

	vector< shared_ptr<TrackProxy> > pip;
	vector< shared_ptr<TrackProxy> > pim;


	void analyzePairs( StPicoEvent * _event ){

		float bField = _event->bField() * kilogauss;

		TLorentzVector lvp, lvn, lv;
		for ( auto pos : pip ){

			StPhysicalHelixD h1 = pos->_track->helix();
			StThreeVectorD p1 = h1.momentum( bField );
			lvp.SetPtEtaPhiM( p1.perp(), p1.pseudoRapidity(), p1.phi(), MASS_PI );
			

			for ( auto neg : pim ){

				StPhysicalHelixD h2 = neg->_track->helix();
				StThreeVectorD p2 = h2.momentum( bField );
				lvn.SetPtEtaPhiM( p2.perp(), p2.pseudoRapidity(), p2.phi(), MASS_PI );

				lv = lvp + lvn;
				// book->fill( "pt_mass", lv.M(), lv.perp() );


				// continue;
				pair<double,double> pathLengths = h1.pathLengths(h2);

				StThreeVectorD pNegPosAtDca = h1.at(pathLengths.first);
				StThreeVectorD pPosPosAtDca = h2.at(pathLengths.second);
				
				StThreeVectorD dcaVector = pNegPosAtDca-pPosPosAtDca;
				StThreeVectorD primVtxPos = _event->primaryVertex();
				StThreeVectorD secVtxPos = pPosPosAtDca + ( dcaVector * 0.5 );
				StThreeVectorD decLenVec = secVtxPos - primVtxPos;
				

				if (decLenVec.mag() < 2.7) continue; // cut on path length
				if (dcaVector.mag() > 1.5) continue; // cut on track mutual dca

				StThreeVectorD pNegMomAtDca = h1.momentumAt(pathLengths.first, bField );
				StThreeVectorD pPosMomAtDca = h2.momentumAt(pathLengths.second, bField );
				StThreeVectorD K0sMomAtDCA = pNegMomAtDca + pPosMomAtDca;
				lvp.SetPtEtaPhiM( pNegMomAtDca.perp(), pNegMomAtDca.pseudoRapidity(), pNegMomAtDca.phi(), MASS_PI );
				lvn.SetPtEtaPhiM( pPosMomAtDca.perp(), pPosMomAtDca.pseudoRapidity(), pPosMomAtDca.phi(), MASS_PI );
				

				double K0pT = lv.Pt();
				double pointingAngle = K0sMomAtDCA.angle(decLenVec);
				if (TMath::Abs(pointingAngle) > 0.106+0.056-0.1123*K0pT+0.025*K0pT*K0pT) continue;


				// book->fill( "pointingAngle_pt", K0pT, pointingAngle );
				double openingAngle = pNegMomAtDca.angle(pPosMomAtDca);
				// book->fill( "openingAngle", openingAngle );
				book->fill( "pt_mass", lv.M(), K0pT );

				if ( abs(lv.M() - 0.497) < 0.2 ){
					if ( nullptr != pos->_mtdPid ){
						book->fill( "mtd_pt_mass", lv.M(), K0pT );
						book->fill( "posDeltaY", pos->_mtdPid->deltaY() );
						book->fill( "posDeltaZ", pos->_mtdPid->deltaZ() );
						book->fill( "posDeltaTOF", pos->_mtdPid->deltaTimeOfFlight() );
						book->fill( "posCell", pos->_mtdPid->cell() );
					}

					if ( nullptr != neg->_mtdPid ){
						book->fill( "mtd_pt_mass", lv.M(), K0pT );
						book->fill( "negDeltaZ", neg->_mtdPid->deltaZ() );
						book->fill( "negDeltaY", neg->_mtdPid->deltaY() );
						book->fill( "negDeltaTOF", neg->_mtdPid->deltaTimeOfFlight() );
						book->fill( "negCell", neg->_mtdPid->cell() );
					}

					
					

					
					

					// book->fill( "pointingAngle_pt", K0pT, pointingAngle );
					// book->fill( "openingAngle", openingAngle );
					// book->fill( "secVtx", secVtxPos.mag() );
					// book->fill( "dcaLenZ", decLenVec.z() ); 
				}


			}
		}

	}


	virtual void analyzeEvent() {
		StPicoEvent *event = _rEvent.get( 0 );

		if ( nullptr == event ){
			return;
		}

		// if ( event->refMult() > 100 ) return;

		auto pVtx = event->primaryVertex();

		size_t nmtd = 0;
		pip.clear();
		pim.clear();
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
			float gP = track->helix().momentum( event->bField() * kilogauss ).magnitude();
			// book->fill( "gMom", gP );
			if ( gP <= 0 ) continue;
			float gDCA = track->globalDCA( event->bField() * kilogauss, pVtx );
			// book->fill( "gDCA", gDCA );
			if ( gDCA < 1.3 || gDCA > 9 ) continue;
			if ( abs( track->nSigmaPion() ) > 3.0 ) continue;
			
			auto tp = shared_ptr<TrackProxy>( new TrackProxy() );
			tp->_track = track;
			tp->_mtdPid = mtdPid;

			if ( track->charge() > 0 )
				pip.push_back( tp );
			else 
				pim.push_back( tp );
			if ( nullptr != mtdPid )
				nmtd++;
			
		}

		if ( pip.size() < 1 || pim.size() < 1 )
			return;
		if ( nmtd <= 0 )
			return;

		analyzePairs( event );

		// LOG_IF_F( INFO, DEBUG, "#pi+=%lu, #pi-=%lu, #mtd=%lu", pip.size(), pim.size(), nmtd );
		book->fill( "n_pip", pip.size() );
		book->fill( "n_pim", pim.size() );

		// LOG_IF_F( INFO, DEBUG, "RunId: %d", event->runId() );
		// LOG_IF_F( INFO, DEBUG, "#Tracks: %u", _rTrack.N() );
		// LOG_IF_F( INFO, DEBUG, "#MtdHits: %u", _rMtdHit.N() );
		// LOG_IF_F( INFO, DEBUG, "#MtdPids: %u", _rMtdPid.N() );

		// size_t n = event->triggerIds().size();

		// LOG_IF_F( INFO, DEBUG, "# of trigger Ids = %lu", n );
		// for ( unsigned int id : event->triggerIds() ){
		// 	LOG_IF_F( INFO, DEBUG, "id: %lu", id );
		// } 


	}
	
};


#endif