#include "QAG4SimulationD0.h"
#include "QAHistManagerDef.h"

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxVertexEval.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h> // for cluskey

#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenRanges.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/getClass.h>

#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

QAG4SimulationD0::QAG4SimulationD0(const std::string &name)
  : SubsysReco(name)
{
  cout << "QAG4SimulationD0::QAG4SimulationD0(const std::string &name) Calling ctor" << endl;
}

int QAG4SimulationD0::InitRun(PHCompositeNode *topNode)
{
  if (!m_svtxEvalStack)
    {
      m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
      m_svtxEvalStack->set_strict(false);
      m_svtxEvalStack->set_verbosity(Verbosity() + 1);
    }
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  m_Geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_Geneventmap)
    {
      std::cout << PHWHERE << " - Fatal error - missing node PHHepMCGenEventMap" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationD0::Init(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  // Invariant mass assuming K+,Pi-
  h = new TH1F(TString(get_histo_prefix()) + "InvM_KpPm",
               "Invariant Mass Distribution (K+ #pi -);Invariant Mass (GeV/c);Counts", 50, 1.7, 2.);
  hm->registerHisto(h);
  // Invariant mass assuming K-,Pi+
  h = new TH1F(TString(get_histo_prefix()) + "InvM_KmPp",
               "Invariant Mass Distribution (K- #pi +);Invariant Mass (GeV/c);Counts", 50, 1.7, 2.);
  hm->registerHisto(h);

  // Invariant mass from HepMC event true
  h = new TH1F(TString(get_histo_prefix()) + "HepMC_InvM_true",
               "Invariant Mass Distribution from HepMC Events (true);Invariant Mass (GeV/c);Counts", 50, 1.7, 2.);
  hm->registerHisto(h);
  // Invariant mass from HepMC event false
  h = new TH1F(TString(get_histo_prefix()) + "HepMC_InvM_false",
               "Invariant Mass Distribution from HepMC Events (false);Invariant Mass (GeV/c);Counts", 50, 1.7, 2.);
  hm->registerHisto(h);
  

  return Fun4AllReturnCodes::EVENT_OK;
}

void QAG4SimulationD0::addEmbeddingID(int embeddingID)
{
  m_embeddingIDs.insert(embeddingID);
}

int QAG4SimulationD0::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    cout << "QAG4SimulationD0::process_event() entered" << endl;

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  
  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);

  SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();
  assert(trackeval);
  /*SvtxTruthEval *trutheval = m_svtxEvalStack->get_truth_eval();
  assert(trutheval);
  */
  // SvtxVertexEval *vertexeval = m_svtxEvalStack->get_vertex_eval();
  // SvtxClusterEval *clustereval = m_svtxEvalStack->get_cluster_eval();

  // K+Pi- Inv mass hist
  TH1 *h_InvM_KpPm = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "InvM_KpPm"));
  assert(h_InvM_KpPm);
  // K-Pi+ Inv mass hist
  TH1 *h_InvM_KmPp = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "InvM_KmPp"));
  assert(h_InvM_KmPp);

  // K+Pi- Inv mass hist from HepMC
  TH1 *h_HepMC_InvM_true = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "HepMC_InvM_true"));
  assert(h_HepMC_InvM_true);
  // K-Pi+ Inv mass hist from HepMC
  TH1 *h_HepMC_InvM_false = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "HepMC_InvM_false"));
  assert(h_HepMC_InvM_false);

  SvtxVertexMap *vertexmap = nullptr;


  assert(m_Geneventmap);

  PHHepMCGenEvent* genevt = m_Geneventmap->get(_embedding_id);
  if (!genevt)
    {
      std::cout << PHWHERE << " - Fatal error - node PHHepMCGenEventMap missing subevent with embedding ID " << _embedding_id;
      std::cout << ". Print PHHepMCGenEventMap:";
      m_Geneventmap->identify();
      return Fun4AllReturnCodes::ABORTRUN;
    }

  HepMC::GenEvent* theEvent = genevt->getEvent();
  assert(theEvent);
  if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "HFMLTriggerHepMCTrigger::process_event - process HepMC::GenEvent with signal_process_id = "
	   << theEvent->signal_process_id();
      if (theEvent->signal_process_vertex())
	{
	  cout << " and signal_process_vertex : ";
	  theEvent->signal_process_vertex()->print();
	}
      cout << "  - Event record:" << endl;
      theEvent->print();
    }

  TDatabasePDG* pdg = TDatabasePDG::Instance();

  int targetPID = std::abs(pdg->GetParticle("D0")->PdgCode());
  int daughter1PID = std::abs(pdg->GetParticle("pi+")->PdgCode());
  // int daughter1PID_2 = std::abs(pdg->GetParticle("pi-")->PdgCode());
  int daughter2PID = std::abs(pdg->GetParticle("K+")->PdgCode());
  // int daughter2PID_2 = std::abs(pdg->GetParticle("K-")->PdgCode());

  float K_mass = .493677; // GeV/c^2
  float Pi_mass = .139570; // GeV/c^2 

  vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (vertexmap && truthinfo)
    {
      vector<SvtxTrack*> positive_tracks;
      vector<SvtxTrack*> negative_tracks;
      double charge;
      for (SvtxVertexMap::Iter iter = vertexmap->begin();
	   iter != vertexmap->end();
	   ++iter)
	{
	  SvtxVertex *vertex = iter->second;

	  // PHG4VtxPoint *point = vertexeval->max_truth_point_by_ntracks(vertex);
	  // float ntracks = vertex->size_tracks();
	  // float gntracks = truthinfo->GetNumPrimaryVertexParticles();
	  for (SvtxVertex::TrackIter iter2 = vertex->begin_tracks();
	       iter2 != vertex->end_tracks();
	       ++iter2)
	    {
	      SvtxTrack* track = m_trackMap->get(*iter2);
	      if (false)
		{
		  assert(track);
		}
	      else if (!track)
		{
		  continue;
		} 
	      int MVTX_hits = 0;
	      int INTT_hits = 0;
	      int TPC_hits = 0;

	      for (auto cluster_iter = track->begin_cluster_keys(); cluster_iter != track->end_cluster_keys(); ++cluster_iter)
		{
		  const auto &cluster_key = *cluster_iter;
		  const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

		  if (trackerID == TrkrDefs::mvtxId)
		    ++MVTX_hits;
		  else if (trackerID == TrkrDefs::inttId)
		    ++INTT_hits;
		  else if (trackerID == TrkrDefs::tpcId)
		    ++TPC_hits;
		  else
		    {
		      if (Verbosity())
			cout << "QAG4SimulationD0::process_event - unkown tracker ID = " << trackerID << " from cluster " << cluster_key << endl;
		    }
		}
	      if (MVTX_hits >= 2 && INTT_hits >= 1 && TPC_hits >=30)
		{
		  charge = track->get_charge();
		  if (charge > 0)
		    {
		      positive_tracks.push_back(track);
		    }
		  else if (charge < 0)
		    {
		      negative_tracks.push_back(track);
		    }
		}
	    }
	  float p_energy_K;
	  float n_energy_K;
	  float p_energy_P;
	  float n_energy_P;
	  // TVector3 pp;
	  // TVector3 np;
	  float px;
	  float py;
	  float pz;
	  float inv_mass;
	  for(SvtxTrack* tr_p : positive_tracks)
	    {
	      px = tr_p->get_px();
	      py = tr_p->get_py();
	      pz = tr_p->get_pz();
	      TVector3 pp(px,py,pz);
	      p_energy_K = sqrt(pp.Mag2() + pow(K_mass,2)); // energy-momentum relation, assuming K+
	      p_energy_P = sqrt(pp.Mag2() + pow(Pi_mass,2)); // assuming Pi+
	           
	      for(SvtxTrack* tr_n : negative_tracks)
		{
		  px = tr_n->get_px();
		  py = tr_n->get_py();
		  pz = tr_n->get_pz();
		  TVector3 np(px,py,pz);
		  n_energy_K = sqrt(np.Mag2() + pow(K_mass,2)); // assuming K-
		  n_energy_P = sqrt(np.Mag2() + pow(Pi_mass,2)); // assuming Pi-

		  inv_mass = sqrt(pow((p_energy_K + n_energy_P),2) - ((pp+np).Mag2())); // K+Pi-
		  h_InvM_KpPm->Fill(inv_mass);
		  inv_mass = sqrt(pow((p_energy_P + n_energy_K),2) - ((pp+np).Mag2())); // K-Pi+
		  h_InvM_KmPp->Fill(inv_mass);
		}
	    }
	}
    }

  auto range = theEvent->particle_range();
  for (HepMC::GenEvent::particle_const_iterator piter = range.begin(); piter != range.end(); ++piter)
    {
      const HepMC::GenParticle* p = *piter;
      assert(p);

      if (std::abs(p->pdg_id()) == targetPID)
	{
	  if (Verbosity())
	    {
	      cout << "HFMLTriggerHepMCTrigger::process_event - Accept signal particle : ";
	      p->print();
	      cout << endl;
	    }

	  const HepMC::GenVertex* decayVertex = p->end_vertex();

	  int hasDecay1(0);
	  int hasDecay2(0);
	  int hasDecayOther(0);

	  if (decayVertex)
	    {
	      for (auto diter = decayVertex->particles_out_const_begin();
		   diter != decayVertex->particles_out_const_end();
		   ++diter)

		{
		  const HepMC::GenParticle* pd = *diter;
		  assert(pd);

		  if (Verbosity())
		    {
		      cout << "HFMLTriggerHepMCTrigger::process_event - Testing daughter particle: ";
		      pd->print();
		      cout << endl;
		    }

		  if (std::abs(pd->pdg_id()) == daughter1PID)
		    {
		      const double eta = pd->momentum().eta();

		      if (eta > -1. and eta < 1.)
			++hasDecay1;
		    }
		  else if (std::abs(pd->pdg_id()) == daughter2PID)
		    {
		      const double eta = pd->momentum().eta();

		      if (eta > -1. and eta < 1.)
			++hasDecay2;
		    }
		  else
		    ++hasDecayOther;
		}

	      if (hasDecay1 == 1 and hasDecay2 == 1 and hasDecayOther == 0)
		{
		  int hepmc_barcode1;
		  int hepmc_barcode2;
		  
		  int p1ID;
		  int p2ID;

		  for (auto diter = decayVertex->particles_out_const_begin();
		       diter != decayVertex->particles_out_const_end();
		       ++diter)
		    {
		      const HepMC::GenParticle* pd1 = *diter;
		      assert(pd1);
		      ++diter;
		      const HepMC::GenParticle* pd2 = *diter;
		      assert(pd2);
		            
		      hepmc_barcode1 = pd1->barcode();
		      hepmc_barcode2 = pd2->barcode();

		      p1ID = std::abs(pd1->pdg_id());
		      p2ID = std::abs(pd2->pdg_id());
 
		      float energy_K1 = 0;
		      float energy_K2 = 0;
		      float energy_P1 = 0;
		      float energy_P2 = 0;
		      float px1;
		      float py1;
		      float pz1;
		      float px2;
		      float py2;
		      float pz2;
		      float inv_mass;
		      TVector3 p1;
		      TVector3 p2;

		      const auto prange = truthinfo->GetPrimaryParticleRange();
		      for (auto iter = prange.first; iter != prange.second; ++iter)  // process all primary paricle
			{
			  PHG4Particle* g4particle = iter->second;
			  if (g4particle->get_barcode() == hepmc_barcode1)
			    {
			      SvtxTrack* track = trackeval->best_track_from(g4particle);
			      if (track)
			      {
				px1 = track->get_px();
				py1 = track->get_py();
				pz1 = track->get_pz();
				p1.SetXYZ(px1,py1,pz1);
				energy_K1 = sqrt(p1.Mag2() + pow(K_mass,2)); // energy-momentum relation, assuming K
				energy_P1 = sqrt(p1.Mag2() + pow(Pi_mass,2)); // assuming Pi         
			      }
			    }
			  else if (g4particle->get_barcode() == hepmc_barcode2)
			    {
			      SvtxTrack* track = trackeval->best_track_from(g4particle);
			      if (track)
			      {
				px2 = track->get_px();
				py2 = track->get_py();
				pz2 = track->get_pz();
				p2.SetXYZ(px2,py2,pz2);
				energy_K2 = sqrt(p2.Mag2() + pow(K_mass,2)); // energy-momentum relation, assuming 
				energy_P2 = sqrt(p2.Mag2() + pow(Pi_mass,2)); // assuming Pi 
			      }
			    }
			}
		      if (p1ID == daughter1PID && p2ID == daughter2PID)
			{
			  inv_mass = sqrt(pow((energy_P1 + energy_K2),2) - (p1+p2).Mag2());
			  h_HepMC_InvM_true->Fill(inv_mass);
			  inv_mass = sqrt(pow((energy_P2 + energy_K1),2) - (p1+p2).Mag2());
                          h_HepMC_InvM_false->Fill(inv_mass);
			}
		      else if (p1ID == daughter2PID && p2ID == daughter1PID)
			{
			  inv_mass = sqrt(pow((energy_P2 + energy_K1),2) - (p1+p2).Mag2());
			  h_HepMC_InvM_true->Fill(inv_mass);
			  inv_mass = sqrt(pow((energy_P1 + energy_K2),2) - (p1+p2).Mag2());
                          h_HepMC_InvM_false->Fill(inv_mass);
			}
		    }
		}

	    }  //      if (decayVertex)
	  else
	    {
	      cout << "HFMLTriggerHepMCTrigger::process_event - Warning - target particle did not have decay vertex : ";
	      p->print();
	      cout << endl;
	    }

	}  //    if (std::abs(p-> pdg_id()) == targetPID)
    }    //  for (HepMC::GenEvent::particle_const_iterator piter = range.begin(); piter != range.end(); ++piter)


  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationD0::load_nodes(PHCompositeNode *topNode)
{
  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthContainer)
    {
      cout << "QAG4SimulationD0::load_nodes - Fatal Error - "
	   << "unable to find DST node "
	   << "G4TruthInfo" << endl;
      assert(m_truthContainer);
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

string
QAG4SimulationD0::get_histo_prefix()
{
  return string("h_") + Name() + string("_");
}
