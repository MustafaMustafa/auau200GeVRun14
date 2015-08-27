#include <iostream>
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TParticle.h"

#include "StFastSimConstants.h"
#include "StFastSimUtilities.h"

using namespace FastSimUtilitiesConstants;
using namespace std;

StFastSimUtilities::StFastSimUtilities()
{
   gRandom->SetSeed();

   // Note that this destructor can fail if any of the input files do not exist.
   // RAII

   TFile f(momResolutionFileName.c_str());
   mf1PionMomResolution = (TF1*)f.Get("fPion")->Clone("fPionMomResolution");
   mf1KaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaonMomResolution");
   f.Close();

   TFile fVertex(vertexFileName.c_str());
   for (int ii = 0; ii < nCent; ++ii)
   {
      mh1Vz[ii]      = (TH1D*)(fVertex.Get(Form("mh1Vz_%i", ii))->Clone(Form("mh1Vz_%i", ii)));
   }
   fVertex.Close();

   TFile fHftRatio(hftRatioFileName.c_str());
   TFile fDca(dcaFileName.c_str());

   for (int iParticle = 0; iParticle < nParticles; ++iParticle)
   {
      for (int iEta = 0; iEta < nEtas; ++iEta)
      {
         for (int iVz = 0; iVz < nVzs; ++iVz)
         {
            for (int ii = 0; ii < nCent; ++ii)
            {
               mh1HftRatio[iParticle][iEta][iVz][ii] =
                  (TH1D*)(fHftRatio.Get(Form("mh1HFT1PtCentPartEtaVzRatio_%i_%i_%i_%i", iParticle, iEta, iVz, ii))->Clone(Form("mh1HFT1PtCentPartEtaVzRatio_%i_%i_%i_%i", iParticle, iEta, iVz, ii)));
               for (int jj = 0; jj < nPtBins; ++jj)
               {
                  mh1DcaXY[iParticle][iEta][iVz][ii][jj] =
                     (TH1D*)((fDca.Get(Form("mh1DcaXyPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)))->Clone(Form("mh1DcaXyPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)));
                  mh1DcaZ[iParticle][iEta][iVz][ii][jj] =
                     (TH1D*)((fDca.Get(Form("mh1DcaZPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)))->Clone(Form("mh1DcaZPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)));
               }
            }
         }
      }
   }

   fHftRatio.Close();
   fDca.Close();
}

StFastSimUtilities::~StFastSimUtilities()
{
   // needs to be filled
}

TLorentzVector StFastSimUtilities::smearPos(TParticle const* const mcParticle,TVector3 const& vertex,int centrality) const
{
  int iParticleIndex = 0;

  switch(abs(mcParticle->GetPdgCode()))
   {
     case 211:
       iParticleIndex = 0;
       break;
     case 321:
       iParticleIndex = 1;
       break;
     default:
       cout << "There are no specific DCA resolution distributions available for PDG code = " << mcParticle->GetPdgCode() << "\n";
       cout << "Using Pions DCA resolution" << endl;
       iParticleIndex = 0;
   }

  TLorentzVector mom;
  TVector3 pos;
  mcParticle->Momentum(mom);
  pos.SetXYZ(mcParticle->Vx(), mcParticle->Vy(), mcParticle->Vz());

  return TLorentzVector(smearPos(iParticleIndex,vertex.z(),centrality,mom,pos),mcParticle->T());
}

TVector3 StFastSimUtilities::smearPos(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& rMom, TVector3 const& pos) const
{
   int const iEtaIndex = getEtaIndex(rMom.PseudoRapidity());
   int const iVzIndex = getVzIndex(vz);
   int const iPtIndex = getPtIndex(rMom.Perp());

   float sigmaPosZ = 0;
   float sigmaPosXY = 0;

   if (mh1DcaZ[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetEntries())
      sigmaPosZ = mh1DcaZ[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom();

   if (mh1DcaXY[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetEntries())
      sigmaPosXY = mh1DcaXY[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom();

   TVector3 newPos(pos);
   newPos.SetZ(0);
   TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
   newPos += momPerp.Unit() * sigmaPosXY;

   return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TVector3 StFastSimUtilities::getVertex(int const centrality) const
{
   double rdmVz;
   if (mh1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
   else rdmVz = mh1Vz[centrality]->GetRandom();
   return TVector3(0., 0., rdmVz);
}

bool StFastSimUtilities::matchHft(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& mom) const
{
   int const iEtaIndex = getEtaIndex(mom.PseudoRapidity());
   int const iVzIndex = getVzIndex(vz);
   // int const iPhiIndex = getPhiIndex(mom.Phi());
   int const bin = mh1HftRatio[iParticleIndex][iEtaIndex][iVzIndex][cent]->FindBin(mom.Perp());
   return gRandom->Rndm() < mh1HftRatio[iParticleIndex][iEtaIndex][iVzIndex][cent]->GetBinContent(bin);
}

float StFastSimUtilities::dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex) const
{
   TVector3 posDiff = pos - vertex;
   return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float StFastSimUtilities::dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0) const
{
   TVector3 posDiff = pos2 - pos1;
   TVector3 pu1 = p1.Unit();
   TVector3 pu2 = p2.Unit();
   double pu1Pu2 = pu1.Dot(pu2);
   double g = posDiff.Dot(pu1);
   double k = posDiff.Dot(pu2);
   double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
   double s1 = g + s2 * pu1Pu2;
   TVector3 posDca1 = pos1 + pu1 * s1;
   TVector3 posDca2 = pos2 + pu2 * s2;
   v0 = 0.5 * (posDca1 + posDca2);
   return (posDca1 - posDca2).Mag();
}

TParticle StFastSimUtilities::smear(TParticle const* mcParticle,TVector3 const& vertex,int const centrality) const
{
  return TParticle(mcParticle->GetPdgCode(),mcParticle->GetStatusCode(),
                   mcParticle->GetMother(0),mcParticle->GetMother(1),
                   mcParticle->GetDaughter(0),mcParticle->GetDaughter(1),
                   smearMom(mcParticle),smearPos(mcParticle,vertex,centrality));
}

TLorentzVector StFastSimUtilities::smearMom(TParticle const* const particle) const
{
   TF1 const* fMomResolution = NULL;

   switch(abs(particle->GetPdgCode()))
   {
     case 211:
       fMomResolution = mf1PionMomResolution;
       break;
     case 321:
       fMomResolution = mf1KaonMomResolution;
       break;
     default:
       cout << "There is no specific momentum resolution parametrization available for PDG code = " << particle->GetPdgCode() << "\n";
       cout << "Using Pions momentum resolution" << endl;
       fMomResolution = mf1PionMomResolution;
   }

   float const pt = particle->Pt();
   float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));

   TLorentzVector sMom;
   sMom.SetPtEtaPhiM(sPt,particle->Eta(),particle->Phi(),particle->GetMass());
   return sMom;
}

int StFastSimUtilities::getPtIndex(double const pT) const
{
   for (int i = 0; i < nPtBins; i++)
   {
      if ((pT >= ptEdge[i]) && (pT < ptEdge[i + 1]))
         return i;
   }
   return nPtBins - 1 ;
}

int StFastSimUtilities::getEtaIndex(double const Eta) const
{
   for (int i = 0; i < nEtas; i++)
   {
      if ((Eta >= EtaEdge[i]) && (Eta < EtaEdge[i + 1]))
         return i;
   }
   return nEtas - 1 ;
}

int StFastSimUtilities::getVzIndex(double const Vz) const
{
   for (int i = 0; i < nVzs; i++)
   {
      if ((Vz >= VzEdge[i]) && (Vz < VzEdge[i + 1]))
         return i;
   }
   return nVzs - 1 ;
}

int StFastSimUtilities::getPhiIndex(double const Phi) const
{
   for (int i = 0; i < nPtBins; i++)
   {
      if ((Phi >= PhiEdge[i]) && (Phi < PhiEdge[i + 1]))
         return i;
   }
   return nPhis - 1 ;
}
