#ifndef StFastSimUtilities_hh
#define StFastSimUtilities_hh

/* **************************************************
 *
 *  Authors:  Xin Dong          (xdong@lbl.gov),
 *            Hao Qiu           (hqiu@lbl.gov),
 *            **Mustafa Mustafa (mmustafa@lbl.gov),
 *            **Guannan Xie     (guannanxie@lbl.gov)
 *
 *  ** Code Maintainers
 *
 * **************************************************
 */

#include "TVector3.h"
#include "TLorentzVector.h"
#include "StFastSimConstants.h"

class TF1;
class TH1D;
class TParticle;

using namespace FastSimUtilitiesConstants;

class StFastSimUtilities
{
public:
   StFastSimUtilities();
   ~StFastSimUtilities();

   TLorentzVector smearMom(TParticle const*) const;
   TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos) const;
   TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos) const;
   float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex) const;
   float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0) const;
   TVector3 getVertex(int centrality) const;
   bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom) const;

private:
   int getPtIndex(double) const;
   int getEtaIndex(double) const;
   int getVzIndex(double) const;
   int getPhiIndex(double) const;

   // vertex distributions
   TH1D* mh1Vz[nCent];

   // Momentum resolution
   TF1* mf1KaonMomResolution = nullptr;
   TF1* mf1PionMomResolution = nullptr;

   // HFT and DCA distributions
   TH1D* mh1HftRatio[nParticles][nEtas][nVzs][nCent];
   TH1D* mh1DcaZ[nParticles][nEtas][nVzs][nCent][nPtBins];
   TH1D* mh1DcaXY[nParticles][nEtas][nVzs][nCent][nPtBins];
};
#endif
