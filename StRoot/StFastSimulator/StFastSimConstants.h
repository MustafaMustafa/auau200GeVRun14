#ifndef CUTS_H
#define CUTS_H

#include <string>

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

namespace FastSimUtilitiesConstants
{
   // input files
   std::string const momResolutionFileName = "momentum_resolution.root";
   std::string const hftRatioFileName      = "HFT_Ratio_VsPt_Centrality_Eta_Phi_Vz_Zdcx.root";
   std::string const dcaFileName           = "Dca_VsPt_Centrality_Eta_Phi_Vz_Zdcx.root";
   std::string const vertexFileName        = "Run14_After107_Vz_Cent.root";

   enum particleName {Pion, Kaon};

   // array constans and binning
   int const nParticles = 2;
   int const nEtas = 10;
   int const nVzs = 6;
   int const nPhis = 30;
   int const nCent = 9;
   int const nPtBins = 35;
   double const EtaEdge[nEtas + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
   double const VzEdge[nVzs + 1]   = { -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0};
   double const PhiEdge[nPhis + 1] = { -3.14159, -2.93215, -2.72271, -2.51327, -2.30383, -2.0944,
                                       -1.88496, -1.67552, -1.46608, -1.25664, -1.0472, -0.837758,
                                       -0.628319, -0.418879, -0.20944, 0.0, 0.20944, 0.418879,
                                       0.628319, 0.837758, 1.0472, 1.25664, 1.46608, 1.67552, 1.88496,
                                       2.0944, 2.30383, 2.51327, 2.72271, 2.93215, 3.14159
                                     };

   double const ptEdge[nPtBins + 1]= { 0.0, 0.2, 0.4,  0.6,  0.8, 1.0, 1.2, 1.4,  1.6,  1.8, 2.0, 2.2,
                                        2.4,  2.6,  2.8, 3.0, 3.2, 3.4,  3.6,  3.8, 4.0, 4.2, 4.4,  4.6,
                                        4.8, 5.0, 5.4, 5.8,  6.2,  6.6, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
}
#endif
