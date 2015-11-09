//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4VGammaDeexcitation.cc 87376 2014-12-02 08:25:05Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4VGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//
//      Creation date: 23 October 1998
//
//      Modifications:
//
//        21 Nov 2001, Fan Lei (flei@space.qinetiq.com)
//           Modified GenerateGamma() and UpdateUncleus() for implementation
//           of Internal Conversion processs
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//
//        19 April 2010, J. M. Quesada calculations in CM system
//              pending final boost to lab system  (not critical)
//
//        23 April 2010, V.Ivanchenko rewite kinematic part using PDG formula
//                                    for 2-body decay
//
//        07 May   2011, V.Ivanchenko implement check ICM flag - produce or not e-
//
// -------------------------------------------------------------------

#include "G4VGammaDeexcitation.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LorentzVector.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4DiscreteGammaTransition.hh"

#include "G4NuclearLevelStore.hh" // Evan Rand

G4VGammaDeexcitation::G4VGammaDeexcitation(): _transition(0), _verbose(0),
                          _electronO (0), _vSN(-1)
{
  _tolerance = 2*CLHEP::keV;
  _timeLimit = DBL_MAX;
}

G4VGammaDeexcitation::~G4VGammaDeexcitation()
{
  delete _transition;
}

void G4VGammaDeexcitation::DoChain(G4FragmentVector* products,
                   G4Fragment* nucleus)
{
  if (_verbose > 1) { G4cout << "G4VGammaDeexcitation::DoChain" << G4endl; }

  if(CanDoTransition(nucleus)) {
    for(size_t i=0; i<100; ++i) {
      _transition->SetEnergyFrom(nucleus->GetExcitationEnergy());
      G4Fragment* gamma = GenerateGamma(nucleus);
      if (gamma) { products->push_back(gamma); }
      else { break; }
      //G4cout << i << ".  Egamma(MeV)= " << gamma->GetMomentum().e()
      //	     << "; new Eex(MeV)= " << nucleus->GetExcitationEnergy()
      //       << G4endl;
      if(nucleus->GetExcitationEnergy() <= _tolerance) { break; }
    }
  }
  if (_verbose > 1) {
    G4cout << "G4VGammaDeexcitation::DoChain - end" << G4endl;
  }
}

G4Fragment* G4VGammaDeexcitation::GenerateGamma(G4Fragment* aNucleus)
{
  G4Fragment * thePhoton = 0;
  if(!CanDoTransition(aNucleus)) { return thePhoton; }

  _transition->SelectGamma();  // it can be conversion electron too
  G4double etrans = _transition->GetGammaEnergy();
  //G4cout << "G4VGammaDeexcitation::GenerateGamma - Etrans(MeV)= "
  //	 << etrans << G4endl;

  if(etrans <= 0.0) { return thePhoton; }

  // final excitation
  G4double excitation = aNucleus->GetExcitationEnergy() - etrans;
  if(excitation <= _tolerance) { excitation = 0.0; }
  if (_verbose > 1) {
    G4cout << "G4VGammaDeexcitation::GenerateGamma - Edeexc(MeV)= " << etrans
       << " ** left Eexc(MeV)= " << excitation
       << G4endl;
  }

  G4double gammaTime = _transition->GetGammaCreationTime();

  // Do complete Lorentz computation
  G4LorentzVector lv = aNucleus->GetMomentum();
  G4double Mass = aNucleus->GetGroundStateMass() + excitation;

  // select secondary
  G4ParticleDefinition* gamma = G4Gamma::Gamma();

  G4DiscreteGammaTransition* dtransition =
    dynamic_cast <G4DiscreteGammaTransition*> (_transition);

  if (dtransition && !( dtransition->IsAGamma()) ) {
    gamma = G4Electron::Electron();
    _vSN = dtransition->GetOrbitNumber();
    _electronO.RemoveElectron(_vSN);
    lv += G4LorentzVector(0.0,0.0,0.0,CLHEP::electron_mass_c2);
  }

//  G4double cosTheta = 1. - 2. * G4UniformRand();
//  G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
//  G4double phi = twopi * G4UniformRand();

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  G4ThreeVector gammaWThetaVector;
  G4ThreeVector polarizationVector;
  G4double theta;
  G4double phi;
  G4double cosTheta;
  G4double sinTheta;
  G4double higherLevelEnergy;

  G4double levelEnergy      = dtransition->GetLevelEnergy();
  G4bool firstGammaDecay    = G4NuclearLevelStore::GetInstance()->GetFirstGammaDecay();
  if(dtransition->IsAGamma() && G4NuclearLevelStore::GetInstance()->GetUserFilesMultipole(aNucleus->GetZ(), aNucleus->GetA())) { // if we have given the multipole file
      firstGammaDecay = G4NuclearLevelStore::GetInstance()->GetFirstGammaDecay();

      if(firstGammaDecay) { // is this the first gamma decay?
          // This level was beta fed (or alpha, etc.).
          // There is no higher level that fed this level, set higherLevelEnergy = -1.0
          higherLevelEnergy = -1.0;
          polarizationVector = G4ThreeVector(0,0,1);
      }
      else {
          // This level was fed by a gamma (or IC) from a higher level.
          higherLevelEnergy     = G4NuclearLevelStore::GetInstance()->GetHigherLevelEnergy();
          polarizationVector    = G4NuclearLevelStore::GetInstance()->GetPolarizationVector();
      }
      theta     = dtransition->GetThetaFromWTheta(higherLevelEnergy);
      sinTheta  = std::sin(theta);
      cosTheta  = std::cos(theta);
      phi       = twopi * G4UniformRand();

      gammaWThetaVector = G4ThreeVector(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);

      // Apply polarization vector
      gammaWThetaVector.rotateY(polarizationVector.getTheta());
      gammaWThetaVector.rotateZ(polarizationVector.getPhi());

      sinTheta  = std::sin(gammaWThetaVector.getTheta());
      cosTheta  = std::cos(gammaWThetaVector.getTheta());
      phi       = gammaWThetaVector.getPhi();

      // Set SetFirstGammaDecay to false for next gamma decay
      // This can also be set by the DoDecay method G4RadioactiveDecay.cc
      if(firstGammaDecay) G4NuclearLevelStore::GetInstance()->SetFirstGammaDecay(false);
      G4NuclearLevelStore::GetInstance()->SetHigherLevelEnergy(levelEnergy);
      G4NuclearLevelStore::GetInstance()->SetPolarizationVector(gammaWThetaVector);
  }
  else { // isotropic
      cosTheta = 1. - 2. * G4UniformRand();
      sinTheta = std::sqrt(1. - cosTheta * cosTheta);
      phi = twopi * G4UniformRand();

      // If an internal conversion occurs, this will ruin the correlation
      // More physics is needed if one wishes to include electron-gamma correlations
      gammaWThetaVector = G4ThreeVector(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);
      G4NuclearLevelStore::GetInstance()->SetPolarizationVector(gammaWThetaVector);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  G4double eMass = gamma->GetPDGMass();
  G4LorentzVector Gamma4P;
  /*
  G4cout << " Mass= " << eMass << " t= " << gammaTime
     << " tlim= " << _timeLimit << G4endl;
  */
  if(gammaTime > _timeLimit) {
    // shortcut for long lived levels
    // not correct position of stopping ion gamma emission
    // 4-momentum balance is breaked
    G4double eGamma = aNucleus->GetExcitationEnergy() - excitation;
    G4double e = eGamma + eMass;
    G4double mom = std::sqrt(eGamma*(eGamma + 2*eMass));
    Gamma4P.set(mom * sinTheta * std::cos(phi),
        mom * sinTheta * std::sin(phi),
        mom * cosTheta, e);
    lv -= Gamma4P;
    e = lv.e();
    if(e < Mass) { e = Mass; }
    mom = std::sqrt((e - Mass)*(e + Mass));
    G4ThreeVector v = lv.vect().unit();
    lv.set(mom*v.x(), mom*v.y(), mom*v.z(), e);

  } else {
    // 2-body decay in rest frame
    G4double Ecm       = lv.mag();
    G4ThreeVector bst  = lv.boostVector();

    G4double GammaEnergy = 0.5*((Ecm - Mass)*(Ecm + Mass) + eMass*eMass)/Ecm;
    if(GammaEnergy < eMass) { GammaEnergy = eMass; }

    G4double mom = std::sqrt((GammaEnergy - eMass)*(GammaEnergy + eMass));
    Gamma4P.set(mom * sinTheta * std::cos(phi),
        mom * sinTheta * std::sin(phi),
        mom * cosTheta,
        GammaEnergy);

    Gamma4P.boost(bst);
    lv -= Gamma4P;
  }

  // modified primary fragment
  gammaTime += aNucleus->GetCreationTime();

  aNucleus->SetMomentum(lv);
  aNucleus->SetCreationTime(gammaTime);

  // gamma or e- are produced
  thePhoton = new G4Fragment(Gamma4P,gamma);
  thePhoton->SetCreationTime(gammaTime);

  //G4cout << "G4VGammaDeexcitation::GenerateGamma : " << thePhoton << G4endl;
  //G4cout << "       Left nucleus: " << aNucleus << G4endl;
  return thePhoton;
}
