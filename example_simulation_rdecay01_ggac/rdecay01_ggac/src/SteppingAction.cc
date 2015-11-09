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
// $Id: SteppingAction.cc 68030 2013-03-13 13:51:27Z gcosmo $
//
/// \file radioactivedecay/rdecay02/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
#include "G4ios.hh"

#include "HistoManager.hh"
#include "SteppingAction.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction()
 : G4UserSteppingAction()
{ 
  foundgamma1 = false;
  foundgamma2 = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* aStep) 
{
  G4double dotprod, normgamma1, normgamma2, value;

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

  G4Track* aTrack = aStep->GetTrack();
  G4int StepNo = aTrack->GetCurrentStepNumber();

  // Get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4String volname = volume->GetName();

  G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();


  G4StepPoint* point2 	= aStep->GetPostStepPoint();
  G4ThreeVector pos2   	= point2->GetPosition();

  size_t found;
  found = volname.find("World");
  if (ekin > 1171*keV && ekin < 1175*keV && found!=G4String::npos) {
    gamma1 = pos2;
    foundgamma1 = true;
    foundgamma2 = false; // gamma 1173 is the 1st gammma
  }
  else if (ekin > 1330*keV && ekin < 1334*keV && found!=G4String::npos) {
    gamma2 = pos2;
    foundgamma2 = true;
  }

  if(foundgamma1 == true && foundgamma2 == true) {
    dotprod   	= gamma1.x()*gamma2.x() + gamma1.y()*gamma2.y() + gamma1.z()*gamma2.z();
    normgamma1	= sqrt( gamma1.x()*gamma1.x() + gamma1.y()*gamma1.y() + gamma1.z()*gamma1.z() );
    normgamma2  = sqrt( gamma2.x()*gamma2.x() + gamma2.y()*gamma2.y() + gamma2.z()*gamma2.z() );    
    value       = acos( dotprod/(normgamma1*normgamma2) );

    analysis->FillH1(0,cos(value));

    foundgamma1 = false;    
    foundgamma2 = false;    
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


