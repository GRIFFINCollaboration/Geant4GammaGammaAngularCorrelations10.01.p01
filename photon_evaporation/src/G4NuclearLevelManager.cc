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
// $Id: G4NuclearLevelManager.cc 87163 2014-11-26 08:46:54Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4NuclearLevelManager
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//
//      Creation date: 24 October 1998
//
//      Modifications:
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//        02 May 2003,   Vladimir Ivanchenko remove rublic copy constructor
//        06 Oct 2010, M. Kelsey -- Use object storage, not pointers, drop
//		public access to list, simplify list construction
// -------------------------------------------------------------------
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "G4NuclearLevelManager.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclearLevel.hh"
#include "G4ios.hh"
#include "G4HadronicException.hh"
#include "G4HadTmpUtil.hh"
#include "G4NucleiProperties.hh"
#include "G4PhysicalConstants.hh"

G4double G4NuclearLevelManager::_levelEnergy=0.;
G4double G4NuclearLevelManager::_gammaEnergy=0.;
G4double G4NuclearLevelManager::_probability=0.;
G4double G4NuclearLevelManager::_polarity=0.;
G4double G4NuclearLevelManager::_halfLife=0.;
G4double G4NuclearLevelManager::_angularMomentum=0.;
G4double G4NuclearLevelManager::_kCC=0.;
G4double G4NuclearLevelManager::_l1CC=0.;
G4double G4NuclearLevelManager::_l2CC=0.;
G4double G4NuclearLevelManager::_l3CC=0.;
G4double G4NuclearLevelManager::_m1CC=0.;
G4double G4NuclearLevelManager::_m2CC=0.;
G4double G4NuclearLevelManager::_m3CC=0.;
G4double G4NuclearLevelManager::_m4CC=0.;
G4double G4NuclearLevelManager::_m5CC=0.;
G4double G4NuclearLevelManager::_nPlusCC=0.;
G4double G4NuclearLevelManager::_totalCC=0.;

G4NuclearLevelManager::G4NuclearLevelManager(G4int Z, G4int A,
                         const G4String& filename) :
    _nucleusA(A), _nucleusZ(Z), _validity(false),
    _levels(0)
{
  if (A <= 0 || Z <= 0 || Z > A ) {
    throw G4HadronicException(__FILE__, __LINE__,
                  "==== G4NuclearLevelManager ==== (Z,A) <0, or Z>A");
  }
  MakeLevels(filename);
}

G4NuclearLevelManager::~G4NuclearLevelManager()
{
  ClearLevels();
}

void G4NuclearLevelManager::SetNucleus(G4int Z, G4int A, const G4String& filename)
{
  if (A <= 0 || Z <= 0 || Z > A ) {
    throw G4HadronicException(__FILE__, __LINE__,
                  "==== G4NuclearLevelManager ==== (Z,A) <0, or Z>A");
  }
  if (_nucleusZ != Z || _nucleusA != A)
    {
      _nucleusA = A;
      _nucleusZ = Z;
      MakeLevels(filename);
    }
}

const G4NuclearLevel* G4NuclearLevelManager::GetLevel(G4int i) const {
  return (i>=0 && i<NumberOfLevels()) ? (*_levels)[i] : 0;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Get level to add multipole data
G4NuclearLevel* G4NuclearLevelManager::GetLevelIndexedMultipole(G4int i) {
  return (i>=0 && i<NumberOfLevels()) ? (*_levels)[i] : 0;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

const G4NuclearLevel*
G4NuclearLevelManager::NearestLevel(G4double energy, G4double) const
{
  if (NumberOfLevels() <= 0) return 0;

  G4int iNear = -1;

  //G4cout << "G4NuclearLevelManager::NearestLevel E(MeV)= "
  //	 << energy/MeV << " dEmax(MeV)= " << eDiffMax/MeV << G4endl;

  G4double diff = 1.e+10;
  for (unsigned int i=0; i<_levels->size(); ++i)
    {
      G4double e = GetLevel(i)->Energy();
      G4double eDiff = std::fabs(e - energy);
      //G4cout << i << ".   eDiff(MeV)= " << eDiff/MeV << G4endl;
      if (eDiff <= diff)
    {
      diff = eDiff;
      iNear = i;
    }
    }

  return GetLevel(iNear);	// Includes range checking on iNear
}

G4double G4NuclearLevelManager::MinLevelEnergy() const
{
  return (NumberOfLevels() > 0) ? _levels->front()->Energy() : 9999.*GeV;
}

G4double G4NuclearLevelManager::MaxLevelEnergy() const
{
  return (NumberOfLevels() > 0) ? _levels->back()->Energy() : 0.*GeV;
}

const G4NuclearLevel* G4NuclearLevelManager::HighestLevel() const
{
  return (NumberOfLevels() > 0) ? _levels->front() : 0;
}

const G4NuclearLevel* G4NuclearLevelManager::LowestLevel() const
{
  return (NumberOfLevels() > 0) ? _levels->back() : 0;
}

G4bool G4NuclearLevelManager::Read(std::ifstream& dataFile)
{
  G4bool goodRead = ReadDataLine(dataFile);

  if (goodRead) ProcessDataLine();
  return goodRead;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G4bool G4NuclearLevelManager::ReadMultipole(std::ifstream& dataFile)
{
  G4bool goodRead = ReadDataLineMultipole(dataFile);

  if (goodRead) ProcessDataLineMultipole();
  return goodRead;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// NOTE:  Standard stream I/O generates a 45 byte std::string per item!

G4bool G4NuclearLevelManager::ReadDataLine(std::ifstream& dataFile) {
  /***** DO NOT USE REGULAR STREAM I/O
  G4bool result = true;
  if (dataFile >> _levelEnergy)
    {
      dataFile >> _gammaEnergy >> _probability >> _polarity >> _halfLife
           >> _angularMomentum  >> _totalCC >> _kCC >> _l1CC >> _l2CC
           >> _l3CC >> _m1CC >> _m2CC >> _m3CC >> _m4CC >> _m5CC
           >> _nPlusCC;
    }
  else result = false;
  *****/

  // Each item will return iostream status
  return (ReadDataItem(dataFile, _levelEnergy) &&
      ReadDataItem(dataFile, _gammaEnergy) &&
      ReadDataItem(dataFile, _probability) &&
      ReadDataItem(dataFile, _polarity) &&
      ReadDataItem(dataFile, _halfLife) &&
      ReadDataItem(dataFile, _angularMomentum) &&
      ReadDataItem(dataFile, _totalCC) &&
      ReadDataItem(dataFile, _kCC) &&
      ReadDataItem(dataFile, _l1CC) &&
      ReadDataItem(dataFile, _l2CC) &&
      ReadDataItem(dataFile, _l3CC) &&
      ReadDataItem(dataFile, _m1CC) &&
      ReadDataItem(dataFile, _m2CC) &&
      ReadDataItem(dataFile, _m3CC) &&
      ReadDataItem(dataFile, _m4CC) &&
      ReadDataItem(dataFile, _m5CC) &&
      ReadDataItem(dataFile, _nPlusCC) );
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G4bool G4NuclearLevelManager::ReadDataLineMultipole(std::ifstream& dataFile) {
  // Each item will return iostream status
  return (ReadDataItem(dataFile, _levelEnergyMultipole) &&
      ReadDataItem(dataFile, _gammaEnergyMultipole) &&
      ReadDataItem(dataFile, _L1) &&
      ReadDataItem(dataFile, _L2) &&
      ReadDataItem(dataFile, _mixingRatio) );
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

G4bool
G4NuclearLevelManager::ReadDataItem(std::istream& dataFile, G4double& x)
{
  // G4bool okay = (dataFile >> buffer) != 0;		// Get next token
  // if (okay) x = strtod(buffer, NULL);
  char buffer[30];
  for(G4int i=0; i<30; ++i) { buffer[i] = 0; }
  G4bool okay = true;
  dataFile >> buffer;
  if(dataFile.fail()) { okay = false; }
  else { x = strtod(buffer, NULL); }

  return okay;
}

void G4NuclearLevelManager::ProcessDataLine()
{
  const G4double minProbability = 1e-8;

  // Assign units for dimensional quantities
  _levelEnergy *= keV;
  _gammaEnergy *= keV;
  _halfLife *= second;

  // The following adjustment is needed to take care of anomalies in
  // data files, where some transitions show up with relative probability
  // zero
  if (_probability < minProbability) _probability = minProbability;
  // the folowwing is to convert icc probability to accumulative ones
  _l1CC += _kCC;
  _l2CC += _l1CC;
  _l3CC += _l2CC;
  _m1CC += _l3CC;
  _m2CC += _m1CC;
  _m3CC += _m2CC;
  _m4CC += _m3CC;
  _m5CC += _m4CC;
  _nPlusCC += _m5CC;

  if (_nPlusCC!=0) {	// Normalize to probabilities
    _kCC /= _nPlusCC;
    _l1CC /= _nPlusCC;
    _l2CC /= _nPlusCC;
    _l3CC /= _nPlusCC;
    _m1CC /= _nPlusCC;
    _m2CC /= _nPlusCC;
    _m3CC /= _nPlusCC;
    _m4CC /= _nPlusCC;
    _m5CC /= _nPlusCC;
    _nPlusCC /= _nPlusCC;
  } else {		// Total was zero, reset to unity
    _kCC = 1;
    _l1CC = 1;
    _l2CC = 1;
    _l3CC = 1;
    _m1CC = 1;
    _m2CC = 1;
    _m3CC = 1;
    _m4CC = 1;
    _m5CC = 1;
    _nPlusCC = 1;
  }

  // G4cout << "Read " << _levelEnergy << " " << _gammaEnergy << " " << _probability << G4endl;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void G4NuclearLevelManager::ProcessDataLineMultipole()
{
  // Assign units for dimensional quantities
  _levelEnergyMultipole *= keV;
  _gammaEnergyMultipole *= keV;

  // multipole info is fine, ie _L1, etc
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void G4NuclearLevelManager::ClearLevels()
{
  _validity = false;

  if (NumberOfLevels() > 0) {
    std::for_each(_levels->begin(), _levels->end(), DeleteLevel());
    _levels->clear();
  }

  delete _levels;
  _levels = 0;
}

void G4NuclearLevelManager::MakeLevels(const G4String& filename)
{
  _validity = false;
  if (NumberOfLevels() > 0) ClearLevels();	// Discard existing data

  std::ifstream inFile(filename, std::ios::in);
  if (! inFile)
    {
#ifdef GAMMAFILEWARNING
      if (_nucleusZ > 10) G4cout << " G4NuclearLevelManager: nuclide ("
                                 << _nucleusZ << "," << _nucleusA
                                 << ") does not have a gamma levels file" << G4endl;
#endif
      return;
    }

  _levels = new G4PtrLevelVector;

  // Read individual gamma data and fill levels for this nucleus

  G4NuclearLevel* thisLevel = 0;
  G4int nData = 0;

  while (Read(inFile)) {
    thisLevel = UseLevelOrMakeNew(thisLevel);	// May create new pointer
    AddDataToLevel(thisLevel);
    nData++;					// For debugging purposes
  }

  FinishLevel(thisLevel);		// Final  must be completed by hand

  // ---- MGP ---- Don't forget to close the file
  inFile.close();

  //  G4cout << " ==== MakeLevels ===== " << nData << " data read " << G4endl;

  G4PtrSort<G4NuclearLevel>(_levels);

  _validity = true;

  return;
}

G4NuclearLevel*
G4NuclearLevelManager::UseLevelOrMakeNew(G4NuclearLevel* level)
{
  if (level && _levelEnergy == level->Energy()) return level;	// No change

  if (level) FinishLevel(level);	// Save what we have up to now

  //  G4cout << "Making a new level... " << _levelEnergy << G4endl;
  return new G4NuclearLevel(_levelEnergy, _halfLife, _angularMomentum);
}

void G4NuclearLevelManager::AddDataToLevel(G4NuclearLevel* level)
{
  if (!level) return;		// Sanity check

  level->_energies.push_back(_gammaEnergy);
  level->_weights.push_back(_probability);
  level->_polarities.push_back(_polarity);
  level->_kCC.push_back(_kCC);
  level->_l1CC.push_back(_l1CC);
  level->_l2CC.push_back(_l2CC);
  level->_l3CC.push_back(_l3CC);
  level->_m1CC.push_back(_m1CC);
  level->_m2CC.push_back(_m2CC);
  level->_m3CC.push_back(_m3CC);
  level->_m4CC.push_back(_m4CC);
  level->_m5CC.push_back(_m5CC);
  level->_nPlusCC.push_back(_nPlusCC);
  level->_totalCC.push_back(_totalCC);
}

void G4NuclearLevelManager::FinishLevel(G4NuclearLevel* level)
{
  if (!level || !_levels) return;		// Sanity check

  level->Finalize();
  _levels->push_back(level);
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G4NuclearLevel*
G4NuclearLevelManager::GetLevelAtThisEnergy(G4double levelEnergy)
{
  G4NuclearLevel* mylevel;
  G4int nLevels = NumberOfLevels();

  for (G4int i=0; i<nLevels; i++) {
    mylevel = GetLevelIndexedMultipole(i);
    if (mylevel && levelEnergy == mylevel->Energy() ) return mylevel;
  }
}

G4NuclearLevel*
G4NuclearLevelManager::GetLevelClosestToThisEnergy(G4double levelEnergy)
{
  G4NuclearLevel* thislevel;
  G4NuclearLevel* testlevel;
  G4int nLevels = NumberOfLevels();
  G4bool foundLevel = false;

  for (G4int i=0; i<nLevels; i++) {
    testlevel = GetLevelIndexedMultipole(i);
    if (testlevel) {
      if(foundLevel) { // find closest level
        if( std::fabs(testlevel->Energy()-levelEnergy) < std::fabs(thislevel->Energy()-levelEnergy) ) {
          thislevel = testlevel;
        }
      }
      else {
        thislevel = testlevel;
        foundLevel = true;
      }
    }
  }
  return thislevel;
}

G4NuclearLevel*
G4NuclearLevelManager::FindNextHighestEnergyLevel(G4double levelEnergy)
{
  G4NuclearLevel* thislevel = GetLevelAtThisEnergy(levelEnergy);
  G4NuclearLevel* higherlevel;
  G4NuclearLevel* testlevel;
  G4int nLevels = NumberOfLevels();
  G4bool foundHigherLevel = false;

  for (G4int i=0; i<nLevels; i++) {
    testlevel = GetLevelIndexedMultipole(i);
    if (testlevel && thislevel->Energy() < testlevel->Energy() ) {
      if(!foundHigherLevel) { // this is the first and maybe only higher level
        higherlevel = testlevel;
      }
      else { // if we found another higher level, see if it's closer in energy to thislevel
        if(testlevel->Energy() < higherlevel->Energy() ) {
          higherlevel = testlevel;
        }
      }
      foundHigherLevel = true;
    }
  }

  if(foundHigherLevel) {
    return higherlevel;
  }
  else {
    return 0;
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void G4NuclearLevelManager::PrintAll()
{
  G4int nLevels = NumberOfLevels();

  G4cout << " ==== G4NuclearLevelManager ==== (" << _nucleusZ << ", " << _nucleusA
     << ") has " << nLevels << " levels" << G4endl
     << "Highest level is at energy " << MaxLevelEnergy() << " MeV "
     << G4endl << "Lowest level is at energy " << MinLevelEnergy()
     << " MeV " << G4endl;

  for (G4int i=0; i<nLevels; ++i) {
    GetLevel(i)->PrintAll();
  }
}


// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void G4NuclearLevelManager::AddMulipole(G4int Z, G4int A, G4String myfilename)
{
  G4double  ji, jo, jf, L1, L1p, L2, L2p, delta1, delta2;
  G4double  gammaEnergy, higherGammaEnergy, lowerGammaEnergy, levelEnergy, higherLevelEnergy;
  G4int     nHigherGammas, indexHigher, nGammas;
  G4bool    boolTransitionToLevel;
  G4bool    boolGoodLevelToOutputToScreen;

  if (A <= 0 || Z <= 0 || Z > A ) {
    throw G4HadronicException(__FILE__, __LINE__,
          "==== G4NuclearLevelManager ==== (Z,A) <0, or Z>A");
  }

  std::ifstream inFile(myfilename, std::ios::in);
  if (! inFile)
    {
    return;
    }

  G4NuclearLevel* thisLevel = 0;

  while (ReadMultipole(inFile)) {

    thisLevel   = GetLevelAtThisEnergy(_levelEnergyMultipole);
    nGammas     = thisLevel->NumberOfGammas();

    for (G4int i=0; i<nGammas; i++) {
      if((_levelEnergyMultipole == thisLevel->Energy()) && (_gammaEnergyMultipole == thisLevel->GammaEnergies()[i]) ) {
        thisLevel->FillL1(i, _L1);
        thisLevel->FillL2(i, _L2);
        thisLevel->FillMixingRatio(i, _mixingRatio);
      }
    }
  }
  inFile.close();
  // Multipole information now added to level information

  // Compute WTheta Parameters
  G4NuclearLevel* higherLevel = 0;
  G4NuclearLevel* lowerLevel  = 0;

  G4cout << "------------------------------------------------------" << G4endl;
  G4cout << "The gamma-gamma angular coefficients are calculated " << G4endl;
  G4cout << "for all combinations of transitions provided in the " << G4endl;
  G4cout << "input files. Those transitions which are \"good\", " << G4endl;
  G4cout << "ie. the gamma energies which connect the levels are " << G4endl;
  G4cout << "within 2.0 keV, are output to the screen below for  " << G4endl;
  G4cout << "reference. " << G4endl;
  G4cout << "------------------------------------------------------" << G4endl;

  G4int nLevels = NumberOfLevels();
  for (G4int allLevels = 0; allLevels < nLevels; allLevels++) { // loop over all possible levels

    thisLevel   = GetLevelIndexedMultipole(allLevels);
    nGammas     = thisLevel->NumberOfGammas();
    levelEnergy = thisLevel->Energy();
    jo          = thisLevel->AngularMomentum();
    for (G4int i = 0; i < nGammas; i++) { // loop over the number of gammas in this level

      gammaEnergy   = thisLevel->GammaEnergies()[i];
      L2            = thisLevel->L1()[i];
      L2p           = thisLevel->L2()[i];
      delta2        = thisLevel->MixingRatio()[i];
      lowerGammaEnergy = gammaEnergy;

      // Find jf AngularMomentum
      if( (levelEnergy-gammaEnergy) >= -2.0*keV && (levelEnergy-gammaEnergy) <= 2.0*keV) { // ground state
        jf = std::fabs(_groundStateSpinAngularMomentum);
      }
      else {
        lowerLevel = GetLevelClosestToThisEnergy(levelEnergy-gammaEnergy);
        jf = lowerLevel->AngularMomentum();
      }

      // Now we need to find if a higher level feeds this level
      indexHigher = 0;
      higherLevelEnergy = thisLevel->Energy();
      for (G4int loopLevels = 0; loopLevels < nLevels; loopLevels++) { // loop over all possible levels
        higherLevel     = FindNextHighestEnergyLevel(higherLevelEnergy);
        if(higherLevel) { // found the next highest energy level
          higherLevelEnergy = higherLevel->Energy();
          nHigherGammas     = higherLevel->NumberOfGammas();

          boolTransitionToLevel = false;
          for (G4int j = 0; j < nHigherGammas; j++) { // loop over the all gammas in higher level
            // Find the gamma that feeds this level. Note: It doesn't need to be physically true!
            // The below code will find the gamma from higherLevel that comes closest to thisLevel.
            if(!boolTransitionToLevel) { // first transition found
              boolTransitionToLevel = true;
              higherGammaEnergy = higherLevel->GammaEnergies()[j];
              ji    = higherLevel->AngularMomentum();
              L1    = higherLevel->L1()[j];
              L1p   = higherLevel->L2()[j];
              delta1= higherLevel->MixingRatio()[j];
            }
            else{ // if another transition is found that is closer in energy, choose this level
              if(std::fabs( (higherLevelEnergy-levelEnergy)-higherLevel->GammaEnergies()[j] )  < std::fabs( (higherLevelEnergy-levelEnergy)-higherGammaEnergy ) ) {
                higherGammaEnergy = higherLevel->GammaEnergies()[j];
                ji    = higherLevel->AngularMomentum();
                L1    = higherLevel->L1()[j];
                L1p   = higherLevel->L2()[j];
                delta1= higherLevel->MixingRatio()[j];
              }
            }
          } // end for loop over the all gammas in higher level

          if(boolTransitionToLevel) {
            boolGoodLevelToOutputToScreen = false;
            if(std::fabs( (higherLevelEnergy-levelEnergy)-higherGammaEnergy ) < 2.0*keV )
              boolGoodLevelToOutputToScreen = true;
            thisLevel->GenerateWThetaParameters(i, indexHigher, higherLevelEnergy, lowerGammaEnergy, higherGammaEnergy, ji, jo, jf, L1, L1p, L2, L2p, delta1, delta2, boolGoodLevelToOutputToScreen);
            indexHigher++;
          }
        }
        else { // no more higher levels
          break;
        }
      } // end for loop over all possible level
    } // end for loop over the number of gammas in this level
  } // end for loop over all possible levels

}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void G4NuclearLevelManager::PrintLevels()
{
  G4int nLevels = NumberOfLevels();
  G4double efermi = G4NucleiProperties::GetNuclearMass(_nucleusA-1, _nucleusZ)
    + neutron_mass_c2 -
    G4NucleiProperties::GetNuclearMass(_nucleusA, _nucleusZ);

  G4cout << "Z= " << _nucleusZ << " A= " << _nucleusA
     << "  " << nLevels << " levels"
     << "  Efermi(MeV)= " << efermi << G4endl;

  for (G4int i=0; i<nLevels; ++i) {
    GetLevel(i)->PrintLevels();
  }
}

G4NuclearLevelManager::G4NuclearLevelManager(const G4NuclearLevelManager &right)
{
  _nucleusA = right._nucleusA;
  _nucleusZ = right._nucleusZ;
  _validity = right._validity;

  if (right._levels != 0)
    {
      _levels = new G4PtrLevelVector;
      G4int n = right._levels->size();
      G4int i;
      for (i=0; i<n; ++i)
    {
      _levels->push_back(new G4NuclearLevel(*(right.GetLevel(i))));
    }
      G4PtrSort<G4NuclearLevel>(_levels);
    }
  else
    {
      _levels = 0;
    }
}
