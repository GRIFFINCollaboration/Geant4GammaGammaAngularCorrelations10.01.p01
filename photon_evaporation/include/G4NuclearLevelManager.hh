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
// $Id: G4NuclearLevelManager.hh 87163 2014-11-26 08:46:54Z gcosmo $
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
//      Creation date: 25 October 1998
//
//      Modifications:
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions
//
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//        02 May 2003,   Vladimir Ivanchenko remove rublic copy contructor
//        06 Oct 2010, M. Kelsey -- Use object storage, not pointers, drop
//		public access to list
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVELMANAGER_HH
#define G4NUCLEARLEVELMANAGER_HH 1

#include <iosfwd>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4PtrLevelVector.hh"

class G4NuclearLevel;


class G4NuclearLevelManager
{

public:

  G4NuclearLevelManager(G4int Z, G4int A, const G4String& filename);
  G4NuclearLevelManager(const G4NuclearLevelManager & right);

  ~G4NuclearLevelManager();

  void SetNucleus(G4int Z, G4int A, const G4String& filename);

  inline G4bool IsValid() const { return _validity; }

  inline G4int NumberOfLevels() const { return (_levels ? _levels->size() : 0); }

  const G4NuclearLevel* GetLevel(G4int i) const;

  G4NuclearLevel* GetLevelIndexedMultipole(G4int i); // Evan Rand

  const G4NuclearLevel* NearestLevel(G4double energy,
                     G4double eDiffMax=1.e+8) const;

  const G4NuclearLevel* LowestLevel() const;
  const G4NuclearLevel* HighestLevel() const;

  G4double MinLevelEnergy() const;
  G4double MaxLevelEnergy() const;

  void PrintLevels();

  void PrintAll();

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void AddMulipole(G4int Z, G4int A, G4String myfilename);
  inline void GroundStateSpinAngularMomentum(G4double spin) {_groundStateSpinAngularMomentum = spin;};
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

private:
  const G4NuclearLevelManager& operator=(const G4NuclearLevelManager &right);
  G4bool operator==(const G4NuclearLevelManager &right) const;
  G4bool operator!=(const G4NuclearLevelManager &right) const;

  G4bool Read(std::ifstream& aDataFile);
  G4bool ReadDataLine(std::ifstream& dataFile);
  G4bool ReadDataItem(std::istream& dataFile, G4double& x);
  void ProcessDataLine();

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  G4bool ReadMultipole(std::ifstream& aDataFile);
  G4bool ReadDataLineMultipole(std::ifstream& dataFile);
  void ProcessDataLineMultipole();
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void MakeLevels(const G4String& filename);
  void ClearLevels();

  G4NuclearLevel* UseLevelOrMakeNew(G4NuclearLevel* level);
  void AddDataToLevel(G4NuclearLevel* level);
  void FinishLevel(G4NuclearLevel* level);

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  G4NuclearLevel* GetLevelAtThisEnergy(G4double levelEnergy);
  G4NuclearLevel* GetLevelClosestToThisEnergy(G4double levelEnergy);
  G4NuclearLevel* FindNextHighestEnergyLevel(G4double levelEnergy);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  G4int _nucleusA;
  G4int _nucleusZ;
  //G4String _fileName;
  G4bool _validity;
  G4PtrLevelVector* _levels;

  // Buffers for reading data file
  //char buffer[30];		// For doubles in scientific notation
  static G4double _levelEnergy;
  static G4double _gammaEnergy;
  static G4double _probability;
  static G4double _polarity;
  static G4double _halfLife;
  static G4double _angularMomentum;
  static G4double _kCC;
  static G4double _l1CC;
  static G4double _l2CC;
  static G4double _l3CC;
  static G4double _m1CC;
  static G4double _m2CC;
  static G4double _m3CC;
  static G4double _m4CC;
  static G4double _m5CC;
  static G4double _nPlusCC;
  static G4double _totalCC;

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  G4double _groundStateSpinAngularMomentum;
  G4double _levelEnergyMultipole;
  G4double _gammaEnergyMultipole;
  G4double _L1;
  G4double _L2;
  G4double _mixingRatio;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
};

#endif
