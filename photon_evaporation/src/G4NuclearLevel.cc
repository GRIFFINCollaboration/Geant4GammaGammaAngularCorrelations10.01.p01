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
// $Id: G4NuclearLevel.cc 86986 2014-11-21 13:00:05Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4NuclearLevel
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//
//      Creation date: 24 October 1998
//
//      Modifications:
//	  06 Oct 2010, M. Kelsey (kelsey@slac.stanford.edu)
//		Add friendship for G4NuclearLevelManager; define private
//		constructors without vectors.
//
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added IC probability when calculate the channel probabilities in
//              MakeProbabilities().
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions.
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//
//        28 October 2010, V.Ivanchenko moved copy constructor to source, cleanup
//
// -------------------------------------------------------------------

#include "G4NuclearLevel.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

G4int G4NuclearLevel::Increment(G4int aF)
{
  static G4ThreadLocal G4int instanceCount = 0;
  instanceCount+=aF;
  return instanceCount;
}

G4NuclearLevel::G4NuclearLevel()
  : _energy(0.), _halfLife(0.), _angularMomentum(0.), _nGammas(0) {
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::G4NuclearLevel(G4double energy, G4double halfLife,
                   G4double angularMomentum)
  : _energy(energy), _halfLife(halfLife), _angularMomentum(angularMomentum),
    _nGammas(0) {
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::G4NuclearLevel(G4double energy, G4double halfLife,
                   G4double angularMomentum,
                   const std::vector<G4double>& eGamma,
                   const std::vector<G4double>& wGamma,
                   const std::vector<G4double>& polarities,
                   const std::vector<G4double>& kCC, const std::vector<G4double>& l1CC,
                   const std::vector<G4double>& l2CC, const std::vector<G4double>& l3CC,
                   const std::vector<G4double>& m1CC, const std::vector<G4double>& m2CC,
                   const std::vector<G4double>& m3CC, const std::vector<G4double>& m4CC,
                   const std::vector<G4double>& m5CC, const std::vector<G4double>& nPlusCC,
                   const std::vector<G4double>& totalCC)

  : _energies(eGamma), _weights(wGamma), _polarities(polarities),
     _kCC(kCC), _l1CC(l1CC), _l2CC(l2CC), _l3CC(l3CC),
    _m1CC(m1CC), _m2CC(m2CC), _m3CC(m3CC), _m4CC(m4CC), _m5CC(m5CC),
    _nPlusCC(nPlusCC), _totalCC(totalCC),
    _energy(energy), _halfLife(halfLife), _angularMomentum(angularMomentum)
{
  Finalize();
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::~G4NuclearLevel()
{
 // G4cout << "####### Decrementing "<<Increment(-1)<<G4endl;
}

G4bool G4NuclearLevel::operator==(const G4NuclearLevel &right) const
{
  return (this == (G4NuclearLevel *) &right);
}


G4bool G4NuclearLevel::operator!=(const G4NuclearLevel &right) const
{
  return (this != (G4NuclearLevel *) &right);
}

G4bool G4NuclearLevel::operator<(const G4NuclearLevel &right) const
{
  if (_energy < right.Energy()) return true;
  else return false;
}

const std::vector<G4double>& G4NuclearLevel::GammaEnergies() const
{
  return _energies;
}

const std::vector<G4double>& G4NuclearLevel::GammaWeights() const
{
  return _weights;
}


const std::vector<G4double>& G4NuclearLevel::GammaProbabilities() const
{
  return _prob;
}


const std::vector<G4double>& G4NuclearLevel::GammaCumulativeProbabilities() const
{
  return _cumProb;
}


const std::vector<G4double>& G4NuclearLevel::GammaPolarities() const
{
  return _polarities;
}

const std::vector<G4double>& G4NuclearLevel::KConvertionProbabilities() const
{
  return _kCC;
}

const std::vector<G4double>& G4NuclearLevel::L1ConvertionProbabilities() const
{
  return _l1CC;
}

const std::vector<G4double>& G4NuclearLevel::L2ConvertionProbabilities() const
{
  return _l2CC;
}

const std::vector<G4double>& G4NuclearLevel::L3ConvertionProbabilities() const
{
  return _l3CC;
}

const std::vector<G4double>& G4NuclearLevel::M1ConvertionProbabilities() const
{
  return _m1CC;
}

const std::vector<G4double>& G4NuclearLevel::M2ConvertionProbabilities() const
{
  return _m2CC;
}

const std::vector<G4double>& G4NuclearLevel::M3ConvertionProbabilities() const
{
  return _m3CC;
}

const std::vector<G4double>& G4NuclearLevel::M4ConvertionProbabilities() const
{
  return _m4CC;
}

const std::vector<G4double>& G4NuclearLevel::M5ConvertionProbabilities() const
{
  return _m5CC;
}

const std::vector<G4double>& G4NuclearLevel::NPlusConvertionProbabilities() const
{
  return _nPlusCC;
}

const std::vector<G4double>& G4NuclearLevel::TotalConvertionProbabilities() const
{
  return _totalCC;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

const std::vector<G4double>& G4NuclearLevel::L1() const
{
  return _L1;
}

const std::vector<G4double>& G4NuclearLevel::L2() const
{
  return _L2;
}

const std::vector<G4double>& G4NuclearLevel::MixingRatio() const
{
  return _mixingRatio;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A2() const
{
  return _a2;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A4() const
{
  return _a4;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A6() const
{
  return _a6;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A8() const
{
  return _a8;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A10() const
{
  return _a10;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::MaxWTheta() const
{
  return _maxWTheta;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::HigherLevelEnergy() const
{
  return _higherLevelEnergy;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

G4double G4NuclearLevel::Energy() const
{
  return _energy;
}

G4double G4NuclearLevel::AngularMomentum() const
{
  return _angularMomentum;
}

G4double G4NuclearLevel::HalfLife() const
{
  return _halfLife;
}

G4int G4NuclearLevel::NumberOfGammas() const
{
  return _nGammas;
}

void G4NuclearLevel::PrintAll() const
{
  G4cout << "---- Level energy = " << _energy << ", angular momentum = "
     << _angularMomentum << ", half life " << _halfLife
     << ", " << _nGammas << " photons" << G4endl;
  G4int i;
  G4cout << "     Gammas: ";
  for (i=0; i<_nGammas; i++) { G4cout << _energies[i] << " "; }
  G4cout << G4endl << "     Weights: ";
  for (i=0; i<_nGammas; i++) { G4cout << _weights[i] << " "; }
  G4cout << G4endl << "     Relative transition probabilities ";
  for (i=0; i<_nGammas; i++) { G4cout << _prob[i] << " "; }
  G4cout << G4endl << "     Cumulative probabilities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _cumProb[i] << " "; }
  G4cout << G4endl << "     Polarities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _polarities[i] << " "; }
  G4cout << G4endl;
}

void G4NuclearLevel::PrintLevels() const
{
  G4cout << "   Eexc(MeV)= " << _energy
     << " Time(ns)= " << _halfLife/ns << "  Ntrans= " << _nGammas
     << G4endl;
}

void G4NuclearLevel::Finalize() {
  _nGammas = _energies.size();
  MakeProbabilities();
  MakeCumProb();
  MakeMultipoleData(); // Evan Rand
}

void G4NuclearLevel::MakeProbabilities()
{
  G4double sum = 0.;
  G4int i = 0;
  for (i=0; i<_nGammas; i++) {
    sum += _weights[i]*(1.+_totalCC[i]);
  }

  if (sum <= 0.) _prob.resize(_nGammas, 1./_nGammas);	// Fast fill
  else {
    _prob.reserve(_nGammas);
    for (i=0; i<_nGammas; i++) {
      _prob.push_back(_weights[i]*(1.+_totalCC[i])/sum);
    }
  }
}


void G4NuclearLevel::MakeCumProb()
{
  if (_nGammas <= 0) return;

  _cumProb.reserve(_nGammas);

  G4double sum = _prob[0];
  _cumProb.push_back(sum);

  for (G4int i=1; i<_nGammas; i++) {
    sum += _prob[i];
    _cumProb.push_back(sum);
  }
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void G4NuclearLevel::MakeMultipoleData()
{
  if (_nGammas <= 0) return;

  for (G4int i=0; i<_nGammas; i++) {
    _L1.push_back(0);
    _L2.push_back(0);
    _mixingRatio.push_back(0);

    _higherLevelEnergy.push_back(std::vector< G4double >());
    _a2.push_back(std::vector< G4double >());
    _a4.push_back(std::vector< G4double >());
    _a6.push_back(std::vector< G4double >());
    _a8.push_back(std::vector< G4double >());
    _a10.push_back(std::vector< G4double >());
    _maxWTheta.push_back(std::vector< G4double >());
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

G4NuclearLevel& G4NuclearLevel::operator=(const G4NuclearLevel &right)
{
  if(this != &right)
    {
      _energies = right._energies;
      _weights =right._weights;
      _prob =right._prob;
      _cumProb =right._cumProb;
      _polarities =right._polarities;
      _kCC = right._kCC;
      _l1CC =right._l1CC;
      _l2CC =right._l2CC;
      _l3CC =right._l3CC;
      _m1CC = right._m1CC;
      _m2CC = right._m2CC;
      _m3CC = right._m3CC;
      _m4CC = right._m4CC;
      _m5CC = right._m5CC;
      _nPlusCC = right._nPlusCC;
      _totalCC = right._totalCC;
      _energy = right._energy;
      _halfLife = right._halfLife;
      _angularMomentum = right._angularMomentum;
      _nGammas = right._nGammas;
    }
  return *this;
}

G4NuclearLevel::G4NuclearLevel(const G4NuclearLevel &right)
{
  _energies = right._energies;
  _weights =right._weights;
  _prob =right._prob;
  _cumProb =right._cumProb;
  _polarities =right._polarities;
  _kCC = right._kCC;
  _l1CC =right._l1CC;
  _l2CC =right._l2CC;
  _l3CC =right._l3CC;
  _m1CC = right._m1CC;
  _m2CC = right._m2CC;
  _m3CC = right._m3CC;
  _m4CC = right._m4CC;
  _m5CC = right._m5CC;
  _nPlusCC = right._nPlusCC;
  _totalCC = right._totalCC;
  _energy = right._energy;
  _halfLife = right._halfLife;
  _angularMomentum = right._angularMomentum;
  _nGammas = right._nGammas;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void G4NuclearLevel::FillHigherLevelEnergy(G4int gamma_i, G4int hlevel_i, G4double hlevel_energy)
{
  G4int nLevels = _higherLevelEnergy[gamma_i].size();
  if (hlevel_i > nLevels){
    G4cout << "Error! Skipped a level?" << G4endl;
    exit(1);
  }
  _higherLevelEnergy[gamma_i].push_back(hlevel_energy);
}

void G4NuclearLevel::GenerateWThetaParameters(G4int gamma_i, G4int hlevel_i, G4double hlevel_energy, G4double lowerGammaEnergy, G4double higherGammaEnergy, G4double ji, G4double jo, G4double jf, G4int L1, G4int L1p, G4int L2, G4int L2p, G4double delta1, G4double delta2, G4bool boolGoodLevelToOutputToScreen) {
  G4double a2, a4, a6, a8, a10, P2, P4, wTheta, thetaRad, thisMaxWTheta;

  // Save higher level energy
  FillHigherLevelEnergy(gamma_i, hlevel_i, hlevel_energy);

  // An additional check that these variables are positive
  G4double Ji = std::fabs(ji);
  G4double Jo = std::fabs(jo);
  G4double Jf = std::fabs(jf);

  a2 = B(2,Jo,Ji,L1,L1p,delta1)*A(2,Jf,Jo,L2,L2p,delta2);
  a4 = B(4,Jo,Ji,L1,L1p,delta1)*A(4,Jf,Jo,L2,L2p,delta2);
  a6 = B(6,Jo,Ji,L1,L1p,delta1)*A(6,Jf,Jo,L2,L2p,delta2);
  a8 = B(8,Jo,Ji,L1,L1p,delta1)*A(8,Jf,Jo,L2,L2p,delta2);
  a10= B(10,Jo,Ji,L1,L1p,delta1)*A(10,Jf,Jo,L2,L2p,delta2);

  if(boolGoodLevelToOutputToScreen) {
    G4cout << "---------- gamma-gamma angular coefficients ----------" << G4endl;
    G4cout << "highest level in cascade is " << hlevel_energy/keV << " keV" << G4endl;
    G4cout << "first gamma in cascade is " << higherGammaEnergy/keV << " keV" << G4endl;
    G4cout << "second gamma in cascade is " << lowerGammaEnergy/keV << " keV" << G4endl;
    G4cout << "ji = " << ji << " --> jo = " << jo <<  " --> jf = " << jf << G4endl;
    G4cout << "L1 = " << L1 << ", L1p = " << L1p <<  " - " << "L2 = " << L2 << ", L2p = " << L2p << G4endl;
    G4cout << "delta1 = " << delta1 << " - " << "delta2 = " << delta2 << G4endl;
    G4cout << "a2  = " << a2 << G4endl;
    G4cout << "a4  = " << a4 << G4endl;
    G4cout << "a6  = " << a6 << G4endl;
    G4cout << "a8  = " << a8 << G4endl;
    G4cout << "a10 = " << a10 << G4endl;
    G4cout << "------------------------------------------------------" << G4endl;
  }
  _a2[gamma_i].push_back(a2);
  _a4[gamma_i].push_back(a4);
  _a6[gamma_i].push_back(a6);
  _a8[gamma_i].push_back(a8);
  _a10[gamma_i].push_back(a10);

  // find max wTheta to normalize to later.
  thisMaxWTheta = 0;
  for(G4int i = 0; i < 180000 ; i++ ) {
    thetaRad = (i/1000)*(M_PI/180.0);
    P2 = LegendreP(2,cos(thetaRad));
    P4 = LegendreP(4,cos(thetaRad));
    wTheta = (1+a2*P2+a4*P4)*std::sin(thetaRad);
    if (wTheta > thisMaxWTheta)
      thisMaxWTheta = wTheta;
  }
  // Add a positive 1% error on maxWTheta, this ensures we monte carlo over the entire y-range later.
  // This may be an unnecessary check, but better to be safe than sorry.
  // It will be negligibly slower to add this insurance.
  _maxWTheta[gamma_i].push_back(thisMaxWTheta+(0.01*thisMaxWTheta));
}


G4double G4NuclearLevel::ClebschGordan(G4double j1, G4double m1, G4double j2, G4double m2, G4double j, G4double m)
{
  // Conditions check
  if( 2*j1 !=   std::floor(2*j1) || 2*j2 !=   std::floor(2*j2) || 2*j !=   std::floor(2*j) || 2*m1 !=   std::floor(2*m1) || 2*m2 !=   std::floor(2*m2) || 2*m !=   std::floor(2*m) )
  {
    G4cout << "All arguments must be integers or half-integers." << G4endl;
    return 0;
  }
  if(m1 + m2 != m)
  {
    //G4cout << "m1 + m2 must equal m." << G4endl;
    return 0;
  }
  if( j1 - m1 !=   std::floor ( j1 - m1 ) )
  {
    //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
    return 0;
  }
  if( j2 - m2 !=   std::floor ( j2 - m2 ) )
  {
    //G4cout << "2*j2 and 2*m2 must have the same parity" << G4endl;
    return 0;
  }
  if( j - m !=   std::floor ( j - m ) )
  {
    //G4cout << "2*j and 2*m must have the same parity" << G4endl;
    return 0;
  }
  if(j > j1 + j2 || j < std::fabs(j1 - j2))
  {
    //G4cout << "j is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m1) > j1)
  {
    //G4cout << "m1 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m2) > j2)
  {
    //G4cout << "m2 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m) > j)
  {
    //warning('m is out of bounds." << G4endl;
    return 0 ;
  }
  G4double term, cg;
  G4double term1 = std::pow((((2*j+1)/Factorial(j1+j2+j+1))*Factorial(j2+j-j1)*Factorial(j+j1-j2)*Factorial(j1+j2-j)*Factorial(j1+m1)*Factorial(j1-m1)*Factorial(j2+m2)*Factorial(j2-m2)*Factorial(j+m)*Factorial(j-m)),(0.5));
  G4double sum = 0;

  for(G4int k = 0 ; k <= 99 ; k++ )
  {
    if( (j1+j2-j-k < 0) || (j-j1-m2+k < 0) || (j-j2+m1+k < 0) || (j1-m1-k < 0) || (j2+m2-k < 0) )
    {
    }
    else
    {
      term = Factorial(j1+j2-j-k)*Factorial(j-j1-m2+k)*Factorial(j-j2+m1+k)*Factorial(j1-m1-k)*Factorial(j2+m2-k)*Factorial(k);
      if((k%2) == 1)
      {
        term = -1*term;
      }
      sum = sum + 1.0/term;
    }
  }

  cg = term1*sum;
  return cg;
  // Reference: An Effective Algorithm for Calculation of the C.G.
  // Coefficients Liang Zuo, et. al.
  // J. Appl. Cryst. (1993). 26, 302-304
}

G4double G4NuclearLevel::Factorial(G4double value)
{
  G4double fac;
  if(value > 1)
  {
    fac = value*Factorial(value-1);
  }
  else
  {
    fac = 1;
  }
  return fac;
}


G4double G4NuclearLevel::LegendreP(G4int n, G4double x)
{
  G4double PP;
  if(n == 0) {
    PP = 1;
  }
  else if(n == 2) {
    PP = (1.0/2.0)*(3*(std::pow(x,2))-1);
  }
  else if(n == 4) {
    PP = (1.0/8.0)*(35*(std::pow(x,4))-30*(std::pow(x,2))+3);
  }
  else if(n == 6) {
    PP = (1.0/16.0)*(231*(std::pow(x,6))-315*(std::pow(x,4))+105*(std::pow(x,2))-5);
  }
  else if(n == 8) {
    PP = (1.0/128.0)*(6435*(std::pow(x,8))-12012*(std::pow(x,6))+6930*(std::pow(x,4))-1260*(std::pow(x,2))+35);
  }
  else if(n == 10) {
    PP = (1.0/256.0)*(46189*(std::pow(x,10))-109395*(std::pow(x,8))+90090*(std::pow(x,6))-30030*(std::pow(x,4))+3465*(std::pow(x,2))-63);
  }
  else {
    G4cout << "Legendre Polynomial Not Found" << G4endl;
    exit(1);
  }
  return PP;
}

G4double G4NuclearLevel::F(G4int k, G4double jf, G4int L1, G4int L2, G4double ji)
{
  G4double out;
  G4double CG = ClebschGordan(L1,1,L2,-1,k,0);
  if(CG == 0)
  {
    return 0;
  }
  G4double W = RacahW(ji,ji,L1,L2,k,jf);
  if(W == 0)
  {
    return 0;
  }
  out = std::pow((-1),(jf-ji-1))*(std::pow((2*L1+1)*(2*L2+1)*(2*ji+1),(1.0/2.0)))*CG*W;
  return out;
  // Reference: Tables of coefficients for angular distribution of gamma rays from aligned nuclei
  // T. Yamazaki. Nuclear Data A, 3(1):1?23, 1967.
}

G4double G4NuclearLevel::B(G4int k, G4double ji, G4double jf, G4int L1, G4int L2, G4double delta)
{
  G4double out;
  out = (1/(1+std::pow(delta,2)))*(  F(k,jf,L1,L1,ji)+(std::pow((-1),((G4double)(L1+L2))))*2*delta*F(k,jf,L1,L2,ji)+delta*delta*F(k,jf,L2,L2,ji) );
  return out;
}


G4double G4NuclearLevel::RacahW(G4double a, G4double b, G4double c, G4double d, G4double e, G4double f)
{
  G4double out = std::pow((-1),(a+b+d+c))*Wigner6j(a,b,e,d,c,f);
  return out;
}

G4double G4NuclearLevel::Wigner6j(G4double J1, G4double J2, G4double J3, G4double J4, G4double J5, G4double J6)
{
  // Conditions check
  if(J3 > J1 + J2 || J3 < std::fabs(J1 - J2))
  {
    //G4cout << "first J3 triange condition not satisfied. J3 > J1 + J2 || J3 < std::fabs(J1 - J2)" << G4endl;
    return 0;
  }
  if(J3 > J4 + J5 || J3 < std::fabs(J4 - J5))
  {
    //G4cout << "second J3 triange condition not satisfied. J3 > J4 + J5 || J3 < std::fabs(J4 - J5)" << G4endl;
    return 0;
  }
  if(J6 > J2 + J4 || J6 < std::fabs(J2 - J4))
  {
    //G4cout << "first J6 triange condition not satisfied. J6 > J2 + J4 || J6 < std::fabs(J2 - J4)" << G4endl;
    return 0;
  }
  if(J6 > J1 + J5 || J6 < std::fabs(J1 - J5))
  {
    //G4cout << "second J6 triange condition not satisfied. J6 > J1 + J5 || J6 < std::fabs(J1 - J5)" << G4endl;
    return 0;
  }

  G4double j1 = J1;
  G4double j2 = J2;
  G4double j12 = J3;
  G4double j3 = J4;
  G4double j = J5;
  G4double j23 = J6;
  G4double sum = 0;

  for(G4double m1 = -j1 ; m1 <= j1 ; m1++ )
  {
    for(G4double m2 = -j2 ; m2 <= j2 ; m2++ )
    {
      for(G4double m3 = -j3 ; m3 <= j3 ; m3++ )
      {
        for(G4double m12 = -j12 ; m12 <= j12 ; m12++ )
        {
          for(G4double m23 = -j23 ; m23 <= j23 ; m23++ )
          {
            for(G4double m = -j ; m <= j ; m++ )
            {
              sum = sum + std::pow((-1),(j3+j+j23-m3-m-m23))*Wigner3j(j1,j2,j12,m1,m2,m12)*Wigner3j(j1,j,j23,m1,-m,m23)*Wigner3j(j3,j2,j23,m3,m2,-m23)*Wigner3j(j3,j,j12,-m3,m,m12);
            }
          }
        }
      }
    }
  }
  return sum;
}

G4double G4NuclearLevel::Wigner3j(G4double j1, G4double j2, G4double j3, G4double m1, G4double m2, G4double m3)
{
  // Conditions check
  if( 2*j1 !=   std::floor(2*j1) || 2*j2 !=   std::floor(2*j2) || 2*j3 !=   std::floor(2*j3) || 2*m1 !=   std::floor(2*m1) || 2*m2 !=   std::floor(2*m2) || 2*m3 !=   std::floor(2*m3) )
  {
    G4cout << "All arguments must be integers or half-integers." << G4endl;
    return 0;
  }
  if(m1 + m2 + m3 != 0)
  {
    //G4cout << "m1 + m2 + m3 must equal zero." << G4endl;
    return 0;
  }
  if( j1 + j2 + j3 !=   std::floor(j1 + j2 + j3) )
  {
    //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
    return 0;
  }
  if(j3 > j1 + j2 || j3 < std::fabs(j1 - j2))
  {
    //G4cout << "j3 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m1) > j1)
  {
    //G4cout << "m1 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m2) > j2)
  {
    //G4cout << "m2 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m3) > j3)
  {
    //G4cout << "m3 is out of bounds." << G4endl;
    return 0;
  }

  G4double out = (std::pow((-1),(j1-j2-m3)))/(std::pow((2*j3+1),(1.0/2.0)))*ClebschGordan(j1,m1,j2,m2,j3,-1*m3);
  return out;
}

G4double G4NuclearLevel::A(G4int k, G4double ji, G4double jf, G4int L1, G4int L2, G4double delta)
{
  G4double out;
  out = (1/(1+std::pow(delta,2)))*(F(k,ji,L1,L1,jf)+2*delta*F(k,ji,L1,L2,jf)+delta*delta*F(k,ji,L2,L2,jf));
  return out;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
