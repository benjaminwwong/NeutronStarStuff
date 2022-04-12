#ifndef DarkMatterNREFT_hh
#define DarkMatterNREFT_hh

#include "ModelSpace.hh"
#include "Operator.hh"



// Operators for Dark matter scattering as defined in
// Gazda, Catena, and Forssen, Phys Rev D 95 103011 (2017)
// The implementation is essentially copied from Daniel Gazda's
// fortran implementation and translated to C++. Any errors
// are almost certainly Ragnar's fault.
//                                         - SRS  June 21 2018
namespace DM_NREFT
{

  double jho( int a, int la, int nb, int lb, int L, double y);
  double jdmho( int a, int la, int nb, int lb, int L, double y);
  double jdpho( int a, int la, int nb, int lb, int L, double y);

  double PhiF(  int la, int j2a, int j2b );
  double PhiS1( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );
  double PhiS2( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );
  double PhiS3( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );
  double PhiS4( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );

  Operator Ms(      ModelSpace&, std::string IsoSV, int J, int L, double q ); // helper function
  Operator Mg(      ModelSpace&, std::string IsoSV, int J, int L, double q ); // helper function

  // These are the operators
  Operator M(       ModelSpace&, std::string IsoSV, int J, double q );
  Operator Sigma(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Sigmap(  ModelSpace&, std::string IsoSV, int J, double q );
  Operator Sigmapp( ModelSpace&, std::string IsoSV, int J, double q );
  Operator Delta(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Deltap(  ModelSpace&, std::string IsoSV, int J, double q );
  Operator Phip(    ModelSpace&, std::string IsoSV, int J, double q );
  Operator Phitp(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Phipp(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Omega(   ModelSpace&, std::string IsoSV, int J, double q ); 
  Operator Omegat(  ModelSpace&, std::string IsoSV, int J, double q ); // Not implemented

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
// Operators for spin-dependent WIMP scattering as defined in
// P. Klos, J. Menendez, D. Gazit and A. Schwenk,
// Phys Rev D 88, 083516 (2013) and Phys Rev D 89, 029901(E) (2014).
// The implementation is based on José de Jesús Padua Arguelles's code.
//	Copyright (C) 2020 Baishan Hu
//	E-mail: bhu@triumf.ca
//	Last upgrade: July 2020
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

 //const double gA = 1.2670; // From Park et al. (2003) "Parameter Free Effective Field Theory ..."
 const double gA = 1.29; // From N3LO(EM), NNLOsat
 // Values of constants from Men\'edez, Gazit, Schwenk (2012)
 const double gpipn = 13.05;
 const double Fpi = 92.4;
 const double mpi = 138.04; // This is the value of Men\'edez et al. (isospin average)
 const double LambdaA = 1040.0; // in MeV
 const double LambdaChi = 700;  // in MeV

// These are helper functions
 double Bessel(ModelSpace&, int np, int lp, int Lp, int n, int l, double p);
 double Binomial_Coeff(double n, double k);
 double delta_a1(double rho, double c3, double c4, double cD, double p, double Fpi, double mpi);
 double delta_a1_P(double rho, double c3, double c4, double p, double Fpi, double mpi);
 double Integral_c6(double p, double kF, double mpi);
 Operator Y_tensor_Sigma(ModelSpace&, int L, int Lp, double p);
 Operator Bessel_Op(ModelSpace&, int Lp, double p);

// These are the operators
 Operator Longitudinal_5(ModelSpace&, double rho, double c3, double c4, double cD, std::string sfact_type, int L, double p);
 Operator Electric_Transverse_5(ModelSpace&, double rho, double c3, double c4, double cD, std::string sfact_type, int L, double p);
 Operator Magnetic_Transverse_5(ModelSpace&, double rho, double c3, double c4, double cD, std::string sfact_type, int L, double p);

}


#endif
