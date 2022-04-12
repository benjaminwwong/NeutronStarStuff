
#include "Commutator.hh"
#include "Operator.hh"
//#include "DaggerOperator.hh"
#include "TwoBodyME.hh"
#include "ThreeBodyME.hh"
#include "armadillo"
#include "PhysicalConstants.hh" // for SQRT2
#include "AngMom.hh"
#include <map>
#include <deque>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>



/// Commutator expressions for second-quantized operators
namespace Commutator {

bool use_goose_tank_correction = false;
bool use_brueckner_bch = false;
bool use_imsrg3 = false;
bool use_imsrg3_n7 = false;
bool only_2b_omega = false;
bool perturbative_triples = false;
double bch_transform_threshold = 1e-9;
double bch_product_threshold = 1e-4;


void Set_BCH_Transform_Threshold(double x)
{bch_transform_threshold=x;}

void Set_BCH_Product_Threshold(double x)
{bch_product_threshold=x;}

void SetUseBruecknerBCH(bool tf)
{use_brueckner_bch = tf;}

void SetUseGooseTank(bool tf)
{use_goose_tank_correction = tf;}

void SetUseIMSRG3(bool tf)
{use_imsrg3 = tf;}

void SetUseIMSRG3N7(bool tf)
{use_imsrg3_n7 = tf;}

void SetOnly2bOmega(bool tf)
{only_2b_omega = tf;}

//Operator Operator::Commutator( Operator& opright)
/// Returns \f$ Z = [X,Y] \f$
Operator Commutator( const Operator& X, const Operator& Y)
{
//  int jrank = std::max(X.rank_J,Y.rank_J);
//  int trank = std::max(X.rank_T,Y.rank_T);
//  int parity = (X.parity+Y.parity)%2;
//  int particlerank = std::max(X.particle_rank,Y.particle_rank);
  int xrank = X.rank_J + X.rank_T + X.parity;
  int yrank = Y.rank_J + Y.rank_T + Y.parity;
  int xlegs = X.GetNumberLegs();
  int ylegs = Y.GetNumberLegs();

  X.modelspace->PreCalculateSixJ();

  if (xrank==0)
   {
      if ( (xlegs%2==0) and (ylegs%2==1) )
      {
        return CommutatorScalarDagger(X,Y);
      }
      if ( (xlegs%2==1) and (ylegs%2==0) )
      {
        return -CommutatorScalarDagger(Y,X);
      }
      if ( (xlegs%2==1) and (ylegs%2==1) )
      {
        std::cout << "TROUBLE!!!! Called commutator with two Dagger operators. This isn't implemented. Dying..." << std::endl;
	std::exit(EXIT_FAILURE);
      }
      if (yrank==0)
      {
         return CommutatorScalarScalar(X,Y); // [S,S]
      }
      else
      {
         return CommutatorScalarTensor(X,Y); // [S,T]
      }
   }
   else if(yrank==0)
   {
      return -CommutatorScalarTensor(Y,X); // [T,S]
   }
   else
   {
      std::cout << "In Tensor-Tensor because X.rank_J = " << X.rank_J << "  X.rank_T = " << X.rank_T << "  X.parity = " << X.parity << "   ";
      std::cout <<                        "  Y.rank_J = " << Y.rank_J << "  Y.rank_T = " << Y.rank_T << "  Y.parity = " << Y.parity << std::endl;
      std::cout << " Tensor-Tensor commutator not yet implemented." << std::endl;
   }
   return  0*Y;
}




/// Commutator where \f$ X \f$ and \f$Y\f$ are scalar operators.
/// Should be called through Commutator()
Operator CommutatorScalarScalar( const Operator& X, const Operator& Y) 
{
//   std::cout << "Calling CommutatorScalarScalar..." << std::endl;
   X.profiler.counter["N_ScalarCommutators"] += 1;
   double t_css = omp_get_wtime();
   int z_Jrank = std::max( X.GetJRank(),Y.GetJRank());
   int z_Trank = std::max( X.GetTRank(),Y.GetTRank());
   int z_parity = (X.GetParity()+Y.GetParity())%2;
   int z_particlerank = std::max(X.GetParticleRank(),Y.GetParticleRank());
   if ( use_imsrg3 )  z_particlerank = std::max(z_particlerank, 3);
   ModelSpace& ms = *(Y.GetModelSpace());
   Operator Z( ms, z_Jrank, z_Trank, z_parity, z_particlerank );
//   Operator Z( *(Y.GetModelSpace()), std::max(X.GetJRank(),Y.GetJRank()), std::max(X.GetTRank(),Y.GetTRank()), (X.GetParity()+Y.GetParity())%2, std::max(X.GetParticleRank(),Y.GetParticleRank()) );

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();


   if ( Z.GetParticleRank() > 2 )
   {
     Z.ThreeBody.SwitchToPN_and_discard();
   }

   if ( not Z.IsAntiHermitian() )
   {
      comm110ss(X, Y, Z);
      if (X.particle_rank>1 and Y.particle_rank>1)
        comm220ss(X, Y, Z) ;
   }

   double t_start = omp_get_wtime();
   comm111ss(X, Y, Z);
   X.profiler.timer["comm111ss"] += omp_get_wtime() - t_start;

    t_start = omp_get_wtime();
   comm121ss(X, Y, Z);
   X.profiler.timer["comm121ss"] += omp_get_wtime() - t_start;

    t_start = omp_get_wtime();
   comm122ss(X, Y, Z); 
   X.profiler.timer["comm122ss"] += omp_get_wtime() - t_start;


   if (X.particle_rank>1 and Y.particle_rank>1)
   {
     t_start = omp_get_wtime();
     comm222_pp_hh_221ss(X, Y, Z);
     X.profiler.timer["comm222_pp_hh_221ss"] += omp_get_wtime() - t_start;
      
     t_start = omp_get_wtime();
     comm222_phss(X, Y, Z);
     X.profiler.timer["comm222_phss"] += omp_get_wtime() - t_start;
   }


   if (use_imsrg3)
   {
       X.profiler.counter["N_ScalarCommutators_3b"] += 1;

       // This gets the perturbative energy from the induced 3 body
       std::cout << " comm330 " << std::endl;
       comm330ss(X, Y, Z); // scales as n^6

//     Maybe not so important, but I think relatively cheap
       std::cout << " comm331 " << std::endl;
       comm331ss(X, Y, Z); // scales as n^7

//     This one is essential. If it's not here, then there are no induced 3 body terms
       std::cout << " comm223 " << std::endl;
       comm223ss(X, Y, Z); // scales as n^7
//       comm223ss_debug(X, Y, Z); // scales as n^7

//     if (X.GetParticleRank()>2 or Y.GetParticleRank()>2)
//     {
       // Demonstrated that this can have some effect
       std::cout << " comm231 " << std::endl;
       comm231ss(X, Y, Z);  // scales as n^6

//     no demonstrated effect yet, but it's cheap
       std::cout << " comm132 " << std::endl;
       comm132ss(X, Y, Z); // scales as n^6

//     one of the two most important IMSRG(3) terms
       std::cout << " comm232 " << std::endl;
       comm232ss(X, Y, Z);   // this is the slowest n^7 term
//       comm232ss_debug(X, Y, Z);   // this is the slowest n^7 term

//     important for suppressing off-diagonal H3
       std::cout << " comm133 " << std::endl;
       comm133ss(X, Y, Z);  // scales as n^7, but really more like n^6

      if ( not use_imsrg3_n7 )
      {
////    Not too bad, though naively n^8
       std::cout << " comm233_pp_hh " << std::endl;
       comm233_pp_hhss(X, Y, Z);
//       comm233_pp_hhss_debug(X, Y, Z);
//       X.profiler.timer["comm233_pp_hhss"] += omp_get_wtime() - t_start;

//     This one is super slow too. It involves 9js
//     mat mult makes everything better!
       std::cout << " comm233_ph " << std::endl;
       comm233_phss(X, Y, Z);
//       comm233_phss_debug(X, Y, Z);

//       not too bad, though naively n^8
       std::cout << " comm332_ppph_hhhp " << std::endl;
       comm332_ppph_hhhpss(X, Y, Z);

//      naively n^8, but reasonably fast when implemented as a mat mult
       std::cout << " comm332_pphh " << std::endl;
       comm332_pphhss(X, Y, Z);
//       comm332_pphhss_debug(X, Y, Z);

//       naively n^9 but pretty fast as a mat mult
       std::cout << " comm333_ppp_hhhss " << std::endl;
       comm333_ppp_hhhss(X, Y, Z);

//     This one works, but it's incredibly slow.  naively n^9.
//     Much improvement by going to mat mult
       std::cout << " comm333_pph_hhpss " << std::endl;
       comm333_pph_hhpss(X, Y, Z);
//       comm333_pph_hhpss_debug(X, Y, Z);

      } // if not use_imsrg3_n7

     // after going through once, we've stored all the 6js (and maybe 9js), so we can run in OMP loops from now on
     X.modelspace->scalar3b_transform_first_pass = false;
   }


   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();

   X.profiler.timer["CommutatorScalarScalar"] += omp_get_wtime() - t_css;
   return Z;

}


/// Commutator \f$[X,Y]\f$ where \f$ X \f$ is a scalar operator and \f$Y\f$ is a tensor operator.
/// Should be called through Commutator()
Operator CommutatorScalarTensor( const Operator& X, const Operator& Y) 
{
   X.profiler.counter["N_TensorCommutators"] += 1;
   double t_cst = omp_get_wtime();
   Operator Z = Y; // This ensures the commutator has the same tensor rank as Y
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

   double t_start = omp_get_wtime();
   comm111st(X, Y, Z);
   X.profiler.timer["comm111st"] += omp_get_wtime() - t_start;
   t_start = omp_get_wtime();
   comm121st(X, Y, Z);
   X.profiler.timer["comm121st"] += omp_get_wtime() - t_start;

   t_start = omp_get_wtime();
   comm122st(X, Y, Z);
   X.profiler.timer["comm122st"] += omp_get_wtime() - t_start;
   t_start = omp_get_wtime();
   comm222_pp_hh_221st(X, Y, Z);
   X.profiler.timer["comm222_pp_hh_221st"] += omp_get_wtime() - t_start;
   t_start = omp_get_wtime();
   comm222_phst(X, Y, Z);
   X.profiler.timer["comm222_phst"] += omp_get_wtime() - t_start;

  if (use_imsrg3 and X.GetJRank()==0 and Y.GetJRank()==0 and Z.GetJRank()==0 )
  {
    if ( Z.GetParticleRank()<3 )
    {
      Z.ThreeBody.SwitchToPN_and_discard();
    }
    std::cout << "tensor comm223ss" << std::endl;
    comm223ss(X,Y,Z);
    std::cout << "tensor comm232ss" << std::endl;
    comm232ss_slow(X,Y,Z);
    std::cout << "tensor comm231ss" << std::endl;
    comm231ss_slow(X,Y,Z);
  }


   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();
   X.profiler.timer["CommutatorScalarTensor"] += omp_get_wtime() - t_cst;
   return Z;
}




/// Commutator where \f$ X \f$ is a scalar and \f$Y\f$ is a dagger (i.e. it creates an additional particle).
/// Should be called through Commutator()
Operator CommutatorScalarDagger( const Operator& X, const Operator& Y) 
{
//   std::cout << "Call " << __func__ << std::endl;
   X.profiler.counter["N_DaggerCommutators"] += 1;
   double t_css = omp_get_wtime();
   Operator Z = Y;
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseThreeLeg();

//   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
//   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
//   else Z.SetNonHermitian();

   comm211sd( X, Y, Z ) ; 
   comm231sd( X, Y, Z ) ;
   comm413_233sd( X, Y, Z ) ; 
   comm433_pp_hh_431sd( X, Y, Z ) ; 
   comm433sd_ph( X, Y, Z ) ; 

   X.profiler.timer["CommutatorScalarDagger"] += omp_get_wtime() - t_css;
//   std::cout << __func__ <<  " done." << std::endl;
   return Z;

}







//*****************************************************************************************
/// BCH_Transform(X,Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator BCH_Transform(  const Operator& OpIn, const Operator& Omega)
{
   return use_brueckner_bch ? Brueckner_BCH_Transform( OpIn, Omega ) :  Standard_BCH_Transform( OpIn, Omega );
}

/// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator Standard_BCH_Transform( const Operator& OpIn, const Operator &Omega)
{
//   std::cout << "!!! " << __func__ << " !!!   particles ranks are " << OpIn.GetParticleRank() << "  and  " << Omega.GetParticleRank() 
//             << "  PN mode is " << OpIn.ThreeBody.PN_mode << "   and  " << Omega.ThreeBody.PN_mode  << std::endl;
   double t_start = omp_get_wtime();
   int max_iter = 40;
   int warn_iter = 12;
   double nx = OpIn.Norm();
   double ny = Omega.Norm();
   Operator OpOut = OpIn;
//   if (use_imsrg3 and not OpOut.ThreeBody.is_allocated )
   if (use_imsrg3 and not OpOut.ThreeBody.IsAllocated() )
   {
     OpOut.ThreeBody.Allocate();
//     OpOut.ThreeBody.Allocate_PN();
   }
   double factorial_denom = 1.0;
   Operator goosetank_chi;  // auxiliary one-body operator used to recover 4th-order quadruples.
   if (use_goose_tank_correction)
   {
     goosetank_chi = OpIn;
     goosetank_chi.SetParticleRank(1);
     goosetank_chi.Erase();
   }
   if (nx>bch_transform_threshold)
   {
     Operator OpNested = OpIn;
     double epsilon = nx * exp(-2*ny) * bch_transform_threshold / (2*ny);
     for (int i=1; i<=max_iter; ++i)
     {

        if (use_goose_tank_correction  )
        {
          auto chi_last = goosetank_chi.OneBody;
          goosetank_chi = GooseTankUpdate( Omega, OpNested);
          OpNested.OneBody += chi_last;  // add the chi from the previous step to OpNested.
        }
        

        OpNested = Commutator(Omega,OpNested); // the ith nested commutator
        factorial_denom /= i;
        OpOut += factorial_denom * OpNested;
  
        if (OpOut.rank_J > 0)
        {
            std::cout << "Tensor BCH, i=" << i << "  Norm = " << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.OneBodyNorm() << " " 
                                                              << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.TwoBodyNorm() << " "
                                                              << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.Norm() << std::endl;
        }
        epsilon *= i+1;
        if (OpNested.Norm() < epsilon)  break;
        if (i == warn_iter)  std::cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << std::endl;
        else if (i == max_iter)   std::cout << "Warning: BCH_Transform didn't coverge after "<< max_iter << " nested commutators" << std::endl;
     }
   }
//   std::cout << "Done with BCH_Transform, 3-body norm of OpOut = " << OpOut.ThreeBodyNorm() << std::endl;
   OpIn.profiler.timer["BCH_Transform"] += omp_get_wtime() - t_start;
   return OpOut;
}


//  Update the auxiliary one-body operator chi, using Omega and the ith nested commutator
//  This has not been tested for tensor commutators, but it *should* work.
// 
Operator GooseTankUpdate( const Operator& Omega, const Operator& OpNested)
{
   double t_start = omp_get_wtime();
   Operator goosetank_chi = Operator( *(OpNested.modelspace), OpNested.rank_J, OpNested.rank_T, OpNested.parity, 1) ;
   goosetank_chi.EraseOneBody();
   if (goosetank_chi.rank_J==0 )
   {
     comm221ss( Omega, OpNested, goosetank_chi );  // update chi.
   }
   else
   {
     comm222_pp_hh_221st( Omega, OpNested, goosetank_chi );  // update chi.
   }
   goosetank_chi.Symmetrize(); // the commutator call only does half the matrix, so we symmetrize
//   int norbits = OpNested.modelspace->GetNumberOrbits();
//   for (int i=0;i<norbits;++i)  // enforce n_in_j + nbar_i nbar_j
   for (auto i : OpNested.modelspace->all_orbits )  // enforce n_in_j + nbar_i nbar_j
   {
     Orbit &oi = OpNested.modelspace->GetOrbit(i);
//     for (int j=0;j<norbits;++j)
     for ( auto j : OpNested.modelspace->all_orbits )
     {
      Orbit &oj = OpNested.modelspace->GetOrbit(j);
      goosetank_chi.OneBody(i,j) *=  oi.occ*oj.occ + (1.0-oi.occ)*(1.0-oj.occ) ;
      }
   }
   OpNested.profiler.timer["GooseTankUpdate"] += omp_get_wtime() - t_start;
   return goosetank_chi;
}




/// Variation of the BCH transformation procedure
/// requested by a one Dr. T.D. Morris
/// \f[ e^{\Omega_1 + \Omega_2} X e^{-\Omega_1 - \Omega_2}
///    \rightarrow 
///  e^{\Omega_2} e^{\Omega_1}  X e^{-\Omega_1} e^{-\Omega_2} \f]
Operator Brueckner_BCH_Transform( const Operator& OpIn, const Operator& Omega)
{
   Operator Omega1 = Omega;
   Operator Omega2 = Omega;
   Omega1.SetParticleRank(1);
   Omega1.EraseTwoBody();
   Omega2.EraseOneBody();
   Operator OpOut = Standard_BCH_Transform(OpIn, Omega1);
   OpOut = Standard_BCH_Transform(OpOut, Omega2);
   return OpOut;
}


//*****************************************************************************************
// Baker-Campbell-Hausdorff formula
//  returns Z, where
//  exp(Z) = exp(X) * exp(Y).
//  Z = X + Y + 1/2[X, Y]
//     + 1/12 [X,[X,Y]] + 1/12 [Y,[Y,X]]
//     - 1/24 [Y,[X,[X,Y]]]
//     - 1/720 [Y,[Y,[Y,[Y,X]]]] - 1/720 [X,[X,[X,[X,Y]]]]
//     + ...
//*****************************************************************************************
/// X.BCH_Product(Y) returns \f$Z\f$ such that \f$ e^{Z} = e^{X}e^{Y}\f$
/// by employing the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + Y + \frac{1}{2}[X,Y] + \frac{1}{12}([X,[X,Y]]+[Y,[Y,X]]) + \ldots \f]
//*****************************************************************************************
Operator BCH_Product(  Operator& X, Operator& Y)
{
//   std::cout << "!!! " << __func__ << " !!! " << std::endl;
   double tstart = omp_get_wtime();
   double nx = X.Norm();
   std::vector<double> bernoulli = {1.0, -0.5, 1./6, 0.0, -1./30,  0.0 ,  1./42,     0,  -1./30};
   std::vector<double> factorial = {1.0,  1.0,  2.0, 6.0,    24.,  120.,   720., 5040.,  40320.};


   Operator Z = X + Y;

   // if we set the option only_2b_omega
   // then temporarily switch to imsrg2 for updating Omega
   bool _save_imsrg3 = use_imsrg3;
   if (only_2b_omega)
   {
     Z.ThreeBody.Erase();
     use_imsrg3 = false;
   }

   Operator Nested = Commutator(Y,X);  // [Y,X]


   double nxy = Nested.Norm();
   // We assume X is small, but just in case, we check if we should include the [X,[X,Y]] term.
   if ( nxy*nx > bch_product_threshold)
   {
     Z += (1./12) * Commutator(Nested,X);
   }
   
   int k = 1;
   // k=1 adds 1/2[X,Y],  k=2 adds 1/12 [Y,[Y,X]], k=4 adds -1/720 [Y,[Y,[Y,[Y,X]]]], and so on.
   while( Nested.Norm() > bch_product_threshold and k<9)
   {
     if ((k<2) or (k%2==0))
        Z += (bernoulli[k]/factorial[k]) * Nested;

     Nested = Commutator(Y,Nested);
     k++;
   }

   use_imsrg3 = _save_imsrg3;// set it back to how it was.

   X.profiler.timer["BCH_Product"] += omp_get_wtime() - tstart;
   return Z;
}







double EstimateBCHError( Operator& Omega, Operator H)
{

  double normOmega = Omega.Norm();
  double normOmega2 = Omega.TwoBodyNorm();
  double normH2 = H.TwoBodyNorm();
  Operator H223 = H;
  H223.Erase();
  comm222_phss( Omega, H, H223 ) ;
  comm222_pp_hh_221ss( Omega, H, H223 ) ;
  if ( H223.GetParticleRank()<3) H223.SetParticleRank(3);
//  if ( not H223.ThreeBody.IsAllocated() )  H223.ThreeBody.SwitchToPN_and_discard();
  if ( not H223.ThreeBody.IsAllocated() )  H223.ThreeBody.Allocate();
  
  comm223ss(Omega, H, H223 );
  comm220ss(Omega,H, H223);

  double Norm220 = H223.ZeroBody;
  

  double Norm223 = H223.ThreeBodyNorm();
  Operator H2223 = H223;
  H2223.Erase();
  comm223ss(Omega, H223,H2223);
  double Norm2223 = H2223.ThreeBodyNorm();

  // now 4 nested commutators
  H223.Erase();
  H2223.Erase();
  H223 = Commutator(Omega,H); // 1 nested
  H2223 = Commutator(Omega,H223); // 2 nested
  H223 = Commutator(Omega,H2223); // 3 nested
  H2223 = Commutator(Omega,H223);
  double Norm4nested = H2223.Norm();

  double est_err = 2./3*normOmega*normOmega * Norm223  +  1./6 * normOmega*normOmega * Norm2223 + exp(2.*normOmega)/24 * Norm4nested;
  std::cout << "Contributions to err " << 2./3*normOmega*normOmega * Norm223 << "  " << 1./6 * normOmega*normOmega * Norm2223
            << "  " << exp(2.*normOmega)/24 * Norm4nested << "    220 = " << Norm220 << "  ||H2|| ,||Omega2|| = " << normH2 << " " << normOmega2 << std::endl;

  return est_err;

}









///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// Below is the implementation of the commutators in the various channels
///////////////////////////////////////////////////////////////////////////////////////////

//*****************************************************************************************
//                ____Y    __         ____X
//          X ___(_)             Y___(_) 
//
//  [X1,Y1](0) = Sum_ab (2j_a+1) x_ab y_ba  (n_a-n_b) 
//             = Sum_a  (2j_a+1)  (xy-yx)_aa n_a
//
// -- AGREES WITH NATHAN'S RESULTS
/// \f[
///  [X_{1)},Y_{(1)}]_{(0)} = \sum_{a} n_a (2j_a+1) \left(X_{(1)}Y_{(1)}-Y_{(1)}X_{(1)}\right)_{aa}
/// \f]
//void Operator::comm110ss( const Operator& X, const Operator& Y) 
void comm110ss( const Operator& X, const Operator& Y, Operator& Z) 
{
  if (X.IsHermitian() and Y.IsHermitian()) return ; // I think this is the case
  if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return ; // I think this is the case
  if (Z.GetJRank()>0 or Z.GetTRank()>0 or Z.GetParity()!=0) return;

   arma::mat xyyx = X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
   for ( auto& a : Z.modelspace->holes) 
   {
      Orbit& oa = Z.modelspace->GetOrbit(a);
      Z.ZeroBody += (oa.j2+1) * oa.occ * xyyx(a,a);
   }
}


//*****************************************************************************************
//         __Y__       __X__
//        ()_ _()  -  ()_ _()
//           X           Y
//
//  [ X^(2), Y^(2) ]^(0) = 1/2 Sum_abcd  Sum_J (2J+1) x_abcd y_cdab (n_a n_b nbar_c nbar_d)
//                       = 1/2 Sum_J (2J+1) Sum_abcd x_abcd y_cdab (n_a n_b nbar_c nbar_d)  
//                       = 1/2 Sum_J (2J+1) Sum_ab  (X*P_pp*Y)_abab  P_hh
//
//  -- AGREES WITH NATHAN'S RESULTS (within < 1%)
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{abcd} (n_a n_b \bar{n}_c \bar{n}_d) \tilde{X}_{abcd}^{J} \tilde{Y}_{cdab}^{J}
/// \f]
/// may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = 2 \sum_{J} (2J+1) Tr(X_{hh'pp'}^{J} Y_{pp'hh'}^{J})
/// \f] where we obtain a factor of four from converting two unrestricted sums to restricted sums, i.e. \f$\sum_{ab} \rightarrow \sum_{a\leq b} \f$,
/// and using the normalized TBME.
//void Operator::comm220ss( const Operator& X, const Operator& Y) 
void comm220ss( const Operator& X, const Operator& Y, Operator& Z) 
{
   if (X.IsHermitian() and Y.IsHermitian()) return; // I think this is the case
   if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return; // I think this is the case
   if (Z.GetJRank()>0 or Z.GetTRank()>0 or Z.GetParity()!=0) return;

   for (int ch=0;ch<Y.nChannels;++ch)
   {
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      auto hh = tbc.GetKetIndex_hh();
      auto ph = tbc.GetKetIndex_ph();
      auto pp = tbc.GetKetIndex_pp();
      arma::uvec nbar_indices = arma::join_cols(hh,ph);
      nbar_indices = arma::join_cols(nbar_indices,pp);
      if (hh.size()==0 ) continue;
      auto nn = tbc.Ket_occ_hh;
      arma::vec nbarnbar = arma::join_cols(tbc.Ket_unocc_hh, tbc.Ket_unocc_ph);
      auto & X2 = X.TwoBody.GetMatrix(ch).submat(hh,nbar_indices);
      arma::mat Y2 = Y.TwoBody.GetMatrix(ch).submat(nbar_indices,hh);
      Y2.head_rows(nbarnbar.size()).each_col() %= nbarnbar;
      Z.ZeroBody += 2 * (2*tbc.J+1) * arma::sum(arma::diagvec( X2 * Y2 ) % nn); // This could be made more efficient, but who cares?
   }
}

//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// -- AGREES WITH NATHAN'S RESULTS
/// \f[
/// [X_{(1)},Y_{(1)}]_{(1)} = X_{(1)}Y_{(1)} - Y_{(1)}X_{(1)}
/// \f]
//void Operator::comm111ss( Operator & Y, Operator& Z) 
//void Operator::comm111ss( const Operator & X, const Operator& Y) 
void comm111ss( const Operator & X, const Operator& Y, Operator& Z) 
{
   Z.OneBody += X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
}

//*****************************************************************************************
//                                       |
//      i |              i |             |
//        |    ___.Y       |__X__        |
//        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)  
//                                       |      
//---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
//                     (note: I think this should actually be)
//                                                = sum_ab (n_a nbar_b) sum_J (2J+1)/(2j_i+1)
//                                                      * y_ab xbiaj - yba x_aibj
//
// -- AGREES WITH NATHAN'S RESULTS 
/// Returns \f$ [X_{(1)},Y_{(2)}] - [Y_{(1)},X_{(2)}] \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1}\sum_{ab} (n_a \bar{n}_b) \sum_{J} (2J+1) (X_{ab} Y^J_{biaj} - X_{ba} Y^J_{aibj})
/// \f]
//void Operator::comm121ss( const Operator& X, const Operator& Y) 
void comm121ss( const Operator& X, const Operator& Y, Operator& Z) 
{
//   index_t norbits = Z.modelspace->GetNumberOrbits();
   index_t norbits = Z.modelspace->all_orbits.size();
  
   #pragma omp parallel for 
   for (index_t indexi=0;indexi<norbits;++indexi)
//   for (index_t i=0;i<norbits;++i)
   {
//      auto i = Z.modelspace->all_orbits[indexi];
      auto i = indexi;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      index_t jmin = Z.IsNonHermitian() ? 0 : i;
      for (auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : Z.modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = Z.modelspace->GetOrbit(a);
//             for (index_t b=0; b<norbits; ++b)
//             for (auto b : Z.modelspace->all_orbits)
             if (Y.particle_rank>1)
             {
               for (auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2} ))
               {
                  Orbit &ob = Z.modelspace->GetOrbit(b);
                  double nanb = oa.occ * (1-ob.occ);
                  if (std::abs(nanb)<ModelSpace::OCC_CUT) continue;
                  Z.OneBody(i,j) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,j) ;
                }
             }
//             for (index_t b=0; b<norbits; ++b)
//             for (auto b : Z.modelspace->all_orbits)
             if (X.particle_rank>1)
             {
                for (auto b : Y.OneBodyChannels.at({oa.l,oa.j2,oa.tz2} ))
                {
                  // comm211 part
                  Orbit &ob = Z.modelspace->GetOrbit(b);
                  double nanb = oa.occ * (1-ob.occ);
                  if (std::abs(nanb)<ModelSpace::OCC_CUT) continue;
                  Z.OneBody(i,j) -= (ob.j2+1) * nanb * Y.OneBody(a,b) * X.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) += (oa.j2+1) * nanb * Y.OneBody(b,a) * X.TwoBody.GetTBMEmonopole(a,i,b,j) ;
                }
             }
          }
      }
   }
}



//*****************************************************************************************
//
//      i |              i |            [X2,Y2](1)  =  1/(2(2j_i+1)) sum_J (2J+1) 
//        |__Y__           |__X__           * sum_abc (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
//        |    /\          |    /\          * (x_ciab y_abcj - y_ciab xabcj)
//        |   (  )   _     |   (  )                                                                                      
//        |____\/          |____\/       = 1/(2(2j+1)) sum_J (2J+1)
//      j | X            j |  Y            *  sum_c ( Pp*X*Phh*Y*Pp - Pp*Y*Phh*X*Pp)  - (Ph*X*Ppp*Y*Ph - Ph*Y*Ppp*X*Ph)_cicj
//                                     
//
// -- AGREES WITH NATHAN'S RESULTS 
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
// \f[
// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{c}
// \left( \mathcal{P}_{pp} (X \mathcal{P}_{hh} Y^{J} 
// - Y^{J} \mathcal{P}_{hh} X^{J}) \mathcal{P}_{pp}
//  - \mathcal{P}_{hh} (X^{J} \mathcal{P}_{pp} Y^{J} 
//  -  Y^{J} \mathcal{P}_{pp} X^{J}) \mathcal{P}_{hh} \right)_{cicj}
// \f]
/// \f[
/// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{abc} (\bar{n}_a\bar{n}_bn_c + n_an_b\bar{n}_c)
///  (X^{J}_{ciab} Y^{J}_{abcj} - Y^{J}_{ciab}X^{J}_{abcj})
/// \f]
/// This may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1} \sum_{c} \sum_{J} (2J+1) \left( n_c \mathcal{M}^{J}_{pp,icjc} + \bar{n}_c\mathcal{M}^{J}_{hh,icjc} \right)
/// \f]
/// With the intermediate matrix \f[ \mathcal{M}^{J}_{pp} \equiv \frac{1}{2} (X^{J}\mathcal{P}_{pp} Y^{J} - Y^{J}\mathcal{P}_{pp}X^{J}) \f]
/// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$
//void Operator::comm221ss( const Operator& X, const Operator& Y) 
void comm221ss( const Operator& X, const Operator& Y, Operator& Z) 
{

   double t_start = omp_get_wtime();

//   static TwoBodyME Mpp = Y.TwoBody; // SRS: Is there a good reason to make these static?
//   static TwoBodyME Mhh = Y.TwoBody;

   TwoBodyME Mpp = Z.TwoBody;
   TwoBodyME Mhh = Z.TwoBody;
   Mpp.Erase();
   Mhh.Erase();
   ConstructScalarMpp_Mhh( X, Y, Z, Mpp, Mhh);


//   int norbits = Z.modelspace->GetNumberOrbits();
   int norbits = Z.modelspace->all_orbits.size();
   std::vector<index_t> allorb_vec( Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
   #pragma omp parallel for schedule(dynamic,1)
//   for (int i=0;i<norbits;++i)
   for (int indexi=0;indexi<norbits;++indexi)
   {
//      auto i = Z.modelspace->all_orbits[indexi];
      auto i = allorb_vec[indexi];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         if (j<jmin) continue;
         Orbit &oj = Z.modelspace->GetOrbit(j);
         double cijJ = 0;
         for ( auto& c : Z.modelspace->all_orbits)
         {
           Orbit& oc = Z.modelspace->GetOrbit(c);
           double nc = oc.occ;
           double nbarc = 1.0 - nc;
           int Jmin = std::max( std::abs( oc.j2-oi.j2), std::abs( oc.j2-oj.j2))/2;
           int Jmax = ( oc.j2 + std::min(oi.j2,oj.j2) ) /2;

           if ( std::abs(nc)>1e-9 )
           {
             for (int J=Jmin; J<=Jmax; J++ )
             {
               cijJ += (2*J+1) * nc   * Mpp.GetTBME_J(J,c,i,c,j) ;
             }
           }
           if ( std::abs(nbarc)>1e-9 )
           {
             for (int J=Jmin; J<=Jmax; J++ )
             {
               cijJ += (2*J+1) * nbarc  * Mhh.GetTBME_J(J,c,i,c,j);
             }
           }
         }


         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   }

   X.profiler.timer["comm221ss"] += omp_get_wtime() - t_start;

}





//*****************************************************************************************
//
//    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
//    |     |___.Y          |__X___|         
//    |     |         _     |      |          
//    |_____|               |      |_____.Y        
//    |  X  |               |      |            
//
// -- AGREES WITH NATHAN'S RESULTS
/// Returns \f$ [X_{(1)},Y_{(2)}]_{(2)} - [Y_{(1)},X_{(2)}]_{(2)} \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]^{J}_{ijkl} = \sum_{a} ( X_{ia}Y^{J}_{ajkl} + X_{ja}Y^{J}_{iakl} - X_{ak} Y^{J}_{ijal} - X_{al} Y^{J}_{ijka} )
/// \f]
/// here, all TBME are unnormalized, i.e. they should have a tilde.
// This is still too slow...
//void Operator::comm122ss( Operator& Y, Operator& Z ) 
//void Operator::comm122ss( const Operator& X, const Operator& Y ) 
// TODO: MODIFY THIS TO ACCOMMODATE J=0, T!=0 operators
void comm122ss( const Operator& X, const Operator& Y, Operator& Z ) 
{
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
   int hZ = Z.IsHermitian() ? 1 : -1;

   int n_nonzero = Z.modelspace->SortedTwoBodyChannels.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      auto& X2 = X.TwoBody.GetMatrix(ch,ch);
      auto& Y2 = Y.TwoBody.GetMatrix(ch,ch);
      auto& Z2 = Z.TwoBody.GetMatrix(ch,ch);
      arma::mat W2(size(Z2),arma::fill::zeros); // temporary intermediate matrix

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         int flipphaseij = - Z.modelspace->phase((oi.j2+oj.j2)/2-tbc.J);

         // make lists of the indices we want, then do matrix multiplication.
         // there may be a more efficient way to find these
         std::vector<index_t> ind1_ia, ind1_ja,ind2_aj,ind2_ai;
         std::vector<double> factor_ia,factor_ja;
         for (int a : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
         {
            size_t ind2 = tbc.GetLocalIndex( std::min(a,j), std::max(a,j) );
            if (ind2<0 or ind2>=tbc.GetNumberKets()) continue;
            ind1_ia.push_back(a);
            ind2_aj.push_back(ind2);
            factor_ia.push_back( a>j ? flipphaseij : (a==j ? PhysConst::SQRT2 : 1));
         }
         if (i!=j)
         {
           for (int a : Z.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
           {
              size_t ind2 = tbc.GetLocalIndex( std::min(a,i), std::max(a,i) );
              if (ind2<0 or ind2>=tbc.GetNumberKets()) continue;
              ind1_ja.push_back(a);
              ind2_ai.push_back(ind2);
              factor_ja.push_back( i>a ? flipphaseij : (i==a ? PhysConst::SQRT2 : 1));
           }

         }

         arma::uvec u_ind1_ia(ind1_ia);
         arma::uvec u_ind1_ja(ind1_ja);
         arma::uvec u_ind2_aj(ind2_aj);
         arma::uvec u_ind2_ai(ind2_ai);
         arma::vec  v_factor_ia(factor_ia);
         arma::vec  v_factor_ja(factor_ja);
         if (i==j)
         {
           v_factor_ia /= PhysConst::SQRT2;
           v_factor_ja /= PhysConst::SQRT2;
         }

         // This is fairly obfuscated, but hopefully faster for bigger calculations
         // SRS THERE MUST BE A BETTER WAY...
         if (X.particle_rank>1 and Y.particle_rank>1)
         {
            W2.col(indx_ij) = join_horiz(    Y2.cols(join_vert( u_ind2_aj,u_ind2_ai))  , X2.cols(join_vert(u_ind2_aj,u_ind2_ai) ) )
                             * join_vert(  join_vert( X1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja ),
                                          -join_vert( Y1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja ));
         }
         else if (X.particle_rank<2 and Y.particle_rank>1)
         {
            W2.col(indx_ij) =     Y2.cols(join_vert( u_ind2_aj,u_ind2_ai))    
                             *   join_vert( X1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja );
         }
         else if (X.particle_rank>1 and Y.particle_rank<2)
         {
            W2.col(indx_ij) =      -X2.cols(join_vert(u_ind2_aj,u_ind2_ai) ) 
                                      *    join_vert( Y1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja );
         }

      if (i==j) W2.col(indx_ij) *= 2;

      }
      Z2 -= W2 + hZ*W2.t();
   }

}





//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_pp(hh) = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |   
//
// -- AGREES WITH NATHAN'S RESULTS
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
/// Calculates the part of the commutator \f$ [X_{(2)},Y_{(2)}]_{(2)} \f$ which involves particle-particle
/// or hole-hole intermediate states.
/// \f[
/// [X_{(2)},Y_{(2)}]^{J}_{ijkl} = \frac{1}{2} \sum_{ab} (\bar{n}_a\bar{n}_b - n_an_b) (X^{J}_{ijab}Y^{J}_{ablk} - Y^{J}_{ijab}X^{J}_{abkl})
/// \f]
/// This may be written as
/// \f[
/// [X_{(2)},Y_{(2)}]^{J} = \mathcal{M}^{J}_{pp} - \mathcal{M}^{J}_{hh}
/// \f]
/// With the intermediate matrices
/// \f[
/// \mathcal{M}^{J}_{pp} \equiv \frac{1}{2}(X^{J} \mathcal{P}_{pp} Y^{J} - Y^{J} \mathcal{P}_{pp} X^{J})
/// \f]
/// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$.
//void Operator::comm222_pp_hhss( Operator& opright, Operator& opout ) 
//void Operator::comm222_pp_hhss( const Operator& X, const Operator& Y ) 
void comm222_pp_hhss( const Operator& X, const Operator& Y, Operator& Z ) 
{
   double t_start = omp_get_wtime();

   TwoBodyME Mpp = Z.TwoBody;
   TwoBodyME Mhh = Z.TwoBody;
   Mpp.Erase();
   Mhh.Erase();

   ConstructScalarMpp_Mhh( X, Y, Z, Mpp, Mhh);

   Z.TwoBody += Mpp;
   Z.TwoBody -= Mhh;
   X.profiler.timer[__func__] += omp_get_wtime() - t_start;
}




void ConstructScalarMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh)
{

   double t_start = omp_get_wtime();

   std::vector<size_t> ch_bra_list,ch_ket_list;
   for ( auto& iter : Z.TwoBody.MatEl )
   {
      ch_bra_list.push_back( iter.first[0] );
      ch_ket_list.push_back( iter.first[1] );
   }
   int nch = ch_bra_list.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch_bra = ch_bra_list[ich];
      int ch_ket = ch_ket_list[ich];
      TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);

      auto& X2_bra = X.TwoBody.GetMatrix(ch_bra,ch_bra);
      auto& X2_ket = X.TwoBody.GetMatrix(ch_ket,ch_ket);
      auto& Y2 = Y.TwoBody.GetMatrix(ch_bra,ch_ket);

      auto& Matrixpp = Mpp.GetMatrix(ch_bra,ch_ket);
      auto& Matrixhh = Mhh.GetMatrix(ch_bra,ch_ket);

      auto& bras_pp = tbc_bra.GetKetIndex_pp();
      auto& bras_hh = tbc_bra.GetKetIndex_hh();
      auto& bras_ph = tbc_bra.GetKetIndex_ph();
      auto& kets_pp = tbc_ket.GetKetIndex_pp();
      auto& kets_hh = tbc_ket.GetKetIndex_hh();
      auto& kets_ph = tbc_ket.GetKetIndex_ph();
      auto& nanb_bra = tbc_bra.Ket_occ_hh;
      auto& nanb_ket = tbc_ket.Ket_occ_hh;
      auto& nbarnbar_hh_bra = tbc_bra.Ket_unocc_hh;
      auto& nbarnbar_ph_bra = tbc_bra.Ket_unocc_ph;
      auto& nbarnbar_hh_ket = tbc_ket.Ket_unocc_hh;
      auto& nbarnbar_ph_ket = tbc_ket.Ket_unocc_ph;
      
      Matrixpp =  X2_bra.cols(bras_pp) * Y2.rows(bras_pp);
      Matrixhh =  X2_bra.cols(bras_hh) * arma::diagmat(nanb_bra) *  Y2.rows(bras_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  X2_bra.cols(bras_hh) * arma::diagmat(nbarnbar_hh_bra) *  Y2.rows(bras_hh); 
      if (kets_ph.size()>0)
        Matrixpp += X2_bra.cols(bras_ph) * arma::diagmat(nbarnbar_ph_bra) *  Y2.rows(bras_ph) ;


      if (Z.IsHermitian() and ch_bra==ch_ket)
        {
           Matrixpp +=  Matrixpp.t();
           Matrixhh +=  Matrixhh.t();
        }
        else if (Z.IsAntiHermitian() and ch_bra==ch_ket) // i.e. LHS and RHS are both hermitian or ant-hermitian
        {
           Matrixpp -=  Matrixpp.t();
           Matrixhh -=  Matrixhh.t();
        }
        else
        {
          Matrixpp -=  Y2.cols(kets_pp) * X2_ket.rows(kets_pp);
          Matrixhh -=  Y2.cols(kets_hh) * arma::diagmat(nanb_ket) *  X2_ket.rows(kets_hh) ;
          Matrixpp -=  Y2.cols(kets_hh) * arma::diagmat(nbarnbar_hh_ket) *  X2_ket.rows(kets_hh) ;
          if (kets_ph.size()>0)
            Matrixpp -=  Y2.cols(kets_ph) * arma::diagmat(nbarnbar_ph_ket) *  X2_ket.rows(kets_ph) ;
        }

   } //for ch

   X.profiler.timer[__func__] += omp_get_wtime() - t_start;
}




/// Since comm222_pp_hhss() and comm221ss() both require the construction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
void comm222_pp_hh_221ss( const Operator& X, const Operator& Y, Operator& Z )  
{

   double t_start = omp_get_wtime();
//   int herm = Z.IsHermitian() ? 1 : -1;
//   Operator& Z = *this;

//   static TwoBodyME Mpp = Z.TwoBody;
//   static TwoBodyME Mhh = Z.TwoBody;
   TwoBodyME Mpp = Z.TwoBody;
   TwoBodyME Mhh = Z.TwoBody;
   Mpp.Erase();
   Mhh.Erase();

   ConstructScalarMpp_Mhh( X, Y, Z, Mpp, Mhh);

//   Z.TwoBody += (Mpp - Mhh);
   Z.TwoBody += Mpp;
   Z.TwoBody -= Mhh;
//   OUT += Matrixpp - Matrixhh;
   X.profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t_start;

   t_start = omp_get_wtime();
//   int norbits = Z.modelspace->GetNumberOrbits();
   int norbits = Z.modelspace->all_orbits.size();
   // The one body part
   std::vector<index_t> allorb_vec( Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
   #pragma omp parallel for schedule(dynamic,1)
//   for (int i=0;i<norbits;++i)
   for (int indexi=0;indexi<norbits;++indexi)
   {
//      auto i = Z.modelspace->all_orbits[indexi];
      auto i = allorb_vec[indexi];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         if (j<jmin) continue;

         Orbit &oj = Z.modelspace->GetOrbit(j);
         double cijJ = 0;
         for ( auto& c : Z.modelspace->all_orbits)
         {
           Orbit& oc = Z.modelspace->GetOrbit(c);
           double nc = oc.occ;
           double nbarc = 1.0 - nc;
           int Jmin = std::max( std::abs( oc.j2-oi.j2), std::abs( oc.j2-oj.j2))/2;
           int Jmax = ( oc.j2 + std::min(oi.j2,oj.j2) ) /2;

           if ( std::abs(nc)>1e-9 )
           {
             for (int J=Jmin; J<=Jmax; J++ )
             {
               cijJ += (2*J+1) * nc   * Mpp.GetTBME_J(J,c,i,c,j) ;
             }
           }
           if ( std::abs(nbarc)>1e-9 )
           {
             for (int J=Jmin; J<=Jmax; J++ )
             {
               cijJ += (2*J+1) * nbarc  * Mhh.GetTBME_J(J,c,i,c,j);
             }
           }
         }

         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   } // for i
   X.profiler.timer["pphh One Body bit"] += omp_get_wtime() - t_start;
}



//**************************************************************************
//
//  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
//                        { k l J'}
// SCALAR VARIETY
/// The scalar Pandya transformation is defined as
/// \f[
///  \bar{X}^{J}_{i\bar{j}k\bar{l}} = - \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  X^{J}_{ilkj}
/// \f]
/// where the overbar indicates time-reversed orbits.
//void Operator::DoPandyaTransformation_SingleChannel(arma::mat& TwoBody_CC_ph, int ch_cc, std::string orientation="normal") const
void DoPandyaTransformation_SingleChannel(const Operator& Z, arma::mat& TwoBody_CC_ph, int ch_cc, std::string orientation="normal")
{
   int herm = Z.IsHermitian() ? 1 : -1;
   TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
   int nKets_cc = tbc_cc.GetNumberKets();
   arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
   int nph_kets = kets_ph.n_rows;
   int J_cc = tbc_cc.J;

   if (orientation=="normal") TwoBody_CC_ph.zeros( 2*nph_kets, nKets_cc);
   else if (orientation=="transpose") TwoBody_CC_ph.zeros( nKets_cc, 2*nph_kets);
   else
   {
     std::cout << __PRETTY_FUNCTION__ << " =>  Unknown orientation input  " << orientation << ". Don't know what to do with this." << std::endl;
     return;
   }

   // loop over cross-coupled ph bras <ab| in this channel
   // (this is the side that gets summed over in the matrix multiplication)
   for (int ibra=0; ibra<nph_kets; ++ibra)
   {
      Ket & bra_cc = tbc_cc.GetKet( kets_ph[ibra] );
      // we want to evaluate a<=b and a>=b, so to avoid code duplication, we turn this into a loop over the two orderings
      std::vector<size_t> ab_switcheroo = { bra_cc.p, bra_cc.q };
      for ( int ab_case=0; ab_case<=1; ab_case++)
      {
        int a = ab_switcheroo[ab_case];   // this little bit gives us a,b if ab_case=0 and b,a if ab_case=1
        int b = ab_switcheroo[1-ab_case];
        size_t bra_shift = ab_case*nph_kets;  // if we switch a<->b, we offset the bra index by nph_kets

        Orbit & oa = Z.modelspace->GetOrbit(a);
        Orbit & ob = Z.modelspace->GetOrbit(b);
        double ja = oa.j2*0.5;
        double jb = ob.j2*0.5;
        double na_nb_factor = oa.occ - ob.occ;

        // loop over cross-coupled kets |cd> in this channel
        for (int iket_cc=0; iket_cc<nKets_cc; ++iket_cc)
        {
           Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
           int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
           int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
           Orbit & oc = Z.modelspace->GetOrbit(c);
           Orbit & od = Z.modelspace->GetOrbit(d);
           double jc = oc.j2*0.5;
           double jd = od.j2*0.5;

           int jmin = std::max(std::abs(ja-jd),std::abs(jc-jb));
           int jmax = std::min(ja+jd,jc+jb);
           double Xbar = 0;
           for (int J_std=jmin; J_std<=jmax; ++J_std)
           {
              double sixj = Z.modelspace->GetSixJ(ja,jb,J_cc,jc,jd,J_std);
              if (std::abs(sixj) < 1e-8) continue;
              double tbme = Z.TwoBody.GetTBME_J(J_std,a,d,c,b);
              Xbar -= (2*J_std+1) * sixj * tbme ;
           }
           if (orientation=="normal")
           {
             TwoBody_CC_ph(ibra+bra_shift, iket_cc) = Xbar;
           }
           else // "transpose"
           {
             TwoBody_CC_ph(iket_cc,ibra+bra_shift) = herm * Xbar * na_nb_factor;  // we slap the (na-nb) on the transposed one.
           }

        }
      }
   }
}

void DoPandyaTransformation_SingleChannel_XandY(const Operator& X, const Operator& Y, arma::mat& X2_CC_ph, arma::mat& Y2_CC_ph, int ch_cc)
{
//   int hX = X.IsHermitian() ? 1 : -1;
//   int hY = X.IsHermitian() ? 1 : -1;
   TwoBodyChannel& tbc_cc = X.modelspace->GetTwoBodyChannel_CC(ch_cc);
   int nKets_cc = tbc_cc.GetNumberKets();
   arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
   int nph_kets = kets_ph.n_rows;
   int J_cc = tbc_cc.J;

   X2_CC_ph.zeros( nKets_cc, 2*nph_kets );
   Y2_CC_ph.zeros( 2*nph_kets, nKets_cc);


   // loop over cross-coupled ph bras <ab| in this channel
   // (this is the side that gets summed over in the matrix multiplication)
   for (int ibra=0; ibra<nph_kets; ++ibra)
   {
      Ket & bra_cc = tbc_cc.GetKet( kets_ph[ibra] );
      // we want to evaluate a<=b and a>=b, so to avoid code duplication, we turn this into a loop over the two orderings
      std::vector<size_t> ab_switcheroo = { bra_cc.p, bra_cc.q };
      for ( int ab_case=0; ab_case<=1; ab_case++)
      {
        int a = ab_switcheroo[ab_case];   // this little bit gives us a,b if ab_case=0 and b,a if ab_case=1
        int b = ab_switcheroo[1-ab_case];
        size_t bra_shift = ab_case*nph_kets;  // if we switch a<->b, we offset the bra index by nph_kets

        Orbit & oa = X.modelspace->GetOrbit(a);
        Orbit & ob = X.modelspace->GetOrbit(b);
        double ja = oa.j2*0.5;
        double jb = ob.j2*0.5;
        double na_nb_factor = oa.occ - ob.occ;

        // loop over cross-coupled kets |cd> in this channel
        for (int iket_cc=0; iket_cc<nKets_cc; ++iket_cc)
        {
           Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
           int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
           int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
           Orbit & oc = X.modelspace->GetOrbit(c);
           Orbit & od = X.modelspace->GetOrbit(d);
           double jc = oc.j2*0.5;
           double jd = od.j2*0.5;

           int jmin = std::max(std::abs(ja-jd),std::abs(jc-jb));
           int jmax = std::min(ja+jd,jc+jb);
           double Xbar = 0;
           double Ybar = 0;
           for (int J_std=jmin; J_std<=jmax; ++J_std)
           {
              double sixj = X.modelspace->GetSixJ(ja,jb,J_cc,jc,jd,J_std);
              if (std::abs(sixj) < 1e-8) continue;
//              double tbme = Z.TwoBody.GetTBME_J(J_std,a,d,c,b);
              double xcbad = X.TwoBody.GetTBME_J(J_std,c,b,a,d);
              double yadcb = Y.TwoBody.GetTBME_J(J_std,a,d,c,b);
              Xbar -= (2*J_std+1) * sixj * xcbad ;
              Ybar -= (2*J_std+1) * sixj * yadcb ;
           }
           X2_CC_ph( iket_cc, ibra+bra_shift ) = Xbar * na_nb_factor;
           Y2_CC_ph( ibra+bra_shift, iket_cc ) = Ybar;

        }// for iket_cc
      }// for ab_case
   }// for ibra
}






//void Operator::DoPandyaTransformation(deque<arma::mat>& TwoBody_CC_ph, std::string orientation="normal") const
void DoPandyaTransformation(const Operator& Z, std::deque<arma::mat>& TwoBody_CC_ph, std::string orientation="normal")
{
   // loop over cross-coupled channels
   int n_nonzero = Z.modelspace->SortedTwoBodyChannels_CC.size();
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      DoPandyaTransformation_SingleChannel(Z, TwoBody_CC_ph[ch_cc] , ch_cc, orientation);
   }
}



void AddInversePandyaTransformation_SingleChannel( Operator& Z,  arma::mat& Zbar, int ch_cc)
{
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
   TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
   int Jcc = tbc_cc.J;
   int Jhat2 = 2*Jcc+1;
   int parity_cc = tbc_cc.parity;
   int Tz_cc = tbc_cc.Tz;
   int nkets_cc = tbc_cc.GetNumberKets();
   int n_nonzeroChannels = Z.modelspace->SortedTwoBodyChannels.size();
//   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nKets = tbc.GetNumberKets();
      auto& Zmat = Z.TwoBody.GetMatrix(ch,ch);

      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = Z.modelspace->GetOrbit(i);
         Orbit & oj = Z.modelspace->GetOrbit(j);
         double ji = 0.5*oi.j2;
         double jj = 0.5*oj.j2;
         int ketmin = Z.IsHermitian() ? ibra : ibra+1;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = Z.modelspace->GetOrbit(k);
            Orbit & ol = Z.modelspace->GetOrbit(l);
            double jk = 0.5*ok.j2;
            double jl = 0.5*ol.j2;

            double commij = 0;
            double commji = 0;

            int jmin = std::max(std::abs(int(ji-jl)),std::abs(int(jk-jj)));
            int jmax = std::min(int(ji+jl),int(jk+jj));
//            if ( ((oi.l+ol.l)%2==parity_cc)  and  (std::abs(oi.tz2+ol.tz2)==Tz_cc*2) and (Jcc>=jmin) and (Jcc<=jmax) )
            if ( ((oi.l+ol.l)%2==parity_cc)  and  (std::abs(oi.tz2-ol.tz2)==Tz_cc*2) and (Jcc>=jmin) and (Jcc<=jmax) )
            {
               double sixj = Z.modelspace->GetSixJ(ji,jj,J,jk,jl,Jcc);
               int indx_il = tbc_cc.GetLocalIndex(i,l) ;
               int indx_kj = tbc_cc.GetLocalIndex( std::min(j,k), std::max(j,k) ) +(k>j?nkets_cc:0);
               commij -= Jhat2 * sixj * Zbar(indx_il,indx_kj) ;

            }

            if (k==l)
            {
              commji = commij;
            }
            else if (i==j)
            {
              commji = Z.modelspace->phase(ji+jj+jk+jl) * commij;
            }
            else
            {
              // now loop over the cross coupled TBME's
              jmin = std::max(std::abs(int(jj-jl)),std::abs(int(jk-ji)));
              jmax = std::min(int(jj+jl),int(jk+ji));
//              if ( (oi.l+ok.l)%2==parity_cc  and  std::abs(oi.tz2+ok.tz2)==Tz_cc*2 and Jcc>=jmin and Jcc<=jmax)
              if ( (oi.l+ok.l)%2==parity_cc  and  std::abs(oi.tz2-ok.tz2)==Tz_cc*2 and Jcc>=jmin and Jcc<=jmax)
              {
                 double sixj = Z.modelspace->GetSixJ(jj,ji,J,jk,jl,Jcc);
                 int indx_ik = tbc_cc.GetLocalIndex(i,k) ;
                 int indx_lj = tbc_cc.GetLocalIndex(std::min(l,j),std::max(l,j)) +(l>j?nkets_cc:0);
                 // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
                 // the phase we get from that flip combines with the phase from Pij, to give the phase included below
                 commji -= Jhat2 *  sixj *  Zbar(indx_ik, indx_lj) ;


              }
            }
            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : PhysConst::SQRT2;
            #pragma omp atomic
            Zmat(ibra,iket) -= (commij - Z.modelspace->phase(jk+jl-J ) * commji) / norm;
         }
      }
   }
}



//void Operator::AddInversePandyaTransformation(const deque<arma::mat>& Zbar)
void AddInversePandyaTransformation(const std::deque<arma::mat>& Zbar, Operator& Z)
{
    // Do the inverse Pandya transform
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
   int n_nonzeroChannels = Z.modelspace->SortedTwoBodyChannels.size();

   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nKets = tbc.GetNumberKets();

      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = Z.modelspace->GetOrbit(i);
         Orbit & oj = Z.modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;
         int ketmin = Z.IsHermitian() ? ibra : ibra+1;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = Z.modelspace->GetOrbit(k);
            Orbit & ol = Z.modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;

            double commij = 0;
            double commji = 0;

            int parity_cc = (oi.l+ol.l)%2;
//            int Tz_cc = std::abs(oi.tz2+ol.tz2)/2;
            int Tz_cc = std::abs(oi.tz2-ol.tz2)/2;
            int Jpmin = std::max(std::abs(int(ji-jl)),std::abs(int(jk-jj)));
            int Jpmax = std::min(int(ji+jl),int(jk+jj));
            for (int Jprime=Jpmin; Jprime<=Jpmax; ++Jprime)
            {
               double sixj = Z.modelspace->GetSixJ(ji,jj,J,jk,jl,Jprime);
               if (std::abs(sixj)<1e-8) continue;
               int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
               int nkets_cc = tbc_cc.GetNumberKets();
               int indx_il = tbc_cc.GetLocalIndex(std::min(i,l),std::max(i,l)) +(i>l?nkets_cc:0);
               int indx_kj = tbc_cc.GetLocalIndex(std::min(j,k),std::max(j,k)) +(k>j?nkets_cc:0);
               double me1 = Zbar.at(ch_cc)(indx_il,indx_kj);
               commij -= (2*Jprime+1) * sixj * me1;
            }

            if (k==l)
            {
              commji = commij;
            }
            else if (i==j)
            {
              commji = Z.modelspace->phase(ji+jj+jk+jl) * commij;
            }
            else
            {
              // now loop over the cross coupled TBME's
              parity_cc = (oi.l+ok.l)%2;
              Tz_cc = std::abs(oi.tz2-ok.tz2)/2;
//              Tz_cc = std::abs(oi.tz2+ok.tz2)/2;
              Jpmin = std::max(std::abs(int(jj-jl)),std::abs(int(jk-ji)));
              Jpmax = std::min(int(jj+jl),int(jk+ji));
              for (int Jprime=Jpmin; Jprime<=Jpmax; ++Jprime)
              {
                 double sixj = Z.modelspace->GetSixJ(jj,ji,J,jk,jl,Jprime);
                 if (std::abs(sixj)<1e-8) continue;
                 int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
                 TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
                 int nkets_cc = tbc_cc.GetNumberKets();
                 int indx_ik = tbc_cc.GetLocalIndex(std::min(i,k),std::max(i,k)) +(i>k?nkets_cc:0);
                 int indx_lj = tbc_cc.GetLocalIndex(std::min(l,j),std::max(l,j)) +(l>j?nkets_cc:0);
                 // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
                 // the phase we get from that flip combines with the phase from Pij, to give the phase included below
                 double me1 = Zbar.at(ch_cc)(indx_ik, indx_lj) ;
                 commji -= (2*Jprime+1) *  sixj * me1;

              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : PhysConst::SQRT2;
            Z.TwoBody.GetMatrix(ch,ch)(ibra,iket) -= (commij - Z.modelspace->phase(jk+jl-J ) * commji) / norm;
         }
      }
   }
}

///*************************************
/// convenience function
/// called by comm222_phss
///*************************************
//deque<arma::mat> Operator::InitializePandya(size_t nch, std::string orientation="normal")
std::deque<arma::mat> InitializePandya(Operator& Z, size_t nch, std::string orientation="normal")
{
   std::deque<arma::mat> X(nch);
   int n_nonzero = Z.modelspace->SortedTwoBodyChannels_CC.size();
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();
      if (orientation=="normal")
         X[ch_cc] = arma::mat(2*nph_kets,   nKets_cc, arma::fill::zeros);
      else if (orientation=="transpose")
         X[ch_cc] = arma::mat(nKets_cc, 2*nph_kets,   arma::fill::zeros);
   }
   return X;
}

//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.     
//                                             
//   |          |      |          |           
//   |     __Y__|      |     __X__|            
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |            
//   |____\/    |      |____\/    |            
//   |  X       |      |  Y       |            
//           
//            
// -- This appears to agree with Nathan's results
//
/// Calculates the part of \f$ [X_{(2)},Y_{(2)}]_{ijkl} \f$ which involves ph intermediate states, here indicated by \f$ Z^{J}_{ijkl} \f$
/// \f[
/// Z^{J}_{ijkl} = \sum_{ab}(n_a\bar{n}_b-\bar{n}_an_b)\sum_{J'} (2J'+1)
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
/// \left( \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} - 
///   \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \right)
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_i  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
/// \left( \bar{X}^{J'}_{j\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{i}} - 
///   \bar{Y}^{J'}_{j\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{i}} \right)
/// \right]
/// \f]
/// This is implemented by defining an intermediate matrix
/// \f[
/// \bar{Z}^{J}_{i\bar{l}k\bar{j}} \equiv \sum_{ab}(n_a\bar{n}_b)
/// \left[ \left( \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} - 
///   \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \right)
/// -\left( \bar{X}^{J'}_{i\bar{l}b\bar{a}}\bar{Y}^{J'}_{b\bar{a}k\bar{j}} - 
///    \bar{Y}^{J'}_{i\bar{l}b\bar{a}}\bar{X}^{J'}_{b\bar{a}k\bar{j}} \right)\right]
/// \f]
/// The Pandya-transformed matrix elements are obtained with DoPandyaTransformation().
/// The matrices \f$ \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} \f$
/// and \f$ \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \f$
/// are related by a Hermitian conjugation, which saves two matrix multiplications.
/// The commutator is then given by
/// \f[
/// Z^{J}_{ijkl} = \sum_{J'} (2J'+1)
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  \bar{Z}^{J'}_{i\bar{l}k\bar{j}}
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_i  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  \bar{Z}^{J'}_{j\bar{l}k\bar{i}}
///  \right]
///  \f]
///
//void Operator::comm222_phss( const Operator& X, const Operator& Y ) 
void comm222_phss( const Operator& X, const Operator& Y, Operator& Z ) 
{

   int hy = Y.IsHermitian() ? 1 : -1;
   // Create Pandya-transformed hp and ph matrix elements
   double t_start = omp_get_wtime();

   // Construct the intermediate matrix Z_bar
//   size_t nch = Z.modelspace->SortedTwoBodyChannels_CC.size();
   size_t nch = Z.modelspace->GetNumberTwoBodyChannels_CC();

//   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup(Z.GetJRank(), Z.GetTRank(), Z.GetParity());
//   t_start = omp_get_wtime();
//   std::deque<arma::mat> Z_bar ( Z.nChannels );
   std::deque<arma::mat> Z_bar ( nch );
   for ( size_t ch=0; ch<nch; ch++)
   {
     size_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
     Z_bar[ch].zeros( nKets_cc, 2*nKets_cc );
   }
//////   std::vector<bool> lookup_empty(Z.nChannels,true);
////   std::vector<size_t> non_empty_channels;
//   for (size_t ich=0;ich<nch;++ich)
//   {
//      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC[ich];
//      index_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
//      Z_bar[ch].zeros( nKets_cc, 2*nKets_cc );
////      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
//      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) non_empty_channels.push_back(ch);
//   }
//
//   size_t nch_nonempty = non_empty_channels.size();

   #ifndef OPENBLAS_NOUSEOMP
//   #pragma omp parallel for schedule(dynamic,1)
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   #endif
//   for (size_t ich=0; ich<nch; ++ich )
//   for (size_t ich=0; ich<nch_nonempty; ++ich )
   for (size_t ch=0; ch<nch; ++ch )
   {
//      if (lookup_empty.at(ich)) continue;
//      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC.at(ich);
//      size_t ch = non_empty_channels[ich];
      const TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      size_t nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

      arma::mat Y_bar_ph;
      arma::mat Xt_bar_ph;

      // Y has dimension (2*nph , nKets_CC)
      // X has dimension (nKets_CC, 2*nph )
//      DoPandyaTransformation_SingleChannel(Y,Y_bar_ph,ch,"normal");
//      DoPandyaTransformation_SingleChannel(X,Xt_bar_ph,ch,"transpose");


      // This function gives  <ij`| X |ab`> *na(1-nb)  and   <ab`| Y |ij`>,  for both orderings ab` and ba`
      DoPandyaTransformation_SingleChannel_XandY( X,  Y, Xt_bar_ph, Y_bar_ph, ch);

      auto& Zbar_ch = Z_bar.at(ch);


      if (Y_bar_ph.size()<1 or Xt_bar_ph.size()<1)
      {
//        Zbar_ch = arma::zeros( Xt_bar_ph.n_rows, Y_bar_ph.n_cols*2);
        continue;
      }

      // get the phases for taking the transpose of a Pandya-transformed operator
      arma::mat PhaseMatZ(nKets_cc, nKets_cc, arma::fill::ones );
      for (index_t iket=0;iket<nKets_cc;iket++)
      {
         const Ket& ket = tbc_cc.GetKet(iket);
         if ( Z.modelspace->phase( (ket.op->j2 + ket.oq->j2)/2 ) < 0) 
         {
            PhaseMatZ.col( iket ) *= -1;
            PhaseMatZ.row( iket ) *= -1;
         }
      }
      arma::uvec phkets = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
      auto PhaseMatY = PhaseMatZ.rows(phkets) * hy;


//                                           [      |     ]
//     create full Y matrix from the half:   [  Yhp | Y'ph]   where the prime indicates multiplication by (-1)^(i+j+k+l) h_y
//                                           [      |     ]   Flipping hp <-> ph and multiplying by the phase is equivalent to
//                                           [  Yph | Y'hp]   having kets |kj> with k>j.
//
//
//      so <il|Zbar|kj> =  <il|Xbar|hp><hp|Ybar|kj> + <il|Xbar|ph><ph|Ybar|kj>
//
      arma::mat Y_bar_ph_flip = arma::join_vert ( Y_bar_ph.tail_rows(nph_kets)%PhaseMatY ,  Y_bar_ph.head_rows(nph_kets)%PhaseMatY ) ;
      Zbar_ch =  Xt_bar_ph * arma::join_horiz( Y_bar_ph ,  Y_bar_ph_flip );




      // If Z is hermitian, then XY is anti-hermitian, and so XY - YX = XY + (XY)^T
      if ( Z.IsHermitian() )
      {
         Zbar_ch.head_cols(nKets_cc) += Zbar_ch.head_cols(nKets_cc).t();
      }
      else // Z is antihermitian, so XY is hermitian, XY - YX = XY - XY^T
      {
         Zbar_ch.head_cols(nKets_cc) -= Zbar_ch.head_cols(nKets_cc).t();
      }

      // By taking the transpose, we get <il|Zbar|kj> with i>l and k<j, and we want the opposite
      // By the symmetries of the Pandya-transformed matrix element, that means we pick
      // up a factor hZ * phase(i+j+k+l). The hZ cancels the hXhY we have for the "head" part of the matrix
      // so we end up adding in either case.
      Zbar_ch.tail_cols(nKets_cc) += Zbar_ch.tail_cols(nKets_cc).t()%PhaseMatZ;

   }

   X.profiler.timer["Build Z_bar"] += omp_get_wtime() - t_start;

   // Perform inverse Pandya transform on Z_bar to get Z
   t_start = omp_get_wtime();

   // Actually, the Pandya transform has a minus sign in the definition,
   // and the ph commutator has an overall minus sign, so we're technically subtracting
   // the inverse Pandya transformation. Also, the inverse Pandya transformation
   // is just the regular Pandya transformation. The distinction in the code
   // is because some other commutator-specific things are done at the same time.
   AddInversePandyaTransformation(Z_bar, Z);

   Z.modelspace->scalar_transform_first_pass = false;
   X.profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;

}


/*
void comm222_phss_alternative( const Operator& X, const Operator& Y, Operator& Z ) 
{
   auto& X2 = X.TwoBody;
   auto& Y2 = Y.TwoBody;
   auto& Z2 = Z.TwoBody;

   size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
   size_t nch_CC = Z.modelspace->GetNumberTwoBodyChannels_CC();
   std::deque<arma::mat> Z_bar ( nch_CC );

   for (size_t ch_CC=0; ch_CC<nch_CC; ch_CC++)
   {
     TwoBodyChannel_CC& tbc_CC = Z.modelspace->GetTwoBodyChanel_CC(ch_CC);
     size_t nkets_CC = tbc_CC.GetNumberKets();
     int Jcc = tbc_CC.J;
     std::vector<size_t> ph_kets;
     std::vector<double> ph_occs;
     // count how many ab states contribute
     for ( size_t iket_CC=0; iket_CC<nkets_CC; iket_CC++)
     {
       Ket& ket_CC = tbc_CC.GetKet(iket_CC);
       double occfactor = ket_CC.op->occ - ket_CC.oq->occ;
       if ( std::abs(occfactor)<1e-8 ) continue;
       ph_kets.push_back(iket_CC);
       ph_occs.push_back(occfactor);
     }
     /// Allocate Xbar and Ybar 
     size_t nkets_ab = ph_kets.size();
     arma::mat Xbar( nkets_CC,  2*nkets_ab );
     arma::mat Ybar( 2*nkets_ab, nkets_CC );
  
     /// Fill Xbar and Ybar    Xbar^Jcc_ijab  =  - sum_J  (2J+1) { i  j  Jcc } X^J_ibaj
     ///                                                         { a  b  J   } 
     for ( size_t ibra_CC=0; ibra_CC<nkets_CC; ibra_CC++)
     {
       Ket& bra_CC = tbc_CC.GetKet(ibra_CC);
       int j2i = bra_CC.op->j2;
       int j2j = bra_CC.oq->j2;
       double ji = 0.5*j2i;
       double jj = 0.5*j2j;
       for ( size_t index_ab=0; index_ab<nkets_ab; index_ab++)
       {
         size_t iket_ab = ph_kets[index_ab];
         double occ_ab = ph_occs[index_ab];
         Ket& ket_ab = tbc_CC.GetKet(iket_ab);
         int j2a = ket_ab.op->j2;
         int j2b = ket_ab.oq->j2;
  
         double xbar_ijab = 0;
         double ybar_abij = 0;
         int Jmin = std::max( std::abs( oi.j2-ob.j2), std::abs( oj.j2-oa.j2) )/2;
         int Jmax = std::max( ( oi.j2+ob.j2), ( oj.j2+oa.j2) )/2;
         for (int J=Jmin; J<=Jmin; J++)
         {
           double sixj = Z.modelspace->GetSixJ(ja,jb,Jcc,ji,jj,J);
           xbar_ijab -= (2*J+1) * sixj * X2.GetTBME_J(J,i,b,a,j);
           ybar_abij -= (2*J+1) * sixj * X2.GetTBME_J(J,a,j,i,b);
         }
         Xbar(ibra_CC, index_ab) = xbar_ijab;
         Ybar(index_ab, ibra_CC) = ybar_abij * occ_ab;

         double xbar_ijba = 0;
         double ybar_baij = 0;
         int Jmin = std::max( std::abs( oi.j2-oa.j2), std::abs( oj.j2-ob.j2) )/2;
         int Jmax = std::max( ( oi.j2+oa.j2), ( oj.j2+ob.j2) )/2;
         for (int J=Jmin; J<=Jmin; J++)
         {
           double sixj = Z.modelspace->GetSixJ(jb,ja,Jcc,ji,jj,J);
           xbar_ijab -= (2*J+1) * sixj * X2.GetTBME_J(J,i,a,b,j);
           ybar_abij -= (2*J+1) * sixj * X2.GetTBME_J(J,b,j,i,a);
         }
         Xbar(ibra_CC, index_ab + nkets_ab) = xbar_ijba;
         Ybar(index_ab + nkets_ab, ibra_CC) = ybar_baij * (-occ_ab_;

       }// for index_ab
     }// for ibra_CC

   /// MatMult so Zbar = Xbar*Ybar
   arma::mat zbar = Xbar*Ybar;
   zbar -= zbar.t();

   }// for ch_CC
   /// Transform Zbar to Z
   // Z^J_ijkl = -sum_Jcc (2*Jcc+1) { i j J   } Z^Jcc_ilkj
   //                               { k l Jcc }
   // and we need to do the four permutations (1-Pij)(1-Pkl)
   for (size_t ch=0; ch<nch; ch++)
   {
     TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
     int J = tbc.J;
     size_t nkets = tbc.GetNumberKets();
     for (size_t ibra=0; ibra<nkets; ibra++)
     {
       Ket& bra = tbc.GetKet(ibra);
       size_t i = bra.p;
       size_t j = bra.q;
       int j2i = bra.op->j2
       int j2j = bra.oq->j2
       double ji = 0.5 * j2i;
       double jj = 0.5 * j2j;
       for (size_t iket=ibra; iket<nkets; iket++)
       {
         Ket& ket = tbc.GetKet(iket);
         size_t k = ket.p;
         size_t l = ket.q;
         int j2k = ket.op->j2
         int j2l = ket.oq->j2
         double jk = 0.5 * j2k;
         double jl = 0.5 * j2l;
         int Jcc_min = std::max( std::abs(j2i-j2l), std::abs(j2j-j2k) )/2;
         int Jcc_max = std::min( j2i+j2l, j2j+j2k )/2;
         for (int Jcc=Jcc_min; Jcc<=Jcc_max; Jcc++)
         {
           double sixj = Z.modelspace->GetSixJ( ji, jj, J, jk, jl, Jcc );
           
           double zbar_ilkj = 
         }//for Jcc
         
       }// for iket
     }// for ibra
   }// for ch
}
*/



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////   BEGIN SCALAR-SCALAR COMMUTATORS WITH 3-body ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//  For the three-body operators, we often encounter formulas involving the J-coupled
//  permutation operator P(ij/k)^{J1,J}, which is defined as
//
//  P(ij/k)^{J1,J} = 1 - sum_J2 sqrt{(2J1+1)(2J2+1)} (-1)^2(ji+jj+jk)   {ji jk J1} P_ik
//                                                                      {jk J  J2}
//
//                     - sum_J2 sqrt{(2J1+1)(2J2+1)} (-1)^{jj+jk+J1+J2) {jk ji J1} P_jk
//                                                                      {jk J  J2}
//



//*****************************************************************************************
//
//    ~~~~~~~~~~~~~~~~        Uncoupled expression: 
//   /\      /\      /\         Z_0 = 1/36 sum_abcdef (nanbnc n`dn`en`f) (X_abcdef Y_defabc - Y_abcdef X_defabc)
// a(  )d  b(  )e  c(  )f
//   \/      \/      \/      Coupled expression:
//    ~~~~~~~~~~~~~~~~        Z_0 = 1/36 sum_abcdef sum_J1,J2,J  (nanbnc n`dn`en`f)  (2J+1)
//                                    (X^{J1J2J}_abcdef Y^{J1J2J}_defabc - Y^{J1J2J}_abcdef X^{J1J2J}_defabc)
//    Verified with UnitTest
//   
/// Three body commutator expression
/// \f{equation}{
/// Z_0 = \frac{1}{36} \sum_{abcdef} \sum_{J_{1}J_{2}J} n_a n_b n_c \bar{n}_d\bar{n}_e\bar{n}_f (2J+1)
///        ( X^{J_{1}J_{2}J}_{abcdef} Y^{J_1J_2J}_{defabc} - Y^{J_{1}J_{2}J}_{abcdef} X^{J_1J_2J}_{defabc})
/// \f}
///
void comm330ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  double z0 = 0;
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  if (X3.Norm()<1e-6 or Y3.Norm()<1e-6 ) return;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) reduction(+:z0)
  for ( size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    int twoJ = Tbc.twoJ;
    size_t nkets = Tbc.GetNumberKets();
    for (size_t ibra=0; ibra<nkets; ibra++)
    {
      Ket3& bra = Tbc.GetKet(ibra);
      int ea = 2*bra.op->n+bra.op->l;
      int eb = 2*bra.oq->n+bra.oq->l;
      int ec = 2*bra.oR->n+bra.oR->l;
      int tza = bra.op->tz2;
      int tzb = bra.oq->tz2;
      int tzc = bra.oR->tz2;
      double na = bra.op->occ;
      double nb = bra.oq->occ;
      double nc = bra.oR->occ;
      double occnat_a = bra.op->occ_nat;
      double occnat_b = bra.oq->occ_nat;
      double occnat_c = bra.oR->occ_nat;
      double abc_symm = 6;
      if (bra.p==bra.q and bra.q==bra.r) abc_symm = 1;
      else if (bra.p==bra.q or bra.q==bra.r) abc_symm = 3;
      if ( (std::abs( ea-e_fermi[tza]) + std::abs(eb-e_fermi[tzb]) + std::abs(ec-e_fermi[tzc])) > Z.modelspace->GetdE3max() ) continue;
      if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_c*(1-occnat_c)) < Z.modelspace->GetOccNat3Cut() ) continue;
//      for (size_t iket=0;iket<nkets; iket++)
//      for (size_t iket=ibra;iket<nkets; iket++)
      for (size_t iket=0;iket<ibra; iket++)  // dont need iket=ibra because the commutator will be zero
      {
        Ket3& ket = Tbc.GetKet(iket);
        double nd = ket.op->occ;
        double ne = ket.oq->occ;
        double nf = ket.oR->occ;
        double occnat_d = ket.op->occ_nat;
        double occnat_e = ket.oq->occ_nat;
        double occnat_f = ket.oR->occ_nat;
        int ed = 2*ket.op->n+ket.op->l;
        int ee = 2*ket.oq->n+ket.oq->l;
        int ef = 2*ket.oR->n+ket.oR->l;
        int tzd = ket.op->tz2;
        int tze = ket.oq->tz2;
        int tzf = ket.oR->tz2;
        if ( (std::abs( ed-e_fermi[tzd]) + std::abs(ee-e_fermi[tze]) + std::abs(ef-e_fermi[tzf])) > Z.modelspace->GetdE3max() ) continue;
        if ( (occnat_d*(1-occnat_d) * occnat_e*(1-occnat_e) * occnat_f*(1-occnat_f)) < Z.modelspace->GetOccNat3Cut() ) continue;
        // account for the iket>ibra case, which we don't do explicitly
        double occfactor = na*nb*nc*(1-nd)*(1-ne)*(1-nf) - (1-na)*(1-nb)*(1-nc)*nd*ne*nf  ;
        if (std::abs(occfactor)<1e-6) continue;
        double def_symm = 6;
        if (ket.p==ket.q and ket.q==ket.r) def_symm = 1;
        else if (ket.p==ket.q or ket.q==ket.r) def_symm = 3;

        double xabcdef = X3.GetME_pn_ch(ch3,ch3,ibra,iket);
        double yabcdef = Y3.GetME_pn_ch(ch3,ch3,ibra,iket);
        double xdefabc = X3.GetME_pn_ch(ch3,ch3,iket,ibra);
        double ydefabc = Y3.GetME_pn_ch(ch3,ch3,iket,ibra);

        z0 += 1./36 * occfactor * abc_symm * def_symm * (twoJ+1) * (xabcdef * ydefabc  -  yabcdef*xdefabc);

      }// for iket
    }// for ibra
  }// for ch3

  std::cout << "Adding " << z0 << "  to Zero Body" << std::endl;
  Z.ZeroBody += z0;
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}



//*****************************************************************************************
//                   |i
//    *~~~~~~~[X]~~~~*        Uncoupled expression: 
//   / \      / \    |          Z_ij = 1/12 sum_abcde (nanb n`c n`dn`e) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
// a(   )c  b(   )d  |e 
//   \ /      \ /    |       Coupled expression:
//    *~~~~~~~[Y]~~~~*        Z_ij = 1/12 sum_abcde sum_J1,J2,J  (nanb n`c n`dn`e)  (2J+1)/(2ji+1)
//                   |j               (X^{J1J2J}_abicde Y^{J1J2J}_cdeabj - Y^{J1J2J}_abicde X^{J1J2J}_cdeabj)
//
//  Verfied with UnitTest
//
void comm331ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z1 = Z.OneBody;
  int herm = Z.IsAntiHermitian() ? -1 : 1 ;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  if (X.GetParticleRank() < 3) return;
  if (Y.GetParticleRank() < 3) return;

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
  size_t norb = Z.modelspace->GetNumberOrbits();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t i=0; i<norb; i++)
  {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    int ei = 2*oi.n + oi.l;
    double occnat_i = oi.occ_nat;
    int tzi = oi.tz2;
    for ( auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      if (j>i) continue;
      double zij=0;
      Orbit& oj = Z.modelspace->GetOrbit(j);
      double occnat_j = oj.occ_nat;
      int ej = 2*oj.n + oj.l;
      int tzj = oj.tz2;

      for (size_t ch_ab=0; ch_ab<nch2; ch_ab++)
      {
        auto& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
        int Jab = tbc_ab.J;
        size_t nkets_ab = tbc_ab.GetNumberKets();
        int twoJ_min = std::abs(2*Jab - oi.j2);
        int twoJ_max = 2*Jab+oi.j2;
        
        for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
        {
           Ket& ket_ab = tbc_ab.GetKet(iket_ab);
           size_t a = ket_ab.p;
           size_t b = ket_ab.q;
           if (std::abs( ket_ab.op->occ * ket_ab.oq->occ )<1e-6 ) continue;
           int ea = 2*ket_ab.op->n + ket_ab.op->l;
           int eb = 2*ket_ab.oq->n + ket_ab.oq->l;
           int tza = ket_ab.op->tz2;
           int tzb = ket_ab.oq->tz2;
           double occnat_a = ket_ab.op->occ_nat;
           double occnat_b = ket_ab.oq->occ_nat;
//      double occnat_c = bra.oR->occ_nat;
           if ( (std::abs( ea-e_fermi[tza]) + std::abs(eb-e_fermi[tzb]) + std::abs(ei-e_fermi[tzi])) > Z.modelspace->GetdE3max() ) continue;
           if ( (std::abs( ea-e_fermi[tza]) + std::abs(eb-e_fermi[tzb]) + std::abs(ej-e_fermi[tzj])) > Z.modelspace->GetdE3max() ) continue;
           if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i)) < Z.modelspace->GetOccNat3Cut() ) continue;
           if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j)) < Z.modelspace->GetOccNat3Cut() ) continue;

           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
             size_t ch_abi = Z.modelspace->GetThreeBodyChannelIndex( twoJ, (tbc_ab.parity +oi.l)%2, tbc_ab.Tz*2 + oi.tz2 );
             if (ch_abi==-1) continue; // maybe that channel doesn't exist
             auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch_abi);
//             size_t index_abi = Tbc.GetLocalIndex(a,b,i,Jab);
//             size_t index_abj = Tbc.GetLocalIndex(a,b,j,Jab);
//             std::cout << "index_abi = " << index_abi << "   index_abj = " << index_abj << std::endl;
//             if ((index_abi==-1) or (index_abj==-1)) continue;
             double Jfactor = (twoJ+1.0)/(oi.j2+1);
             double ab_symmetry_factor = (a==b) ?  1.0 : 2.0;
             size_t nkets3 = Tbc.GetNumberKets();
             for (size_t iket_cde=0; iket_cde<nkets3; iket_cde++)
             {
                Ket3& ket_cde = Tbc.GetKet(iket_cde);
                double occfactor = (ket_ab.op->occ * ket_ab.oq->occ) * (1-ket_cde.op->occ)*(1-ket_cde.oq->occ)*(1-ket_cde.oR->occ);
                if (std::abs(occfactor)<1e-6) continue;
                size_t c = ket_cde.p;
                size_t d = ket_cde.q;
                size_t e = ket_cde.r;
                int ec = 2*ket_cde.op->n + ket_cde.op->l;
                int ed = 2*ket_cde.oq->n + ket_cde.oq->l;
                int ee = 2*ket_cde.oR->n + ket_cde.oR->l;
                double occnat_c = ket_cde.op->occ_nat;
                double occnat_d = ket_cde.oq->occ_nat;
                double occnat_e = ket_cde.oR->occ_nat;
                int tzc = ket_cde.op->tz2;
                int tzd = ket_cde.oq->tz2;
                int tze = ket_cde.oR->tz2;
                if ( (std::abs( ec-e_fermi[tzc]) + std::abs(ed-e_fermi[tzd]) + std::abs(ee-e_fermi[tze])) > Z.modelspace->GetdE3max() ) continue;
                if ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_e*(1-occnat_e)) < Z.modelspace->GetOccNat3Cut() ) continue;
                double cde_symmetry_factor = 6;
                if (c==d and d==e) cde_symmetry_factor = 1;
                else if (c==d or d==e) cde_symmetry_factor = 3;
                int Jcd = ket_cde.Jpq;
                 
                zij += 1./12 * ab_symmetry_factor * cde_symmetry_factor * occfactor * Jfactor
                             * (  X3.GetME_pn( Jab, Jcd, twoJ, a,b,i,c,d,e) * Y3.GetME_pn( Jcd, Jab, twoJ, c,d,e,a,b,j)  
                                - Y3.GetME_pn( Jab, Jcd, twoJ, a,b,i,c,d,e) * X3.GetME_pn( Jcd, Jab, twoJ, c,d,e,a,b,j)  ); 
             }// for iket_cde
           }// for twoJ
        }// for iket_ab
      }// for ch_ab

      Z1(i,j) += zij;
      if ( i!=j ) Z1(j,i) += herm * zij;
    }// for j
  }// for i

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}




//*****************************************************************************************
//                   |i
//    *~~[X]~~*      |        Uncoupled expression: 
//   / \      / \    |          Z_ij = 1/4 sum_abcd (nanb n`cn`d) (X_abcd Y_cdiabj - Y_abicdj X_cdab)
// a(   )c  b(   )d  | 
//   \ /      \ /    |       Coupled expression:
//   *~~~~~~~[Y]~~~~~*        Z_ij = 1/4sum_abcd sum_J1J  (nanb n`c n`d)  (2J+1)/(2ji+1)
//                   |j               (X^{J1}_abcd Y^{J1J1J}_cdiabj - Y^{J1J1J}_abicdj X^{J1}_cdab)
//                   |
//                              We only sum a<=b and c<=d, so we do not explicitly include the factor 1/4,
//                              except for the a==b, or c==d case, where there is no double counting.
//  Verfied with UnitTest
//
void comm231ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y2 = Y.TwoBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z1 = Z.OneBody;
  int x_particle_rank = X.GetParticleRank();
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t norb = Z.modelspace->GetNumberOrbits();
  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  for (size_t i=0; i<norb; i++)
  {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    int ei = 2*oi.n + oi.l;
    double d_ei = std::abs( ei - e_fermi[oi.tz2]);
    double occnat_i = oi.occ_nat;
    for ( auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      if (j>i) continue;
      Orbit& oj = Z.modelspace->GetOrbit(j);
      int ej = 2*oj.n + oj.l;
      double d_ej = std::abs( ej - e_fermi[oj.tz2] );
      double occnat_j = oj.occ_nat;
      double zij=0;
      for (int ch=0; ch<nch; ch++)
      {
        auto tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J = tbc.J;
        size_t nkets = tbc.GetNumberKets();
//        for ( auto ibra : tbc.KetIndex_hh )
        for ( size_t ibra=0; ibra<nkets; ibra++ )
        {
          Ket& bra = tbc.GetKet(ibra);
          int a = bra.p;
          int b = bra.q;
          int ea = 2*bra.op->n + bra.op->l;
          int eb = 2*bra.oq->n + bra.oq->l;
          double occnat_a = bra.op->occ_nat;
          double occnat_b = bra.oq->occ_nat;
          double d_ea = std::abs( ea - e_fermi[bra.op->tz2]);
          double d_eb = std::abs( eb - e_fermi[bra.oq->tz2]);
          if (  (ea+eb+std::min(ei,ej))> Z.modelspace->E3max )  continue;
          if (  (d_ea+d_eb+std::min(d_ei,d_ej))> Z.modelspace->dE3max )  continue;
          if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * std::max( occnat_i*(1-occnat_i), occnat_j*(1-occnat_j)) ) < Z.modelspace->GetOccNat3Cut() ) continue;
          double na = bra.op->occ;
          double nb = bra.oq->occ;
          if ( (std::abs(na*nb) + std::abs( (1-na)*(1-nb))) < 1e-6 ) continue;
//          for ( auto iket : tbc.KetIndex_pp )
//          for ( size_t iket=0; iket<nkets; iket++ )
          for ( size_t iket=0; iket<ibra; iket++ )
          {
            Ket& ket = tbc.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            int ec = 2*ket.op->n + ket.op->l;
            int ed = 2*ket.oq->n + ket.oq->l;
            double d_ec = std::abs( ec - e_fermi[ket.op->tz2]);
            double d_ed = std::abs( ec - e_fermi[ket.oq->tz2]);
            double occnat_c = bra.op->occ_nat;
            double occnat_d = bra.oq->occ_nat;
            if (  (ec+ed+std::min(ei,ej))> Z.modelspace->E3max )  continue;
            if (  (d_ec+d_ed+std::min(d_ei,d_ej))> Z.modelspace->dE3max )  continue;
            if ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * std::max( occnat_i*(1-occnat_i), occnat_j*(1-occnat_j)) ) < Z.modelspace->GetOccNat3Cut() ) continue;
            double nc = ket.op->occ;
            double nd = ket.oq->occ;
//            double prefactor = na*nb*(1-nc)*(1-nd);
            double prefactor =  na*nb*(1-nc)*(1-nd) - (1-na)*(1-nb)*nc*nd;
            if ( std::abs(prefactor)<1e-8) continue;
            double Xabcd = X2.GetTBME(ch,bra,ket);
            double Yabcd = Y2.GetTBME(ch,bra,ket);
            double Xcdab = X2.GetTBME(ch,ket,bra);
            double Ycdab = Y2.GetTBME(ch,ket,bra);
            if (a==b) prefactor /= 2;
            if (c==d) prefactor /= 2;
            int twoJ_min = std::abs( 2*J - oi.j2);
            int twoJ_max = 2*J + oi.j2;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
              double xabicdj = 0;
              double yabicdj = 0;
              double xcdiabj = 0;
              double ycdiabj = 0;
              if ( ( std::max(ea+eb+ej,ec+ed+ei) <= Z.modelspace->E3max )
               and ( std::max(d_ea+d_eb+d_ej,d_ec+d_ed+d_ei) <= Z.modelspace->dE3max )
               and (  (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_i*(1-occnat_i) ) > Z.modelspace->GetOccNat3Cut() )
               and (  (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) > Z.modelspace->GetOccNat3Cut() )   )
              {
//                xcdiabj = X3.GetME_pn(J,J,twoJ,c,d,i,a,b,j);
                xcdiabj = x_particle_rank>2 ? X3.GetME_pn(J,J,twoJ,c,d,i,a,b,j) :  0;
                ycdiabj = Y3.GetME_pn(J,J,twoJ,c,d,i,a,b,j);
              }
              if ( (std::max(ea+eb+ei,ec+ed+ej) <= Z.modelspace->E3max )
               and ( std::max(d_ea+d_eb+d_ei,d_ec+d_ed+d_ej) <= Z.modelspace->dE3max ) 
               and (  (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_j*(1-occnat_j) ) > Z.modelspace->GetOccNat3Cut() )
               and (  (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) > Z.modelspace->GetOccNat3Cut() )   )
              {
//                xabicdj = X3.GetME_pn(J,J,twoJ,a,b,i,c,d,j);
                xabicdj = x_particle_rank>2 ? X3.GetME_pn(J,J,twoJ,a,b,i,c,d,j) : 0 ;
                yabicdj = Y3.GetME_pn(J,J,twoJ,a,b,i,c,d,j);
              }
              zij += prefactor * (twoJ+1) * ( (Xabcd * ycdiabj - yabicdj * Xcdab)
                                           -  (Yabcd * xcdiabj - xabicdj * Ycdab) );
            }
          }
        }
      }
      Z1(i,j) += zij / (oi.j2+1.0);
      if (i!=j)
      {
         Z1(j,i) += zij / (oi.j2+1.0);
      }
    }// for j
  }// for i

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}






void comm231ss_slow( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y2 = Y.TwoBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z1 = Z.OneBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t norb = Z.modelspace->GetNumberOrbits();
//  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  for (size_t i=0; i<norb; i++)
  {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    int ei = 2*oi.n + oi.l;
    double d_ei = std::abs( ei - e_fermi[oi.tz2]);
    double occnat_i = oi.occ_nat;
    for ( auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      if (j>i) continue;
      Orbit& oj = Z.modelspace->GetOrbit(j);
      int ej = 2*oj.n + oj.l;
      double d_ej = std::abs( ej - e_fermi[oj.tz2]);
      double occnat_j = oj.occ_nat;
      double zij=0;

      std::vector<std::array<size_t,2>> Xchannels;
      std::vector<std::array<size_t,2>> Ychannels;
      for ( auto& iter : X.TwoBody.MatEl ) Xchannels.push_back(iter.first);
      for ( auto& iter : Y.TwoBody.MatEl ) Ychannels.push_back(iter.first);
      size_t nxchan = Xchannels.size();
      size_t nychan = Ychannels.size();

//      for (int ch=0; ch<nch; ch++)
      for ( size_t ichan=0; ichan<nxchan; ichan++)
      {
        size_t ch_bra = Xchannels[ichan][0];
        size_t ch_ket = Xchannels[ichan][1];
        auto tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        auto tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int J = tbc_bra.J;
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
//        for ( auto ibra : tbc_bra.KetIndex_hh )
        for ( size_t ibra=0; ibra<nbras; ibra++ )
        {
          Ket& bra = tbc_bra.GetKet(ibra);
          int a = bra.p;
          int b = bra.q;
          int ea = 2*bra.op->n + bra.op->l;
          int eb = 2*bra.oq->n + bra.oq->l;
          double d_ea = std::abs( ea - e_fermi[bra.op->tz2]);
          double d_eb = std::abs( eb - e_fermi[bra.oq->tz2]);
          double occnat_a = bra.op->occ_nat;
          double occnat_b = bra.oq->occ_nat;
          if (  (ea+eb+std::min(ei,ej))> Z.modelspace->E3max )  continue;
          if (  (d_ea+d_eb+std::min(d_ei,d_ej))> Z.modelspace->dE3max )  continue;
          if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * std::max( occnat_i*(1-occnat_i), occnat_j*(1-occnat_j)) ) < Z.modelspace->GetOccNat3Cut() ) continue;
          double na = bra.op->occ;
          double nb = bra.oq->occ;
          if ( std::abs(na*nb)<1e-6) continue;
//          for ( auto iket : tbc.KetIndex_pp )
          for ( size_t iket=0; iket<nkets; iket++ )
          {
            Ket& ket = tbc_ket.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            int ec = 2*ket.op->n + ket.op->l;
            int ed = 2*ket.oq->n + ket.oq->l;
            double d_ec = std::abs( ec - e_fermi[ket.op->tz2]);
            double d_ed = std::abs( ec - e_fermi[ket.oq->tz2]);
            double occnat_c = bra.op->occ_nat;
            double occnat_d = bra.oq->occ_nat;
            if (  (ec+ed+std::min(ei,ej))> Z.modelspace->E3max )  continue;
            if (  (d_ec+d_ed+std::min(d_ei,d_ej))> Z.modelspace->dE3max )  continue;
            if ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * std::max( occnat_i*(1-occnat_i), occnat_j*(1-occnat_j)) ) < Z.modelspace->GetOccNat3Cut() ) continue;
            double nc = ket.op->occ;
            double nd = ket.oq->occ;
            double prefactor = na*nb*(1-nc)*(1-nd);
            if ( std::abs(prefactor)<1e-8) continue;
            double Xabcd = X2.GetTBME(ch_bra,ch_ket,bra,ket);
//            double Yabcd = Y2.GetTBME(ch,bra,ket);
            double Xcdab = X2.GetTBME(ch_ket,ch_bra,ket,bra);
//            double Ycdab = Y2.GetTBME(ch,ket,bra);
            if (a==b) prefactor /= 2;
            if (c==d) prefactor /= 2;
            int twoJ_min = std::abs( 2*J - oi.j2);
            int twoJ_max = 2*J + oi.j2;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
//              double xabicdj = 0;
              double yabicdj = 0;
//              double xcdiabj = 0;
              double ycdiabj = 0;
              if ( ( std::max(ea+eb+ej,ec+ed+ei) <= Z.modelspace->E3max )
               and ( std::max(d_ea+d_eb+d_ej,d_ec+d_ed+d_ei) <= Z.modelspace->dE3max ) 
               and (  (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_i*(1-occnat_i) ) > Z.modelspace->GetOccNat3Cut() )
               and (  (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) > Z.modelspace->GetOccNat3Cut() )   )
              {
//                xcdiabj = X3.GetME_pn(J,J,twoJ,c,d,i,a,b,j);
                ycdiabj = Y3.GetME_pn(J,J,twoJ,c,d,i,a,b,j);
              }
              if ( (std::max(ea+eb+ei,ec+ed+ej) <= Z.modelspace->E3max )
               and ( std::max(d_ea+d_eb+d_ei,d_ec+d_ed+d_ej) <= Z.modelspace->dE3max ) 
               and (  (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_j*(1-occnat_j) ) > Z.modelspace->GetOccNat3Cut() )
               and (  (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) > Z.modelspace->GetOccNat3Cut() )   )
              {
//                xabicdj = X3.GetME_pn(J,J,twoJ,a,b,i,c,d,j);
                yabicdj = Y3.GetME_pn(J,J,twoJ,a,b,i,c,d,j);
              }
              zij += prefactor * (twoJ+1) * ( (Xabcd * ycdiabj - yabicdj * Xcdab) );
//              zij += prefactor * (twoJ+1) * ( (Xabcd * ycdiabj - yabicdj * Xcdab)
//                                           -  (Yabcd * xcdiabj - xabicdj * Ycdab) );
              
            }// for twoJ
          }// for iket
        }// for ibra
      }// for ichan

      for ( size_t ichan=0; ichan<nychan; ichan++)
      {
        size_t ch_bra = Ychannels[ichan][0];
        size_t ch_ket = Ychannels[ichan][1];
        auto tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        auto tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int J = tbc_bra.J;
//        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
        for ( auto ibra : tbc_bra.KetIndex_hh )
        {
          Ket& bra = tbc_bra.GetKet(ibra);
          int a = bra.p;
          int b = bra.q;
          int ea = 2*bra.op->n + bra.op->l;
          int eb = 2*bra.oq->n + bra.oq->l;
          double d_ea = std::abs(ea - e_fermi[bra.op->tz2]);
          double d_eb = std::abs(eb - e_fermi[bra.oq->tz2]);
          double occnat_a = bra.op->occ_nat;
          double occnat_b = bra.oq->occ_nat;
          if (  (ea+eb+std::min(ei,ej))> Z.modelspace->E3max )  continue;
          if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * std::max( occnat_i*(1-occnat_i), occnat_j*(1-occnat_j)) ) < Z.modelspace->GetOccNat3Cut() ) continue;
          double na = bra.op->occ;
          double nb = bra.oq->occ;
//          for ( auto iket : tbc.KetIndex_pp )
          for ( size_t iket=0; iket<nkets; iket++ )
          {
            Ket& ket = tbc_ket.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            int ec = 2*ket.op->n + ket.op->l;
            int ed = 2*ket.oq->n + ket.oq->l;
            double d_ec = std::abs(ec - e_fermi[ket.op->tz2]);
            double d_ed = std::abs(ec - e_fermi[ket.oq->tz2]);
            double occnat_c = bra.op->occ_nat;
            double occnat_d = bra.oq->occ_nat;
            if (  (ec+ed+std::min(ei,ej))> Z.modelspace->E3max )  continue;
            if ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * std::max( occnat_i*(1-occnat_i), occnat_j*(1-occnat_j)) ) < Z.modelspace->GetOccNat3Cut() ) continue;
            double nc = ket.op->occ;
            double nd = ket.oq->occ;
            double prefactor = na*nb*(1-nc)*(1-nd);
            if ( std::abs(prefactor)<1e-8) continue;
//            double Xabcd = X2.GetTBME(ch,bra,ket);
            double Yabcd = Y2.GetTBME(ch_bra,ch_ket,bra,ket);
//            double Xcdab = X2.GetTBME(ch,ket,bra);
            double Ycdab = Y2.GetTBME(ch_ket,ch_bra,ket,bra);
            if (a==b) prefactor /= 2;
            if (c==d) prefactor /= 2;
            int twoJ_min = std::abs( 2*J - oi.j2);
            int twoJ_max = 2*J + oi.j2;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
              double xabicdj = 0;
//              double yabicdj = 0;
              double xcdiabj = 0;
//              double ycdiabj = 0;
              if ( ( std::max(ea+eb+ej,ec+ed+ei) <= Z.modelspace->E3max )
               and ( std::max(d_ea+d_eb+d_ej,d_ec+d_ed+d_ei) <= Z.modelspace->dE3max ) 
               and (  (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_i*(1-occnat_i) ) > Z.modelspace->GetOccNat3Cut() )
               and (  (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) > Z.modelspace->GetOccNat3Cut() )   )
              {
                xcdiabj = X3.GetME_pn(J,J,twoJ,c,d,i,a,b,j);
//                ycdiabj = Y3.GetME_pn(J,J,twoJ,c,d,i,a,b,j);
              }
              if ( (std::max(ea+eb+ei,ec+ed+ej) <= Z.modelspace->E3max )
               and ( std::max(d_ea+d_eb+d_ei,d_ec+d_ed+d_ej) <= Z.modelspace->dE3max ) 
               and (  (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_j*(1-occnat_j) ) > Z.modelspace->GetOccNat3Cut() )
               and (  (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) > Z.modelspace->GetOccNat3Cut() )   )
              {
                xabicdj = X3.GetME_pn(J,J,twoJ,a,b,i,c,d,j);
//                yabicdj = Y3.GetME_pn(J,J,twoJ,a,b,i,c,d,j);
              }
              zij -= prefactor * (twoJ+1) * ( (Yabcd * xcdiabj - xabicdj * Ycdab) );
//              zij += prefactor * (twoJ+1) * ( (Xabcd * ycdiabj - yabicdj * Xcdab)
//                                           -  (Yabcd * xcdiabj - xabicdj * Ycdab) );
              
            }// for twoJ
          }// for iket
        }// for ibra
      }// for ichan



      Z1(i,j) += zij / (oi.j2+1.0);
      if (i!=j)
      {
         Z1(j,i) += zij / (oi.j2+1.0);
      }
    }// for j
  }// for i

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}



//*****************************************************************************************
//
// i|  j|   *~~[X]  Uncoupled expression:
//  |   |  / \          Z_ijkl = sum_ab (nan`b-n`anb) (X_ab * Y_ijbkla)
//  |   | (a  )b
//  |   |  \ /
//  *~~[Y]~~*       Coupled expression:
//  |   |              Z_{ijkl}^{J} = sum_ab (nan`b-n`anb) sum_J' (2J'+1)/(2J+1) ( X_ab * Y_{ijbkla}^{J,J,J'} )
// k|  l|                                          
//                           
//                              
//                 When we call Symmetrize, we copy the upper triangle <ibra|iket> with ibra<=iket onto the lower triangle.
//                 Of course, calling AddToTBME updates both according to the symmetry, so we don't need to worry so much... 
//
//  Verfied with UnitTest
//
void comm132ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X1 = X.OneBody;
  auto& X3 = X.ThreeBody;
  auto& Y1 = Y.OneBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;
  int x_particle_rank = X.GetParticleRank();
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();
  
  int norb = Z.modelspace->GetNumberOrbits();
  size_t nch = Z.modelspace->GetNumberTwoBodyChannels();

  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch=0; ch<nch; ch++)
  {
    auto& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;ibra++ ) // <ij| states
    {
      Ket& bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      int ei = 2*bra.op->n + bra.op->l;
      int ej = 2*bra.oq->n + bra.oq->l;
      double d_ei = std::abs(ei - e_fermi[bra.op->tz2]);
      double d_ej = std::abs(ej - e_fermi[bra.oq->tz2]);
      double occnat_i = bra.op->occ_nat;
      double occnat_j = bra.oq->occ_nat;
      for (int iket=ibra;iket<nkets;iket++ ) // |kl> states
      {
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        int ek = 2*ket.op->n + ket.op->l;
        int el = 2*ket.oq->n + ket.oq->l;
        double d_ek = std::abs( ek - e_fermi[ket.op->tz2]);
        double d_el = std::abs( el - e_fermi[ket.oq->tz2]);
        double occnat_k = ket.op->occ_nat;
        double occnat_l = ket.oq->occ_nat;
        double zijkl = 0;
        for (int a=0;a<norb;a++)
        {
          Orbit& oa = Z.modelspace->GetOrbit(a);
          int ea = 2*oa.n + oa.l;
          double d_ea = std::abs( ea - e_fermi.at(oa.tz2));
          double occnat_a = oa.occ_nat;
          if ( (ek+el+ea)>Z.modelspace->E3max) continue;
          if ( (d_ek+d_el+d_ea)>Z.modelspace->dE3max) continue;
          if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_a*(1-occnat_a)) < Z.modelspace->GetOccNat3Cut() ) continue;
          // this is less efficient than doing the loop twice, but less code
          // and easily lets us treat the case of Y having nonzero Tz or odd parity
          std::set<size_t> blist;
//          std::cout << "HERE AT LINE " << __LINE__ <<  " and a is " << oa.l << " " << oa.j2 << " " << oa.tz2 << std::endl;
          for ( auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2})) blist.insert(b);
          for ( auto b : Y.OneBodyChannels.at({oa.l,oa.j2,oa.tz2})) blist.insert(b);
//          for ( auto b : Z.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) // TODO: We can make this a<=b or a>=b, I think. Just need to mind some factors of 2
          for ( auto b : blist ) 
          {
            Orbit& ob = Z.modelspace->GetOrbit(b);
            int eb = 2*ob.n + ob.l;
            double d_eb = std::abs( eb - e_fermi[ob.tz2]);
            double occnat_b = ob.occ_nat;
            if ( (ei+ej+eb)>Z.modelspace->E3max) continue;
            if ( (d_ei+d_ej+d_eb)>Z.modelspace->dE3max) continue;
            if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_b*(1-occnat_b) ) < Z.modelspace->GetOccNat3Cut() ) continue;
            double occfactor = oa.occ - ob.occ;
            if (std::abs(occfactor)<1e-6) continue;
            int twoJ_min = std::abs( oa.j2 - 2*J );
            int twoJ_max = oa.j2 + 2*J;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
//              double xijbkla = X3.GetME_pn(J,J,twoJ,i,j,b,k,l,a);
              double xijbkla = x_particle_rank > 2 ? X3.GetME_pn(J,J,twoJ,i,j,b,k,l,a) : 0;
              double yijbkla = Y3.GetME_pn(J,J,twoJ,i,j,b,k,l,a);

              zijkl += occfactor * (twoJ+1.)/(2*J+1) * ( X1(a,b) * yijbkla -  Y1(a,b) * xijbkla );

            }
          }
        }
        // normalize the tbme
        zijkl /= sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
        Z2.AddToTBME(ch, ch, ibra, iket, zijkl );
      }
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}






//*****************************************************************************************
//
//  |         |
// i|        j|     Uncoupled expression:
//  *~~[X]~*  |         Z_ijkl = -1/2 sum_abc (nanbn`c+n`an`bnc) ( (1-Pij) X_icab * Y_abjklc - (1-Pkl) Yijcabl * Xabkc )
// a|  b|  c\ |
//  |   |    \|
//  *~~[Y]~~~~*      Coupled expression:
//  |   |              Z_{ijkl}^{J} = -1/2 sum_abc (nanbn`c+n`an`bnc) sum_J'J" (2J'+1)(2J"+1)/sqrt(2J+1) (-1)^{2J"+J'-J}
// k|  l|                           *  [   (1 - (-1)^{i+j-J}Pij) (-1)^{j-c} { j  J" J' } X_icab^{J'} * Y_{abjklc}^{J'JJ"}    
//                                                                          { c  i  J  }
//                           
//                                       -(1 - (-1)^{k+l-J}Pkl)  (-1)^{l-c} { l  J" J' } Y_{ijcabl}^{JJ'J"} * X_{abkc}^{J'}  ]
//                                                                          { c  k  J  }
//
//   The factor 1/2 out front is absorbed by the fact that we only do the ordering a<=b  (need to deal carefully with a==b)
//
//  Verified with UnitTest
//
// This is the time hog of the n^7 scaling terms   (seems to be doing better...)

size_t Hash_comm232_key( std::array<size_t,5>& kljJJ )
{
  return kljJJ[0] + (kljJJ[1] << 8) + (kljJJ[2] << 16) + (kljJJ[3] << 24) + (kljJJ[4]<<32);
}

void comm232ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y2 = Y.TwoBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;

  bool x_has_3 = X3.IsAllocated();
  bool y_has_3 = Y3.IsAllocated();
//  bool x_has_3 = X3.is_allocated;
//  bool y_has_3 = Y3.is_allocated;

  int hermX = X.IsHermitian() ? 1 : -1;
  int hermY = Y.IsHermitian() ? 1 : -1;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  int nch = Z.modelspace->GetNumberTwoBodyChannels();
//  std::cout << __func__ << " begin" << std::endl;

  // first, enumerate the one-body channels |i> => ji,parityi,tzi
  std::map<std::array<int,3>,std::vector<size_t>> local_one_body_channels; //  maps {j,parity,tz} => vector of index_i 
  for ( auto i : Z.modelspace->all_orbits )
  {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    std::array<int,3> obc = {oi.j2,oi.l%2,oi.tz2};
    if ( local_one_body_channels.find(obc) == local_one_body_channels.end() ) local_one_body_channels[obc] = {i};
    else local_one_body_channels[obc].push_back(i);
  }
  
  // next, figure out which three-body states |klj`> and |abc`> exist, make a list, and give them an
  // index for where they'll sit in the matrix
  std::map<std::array<int,3>,std::vector<std::array<size_t,4>>> klj_list; //  maps {j,parity,tz} => {kljJ}
  std::map<std::array<int,3>,arma::mat> ZMAT_list; //  maps {j,parity,tz} => matrix <i|Z|klj`>

  std::vector<std::array<int,3>> obc_keys;
  for ( auto& iter_i : local_one_body_channels) obc_keys.push_back(iter_i.first);
  size_t nkeys = obc_keys.size();

  
  for ( size_t ikey=0; ikey<nkeys; ikey++ )
  {
    auto& obc_key = obc_keys[ikey];
    std::vector<size_t>& obc_orbits = local_one_body_channels[obc_key];
    int j2i = obc_key[0];
    int parityi = obc_key[1];
    int tz2i = obc_key[2];
    klj_list[obc_key] = {};
    auto& klj_list_i = klj_list[obc_key];

    for ( auto& iter_j : local_one_body_channels )
    {
      int j2j = iter_j.first[0];
      int parityj = iter_j.first[1];
      int tz2j = iter_j.first[2];

      int Jkl_min = std::abs(j2i-j2j)/2;
      int Jkl_max = (j2i+j2j)/2;
      int Tzkl = (tz2i + tz2j)/2;
      int paritykl = (parityi+parityj)%2;
      for ( int Jkl = Jkl_min; Jkl<=Jkl_max; Jkl++)
      {

        size_t ch_kl = Z.modelspace->GetTwoBodyChannelIndex(Jkl, paritykl, Tzkl);
        TwoBodyChannel& tbc_kl = Z.modelspace->GetTwoBodyChannel(ch_kl);
        size_t nkets_kl = tbc_kl.GetNumberKets();
        for ( size_t iket_kl=0; iket_kl<nkets_kl; iket_kl++)
        {
          Ket& ket_kl = tbc_kl.GetKet(iket_kl);
          double de_k = std::abs( 2*ket_kl.op->n + ket_kl.op->l - e_fermi[ket_kl.op->tz2]);
          double de_l = std::abs( 2*ket_kl.oq->n + ket_kl.oq->l - e_fermi[ket_kl.oq->tz2]);
          double nk = ket_kl.op->occ;
//          double nl = ket_kl.oq->occ;
          double occnat_k = ket_kl.op->occ_nat;
          double occnat_l = ket_kl.oq->occ_nat;
          for ( size_t j : iter_j.second )
          {
            Orbit& oj = Z.modelspace->GetOrbit(j);
//            double de_j = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
            double occnat_j = oj.occ_nat;
            if ( (de_k + de_l) > Z.modelspace->GetdE3max() ) continue;
            if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_j*(1-occnat_j) ) < Z.modelspace->GetOccNat3Cut() ) continue;
//            if ( (de_k + de_l + de_j) > Z.modelspace->GetdE3max() ) continue;
//            double nj = oj.occ;

            klj_list_i.push_back( { ket_kl.p, ket_kl.q, j, (size_t)tbc_kl.J } );

          }// for j
        }// for iket_kl
      }// for Jkl
    }// for iter_j in local one body channels

    size_t dim_i   = obc_orbits.size(); // how many sp states are in this jpt channel
    size_t dim_klj = klj_list[obc_key].size();
    ZMAT_list[obc_key] =  arma::mat(dim_i,dim_klj) ;
  }// for ikey


  // This loop is what takes all the time.
  // Outer loop over one-body quantum numbers, which also label the three-body pph states |klj`> and |abc`>
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for ( size_t ikey=0; ikey<nkeys; ikey++ )
  {

    auto& obc_key = obc_keys[ikey]; // this is an array {j2,parity,tz2}
    int j2i = obc_key[0];
    std::vector<size_t>& obc_orbits = local_one_body_channels[obc_key]; // the orbits that have quantum numbers {j2,parity,tz2}
    auto& klj_list_i = klj_list[obc_key]; // list of 3-body pph stats |klj`> with quantum numbers {j2,parity,tz2}

    std::vector<size_t> abc_list; // keep track of which klj states should go in the abc loop
    std::vector<double> abc_occ_list;
    for (size_t i_kljJ=0; i_kljJ< klj_list[obc_key].size(); i_kljJ++ )
    {
      auto& kljJ = klj_list[obc_key][i_kljJ];
      size_t a = kljJ[0];
      size_t b = kljJ[1];
      size_t c = kljJ[2];
      Orbit& oa = Z.modelspace->GetOrbit(a);
      Orbit& ob = Z.modelspace->GetOrbit(b);
      Orbit& oc = Z.modelspace->GetOrbit(c);
      double na = oa.occ;
      double nb = ob.occ;
      double nc = oc.occ;
      double occupation_factor = na*nb*(1-nc) + (1-na)*(1-nb)*nc;
      if (std::abs(occupation_factor)<1e-6) continue;
      if (a==b) occupation_factor *=0.5;  // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b
      abc_list.push_back( i_kljJ );
      abc_occ_list.push_back( occupation_factor );
    }

    size_t dim_i   = obc_orbits.size(); // how many sp states are in this jpt channel
    size_t dim_klj = klj_list[obc_key].size(); // how many 3-body pph states in this jpt channel
    size_t dim_abc = abc_list.size(); // how many 3-body states which contribute to the |abc`> sum

    // Now allocate the matrices in this channel
    arma::mat X2MAT(dim_i,   dim_abc, arma::fill::zeros);
    arma::mat Y2MAT(dim_i,   dim_abc, arma::fill::zeros);
    arma::mat X3MAT(dim_abc, dim_klj, arma::fill::zeros);
    arma::mat Y3MAT(dim_abc, dim_klj, arma::fill::zeros);

   //figure out which recouplings we'll need when filling the matrices
    std::set<std::array<size_t,5>> kljJJ_needed; // we use a set to avoid repeats
    for (size_t ind_abc=0; ind_abc<dim_abc; ind_abc++) 
    {
      auto& abcJ = klj_list_i[ abc_list[ind_abc] ];
      size_t a = abcJ[0];
      size_t b = abcJ[1];
      size_t c = abcJ[2];
      int Jab  = abcJ[3];
      Orbit& oa = Z.modelspace->GetOrbit(a);
      Orbit& ob = Z.modelspace->GetOrbit(b);
      Orbit& oc = Z.modelspace->GetOrbit(c);
      double de_a = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
      double de_b = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
      double de_c = std::abs( 2*oc.n + oc.l - e_fermi[oc.tz2]);
      double occnat_a = oa.occ_nat;
      double occnat_b = ob.occ_nat;
      double occnat_c = oc.occ_nat;
      int j2c = oc.j2;
      for (size_t ind_klj=0; ind_klj<dim_klj; ind_klj++)
      {
        auto& kljJ = klj_list_i[ ind_klj ];
        size_t k = kljJ[0];
        size_t l = kljJ[1];
        size_t j = kljJ[2];
        int Jkl = kljJ[3];
        Orbit& oj = Z.modelspace->GetOrbit(j);
        Orbit& ok = Z.modelspace->GetOrbit(k);
        Orbit& ol = Z.modelspace->GetOrbit(l);
        double de_j = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
        double de_k = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
        double de_l = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
//        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        double occnat_l = ol.occ_nat;
        if ( (de_k + de_l + de_c) > Z.modelspace->GetdE3max() ) continue; 
        if ( (de_a + de_b + de_j) > Z.modelspace->GetdE3max() ) continue;
        if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        int twoJp_min = std::max( std::abs(2*Jab - oj.j2), std::abs(2*Jkl-j2c));
        int twoJp_max = std::min( 2*Jab + oj.j2, 2*Jkl+j2c);
        for (int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
        {
           kljJJ_needed.insert({k,l,c,(size_t)Jkl,(size_t)twoJp});
           kljJJ_needed.insert({a,b,j,(size_t)Jab,(size_t)twoJp});
        }// for twoJp
      }// for ind_klj
    }// for ind_abc

   // now compute the couplings once and store them in a hash table
   std::vector<double> recoupling_cache; // a vector storing all the relevant recoupling info
   std::unordered_map<size_t, size_t> recoupling_cache_lookup; // a lookup table so we can find the info easily

   for ( auto kljJJ : kljJJ_needed )
   {
     size_t k = kljJJ[0];
     size_t l = kljJJ[1];
     size_t j = kljJJ[2];
     size_t Jkl = kljJJ[3];
     size_t twoJp = kljJJ[4];
     size_t hash = Hash_comm232_key( kljJJ );

     std::vector<size_t> iketlist;
     std::vector<double> recouplelist;
     size_t ch_check = Y.ThreeBody.GetKetIndex_withRecoupling( Jkl, twoJp, k, l, j, iketlist, recouplelist) ;
     size_t listsize = iketlist.size();
     recoupling_cache_lookup[hash] = recoupling_cache.size();
     recoupling_cache.push_back(ch_check);
     recoupling_cache.push_back(listsize);
     for (size_t ik : iketlist) recoupling_cache.push_back(ik); // Store the size_t (an unsigned int) as a double. I hope this doesn't cause problems...
     for (double rec : recouplelist) recoupling_cache.push_back( rec);// Point to where this element lives in the vector
   }


    // Now fill the matrices
    for (size_t ind_abc=0; ind_abc<dim_abc; ind_abc++) // rows of X2MAT
    {

      auto& abcJ = klj_list_i[ abc_list[ind_abc] ];
      size_t a = abcJ[0];
      size_t b = abcJ[1];
      size_t c = abcJ[2];
      int Jab  = abcJ[3];
      Orbit& oa = Z.modelspace->GetOrbit(a);
      Orbit& ob = Z.modelspace->GetOrbit(b);
      Orbit& oc = Z.modelspace->GetOrbit(c);
      int j2c = oc.j2;
      double jc = 0.5 * j2c;
      double occ_abc = abc_occ_list[ind_abc]; // TODO: we can probably incorporate that hat factor with the occupation factor

      double de_a = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
      double de_b = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
      double de_c = std::abs( 2*oc.n + oc.l - e_fermi[oc.tz2]);
      double occnat_a = oa.occ_nat;
      double occnat_b = ob.occ_nat;
      double occnat_c = oc.occ_nat;

      // fill the 2 body part
      for (size_t ind_i=0; ind_i<dim_i; ind_i++)// ind_i indexes the list of sp states in this one body channel. columns of X2MAT
      {
        size_t i = obc_orbits[ind_i]; // i is the orbit index as per ModelSpace

        X2MAT(ind_i,ind_abc) = -sqrt( (2*Jab+1.)) * occ_abc*  X.TwoBody.GetTBME_J(Jab,c,i,a,b);
        Y2MAT(ind_i,ind_abc) = -sqrt( (2*Jab+1.)) * occ_abc*  Y.TwoBody.GetTBME_J(Jab,c,i,a,b);
      }// for ind_i


      // now the 3 body part
      for (size_t ind_klj=0; ind_klj<dim_klj; ind_klj++)
      {
        // <abi Jab twoJ | X | klc Jkl twoJ >
        auto& kljJ = klj_list_i[ ind_klj ];
        size_t k = kljJ[0];
        size_t l = kljJ[1];
        size_t j = kljJ[2];
        int Jkl = kljJ[3];
        Orbit& oj = Z.modelspace->GetOrbit(j);
        Orbit& ok = Z.modelspace->GetOrbit(k);
        Orbit& ol = Z.modelspace->GetOrbit(l);
        double ji = 0.5 * j2i;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;
        double jl = 0.5 * ol.j2;

        double de_j = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
        double de_k = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
        double de_l = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        double occnat_l = ol.occ_nat;
        // Are these checks redundant?
        if ( (de_a + de_b + de_j) > Z.modelspace->GetdE3max() ) continue;
        if ( (de_k + de_l + de_c) > Z.modelspace->GetdE3max() ) continue;
        if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        
        int twoJp_min = std::max( std::abs(2*Jab - oj.j2), std::abs(2*Jkl-j2c));
        int twoJp_max = std::min( 2*Jab + oj.j2, 2*Jkl+j2c);


        for (int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
        {
           double sixj = X.modelspace->GetSixJ(Jkl,jj,ji,  Jab, jc, 0.5*twoJp );

           std::array<size_t,5> abjJJ = { a,b,j,(size_t)Jab,(size_t)twoJp} ;
           std::array<size_t,5> klcJJ = { k,l,c,(size_t)Jkl,(size_t)twoJp} ;
           size_t hash_abj = Hash_comm232_key( abjJJ );
           size_t hash_klc = Hash_comm232_key( klcJJ );
           size_t pointer_abj = recoupling_cache_lookup.at(hash_abj);
           size_t pointer_klc = recoupling_cache_lookup.at(hash_klc);
           size_t ch_check = recoupling_cache[pointer_abj];
           size_t listsize_abj = recoupling_cache[pointer_abj+1];
           size_t listsize_klc = recoupling_cache[pointer_klc+1];
           
           double xabjklc = 0;
           double yabjklc = 0;
           for ( size_t ilist_abj=0; ilist_abj<listsize_abj; ilist_abj++)
           {
             size_t ibra_abj = (size_t) recoupling_cache[pointer_abj+2+ilist_abj];
             double recouple_bra =      recoupling_cache[pointer_abj+2+listsize_abj+ilist_abj];
             for ( size_t ilist_klc=0; ilist_klc<listsize_klc; ilist_klc++)
             {
               size_t iket_klc = (size_t) recoupling_cache[pointer_klc+2+ilist_klc];
               double recouple_ket =      recoupling_cache[pointer_klc+2+listsize_klc+ilist_klc];

//               xabjklc += recouple_bra*recouple_ket * X3.GetME_pn_ch(ch_check,ch_check, ibra_abj, iket_klc );
               if (x_has_3 ) xabjklc += recouple_bra*recouple_ket * X3.GetME_pn_ch(ch_check,ch_check, ibra_abj, iket_klc );
               yabjklc += recouple_bra*recouple_ket * Y3.GetME_pn_ch(ch_check,ch_check, ibra_abj, iket_klc );
             }
           }

           X3MAT(ind_abc, ind_klj) += (twoJp+1.)/sqrt(2*Jkl+1.) * sixj * xabjklc;
           Y3MAT(ind_abc, ind_klj) += (twoJp+1.)/sqrt(2*Jkl+1.) * sixj * yabjklc;

        }// for twoJp
      }// for ind_klj
    }// for ind_abc

    // now we do the mat mult
    ZMAT_list[obc_key] =  (  X2MAT * Y3MAT - Y2MAT * X3MAT  ) ;

  }// for iter_i in local one body channels

  // now we need to unpack all that mess and store it in Z
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for ( size_t ch=0; ch<nch; ch++)
  {
    TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    size_t nkets = tbc.GetNumberKets();
    size_t J = tbc.J;
    for (size_t ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      size_t i=bra.p;
      size_t j=bra.q;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);

      auto& lobc_i = local_one_body_channels.at({oi.j2,oi.l%2,oi.tz2});
      auto& lobc_j = local_one_body_channels.at({oj.j2,oj.l%2,oj.tz2});
      size_t ind_i = std::distance( lobc_i.begin(),  std::find( lobc_i.begin(),lobc_i.end(), i ) );
      size_t ind_j = std::distance( lobc_j.begin(),  std::find( lobc_j.begin(),lobc_j.end(), j ) );

      auto& list_i = klj_list.at({oi.j2,oi.l%2,oi.tz2});
      auto& list_j = klj_list.at({oj.j2,oj.l%2,oj.tz2});
      auto& ZMat_i = ZMAT_list.at({oi.j2,oi.l%2,oi.tz2});
      auto& ZMat_j = ZMAT_list.at({oj.j2,oj.l%2,oj.tz2});

      for (int iket=ibra; iket<nkets; iket++)
      {
        Ket& ket = tbc.GetKet(iket);
        size_t k = ket.p;
        size_t l = ket.q;
        Orbit& ok = Z.modelspace->GetOrbit(k);
        Orbit& ol = Z.modelspace->GetOrbit(l);
        int phase_ij = X.modelspace->phase( (oi.j2+oj.j2)/2-J);
        int phase_kl = X.modelspace->phase( (ok.j2+ol.j2)/2-J);

        auto& lobc_k = local_one_body_channels.at({ok.j2,ok.l%2,ok.tz2});
        auto& lobc_l = local_one_body_channels.at({ol.j2,ol.l%2,ol.tz2});
        size_t ind_k = std::distance( lobc_k.begin(),  std::find( lobc_k.begin(),lobc_k.end(), k ) );
        size_t ind_l = std::distance( lobc_l.begin(),  std::find( lobc_l.begin(),lobc_l.end(), l ) );

        auto& list_k = klj_list.at({ok.j2,ok.l%2,ok.tz2});
        auto& list_l = klj_list.at({ol.j2,ol.l%2,ol.tz2});
        auto& ZMat_k = ZMAT_list.at({ok.j2,ok.l%2,ok.tz2});
        auto& ZMat_l = ZMAT_list.at({ol.j2,ol.l%2,ol.tz2});

        std::array<size_t,4> key_klj = {k,l,j,J};
        std::array<size_t,4> key_kli = {k,l,i,J};
        std::array<size_t,4> key_ijk = {i,j,k,J};
        std::array<size_t,4> key_ijl = {i,j,l,J};

        size_t ind_klj = std::distance( list_i.begin(),  std::find( list_i.begin(),list_i.end(), key_klj ) );
        size_t ind_kli = std::distance( list_j.begin(),  std::find( list_j.begin(),list_j.end(), key_kli ) );
        size_t ind_ijl = std::distance( list_k.begin(),  std::find( list_k.begin(),list_k.end(), key_ijl ) );
        size_t ind_ijk = std::distance( list_l.begin(),  std::find( list_l.begin(),list_l.end(), key_ijk ) );

        if ( ind_klj == list_i.size()   or ind_kli==list_j.size() or ind_ijl==list_k.size() or ind_ijk==list_l.size() )   continue;


        double zijkl =            ZMat_j(ind_j, ind_kli) - phase_ij * ZMat_i(ind_i, ind_klj) 
                - hermX*hermY * ( ZMat_l(ind_l, ind_ijk) - phase_kl * ZMat_k(ind_k, ind_ijl) ) ;

        // normalize the tbme
        zijkl *= -1.0 / sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
        Z2.AddToTBME(ch,ch,ibra,iket,zijkl);
      }//for iket
    }//for ibra
  }// for ch

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}







void comm232ss_debug( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y2 = Y.TwoBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;

  bool x_has_3 = X3.IsAllocated();
  bool y_has_3 = Y3.IsAllocated();
//  bool x_has_3 = X3.is_allocated;
//  bool y_has_3 = Y3.is_allocated;

  int hermX = X.IsHermitian() ? 1 : -1;
  int hermY = Y.IsHermitian() ? 1 : -1;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  int nch = Z.modelspace->GetNumberTwoBodyChannels();

//  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (int ch=0; ch<nch; ch++)
  {
    auto& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();

//    double tstart = omp_get_wtime();
    // The strategy used here is the following. We reorganize <ic|X|ab> as <i|X|abc'>  and  <abj|Y|klc> as <abc'|Y|klj'>
    // where the prime indicates time reversal. Then we can cast things as a matrix multiplication
    //  <i|Z|klj'> = <i|X|abc'><abc'|Y|klj'>   and then we need to transform Z back to <ij|Z|kl>
    //
    //   i|   c|      i|                a| b| j|       a| b|   /c'
    //    |____|  ==>  |____             |__|__|  ==>   |__|__/
    //    |    |       |    |\           |  |  |        |  |  \
    //   a|   b|      a|   b| \c'       k| l| c|       k| l|   \j'
    //
    // The matrices for X and Z will clearly not be square. The left side is just a single particle orbit, while the right has 3 orbits.
    // Here, I determine which single orbits will be needed in this two-body channel.
    std::set<size_t> ij_orbits_set;
    std::set<int> j_jvals_set;
    for ( int ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      ij_orbits_set.insert( bra.p );
      ij_orbits_set.insert( bra.q );
      j_jvals_set.insert( bra.op->j2);
      j_jvals_set.insert( bra.oq->j2);
    }
    std::vector<size_t> ij_orbits;
    std::vector<int> j_jvals;
    std::map<size_t,size_t> ij_orbits_lookup;
    std::map<int,size_t> jvals_lookup;
    for ( auto i : ij_orbits_set )
    {
      ij_orbits.push_back(i); // list of all the indices an orbit (either i or j ) in bra could have
      ij_orbits_lookup[i] = ij_orbits.size()-1;
    }
    for ( auto j : j_jvals_set )
    {
       j_jvals.push_back(j);  // list of possible j values (angular momentum) an orbit in bra could have
       jvals_lookup[j] = j_jvals.size()-1;
    }
    size_t nij = ij_orbits.size();
    size_t njvals = j_jvals.size();
//    std::cout << "   ch = " << ch << "   njvals = " << njvals << std::endl;

    Z.profiler.timer["comm232_block1"] += omp_get_wtime() - tstart;
    tstart = omp_get_wtime();

    // next, we make a list of the abc' combinations that will inter into the sum
    std::vector<size_t> ch_ab_list;
    std::vector<size_t> iket_ab_list;
    std::vector<size_t> a_list;
    std::vector<size_t> b_list;
    std::vector<size_t> c_list;
    std::vector<double> occ_abc_list;

    for (int ch_ab=0; ch_ab<nch; ch_ab++)
    {
      auto& tbc_ab = X.modelspace->GetTwoBodyChannel(ch_ab);
      int Jab = tbc_ab.J;
      size_t nkets_ab = tbc_ab.GetNumberKets();
      for ( size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++ )
      {
        Ket& ket_ab = tbc_ab.GetKet(iket_ab);
        int a=ket_ab.p;
        int b=ket_ab.q;
        Orbit& oa = Z.modelspace->GetOrbit(a);
        Orbit& ob = Z.modelspace->GetOrbit(b);
        int ea = 2*oa.n + oa.l;
        int eb = 2*ob.n + ob.l;

        for (auto c : Z.modelspace->all_orbits )
        {
          Orbit& oc = Z.modelspace->GetOrbit(c);
          double jc = 0.5*oc.j2;
          int ec = 2*oc.n + oc.l;
//          if ( (std::abs( ea-e_fermi[oa.tz2]) + std::abs(eb-e_fermi[ob.tz2]) + std::abs(ec-e_fermi[oc.tz2])) > Z.modelspace->GetdE3max() ) continue;
          double occfactor = oa.occ * ob.occ * (1-oc.occ) + (1-oa.occ) * (1-ob.occ) * oc.occ;
          if ( std::abs(occfactor) < 1e-6 ) continue;
          if ( std::abs( oa.tz2 + ob.tz2 - oc.tz2)>2) continue;
          if (a==b) occfactor *=0.5;  // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b
//          if (  (std::abs( 2*Jab -oc.j2)>oj.j2)  or  ((2*Jab+oc.j2)<oj.j2) ) continue;
//          if (  (std::abs( 2*Jab -oc.j2)>std::max(oi.j2,oj.j2))  or  ((2*Jab+oc.j2)<std::min(oi.j2,oj.j2)) ) continue;
//          if ( tbc_ab.parity + oj.l + oc.l
          ch_ab_list.push_back(ch_ab);
          iket_ab_list.push_back(iket_ab);
          a_list.push_back(a);
          b_list.push_back(b);
          c_list.push_back(c);
          occ_abc_list.push_back(occfactor);
        }// for c
      }// for iket_ab
    } // for ch_ab

    Z.profiler.timer["comm232_block2"] += omp_get_wtime() - tstart;
    tstart = omp_get_wtime();

    // allocate the matrices
    size_t n_abc = ch_ab_list.size();

    arma::mat X2MAT( nij,   n_abc,             arma::fill::zeros);
    arma::mat Y3MAT( n_abc, nij*njvals*nkets,  arma::fill::zeros);
    arma::mat Y2MAT( nij,   n_abc,             arma::fill::zeros);
    arma::mat X3MAT( n_abc, nij*njvals*nkets,  arma::fill::zeros);



    // Now fill the X2 mat and Y3 mat
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)

      for (size_t ind_abc=0; ind_abc<n_abc; ind_abc++)
      {
        size_t a = a_list[ind_abc];
        size_t b = b_list[ind_abc];
        size_t c = c_list[ind_abc];
        int Jab = X.modelspace->GetTwoBodyChannel( ch_ab_list[ind_abc]).J;
        Orbit& oa = X.modelspace->GetOrbit(a);
        Orbit& ob = X.modelspace->GetOrbit(b);
        Orbit& oc = X.modelspace->GetOrbit(c);
        
        double jc = 0.5 * oc.j2;
        double d_ea = std::abs( 2*oa.n +oa.l - e_fermi[oa.tz2]);
        double d_eb = std::abs( 2*ob.n +ob.l - e_fermi[ob.tz2]);
        double d_ec = std::abs( 2*oc.n +oc.l - e_fermi[oc.tz2]);
        double occnat_a = oa.occ_nat;
        double occnat_b = ob.occ_nat;
        double occnat_c = oc.occ_nat;
    for ( size_t ind_i=0; ind_i<nij; ind_i++ )
    {
      size_t i = ij_orbits[ind_i];
//      if (ch==0) std::cout << "  ind_i = " << ind_i << "   i = " << i << std::endl;
      Orbit& oi = X.modelspace->GetOrbit(i);
      double ji = 0.5*oi.j2;
      double d_ei = std::abs( 2*oi.n +oi.l - e_fermi[oi.tz2]);
      double occnat_i = oi.occ_nat;


        if ( (d_ea+d_eb+d_ei) > Z.modelspace->GetdE3max() ) continue;

        X2MAT(ind_i,ind_abc) = -sqrt( (2*Jab+1.)) * occ_abc_list[ind_abc] *  X.TwoBody.GetTBME_J(Jab,c,i,a,b);
        Y2MAT(ind_i,ind_abc) = -sqrt( (2*Jab+1.)) * occ_abc_list[ind_abc] *  Y.TwoBody.GetTBME_J(Jab,c,i,a,b);
        

        if ( (oa.l+ob.l+oi.l+tbc.parity+oc.l)%2 != Y.parity) continue; //TODO: this will cause problems if we try something with Y.parity =1
        if ( std::abs((oa.tz2+ob.tz2+oi.tz2)-(tbc.Tz*2+oc.tz2)) > 2*Y.rank_T) continue;

        int twoJp_min = std::max( std::abs(2*Jab - oi.j2), std::abs(2*J-oc.j2));
        int twoJp_max = std::min( 2*Jab + oi.j2, 2*J+oc.j2);

        // Why is this part structured like this?
        // In tests, the time for this entire commutator routine was dominated by the time to access 3-body matrix elements.
        // That's because I'm asking for them in an order different from the one in which they're stored, which means recoupling.
        // Fortunately, both the X and Y block need exactly the same recoupling, and we can pull the recoupling for the bra
        // side a few loops out.
        for (int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
        {
           std::vector<size_t> ibra_list;
           std::vector<double> recouple_bra_list;
           size_t ch_check = Y.ThreeBody.GetKetIndex_withRecoupling( Jab, twoJp, a, b, i, ibra_list, recouple_bra_list) ;

            for ( size_t ind_jj=0; ind_jj<njvals; ind_jj++ )
            {
             int jj2 = j_jvals[ind_jj];
             double jj = 0.5*jj2;
             double sixj = X.modelspace->GetSixJ(J,ji,jj,  Jab, jc, 0.5*twoJp );
             if (std::abs(sixj)<1e-7) continue;
          // now loop over |kl> states
             for (int iket=0; iket<nkets; iket++)
             {
               size_t index_kli = (ind_jj+ind_i*njvals)*nkets + iket;
               Ket& ket = tbc.GetKet(iket);
               int k = ket.p;
               int l = ket.q;
//               if (i==0 and k==2 and l==2) std::cout << "i==0,l==2,k==2 => index kli = " <<index_kli << std::endl;
               double d_ek = std::abs( 2*ket.op->n +ket.op->l - e_fermi[ket.op->tz2]);
               double d_el = std::abs( 2*ket.oq->n +ket.oq->l - e_fermi[ket.oq->tz2]);
               double occnat_k = ket.op->occ_nat;
               double occnat_l = ket.oq->occ_nat;
               if ( (d_ek+d_el+d_ec) > Z.modelspace->GetdE3max() ) continue;
               if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
               std::vector<size_t> iket_list;
               std::vector<double> recouple_ket_list;
               ch_check = Y.ThreeBody.GetKetIndex_withRecoupling( J, twoJp, k, l, c, iket_list, recouple_ket_list) ;
      
               double xabiklc = 0;
               double yabiklc = 0;
               for (size_t I=0; I<ibra_list.size(); I++)
               {
                 for (size_t J=0; J<iket_list.size(); J++)
                 {
                   // I explicitly check if x and y have 3-body components because the call GetME_pn_ch goes straight to the data array
                   // without a safety net. If the 3-body structure isn't allocated, then bad things will happen.
                   if (x_has_3) xabiklc += recouple_bra_list[I]*recouple_ket_list[J] * X3.GetME_pn_ch(ch_check,ch_check, ibra_list[I], iket_list[J] );
                   if (y_has_3) yabiklc += recouple_bra_list[I]*recouple_ket_list[J] * Y3.GetME_pn_ch(ch_check,ch_check, ibra_list[I], iket_list[J] );
                 }
               }

//               X3MAT(ind_abc, index_kli) -= (twoJp+1) * sixj * xabiklc;
//               Y3MAT(ind_abc, index_kli) -= (twoJp+1) * sixj * yabiklc;
               X3MAT(ind_abc, index_kli) += (twoJp+1)/sqrt(2*J+1) * sixj * xabiklc;
               Y3MAT(ind_abc, index_kli) += (twoJp+1)/sqrt(2*J+1) * sixj * yabiklc;
//               if (ch==8 and ind_abc ==9 and index_kli==20)
//               {
//          std::cout << "ind_abc = " << ind_abc << "   ind_kli = " << index_kli << "    abiklc = " << a << " " << b << " " << i << " " << k << " " << l << " " << c << "   Jkl Jab twoJp = " << J << " " << Jab << " " << twoJp << " read y " << Y3.GetME_pn(Jab,J,twoJp,a,b,i,k,l,c) << "   yabiklc = " << yabiklc << " , " << Y3MAT(ind_abc, index_kli)  << std::endl;
//               }

             }// for iket
           }
        }// for ind_jj
      }// for ind_abc
          
    } // for ind_i

    Z.profiler.timer["comm232_block3"] += omp_get_wtime() - tstart;
    tstart = omp_get_wtime();
 /// now we're back out to the ch loop level.

    // finally do the matrix multiplication
//    arma::mat ZMat =  -sqrt(1./(2*J+1)) * (  X2MAT * Y3MAT - Y2MAT * X3MAT  ) ;
//    std::cout << "Y3MAT size : " << Y3MAT.n_cols << " x " << Y3MAT.n_rows << std::endl;
    arma::mat ZMat =  (  X2MAT * Y3MAT - Y2MAT * X3MAT  ) ;

//    if (ch==8)
//    {
//      std::cout << "ch = " << ch << "   and X2MAT = " << std::endl << X2MAT.row(5) << std::endl;
//      std::cout << "ch = " << ch << "   and Y3MAT = " << std::endl << Y3MAT.col(20) << std::endl;
//      std::cout << "   product = " << X2MAT.row(5) * Y3MAT.col(20) << std::endl;
//     double xysum =0;
//     for (size_t ii=0; ii<X2MAT.n_cols; ii++)
//     {
//       double z = X2MAT(5,ii) * Y3MAT(ii,20);
//       xysum += z;
//       std::cout << X2MAT(5,ii) << " " << Y3MAT(ii,20) << "  " << z << "  " << xysum << std::endl;
//     }
//      std::cout << "ch = " << ch << "   and Y2MAT = " << std::endl << Y2MAT.row(5) << std::endl;
//      std::cout << "ch = " << ch << "   and X3MAT = " << std::endl << X3MAT.col(20) << std::endl;
//      std::cout << "   product = " << Y2MAT.row(5) * X3MAT.col(20) << std::endl;
//    }

    Z.profiler.timer["comm232_block4"] += omp_get_wtime() - tstart;
    tstart = omp_get_wtime();

    // now convert back from <i|Z|klj'> to <ij|Z|kl>
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (int ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      int i=bra.p;
      int j=bra.q;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      double ji = 0.5*oi.j2;
      double jj = 0.5*oj.j2;
      size_t ind_i = ij_orbits_lookup[i];
      size_t ind_j = ij_orbits_lookup[j];
      size_t ind_ji = jvals_lookup[oi.j2];
      size_t ind_jj = jvals_lookup[oj.j2];

      for (int iket=ibra; iket<nkets; iket++)
      {
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        Orbit& ok = Z.modelspace->GetOrbit(k);
        Orbit& ol = Z.modelspace->GetOrbit(l);
        double jk = 0.5*ok.j2;
        double jl = 0.5*ol.j2;
        size_t ind_k = ij_orbits_lookup[k];
        size_t ind_l = ij_orbits_lookup[l];
        size_t ind_jk = jvals_lookup[ok.j2];
        size_t ind_jl = jvals_lookup[ol.j2];

        int phase_ij = X.modelspace->phase( (oi.j2+oj.j2)/2-J);
        int phase_kl = X.modelspace->phase( (ok.j2+ol.j2)/2-J);
        double zijkl = ZMat(ind_j,  (ind_jj+ind_i*njvals)*nkets + iket) ;
        zijkl -= phase_ij * ZMat(ind_i,  (ind_ji+ind_j*njvals)*nkets + iket) ;

        zijkl -= hermX*hermY* ZMat(ind_l, (ind_jl+ind_k*njvals)*nkets + ibra);
        zijkl += hermX*hermY*phase_kl * ZMat(ind_k, (ind_jk+ind_l*njvals)*nkets + ibra);

////       if (ch==0 and ibra<5 and iket<5)
//       if (ch==8 and i==0 and j==1 and ( (k==2 and l==5) or (k==5 and l==2) ) )
//       {
//         std::cout << " ibra,iket " << ibra << " " << iket << "   ijkl " << i << " " << j << " " << k << " " << l << std::endl;
//         std::cout << "    " << ZMat(ind_j,  (ind_jj+ind_i*njvals)*nkets + iket) << "  -  " <<  phase_ij * ZMat(ind_i,  (ind_ji+ind_j*njvals)*nkets + iket)
//                   << "  -  " << hermX*hermY* ZMat(ind_l, (ind_jl+ind_k*njvals)*nkets + ibra) << "  +   " << hermX*hermY*phase_kl * ZMat(ind_k, (ind_jk+ind_l*njvals)*nkets + ibra)
//                   << "  third term obtained from Z(" << ind_l << " , " << (ind_jl+ind_k*njvals)*nkets + ibra  << std::endl;
//       }


        // normalize the tbme
        zijkl *= -1.0 / sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
        Z2.AddToTBME(ch,ch,ibra,iket,zijkl);

      }//for iket
    }//for ibra

    Z.profiler.timer["comm232_block5"] += omp_get_wtime() - tstart;


  }// for ch

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}








// the old way that also works. It's slower but easier to read.
void comm232ss_slow( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y2 = Y.TwoBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;

  int nch = Z.modelspace->GetNumberTwoBodyChannels();

  std::vector<std::array<size_t,2>> channels;
  for (auto& iter : Z.TwoBody.MatEl ) channels.push_back(iter.first);
  size_t nchans = channels.size();
//  for (int ch=0; ch<nch; ch++)
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ich=0; ich<nchans; ich++)
  {
    size_t ch_bra = channels[ich][0];
    size_t ch_ket = channels[ich][1];
//    auto& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    auto& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    auto& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    int J = tbc_bra.J;
//    int nkets = tbc.GetNumberKets();
    int nbras = tbc_bra.GetNumberKets();
    int nkets = tbc_ket.GetNumberKets();
//    for (int ibra=0; ibra<nkets; ibra++)
    for (int ibra=0; ibra<nbras; ibra++)
    {
      Ket& bra = tbc_bra.GetKet(ibra);
      int i=bra.p;
      int j=bra.q;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double ji = 0.5*oi.j2;
      double jj = 0.5*oj.j2;
      int ket_min = ( ch_bra==ch_ket ) ? ibra : 0;
      for (int iket=ket_min; iket<nkets; iket++)
      {
        double zijkl = 0;
//        Ket& ket = tbc.GetKet(iket);
        Ket& ket = tbc_ket.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        Orbit& ok = Z.modelspace->GetOrbit(k);
        Orbit& ol = Z.modelspace->GetOrbit(l);
        double jk = 0.5*ok.j2;
        double jl = 0.5*ol.j2;
        double occnat_k = ok.occ_nat;
        double occnat_l = ol.occ_nat;
        for (auto c : Z.modelspace->all_orbits )
        {
          Orbit& oc = Z.modelspace->GetOrbit(c);
          double jc = 0.5*oc.j2;
          double occnat_c = oc.occ_nat;

          for (int ch_ab=0; ch_ab<nch; ch_ab++)
          {
            auto& tbc_ab = X.modelspace->GetTwoBodyChannel(ch_ab);
            int Jab = tbc_ab.J;
            size_t nkets_ab = tbc_ab.GetNumberKets();
            for ( size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++ )
            {
              Ket& ket_ab = tbc_ab.GetKet(iket_ab);
              int a=ket_ab.p;
              int b=ket_ab.q;
              Orbit& oa = Z.modelspace->GetOrbit(a);
              Orbit& ob = Z.modelspace->GetOrbit(b);
              double occnat_a = oa.occ_nat;
              double occnat_b = ob.occ_nat;
              double occfactor = oa.occ * ob.occ * (1-oc.occ) + (1-oa.occ) * (1-ob.occ) * oc.occ;
              if ( std::abs(occfactor) < 1e-6 ) continue;
              if (a==b) occfactor *=0.5;  // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b

              // Xicab term
//              if ( ( (oi.l+oc.l+tbc_ab.parity)%2==0) and ((oi.tz2+oc.tz2)==2*tbc_ab.Tz)
//                  and (std::abs(oi.j2-oc.j2)<=2*Jab)  and (oi.j2+oc.j2>=2*Jab) )
              if (  (  ( ((oi.l+oc.l+tbc_ab.parity)%2==X.parity) and (std::abs(oi.tz2+oc.tz2-2*tbc_ab.Tz)==2*X.rank_T ) )
                   or  ( ((oi.l+oc.l+tbc_ab.parity)%2==Y.parity) and (std::abs(oi.tz2+oc.tz2-2*tbc_ab.Tz)==2*Y.rank_T ) ) )
                  and (std::abs(oi.j2-oc.j2)<=2*Jab)  and (oi.j2+oc.j2>=2*Jab) )
              {

                if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                int twoJ_min = std::max( std::abs(oc.j2-2*J), std::abs( oj.j2-2*Jab ) );
                int twoJ_max = std::min( oc.j2+2*J,  oj.j2+2*Jab );
                double xciab = X2.GetTBME_J(Jab,c,i,a,b);
                double yciab = Y2.GetTBME_J(Jab,c,i,a,b);
//                double xciab = X2.GetTBME(ch_ab,c,i,a,b);
//                double yciab = Y2.GetTBME(ch_ab,c,i,a,b);
                int phasefactor = Z.modelspace->phase((oi.j2+oj.j2)/2-J);
                for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                {

                  double Jtot = 0.5 * twoJ;
                  double sixj = Z.modelspace->GetSixJ(jj,ji,J, jc,Jtot,Jab);
                  double hatfactor = (twoJ+1) * sqrt( (2*Jab+1.)/(2*J+1) );
                  double xabjklc = X3.GetME_pn(Jab,J,twoJ,a,b,j,k,l,c);
                  double yabjklc = Y3.GetME_pn(Jab,J,twoJ,a,b,j,k,l,c);
                  zijkl += occfactor * hatfactor * phasefactor * sixj * ( xciab * yabjklc  - yciab * xabjklc);
                }
              }


              // Xjcab term
//              if ( ( (oj.l+oc.l+tbc_ab.parity)%2==0) and ((oj.tz2+oc.tz2)==2*tbc_ab.Tz)
              if (  (  ( ((oj.l+oc.l+tbc_ab.parity)%2==X.parity) and (std::abs(oj.tz2+oc.tz2-2*tbc_ab.Tz)==2*X.rank_T ) )
                   or  ( ((oj.l+oc.l+tbc_ab.parity)%2==Y.parity) and (std::abs(oj.tz2+oc.tz2-2*tbc_ab.Tz)==2*Y.rank_T ) ) )
                  and (std::abs(oj.j2-oc.j2)<=2*Jab)  and (oj.j2+oc.j2>=2*Jab) )
              {
                if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                int twoJ_min = std::max( std::abs(oc.j2-2*J), std::abs( oi.j2-2*Jab ) );
                int twoJ_max = std::min( oc.j2+2*J,  oi.j2+2*Jab );
//                double xcjab = X2.GetTBME(ch_ab,c,j,a,b);
//                double ycjab = Y2.GetTBME(ch_ab,c,j,a,b);
                double xcjab = X2.GetTBME_J(Jab,c,j,a,b);
                double ycjab = Y2.GetTBME_J(Jab,c,j,a,b);
                int phasefactor = 1;
                for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                {
                  double Jtot = 0.5 * twoJ;
                  double sixj = Z.modelspace->GetSixJ(ji,jj,J, jc,Jtot,Jab);
                  double hatfactor = (twoJ+1) * sqrt( (2*Jab+1.)/(2*J+1) );
                  double xabiklc = X3.GetME_pn(Jab,J,twoJ,a,b,i,k,l,c);
                  double yabiklc = Y3.GetME_pn(Jab,J,twoJ,a,b,i,k,l,c);
                  zijkl -= occfactor * hatfactor * phasefactor * sixj * ( xcjab * yabiklc - ycjab * xabiklc);
                }
              }
 

              // Xabkc term
//              if ( ( (ok.l+oc.l+tbc_ab.parity)%2==0) and ((ok.tz2+oc.tz2)==2*tbc_ab.Tz)
              if (  (  ( ((ok.l+oc.l+tbc_ab.parity)%2==X.parity) and (std::abs(ok.tz2+oc.tz2-2*tbc_ab.Tz)==2*X.rank_T ) )
                   or  ( ((ok.l+oc.l+tbc_ab.parity)%2==Y.parity) and (std::abs(ok.tz2+oc.tz2-2*tbc_ab.Tz)==2*Y.rank_T ) ) )
                  and (std::abs(ok.j2-oc.j2)<=2*Jab)  and (ok.j2+oc.j2>=2*Jab) )
              {
                if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_l*(1-occnat_l) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                int twoJ_min = std::max( std::abs(oc.j2-2*J), std::abs( ol.j2-2*Jab ) );
                int twoJ_max = std::min( oc.j2+2*J,  ol.j2+2*Jab );
//                double xabck = X2.GetTBME(ch_ab,a,b,c,k);
//                double yabck = Y2.GetTBME(ch_ab,a,b,c,k);
                double xabck = X2.GetTBME_J(Jab,a,b,c,k);
                double yabck = Y2.GetTBME_J(Jab,a,b,c,k);
                int phasefactor = Z.modelspace->phase((ok.j2+ol.j2)/2-J);
                for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                {
                  double Jtot = 0.5 * twoJ;
                  double sixj = Z.modelspace->GetSixJ(jl,jk,J, jc,Jtot,Jab);
                  double hatfactor = (twoJ+1) * sqrt( (2*Jab+1.)/(2*J+1) );
                  double xijcabl = X3.GetME_pn(J,Jab,twoJ,i,j,c,a,b,l);
                  double yijcabl = Y3.GetME_pn(J,Jab,twoJ,i,j,c,a,b,l);
                  zijkl -= occfactor * hatfactor * phasefactor * sixj * ( yijcabl*xabck - xijcabl*yabck ); 
                }
              }


              // Xablc term
//              if ( ( (ol.l+oc.l+tbc_ab.parity)%2==0) and ((ol.tz2+oc.tz2)==2*tbc_ab.Tz)
              if (  (  ( ((ol.l+oc.l+tbc_ab.parity)%2==X.parity) and (std::abs(ol.tz2+oc.tz2-2*tbc_ab.Tz)==2*X.rank_T ) )
                   or  ( ((ol.l+oc.l+tbc_ab.parity)%2==Y.parity) and (std::abs(ol.tz2+oc.tz2-2*tbc_ab.Tz)==2*Y.rank_T ) ) )
                  and (std::abs(ol.j2-oc.j2)<=2*Jab)  and (ol.j2+oc.j2>=2*Jab) )
              {
                if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
                int twoJ_min = std::max( std::abs(oc.j2-2*J), std::abs( ok.j2-2*Jab ) );
                int twoJ_max = std::min( oc.j2+2*J,  ok.j2+2*Jab );
                double xabcl = X2.GetTBME_J(Jab,a,b,c,l);
                double yabcl = Y2.GetTBME_J(Jab,a,b,c,l);
//                double xabcl = X2.GetTBME(ch_ab,a,b,c,l);
//                double yabcl = Y2.GetTBME(ch_ab,a,b,c,l);
                int phasefactor = 1;
                for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                {
                  double Jtot = 0.5 * twoJ;
                  double sixj = Z.modelspace->GetSixJ(jk,jl,J, jc, Jtot,Jab);
                  double hatfactor = (twoJ+1) * sqrt( (2*Jab+1.)/(2*J+1) );
                  double xijcabk = X3.GetME_pn(J,Jab,twoJ,i,j,c,a,b,k);
                  double yijcabk = Y3.GetME_pn(J,Jab,twoJ,i,j,c,a,b,k);
                  zijkl += occfactor * hatfactor * phasefactor * sixj * ( yijcabk*xabcl - xijcabk*yabcl );
                }
              }

            }// for iket_ab
          }// for ch2
        }// for c

        // normalize the tbme
        zijkl *= -1.0 / sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
//        Z2.AddToTBME(ch,ch,ibra,iket,zijkl);
        Z2.AddToTBME(ch_bra,ch_ket,ibra,iket,zijkl);
      }// for iket
    }// for ibra
  }// for ch

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}







//*****************************************************************************************
//
// i|  j|      
//  |   |           Uncoupled expression:
//  *~~[X]~~~~*        Z_ijkl = 1/6 sum_abcd (nanbncn`d -n`an`bn`cnd)(X_ijdabc Y_abckld - Y_ijdabc Xabckld)
//  |   |     /\                                                                                                    
// a|  b|   c(  )d
//  |   |     \/ 
//  *~~[Y]~~~~*      Coupled expression:
//  |   |              Z_{ijkl}^{J} = 1/6 sum_abcd (nanbncn`d-n`an`bn`cnd) 1/(2J+1) sum_J1J' (2J'+1)(X_{ijdabc}^{J J1 J'} Y_{abckld}^{J1 J J'} - X<->Y )
// k|  l|                                                                                                                                        
//                                                 
//                                                                                                         
//  Checked with UnitTest and passed.
//
void comm332_ppph_hhhpss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;

  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();
  if ( not ( x3_allocated and y3_allocated) ) return;

  int hX = X.IsHermitian() ? 1 : -1;
  int hY = Y.IsHermitian() ? 1 : -1;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();
  
  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (int ch=0; ch<nch; ch++)
  {
    TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();

    std::vector<size_t> ch3_abc_list;
    std::vector<size_t> iket_abc_list;
    std::vector<size_t> d_list;
    std::vector<double> factor_list;
    
    // figure out how bit the abcd side of the matrices should be, and store some stuff for later lookup
    for (size_t ch_abc=0; ch_abc<nch3; ch_abc++)
    {
      ThreeBodyChannel& Tbc_abc = Z.modelspace->GetThreeBodyChannel(ch_abc);
      if ( std::abs( Tbc_abc.twoTz - 2*tbc.Tz) > 1) continue;
      int twoJ = Tbc_abc.twoJ;
      size_t nkets_abc = Tbc_abc.GetNumberKets();
      for (size_t iket_abc=0; iket_abc<nkets_abc; iket_abc++)
      {
        Ket3& ket_abc = Tbc_abc.GetKet(iket_abc);
        size_t a = ket_abc.p;
        size_t b = ket_abc.q;
        size_t c = ket_abc.r;
        double occ_abc = ket_abc.op->occ * ket_abc.oq->occ * ket_abc.oR->occ;
        double occ_abc_bar = (1-ket_abc.op->occ) * (1-ket_abc.oq->occ) * (1-ket_abc.oR->occ);
        double d_ea = std::abs( 2*ket_abc.op->n + ket_abc.op->l - e_fermi[ ket_abc.op->tz2]);
        double d_eb = std::abs( 2*ket_abc.oq->n + ket_abc.oq->l - e_fermi[ ket_abc.oq->tz2]);
        double d_ec = std::abs( 2*ket_abc.oR->n + ket_abc.oR->l - e_fermi[ ket_abc.oR->tz2]);
        double occnat_a = ket_abc.op->occ_nat;
        double occnat_b = ket_abc.oq->occ_nat;
        double occnat_c = ket_abc.oR->occ_nat;
        if ( (d_ea+d_eb+d_ec)>Z.modelspace->dE3max) continue;
        if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( (std::abs(occ_abc)==0) and (std::abs(occ_abc_bar)==0) ) continue;
        int Jab = ket_abc.Jpq;

        double symm_factor = 6;  // 6 possible orderings of abc. If a==b, then only 3 orderings, and if a==b==c, then only 1 ordering.
        if ( a==b and b==c )
           symm_factor = 1;
        else if (a==b or a==c or b==c )
           symm_factor = 3;

        for ( auto d : Z.modelspace->all_orbits ) 
        {
          Orbit& od = Z.modelspace->GetOrbit(d);
          double nd = od.occ;
//          double d_ed = std::abs( 2*od.n + od.l - e_fermi[od.tz2]);
//          double occnat_d = od.occ_nat;
          double occfactor = occ_abc*(1-nd) - occ_abc_bar*nd ;
          if (std::abs(occfactor)<1e-6) continue;
          if ( (std::abs( 2*J-od.j2) > twoJ) or (2*J+od.j2)<twoJ ) continue;
          if ( (od.l+tbc.parity + Tbc_abc.parity)%2>0) continue;
          if ( (Tbc_abc.twoTz)!=(od.tz2 + 2*tbc.Tz) ) continue;


          ch3_abc_list.push_back(ch_abc);
          iket_abc_list.push_back(iket_abc);
          d_list.push_back(d);
          factor_list.push_back( (twoJ+1.)/(2*J+1) * occfactor * symm_factor/6. );
         }// for d
      }//for iket_abc
    }// for ch_abc


    // initialize the matrices
    size_t dim_abcd = ch3_abc_list.size();
    arma::mat XMAT(nkets, dim_abcd, arma::fill::zeros);
    arma::mat YMAT(nkets, dim_abcd, arma::fill::zeros);


    // fill the matrices  TODO: I think we can make this faster by first doing the recoupling lookup for ijd
    for (size_t index_abcd=0; index_abcd<dim_abcd; index_abcd++)
    {
      size_t ch_abc = ch3_abc_list[index_abcd];
      size_t iket_abc = iket_abc_list[index_abcd];
      double factor = factor_list[index_abcd];
      ThreeBodyChannel& Tbc_abc = Z.modelspace->GetThreeBodyChannel(ch_abc);
      Ket3& ket_abc = Tbc_abc.GetKet(iket_abc);
      int Jab = ket_abc.Jpq;
      int twoJ = Tbc_abc.twoJ;
      size_t a = ket_abc.p;
      size_t b = ket_abc.q;
      size_t c = ket_abc.r;
      size_t d = d_list[index_abcd];
      Orbit& od = Z.modelspace->GetOrbit(d);
      double d_ed = std::abs(2*od.n + od.l - e_fermi[od.tz2]);
      double occnat_d = od.occ_nat;

      for (int ibra=0; ibra<nkets; ibra++)
      {
        Ket& bra = tbc.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        double d_ei = std::abs( 2*bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_ej = std::abs( 2*bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double occnat_i = bra.op->occ_nat;
        double occnat_j = bra.oq->occ_nat;
        if ( (d_ei+d_ej+d_ed)>Z.modelspace->dE3max) continue;
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_d*(1-occnat_d) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        double norm2b = (i==j) ? 1./sqrt(2) : 1;

        XMAT(ibra,index_abcd) = norm2b * X3.GetME_pn(J,Jab,twoJ, i,j,d,a,b,c);
        YMAT(ibra,index_abcd) = norm2b * factor * Y3.GetME_pn(J,Jab,twoJ, i,j,d,a,b,c);

        }// for ibra
      }// for index_abcd

      // do the mat mult
      arma::mat ZMAT = hY * XMAT * YMAT.t() - hX*YMAT * XMAT.t();

      // now unpack
      for (int ibra=0; ibra<nkets; ibra++)
      {
        for (int iket=ibra; iket<nkets; iket++)
        {
          Z2.AddToTBME(ch,ch,ibra,iket,ZMAT(ibra,iket) );
        }
      }
  }// for ch
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}


// the old slow way
/*
void comm332_ppph_hhhpss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;
  
  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (int ch=0; ch<nch; ch++)
  {
    TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      for (int iket=ibra; iket<nkets; iket++) 
      {
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;

        double zijkl = 0;

        for (size_t ch_abc=0; ch_abc<nch3; ch_abc++)
        {
          ThreeBodyChannel& Tbc_abc = Z.modelspace->GetThreeBodyChannel(ch_abc);
          if ( std::abs( Tbc_abc.twoTz - 2*tbc.Tz) > 1) continue;
          int twoJ = Tbc_abc.twoJ;
          size_t nkets_abc = Tbc_abc.GetNumberKets();
          for (size_t iket_abc=0; iket_abc<nkets_abc; iket_abc++)
          {
            Ket3& ket_abc = Tbc_abc.GetKet(iket_abc);
            size_t a = ket_abc.p;
            size_t b = ket_abc.q;
            size_t c = ket_abc.r;
            double occ_abc = ket_abc.op->occ * ket_abc.oq->occ * ket_abc.oR->occ;
            double occ_abc_bar = (1-ket_abc.op->occ) * (1-ket_abc.oq->occ) * (1-ket_abc.oR->occ);
            if ( (std::abs(occ_abc)==0) and (std::abs(occ_abc_bar)==0) ) continue;
            int Jab = ket_abc.Jpq;

              double symm_factor = 6;  // 6 possible orderings of abc. If a==b, then only 3 orderings, and if a==b==c, then only 1 ordering.
              if ( a==b and b==c )
                 symm_factor = 1;
              else if (a==b or a==c or b==c )
                 symm_factor = 3;

              for ( auto d : Z.modelspace->all_orbits ) 
              {
                Orbit& od = Z.modelspace->GetOrbit(d);
                double nd = od.occ;
                double occfactor = occ_abc*(1-nd) - occ_abc_bar*nd ;
                if (std::abs(occfactor)<1e-6) continue;
                if ( (std::abs( 2*J-od.j2) > twoJ) or (2*J+od.j2)<twoJ ) continue;
                if ( (od.l+tbc.parity + Tbc_abc.parity)%2>0) continue;
                if ( (Tbc_abc.twoTz)!=(od.tz2 + 2*tbc.Tz) ) continue;
                zijkl += symm_factor * 1./6  * (twoJ+1.)/(2*J+1) * occfactor 
                         *(  X3.GetME_pn(J,Jab,twoJ, i,j,d,a,b,c) * Y3.GetME_pn(Jab, J, twoJ, a,b,c,k,l,d)
                           - Y3.GetME_pn(J,Jab,twoJ, i,j,d,a,b,c) * X3.GetME_pn(Jab, J, twoJ, a,b,c,k,l,d) );

              }// for d
          }// for iket_abc
        }// for ch_abc
        zijkl /=  sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
        Z2.AddToTBME(ch,ch,ibra,iket,zijkl);
      }// for iket
    }// for ibra
  }// for ch
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
*/






//*****************************************************************************************
//
//           i| j|
//            |  |  Uncoupled expression:
//   *~~~[X]~~*  |            Z_ijkl = 1/6 sum_abcd (nanbncn`d -n`an`bn`cnd)(X_ijdabc Y_abckld - Y_ijdabc Xabckld)
//  / \  / \  |  |                  + 1/4 (1-Pij)(1-Pkl) sum_{abcd} nanbn`cn`d ( X_abicdk Y_cdjabl - Y_abicdk X_cdjabl )
// (a c)(b d) |  |
//  \ /  \ /  |  |     
//   *~~~[Y]~~|~~*       
//            |  |   Coupled expression:
//           k| l|     Z_{ijkl}^{J} = 1/6 sum_abcd (nanbncn`d-n`an`bn`cnd) 1/(2J+1) sum_J1J' (2J'+1)(X_{ijdabc}^{J J1 J'} Y_{abckld}^{J1 J J'} - X<->Y )
//                                    + 1/4  (1-Pij^J)(1-Pkl^J) sum_{abcd} nanbn`cn`d sum_{J1 J2 J' J"} (2J'+1)(2J"+1)  (-1)^{2J+J2+J'+J"+2j-i-k}
//                                      { p  q  J }
//                                    * { J1 J" s } ( X_{abicdk}^{J1 J2 J'} Y_{cdjabl}^{J2 J1 J"} - X<->Y )
//                                      { J' J2 r }
//
//  We recouple the expression so that it factorizes and can be cast as matrix multiplication
//
//            i| /k`
//             |/
//    *~~~[X]~~*
//   / \  / \
//  (a c)(b d)
//   \ /  \ /
//    *~~~[Y]~~~*
//              | \
//             l|  \j'
//              
//        Tested with UnitTest and passed.                                                               
//
void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
//  double t_internal = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();
//  std::cout << "  fermi levels " << e_fermi[-1] << " " << e_fermi[+1] << std::endl;
 
  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();
  if ( not ( x3_allocated and y3_allocated) ) return;

   size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
   size_t nch_CC = Z.modelspace->TwoBodyChannels_CC.size();

   std::deque<arma::mat> Z_bar ( nch_CC );
   std::deque<arma::mat> Z_bar_flipbra ( nch_CC );
   std::deque<arma::mat> Z_bar_flipket ( nch_CC );
   std::deque<arma::mat> Z_bar_flipboth ( nch_CC );

   // before entering the parallel loop, lets allocate the Pandya-transformed matrices Zbar.
   // we keep 4 because we want i<j and i>j and k<l and k>l
   for (size_t ch=0;ch<nch_CC;++ch)
   {
      TwoBodyChannel_CC& tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_CC.GetNumberKets();

      Z_bar[ch].zeros( nKets_cc, nKets_cc );
      Z_bar_flipbra[ch].zeros( nKets_cc, nKets_cc );
      Z_bar_flipket[ch].zeros( nKets_cc, nKets_cc );
      Z_bar_flipboth[ch].zeros( nKets_cc, nKets_cc );
   }

   double occnat_factor_max = 0;
   for ( auto i : Z.modelspace->all_orbits)
   {
      double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
      occnat_factor_max = std::max( occnat_factor_max, occnat_i*(1-occnat_i) );
   }

  // Loop through Pandya-transformed channels and compute the matrix Zbar by mat mult
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
   for (size_t ch=0;ch<nch_CC;++ch)
   {
      TwoBodyChannel_CC& tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_CC.GetNumberKets();

      int J_ph = tbc_CC.J;
      int parity_ph = tbc_CC.parity;
      int Tz_ph = tbc_CC.Tz;

      size_t n_abcd_states = 0;

      // count the abcd states that can act in this channel
      for (int ch_ab=0; ch_ab<nch; ch_ab++)
      {
        TwoBodyChannel& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
        int Jab = tbc_ab.J;
        int nkets_ab = tbc_ab.GetNumberKets();

        for (int ch_cd=0; ch_cd<nch; ch_cd++)
        {
          TwoBodyChannel& tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
          int Jcd = tbc_cd.J;
          int nkets_cd = tbc_cd.GetNumberKets();
          if (  std::abs(Jab-Jcd)>J_ph or (Jab+Jcd)<J_ph ) continue;
          if ( (tbc_ab.parity + tbc_cd.parity + tbc_CC.parity)%2 >0 ) continue;
          if ( std::abs( tbc_cd.Tz - tbc_ab.Tz) != tbc_CC.Tz ) continue;
          for (int iket_ab=0; iket_ab<nkets_ab; iket_ab++)
          {
            Ket& ket_ab = tbc_ab.GetKet(iket_ab);
            int a = ket_ab.p;
            int b = ket_ab.q;
            double na = ket_ab.op->occ;
            double nb = ket_ab.oq->occ;
            double occnat_a = ket_ab.op->occ_nat;
            double occnat_b = ket_ab.oq->occ_nat;
            double d_ea = std::abs(2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
            double d_eb = std::abs(2*ket_ab.oq->n + ket_ab.oq->l - e_fermi[ket_ab.oq->tz2]);
            if ( na*nb<1e-8 and (1-na)*(1-nb)<1e-8 ) continue;
            if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue ;
            if ( (d_ea + d_eb) > Z.modelspace->GetdE3max() ) continue;

            for (int iket_cd=0; iket_cd<nkets_cd; iket_cd++)
            {
              Ket& ket_cd = tbc_cd.GetKet(iket_cd);
              int c = ket_cd.p;
              int d = ket_cd.q;
              double nc = ket_cd.op->occ;
              double nd = ket_cd.oq->occ;
              double occnat_c = ket_cd.op->occ_nat;
              double occnat_d = ket_cd.oq->occ_nat;
              double d_ec = std::abs(2*ket_cd.op->n + ket_cd.op->l - e_fermi[ket_cd.op->tz2]);
              double d_ed = std::abs(2*ket_cd.oq->n + ket_cd.oq->l - e_fermi[ket_cd.oq->tz2]);
              double occupation_factor = (1-na)*(1-nb)*nc*nd - na*nb*(1-nc)*(1-nd);
              if (std::abs(occupation_factor)<1e-6) continue;
              if ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue ;
              if ( (d_ec + d_ed) > Z.modelspace->GetdE3max() ) continue;

              n_abcd_states++;

            }// for iket_cd
          }// for iket_ab
        }// for ch_cd
      }// for ch_ab
      

      // Now that we know how big they are, allocate the matrices
      arma::mat X3_ij( nKets_cc, n_abcd_states, arma::fill::zeros );
      arma::mat X3_ji( nKets_cc, n_abcd_states, arma::fill::zeros );
      arma::mat Y3_ij( n_abcd_states, nKets_cc, arma::fill::zeros );
      arma::mat Y3_ji( n_abcd_states, nKets_cc, arma::fill::zeros );


      // Fill the matrices, which are basically <ij`| X | cd a`b`> type things (the ticks ` mean time-reversed state)
      for (size_t ibra_CC=0; ibra_CC<nKets_cc; ibra_CC++)
      {
        Ket& bra_CC = tbc_CC.GetKet(ibra_CC);
        int i = bra_CC.p;
        int j = bra_CC.q;
        Orbit& oi = Z.modelspace->GetOrbit(i);
        Orbit& oj = Z.modelspace->GetOrbit(j);


        int j2i = oi.j2;
        int j2j = oj.j2;
        double ji = 0.5*j2i;
        double jj = 0.5*j2j;
        int phase_ij = Z.modelspace->phase( (j2i+j2j)/2);

        double d_ei = std::abs(2*oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2*oj.n + oj.l - e_fermi[oj.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;

        size_t ind_abcd = -1;
        for (int ch_ab=0; ch_ab<nch; ch_ab++)
        {
          TwoBodyChannel& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
          int Jab = tbc_ab.J;
          int nkets_ab = tbc_ab.GetNumberKets();


          for (int ch_cd=0; ch_cd<nch; ch_cd++)
          {
            TwoBodyChannel& tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
            int Jcd = tbc_cd.J;
            int nkets_cd = tbc_cd.GetNumberKets();
            if (  std::abs(Jab-Jcd)>J_ph or (Jab+Jcd)<J_ph ) continue;
            if ( (tbc_ab.parity + tbc_cd.parity + tbc_CC.parity)%2 >0 ) continue;
            if ( std::abs( tbc_cd.Tz - tbc_ab.Tz) != tbc_CC.Tz ) continue;


            int twoJp_min = std::min( {std::abs(2*Jab-j2i), std::abs(2*Jcd-j2j), std::abs(2*Jab-j2j), std::abs(2*Jcd-j2i) });
            int twoJp_max = std::max( {(2*Jab+j2i), (2*Jcd+j2j), (2*Jab+j2j), (2*Jcd+j2i) });
            std::vector<double> sixj_ij_list;
            std::vector<double> sixj_ji_list;
            for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
            {
              sixj_ij_list.push_back( (twoJp+1) * Z.modelspace->GetSixJ(ji,jj,J_ph, Jcd,Jab,0.5*twoJp) * Z.modelspace->phase((j2i+twoJp)/2) );
              sixj_ji_list.push_back( (twoJp+1) * Z.modelspace->GetSixJ(jj,ji,J_ph, Jcd,Jab,0.5*twoJp) * Z.modelspace->phase((j2j+twoJp)/2) );
            }

            bool abi_cdj_Tz_ok = (  (2*tbc_ab.Tz + oi.tz2) == (2*tbc_cd.Tz + oj.tz2) );
            bool abj_cdi_Tz_ok = (  (2*tbc_ab.Tz + oj.tz2) == (2*tbc_cd.Tz + oi.tz2) );


            for (int iket_ab=0; iket_ab<nkets_ab; iket_ab++)
            {
              Ket& ket_ab = tbc_ab.GetKet(iket_ab);
              int a = ket_ab.p;
              int b = ket_ab.q;
              double na = ket_ab.op->occ;
              double nb = ket_ab.oq->occ;
              if ( na*nb<1e-8 and (1-na)*(1-nb)<1e-8 ) continue;
              double d_ea = std::abs(2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
              double d_eb = std::abs(2*ket_ab.oq->n + ket_ab.oq->l - e_fermi[ket_ab.oq->tz2]);
              double occnat_a = ket_ab.op->occ_nat;
              double occnat_b = ket_ab.oq->occ_nat;
              bool keep_abi = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) >= Z.modelspace->GetOccNat3Cut() ) ;
              bool keep_abj = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) >= Z.modelspace->GetOccNat3Cut() ) ;
              keep_abi = keep_abi and (d_ea+d_eb+d_ei <= Z.modelspace->dE3max );
              keep_abj = keep_abj and (d_ea+d_eb+d_ej <= Z.modelspace->dE3max );

              if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue ;
              if ( (d_ea + d_eb) > Z.modelspace->GetdE3max() ) continue;

              std::vector<size_t> ch_abi_list;
              std::vector<size_t> ch_abj_list;
              std::vector<std::vector<size_t>> kets_abi_list;
              std::vector<std::vector<size_t>> kets_abj_list;
              std::vector<std::vector<double>> recouple_abi_list;
              std::vector<std::vector<double>> recouple_abj_list;

              if ( keep_abi or keep_abj )
              {
                for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                {
                  std::vector<size_t> ketlist;
                  std::vector<double> reclist;
                  size_t ch_check = -1;
                  if ( keep_abi  and abi_cdj_Tz_ok and (std::abs(2*Jab-oi.j2)<=twoJp) and (2*Jab+oi.j2 >=twoJp) )
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling( Jab, twoJp, a, b, i, ketlist, reclist) ;
                  }
                  ch_abi_list.push_back( ch_check );
                  kets_abi_list.push_back( ketlist );
                  recouple_abi_list.push_back( reclist );
                  ketlist.clear();
                  reclist.clear();
                  ch_check = -1;
                  if ( keep_abj  and abj_cdi_Tz_ok  and (std::abs(2*Jab-oj.j2)<=twoJp) and (2*Jab+oj.j2 >=twoJp) )
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling( Jab, twoJp, a, b, j, ketlist, reclist) ;
                  }
                  ch_abj_list.push_back( ch_check );
                  kets_abj_list.push_back( ketlist );
                  recouple_abj_list.push_back( reclist );
                }
              }


              for (int iket_cd=0; iket_cd<nkets_cd; iket_cd++)
              {
                Ket& ket_cd = tbc_cd.GetKet(iket_cd);
                int c = ket_cd.p;
                int d = ket_cd.q;
                double nc = ket_cd.op->occ;
                double nd = ket_cd.oq->occ;
                double d_ec = std::abs(2*ket_cd.op->n + ket_cd.op->l - e_fermi[ket_cd.op->tz2]);
                double d_ed = std::abs(2*ket_cd.oq->n + ket_cd.oq->l - e_fermi[ket_cd.oq->tz2]);
                double occnat_c = ket_cd.op->occ_nat;
                double occnat_d = ket_cd.oq->occ_nat;
                double occupation_factor = (1-na)*(1-nb)*nc*nd - na*nb*(1-nc)*(1-nd);
                if (std::abs(occupation_factor)<1e-6) continue;

                if ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue ;
                if ( (d_ec + d_ed) > Z.modelspace->GetdE3max() ) continue;

                double symmetry_factor = 1;  // we only sum a<=b and c<=d, so we undercount by a factor of 4, canceling the 1/4 in the formula 
                if (a==b) symmetry_factor *= 0.5; // if a==b or c==d, then the permutation doesn't give a new state, so there's less undercounting
                if (c==d) symmetry_factor *= 0.5;

                ind_abcd++; // increment the matrix index.

                if ( not ( keep_abi or keep_abj) ) continue; // need this after the ind_abcd++... maybe be more clever?

                bool keep_cdi = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_i*(1-occnat_i) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                bool keep_cdj = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_j*(1-occnat_j) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                keep_cdi = keep_cdi and (d_ec+d_ed+d_ei <= Z.modelspace->dE3max );
                keep_cdj = keep_cdj and (d_ec+d_ed+d_ej <= Z.modelspace->dE3max );

                if ( not ( keep_cdi or keep_cdj) ) continue;


                std::vector<size_t> ch_cdi_list;
                std::vector<size_t> ch_cdj_list;
                std::vector<std::vector<size_t>> kets_cdi_list;
                std::vector<std::vector<size_t>> kets_cdj_list;
                std::vector<std::vector<double>> recouple_cdi_list;
                std::vector<std::vector<double>> recouple_cdj_list;

                for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                {
                  std::vector<size_t> ketlist;
                  std::vector<double> reclist;
                  size_t ch_check=-1;
                  if ( keep_abj and keep_cdi and abj_cdi_Tz_ok  and (std::abs(2*Jcd-oi.j2)<=twoJp) and (2*Jcd+oi.j2 >=twoJp) )
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling( Jcd, twoJp, c, d, i, ketlist, reclist) ;
                  }
                  ch_cdi_list.push_back( ch_check );
                  kets_cdi_list.push_back( ketlist );
                  recouple_cdi_list.push_back( reclist );
                  ketlist.clear();
                  reclist.clear();
                  ch_check = -1;
                  if ( keep_abi and keep_cdj and abi_cdj_Tz_ok  and (std::abs(2*Jcd-oj.j2)<=twoJp) and (2*Jcd+oj.j2 >=twoJp) )
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling( Jcd, twoJp, c, d, j, ketlist, reclist) ;
                  }
                  ch_cdj_list.push_back( ch_check );
                  kets_cdj_list.push_back( ketlist );
                  recouple_cdj_list.push_back( reclist );
                }// for twoJp

                for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                {
                  double sixj_ij = sixj_ij_list[(twoJp-twoJp_min)/2]; // <-- these also include a phase and a (twoJp+1)
                  double sixj_ji = sixj_ji_list[(twoJp-twoJp_min)/2];

                  if ( keep_abi and keep_cdj and abi_cdj_Tz_ok   )
                  {
                    size_t ch_abi = ch_abi_list[(twoJp-twoJp_min)/2];
                    auto& kets_abi = kets_abi_list[(twoJp-twoJp_min)/2];
                    auto& kets_cdj = kets_cdj_list[(twoJp-twoJp_min)/2];
                    auto& recouple_abi = recouple_abi_list[(twoJp-twoJp_min)/2];
                    auto& recouple_cdj = recouple_cdj_list[(twoJp-twoJp_min)/2];
                    double xabicdj = 0;
                    double ycdjabi = 0;
                    for ( size_t Iabi=0; Iabi<kets_abi.size(); Iabi++)
                    {
                      for ( size_t Icdj=0; Icdj<kets_cdj.size(); Icdj++)
                      {
                         xabicdj += recouple_abi[Iabi] * recouple_cdj[Icdj] * X3.GetME_pn_ch(ch_abi,ch_abi, kets_abi[Iabi], kets_cdj[Icdj] );
                         ycdjabi += recouple_abi[Iabi] * recouple_cdj[Icdj] * Y3.GetME_pn_ch(ch_abi,ch_abi, kets_cdj[Icdj], kets_abi[Iabi] );
                      }
                    }
                    X3_ij(ibra_CC, ind_abcd) +=  sixj_ij * xabicdj;
                    Y3_ij(ind_abcd, ibra_CC) +=  sixj_ij * ycdjabi * occupation_factor * symmetry_factor ;
                  }
                  
                  if (keep_abj and keep_cdi and abj_cdi_Tz_ok)
                  {
                    size_t ch_abj = ch_abj_list[(twoJp-twoJp_min)/2];
                    auto& kets_abj = kets_abj_list[(twoJp-twoJp_min)/2];
                    auto& kets_cdi = kets_cdi_list[(twoJp-twoJp_min)/2];
                    auto& recouple_abj = recouple_abj_list[(twoJp-twoJp_min)/2];
                    auto& recouple_cdi = recouple_cdi_list[(twoJp-twoJp_min)/2];
                    double xabjcdi = 0;
                    double ycdiabj = 0;
                    for ( size_t Iabj=0; Iabj<kets_abj.size(); Iabj++)
                    {
                      for ( size_t Icdi=0; Icdi<kets_cdi.size(); Icdi++)
                      {
                         xabjcdi += recouple_abj[Iabj] * recouple_cdi[Icdi] * X3.GetME_pn_ch(ch_abj,ch_abj, kets_abj[Iabj], kets_cdi[Icdi] );
                         ycdiabj += recouple_abj[Iabj] * recouple_cdi[Icdi] * Y3.GetME_pn_ch(ch_abj,ch_abj, kets_cdi[Icdi], kets_abj[Iabj] );
                      }
                    }
                    X3_ji(ibra_CC, ind_abcd) +=  sixj_ji * xabjcdi;
                    Y3_ji(ind_abcd, ibra_CC) +=  sixj_ji * ycdiabj * occupation_factor * symmetry_factor;
                  }

                }// for twoJp

              }// for iket_cd
            }// for iket_ab
          }// for ch_cd
        }// for ch_ab

      }// for ibra_CC

      Z_bar[ch] = X3_ij * Y3_ij;
      Z_bar_flipbra[ch] = X3_ji * Y3_ij;
      Z_bar_flipket[ch] = X3_ij * Y3_ji;
      Z_bar_flipboth[ch] = X3_ji * Y3_ji;

   }// for ch  (ph-coupled 2-body channels)

   
//  Z.profiler.timer["comm332_CC_loops"] += omp_get_wtime() - t_internal;
//  t_internal = omp_get_wtime();

   // Now transform Zbar to Z

   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
   for (size_t ch=0; ch<nch; ch++)
   {
     TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
     int J = tbc.J;
     size_t nkets = tbc.GetNumberKets();
     for (size_t ibra=0; ibra<nkets; ibra++)
     {
       Ket& bra = tbc.GetKet(ibra);
       size_t i = bra.p;
       size_t j = bra.q;
       Orbit& oi = Z.modelspace->GetOrbit(i);
       Orbit& oj = Z.modelspace->GetOrbit(j);
       double ji = 0.5*oi.j2;
       double jj = 0.5*oj.j2;
       int phase_ij = Z.modelspace->phase( (oi.j2+oj.j2)/2-J);
       for (size_t iket=ibra; iket<nkets; iket++)
       {
         Ket& ket = tbc.GetKet(iket);
         size_t k = ket.p;
         size_t l = ket.q;
         Orbit& ok = Z.modelspace->GetOrbit(k);
         Orbit& ol = Z.modelspace->GetOrbit(l);
         double jk = 0.5*ok.j2;
         double jl = 0.5*ol.j2;
         int phase_kl = Z.modelspace->phase( (ok.j2+ol.j2)/2-J);

         double zijkl = 0;

         int parity_cc = (oi.l+ol.l)%2;
         int Tz_cc = std::abs(oi.tz2-ol.tz2)/2;
         int J_ph_min = std::max( std::abs(oi.j2-ol.j2), std::abs(ok.j2-oj.j2))/2;
         int J_ph_max = std::min( (oi.j2+ol.j2), (ok.j2+oj.j2))/2;
         for (int J_ph=J_ph_min; J_ph<=J_ph_max; J_ph++)
         {
           int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(J_ph,parity_cc,Tz_cc);
           TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
           int indx_il = tbc_cc.GetLocalIndex(std::min(i,l),std::max(i,l)) ;
           int indx_kj = tbc_cc.GetLocalIndex(std::min(k,j),std::max(k,j)) ;
           double sixj = Z.modelspace->GetSixJ(ji,jj,J, jk,jl,J_ph);
           double zbar_ilkj;
           double zbar_jkli;
           if (i<=l and k<=j) // il kj
           {
              zbar_ilkj = Z_bar[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar_flipboth[ch_cc](indx_kj,indx_il);
           }
           else if ( i>l and k<=j) /// li kj
           {
              zbar_ilkj = Z_bar_flipbra[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar_flipbra[ch_cc](indx_kj,indx_il);
           }
           else if ( i<=l and k>j) // il jk
           {
              zbar_ilkj = Z_bar_flipket[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar_flipket[ch_cc](indx_kj,indx_il);
           }
           else  // li jk
           {
              zbar_ilkj = Z_bar_flipboth[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar[ch_cc](indx_kj,indx_il);
           }

           zijkl += (2*J_ph+1)*sixj * ( zbar_ilkj + phase_ij*phase_kl * zbar_jkli);
         }

         parity_cc = (oj.l+ol.l)%2;
         Tz_cc = std::abs(oj.tz2-ol.tz2)/2;
         J_ph_min = std::max( std::abs(oj.j2-ol.j2), std::abs(ok.j2-oi.j2))/2;
         J_ph_max = std::min( (oj.j2+ol.j2), (ok.j2+oi.j2))/2;
         for (int J_ph=J_ph_min; J_ph<=J_ph_max; J_ph++)
         {
           int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(J_ph,parity_cc,Tz_cc);
           TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
           int indx_jl = tbc_cc.GetLocalIndex(std::min(j,l),std::max(j,l)) ;
           int indx_ki = tbc_cc.GetLocalIndex(std::min(k,i),std::max(k,i)) ;
           double sixj = Z.modelspace->GetSixJ(jj,ji,J, jk,jl,J_ph);
           double zbar_jlki;
           double zbar_iklj;

           if ( j<=l and k<=i )
           {
             zbar_jlki = Z_bar[ch_cc](indx_jl, indx_ki); // <-- this one
             zbar_iklj = Z_bar_flipboth[ch_cc](indx_ki, indx_jl);
           }
           else if ( j<=l and k>i)
           {
             zbar_jlki = Z_bar_flipket[ch_cc](indx_jl, indx_ki);
             zbar_iklj = Z_bar_flipket[ch_cc](indx_ki, indx_jl);
           }
           else if ( j>l and k<=i)
           {
             zbar_jlki = Z_bar_flipbra[ch_cc](indx_jl, indx_ki);
             zbar_iklj = Z_bar_flipbra[ch_cc](indx_ki, indx_jl);
           }
           else if ( j>l and k>i)
           {
             zbar_jlki = Z_bar_flipboth[ch_cc](indx_jl, indx_ki);
             zbar_iklj = Z_bar[ch_cc](indx_ki, indx_jl);
           }

           zijkl -= (2*J_ph+1)*sixj * ( phase_ij * zbar_jlki + phase_kl * zbar_iklj);

         }       

        // make it a normalized TBME
        zijkl /=  sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
        // the AddToTBME routine automatically takes care of the hermitian conjugate as well
        Z2.AddToTBME(ch,ch,ibra,iket,zijkl);
 
       }// for iket
     }// for ibra
   }// for ch

//  Z.profiler.timer["comm332_main_loop"] += omp_get_wtime() - t_internal;

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}



void comm332_pphhss_debug( const Operator& X, const Operator& Y, Operator& Z )
//void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z )
{
//  std::cout << "================== BEGIN " << __func__ << " ===================" << std::endl;
  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();


  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();
  if ( not ( x3_allocated and y3_allocated) ) return;
  
  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (int ch=0; ch<nch; ch++)
  {
    TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      int ji2 = bra.op->j2;
      int jj2 = bra.oq->j2;
      double ji = 0.5*ji2;
      double jj = 0.5*jj2;
      double d_ei = std::abs(2*bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
      double d_ej = std::abs(2*bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
      double occnat_i = bra.op->occ_nat;
      double occnat_j = bra.oq->occ_nat;
      for (int iket=ibra; iket<nkets; iket++) 
      {
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        Orbit& ok = Z.modelspace->GetOrbit(k);
        Orbit& ol = Z.modelspace->GetOrbit(l);
        int jk2 = ket.op->j2;
        int jl2 = ket.oq->j2;
        double jk = 0.5*jk2;
        double jl = 0.5*jl2;
        double d_ek = std::abs(2*ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
        double d_el = std::abs(2*ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
        double occnat_k = ket.op->occ_nat;
        double occnat_l = ket.oq->occ_nat;

        int phase_ij = Z.modelspace->phase( (ji2+jj2)/2-J);
        int phase_kl = Z.modelspace->phase( (jk2+jl2)/2-J);
        double zijkl = 0;

        int J_ph_min = std::min( { std::abs(ji2 - jl2), std::abs(jj2-jk2), std::abs(jj2 - jl2), std::abs(ji2 - jk2)} )/2;
        int J_ph_max = std::max( {ji2+jl2, jj2+jk2, jj2+jl2, ji2+jk2} )/2;
//        int J_ph_min = std::max( std::abs(ji2 - jl2), std::abs(jj2-jk2) )/2;
//        int J_ph_max = std::min( ji2+jl2, jj2+jk2 )/2;

        for (int J_ph=J_ph_min; J_ph<=J_ph_max; J_ph++)
        {
          double zbar_ilkj=0;
          double zbar_jlki=0;
          double zbar_iklj=0;
          double zbar_jkli=0;

          for (int ch_ab=0; ch_ab<nch; ch_ab++)
          {
            TwoBodyChannel& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
            int Jab = tbc_ab.J;
            int nkets_ab = tbc_ab.GetNumberKets();

            for (int ch_cd=0; ch_cd<nch; ch_cd++)
            {
              TwoBodyChannel& tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
              int Jcd = tbc_cd.J;
              int nkets_cd = tbc_cd.GetNumberKets();
              if (  std::abs(Jab-Jcd)>J_ph or (Jab+Jcd)<J_ph ) continue;



              for (int iket_ab=0; iket_ab<nkets_ab; iket_ab++)
              {
                Ket& ket_ab = tbc_ab.GetKet(iket_ab);
                int a = ket_ab.p;
                int b = ket_ab.q;
                double na = ket_ab.op->occ;
                double nb = ket_ab.oq->occ;
                if ( na*nb<1e-8 and (1-na)*(1-nb)<1e-8 ) continue;
                double d_ea = std::abs(2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
                double d_eb = std::abs(2*ket_ab.oq->n + ket_ab.oq->l - e_fermi[ket_ab.oq->tz2]);
                double occnat_a = ket_ab.op->occ_nat;
                double occnat_b = ket_ab.oq->occ_nat;
                bool keep_abi = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                bool keep_abj = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                bool keep_abk = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_k*(1-occnat_k) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                bool keep_abl = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_l*(1-occnat_l) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                keep_abi = keep_abi and (d_ea+d_eb+d_ei <= Z.modelspace->dE3max );
                keep_abj = keep_abj and (d_ea+d_eb+d_ej <= Z.modelspace->dE3max );
                keep_abk = keep_abk and (d_ea+d_eb+d_ek <= Z.modelspace->dE3max );
                keep_abl = keep_abl and (d_ea+d_eb+d_el <= Z.modelspace->dE3max );
                if ( not ( keep_abi or keep_abj) ) continue;
                if ( not ( keep_abk or keep_abl) ) continue;


                for (int iket_cd=0; iket_cd<nkets_cd; iket_cd++)
                {
                  Ket& ket_cd = tbc_cd.GetKet(iket_cd);
                  int c = ket_cd.p;
                  int d = ket_cd.q;
                  double nc = ket_cd.op->occ;
                  double nd = ket_cd.oq->occ;
                  double d_ec = std::abs(2*ket_cd.op->n + ket_cd.op->l - e_fermi[ket_cd.op->tz2]);
                  double d_ed = std::abs(2*ket_cd.oq->n + ket_cd.oq->l - e_fermi[ket_cd.oq->tz2]);
                  double occnat_c = ket_cd.op->occ_nat;
                  double occnat_d = ket_cd.oq->occ_nat;
//                  double jc = 0.5*ket_cd.op->j2;
//                  double jd = 0.5*ket_cd.oq->j2;
                  double occupation_factor = (1-na)*(1-nb)*nc*nd - na*nb*(1-nc)*(1-nd);
                  if (std::abs(occupation_factor)<1e-6) continue;

                  bool keep_cdi = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_i*(1-occnat_i) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                  bool keep_cdj = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_j*(1-occnat_j) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                  bool keep_cdk = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_k*(1-occnat_k) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                  bool keep_cdl = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_l*(1-occnat_l) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                  keep_cdi = keep_cdi and (d_ec+d_ed+d_ei <= Z.modelspace->dE3max );
                  keep_cdj = keep_cdj and (d_ec+d_ed+d_ej <= Z.modelspace->dE3max );
                  keep_cdk = keep_cdk and (d_ec+d_ed+d_ek <= Z.modelspace->dE3max );
                  keep_cdl = keep_cdl and (d_ec+d_ed+d_el <= Z.modelspace->dE3max );
                  if ( not ( keep_cdk or keep_cdl) ) continue;
                  if ( not ( keep_cdi or keep_cdj) ) continue;

                  double symmetry_factor = 1;  // we only sum a<=b and c<=d, so we undercount by a factor of 4, canceling the 1/4 in the formula 
                  if (a==b) symmetry_factor *= 0.5; // if a==b or c==d, then the permutation doesn't give a new state, so there's less undercounting
                  if (c==d) symmetry_factor *= 0.5;

                  double X3_il = 0;
                  double X3_jl = 0;
                  double X3_ik = 0;
                  double X3_jk = 0;
                  double Y3_il = 0;
                  double Y3_jl = 0;
                  double Y3_ik = 0;
                  double Y3_jk = 0;


                  if ( keep_abi and keep_cdl and keep_cdj and keep_abk
                     and ( (tbc_ab.parity + tbc_cd.parity + oi.l + ol.l)%2 == X.parity)
                     and  ( std::abs(2*tbc_ab.Tz+oi.tz2 - 2*tbc_cd.Tz - ol.tz2) == 2*X.rank_T ) 
                      and ( std::abs(ji2-jl2)<=2*J_ph ) and ( (ji2+jl2)>=2*J_ph )   
                      and ( std::abs(jj2-jk2)<=2*J_ph ) and ( (jj2+jk2)>=2*J_ph )   )
                  {
                    int twoJp_min = std::max( std::abs(2*Jab-ji2), std::abs(2*Jcd-jl2));
                    int twoJp_max = std::min( (2*Jab+ji2), (2*Jcd+jl2));
//                    if (i==0 and j==1 and k==0 and l==1)
//                    {
////                      std::cout << "Getting a contribution from Xabicdl = " << a << " " << b << " " << i << " " << c << " " << d << " " << l << std::endl;
//                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << i << " " << c << " " << d << " " << l << "   " << na << " " << nb << " " << nc << " " << nd << " -> " << occupation_factor << std::endl;
//                    }
                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(ji,jl,J_ph, Jcd,Jab,0.5*twoJp);
                      X3_il += (twoJp+1) * sixj_x * Z.modelspace->phase((ji2+twoJp)/2) * X3.GetME_pn( Jab, Jcd, twoJp, a,b,i,c,d,l);
                    }

                    twoJp_min = std::max( std::abs(2*Jcd-jj2), std::abs(2*Jab-jk2));
                    twoJp_max = std::min( (2*Jcd+jj2), (2*Jab+jk2));
                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(jj,jk,J_ph, Jab,Jcd,0.5*twoJp);
                      Y3_jk += (twoJp+1) * sixj_x * Z.modelspace->phase((jk2+twoJp)/2) * Y3.GetME_pn( Jcd, Jab, twoJp, c,d,j,a,b,k);
                    }

                  }


                  if ( keep_abj and keep_cdl and keep_cdi and keep_abk
                     and ( (tbc_ab.parity + tbc_cd.parity + oj.l + ol.l)%2 == X.parity)
                     and  ( std::abs(2*tbc_ab.Tz+oj.tz2 - 2*tbc_cd.Tz - ol.tz2) == 2*X.rank_T )  
                      and ( std::abs(jj2-jl2)<=2*J_ph ) and ( (jj2+jl2)>=2*J_ph )   
                      and ( std::abs(ji2-jk2)<=2*J_ph ) and ( (ji2+jk2)>=2*J_ph )   )
                  {
                    int twoJp_min = std::max( std::abs(2*Jab-jj2), std::abs(2*Jcd-jl2));
                    int twoJp_max = std::min( (2*Jab+jj2), (2*Jcd+jl2));


                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(jj,jl,J_ph, Jcd,Jab,0.5*twoJp);
                      X3_jl += (twoJp+1) * sixj_x * Z.modelspace->phase((jj2+twoJp)/2) * X3.GetME_pn( Jab, Jcd, twoJp, a,b,j,c,d,l);
                    }
                    twoJp_min = std::max( std::abs(2*Jcd-ji2), std::abs(2*Jab-jk2));
                    twoJp_max = std::min( (2*Jcd+ji2), (2*Jab+jk2));
                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(ji,jk,J_ph, Jab,Jcd,0.5*twoJp);
                      Y3_ik += (twoJp+1) * sixj_x * Z.modelspace->phase((jk2+twoJp)/2) * Y3.GetME_pn( Jcd, Jab, twoJp, c,d,i,a,b,k);
                    }
//                    if (i==0 and j==1 and k==0 and l==1)
//                    {
////                      std::cout << "Getting a contribution from Xabjcdl = " << a << " " << b << " " << j << " " << c << " " << d << " " << l << std::endl;
//                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << j << " " << c << " " << d << " " << l  << "   "  << na << " " << nb << " " << nc << " " << nd << " -> "<< occupation_factor << " X3_jl = " << X3_jl << "  Y3_ik = " << Y3_ik << std::endl;
//                    }
                  }


                  if ( keep_abi and keep_cdk and keep_cdj and keep_abl
                     and ( (tbc_ab.parity + tbc_cd.parity + oi.l + ok.l)%2 == X.parity)
                     and  ( std::abs(2*tbc_ab.Tz+oi.tz2 - 2*tbc_cd.Tz - ok.tz2) == 2*X.rank_T )  
                      and ( std::abs(ji2-jk2)<=2*J_ph ) and ( (ji2+jk2)>=2*J_ph )   
                      and ( std::abs(jj2-jl2)<=2*J_ph ) and ( (jj2+jl2)>=2*J_ph )   )
                  {
                    int twoJp_min = std::max( std::abs(2*Jab-ji2), std::abs(2*Jcd-jk2));
                    int twoJp_max = std::min( (2*Jab+ji2), (2*Jcd+jk2));
//                    if (i==0 and j==1 and k==0 and l==1)
//                    {
////                      std::cout << "Getting a contribution from Xabicdk = " << a << " " << b << " " << i << " " << c << " " << d << " " << k << std::endl;
//                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << i << " " << c << " " << d << " " << k  << "   "  << na << " " << nb << " " << nc << " " << nd << " -> "<< occupation_factor << std::endl;
//                    }
                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(ji,jk,J_ph, Jcd,Jab,0.5*twoJp);
                      X3_ik += (twoJp+1) * sixj_x * Z.modelspace->phase((ji2+twoJp)/2) * X3.GetME_pn( Jab, Jcd, twoJp, a,b,i,c,d,k);
                    }
                    twoJp_min = std::max( std::abs(2*Jcd-jj2), std::abs(2*Jab-jl2));
                    twoJp_max = std::min( (2*Jcd+jj2), (2*Jab+jl2));
                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(jj,jl,J_ph, Jab,Jcd,0.5*twoJp);
                      Y3_jl += (twoJp+1) * sixj_x * Z.modelspace->phase((jl2+twoJp)/2) * Y3.GetME_pn( Jcd, Jab, twoJp, c,d,j,a,b,l);
                    }
                  }


                  if ( keep_abj and keep_cdk and keep_cdi and keep_abl
                     and ( (tbc_ab.parity + tbc_cd.parity + oj.l + ok.l)%2 == X.parity)
                     and  ( std::abs(2*tbc_ab.Tz+oj.tz2 - 2*tbc_cd.Tz - ok.tz2) == 2*X.rank_T )  
                      and ( std::abs(jj2-jk2)<=2*J_ph ) and ( (jj2+jk2)>=2*J_ph )   
                      and ( std::abs(ji2-jl2)<=2*J_ph ) and ( (ji2+jl2)>=2*J_ph )   )
                  {
                    int twoJp_min = std::max( std::abs(2*Jab-jj2), std::abs(2*Jcd-jk2));
                    int twoJp_max = std::min( (2*Jab+jj2), (2*Jcd+jk2));
//                    if (i==0 and j==1 and k==0 and l==1)
//                    {
////                      std::cout << "Getting a contribution from Xabjcdk = " << a << " " << b << " " << j << " " << c << " " << d << " " << k << std::endl;
//                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << j << " " << c << " " << d << " " << k  << "   "  << na << " " << nb << " " << nc << " " << nd << " -> "<< occupation_factor << std::endl;
//                    }
                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(jj,jk,J_ph, Jcd,Jab,0.5*twoJp);
                      X3_jk += (twoJp+1) * sixj_x * Z.modelspace->phase((jj2+twoJp)/2) * X3.GetME_pn( Jab, Jcd, twoJp, a,b,j,c,d,k);
                    }
                    twoJp_min = std::max( std::abs(2*Jcd-ji2), std::abs(2*Jab-jl2));
                    twoJp_max = std::min( (2*Jcd+ji2), (2*Jab+jl2));
                    for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                    {
                      double sixj_x = Z.modelspace->GetSixJ(ji,jl,J_ph, Jab,Jcd,0.5*twoJp);
                      Y3_il += (twoJp+1) * sixj_x * Z.modelspace->phase((jl2+twoJp)/2) * Y3.GetME_pn( Jcd, Jab, twoJp, c,d,i,a,b,l);
                    }
                  }



                  zbar_ilkj += occupation_factor * symmetry_factor * X3_il * Y3_jk;
                  zbar_jlki += occupation_factor * symmetry_factor * X3_jl * Y3_ik;
                  zbar_iklj += occupation_factor * symmetry_factor * X3_ik * Y3_jl;
                  zbar_jkli += occupation_factor * symmetry_factor * X3_jk * Y3_il;
//                  zbar_ilkj += occupation_factor * symmetry_factor * ( X3_il * Y3_jk + phase_ij*phase_kl*X3_jk * Y3_il);
//                  zbar_jlki += occupation_factor * symmetry_factor * ( X3_jl * Y3_ik + phase_ij*phase_kl*X3_ik * Y3_jl);



                }// for iket_cd
            }// for iket_ab
              }// for ch_cd
          }// for ch_ab


        double sixj1 = Z.modelspace->GetSixJ(ji, jj, J,  jk,jl, J_ph);
        double sixj1_flip = Z.modelspace->GetSixJ(jj, ji, J,  jk,jl, J_ph);
//        int phase_ij = Z.modelspace->phase( (ji2+jj2)/2-J);
//        int phase_kl = Z.modelspace->phase( (jk2+jl2)/2-J);
        zijkl += (2*J_ph+1) * (  sixj1 * (zbar_ilkj + phase_ij*phase_kl * zbar_jkli) - sixj1_flip*( phase_ij*zbar_jlki + phase_kl*zbar_iklj) );
//        zijkl += (2*J_ph+1) * (  sixj1 * zbar_ilkj   - sixj1_flip* phase_ij*zbar_jlki );
//        if ( i==0 and j==1 and k==0 and l==1)
//        {
//           std::cout << " Jph = " << J_ph << " zbar_ilkj " << zbar_ilkj << "   zbar_jkli " << zbar_jkli << "   zbar_jlki " << zbar_jlki << "   zbar_iklj " << zbar_iklj << "  ->  " << zijkl << std::endl;
//        }


        }// for J_ph
        // make it a normalized TBME
        zijkl /=  sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
        // the AddToTBME routine automatically takes care of the hermitian conjugate as well
        Z2.AddToTBME(ch,ch,ibra,iket,zijkl);
      }// for iket
    }// for ibra
  }// for ch
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}




/*
//void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z )
void comm332_pphhss_debug( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();
  
  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (int ch=0; ch<nch; ch++)
  {
    TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      int ji2 = bra.op->j2;
      int jj2 = bra.oq->j2;
      double ji = 0.5*ji2;
      double jj = 0.5*jj2;
      double d_ei = std::abs(2*bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
      double d_ej = std::abs(2*bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
      double occnat_i = bra.op->occ_nat;
      double occnat_j = bra.oq->occ_nat;
      for (int iket=ibra; iket<nkets; iket++) 
      {
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        int jk2 = ket.op->j2;
        int jl2 = ket.oq->j2;
        double jk = 0.5*jk2;
        double jl = 0.5*jl2;
        double d_ek = std::abs(2*ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
        double d_el = std::abs(2*ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
        double occnat_k = ket.op->occ_nat;
        double occnat_l = ket.oq->occ_nat;

        double zijkl = 0;
        // Now the loops on the right hand side
        for (int ch_ab=0; ch_ab<nch; ch_ab++)
        {
          TwoBodyChannel& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
          int Jab = tbc_ab.J;
          int nkets = tbc_ab.GetNumberKets();
          for (int iket_ab=0; iket_ab<nkets; iket_ab++)
          {
            Ket& ket_ab = tbc_ab.GetKet(iket_ab);
            int a = ket_ab.p;
            int b = ket_ab.q;
            double na = ket_ab.op->occ;
            double nb = ket_ab.oq->occ;
            double d_ea = std::abs(2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
            double d_eb = std::abs(2*ket_ab.oq->n + ket_ab.oq->l - e_fermi[ket_ab.oq->tz2]);
            double occnat_a = ket_ab.op->occ_nat;
            double occnat_b = ket_ab.oq->occ_nat;
            bool keep_abi = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) >= Z.modelspace->GetOccNat3Cut() ) ;
            bool keep_abj = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) >= Z.modelspace->GetOccNat3Cut() ) ;
            bool keep_abk = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_k*(1-occnat_k) ) >= Z.modelspace->GetOccNat3Cut() ) ;
            bool keep_abl = ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_l*(1-occnat_l) ) >= Z.modelspace->GetOccNat3Cut() ) ;
            keep_abi = keep_abi and (d_ea+d_eb+d_ei <= Z.modelspace->dE3max );
            keep_abj = keep_abj and (d_ea+d_eb+d_ej <= Z.modelspace->dE3max );
            keep_abk = keep_abk and (d_ea+d_eb+d_ek <= Z.modelspace->dE3max );
            keep_abl = keep_abl and (d_ea+d_eb+d_el <= Z.modelspace->dE3max );
            if ( not ( keep_abi or keep_abj) ) continue;
            if ( not ( keep_abk or keep_abl) ) continue;

            for (int ch_cd=0; ch_cd<nch; ch_cd++)
            {
              TwoBodyChannel& tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
              int Jcd = tbc_cd.J;
              int nkets = tbc_cd.GetNumberKets();
              for (int iket_cd=0; iket_cd<nkets; iket_cd++)
              {
                Ket& ket_cd = tbc_cd.GetKet(iket_cd);
                int c = ket_cd.p;
                int d = ket_cd.q;
                double nc = ket_cd.op->occ;
                double nd = ket_cd.oq->occ;
                double d_ec = std::abs(2*ket_cd.op->n + ket_cd.op->l - e_fermi[ket_cd.op->tz2]);
                double d_ed = std::abs(2*ket_cd.oq->n + ket_cd.oq->l - e_fermi[ket_cd.oq->tz2]);
                double occnat_c = ket_cd.op->occ_nat;
                double occnat_d = ket_cd.oq->occ_nat;
//                double jc = 0.5*ket_cd.op->j2;
//                double jd = 0.5*ket_cd.oq->j2;
                double occupation_factor = (1-na)*(1-nb)*nc*nd - na*nb*(1-nc)*(1-nd);
                if (std::abs(occupation_factor)<1e-6) continue;


                bool keep_cdi = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_i*(1-occnat_i) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                bool keep_cdj = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_j*(1-occnat_j) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                bool keep_cdk = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_k*(1-occnat_k) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                bool keep_cdl = ( (occnat_c*(1-occnat_c) * occnat_d*(1-occnat_d) * occnat_l*(1-occnat_l) ) >= Z.modelspace->GetOccNat3Cut() ) ;
                keep_cdi = keep_cdi and (d_ec+d_ed+d_ei <= Z.modelspace->dE3max );
                keep_cdj = keep_cdj and (d_ec+d_ed+d_ej <= Z.modelspace->dE3max );
                keep_cdk = keep_cdk and (d_ec+d_ed+d_ek <= Z.modelspace->dE3max );
                keep_cdl = keep_cdl and (d_ec+d_ed+d_el <= Z.modelspace->dE3max );
                if ( not ( keep_cdk or keep_cdl) ) continue;
                if ( not ( keep_cdi or keep_cdj) ) continue;

                double symmetry_factor = 1;  // we only sum a<=b and c<=d, so we undercount by a factor of 4, canceling the 1/4 in the formula 
                if (a==b) symmetry_factor *= 0.5; // if a==b or c==d, then the permutation doesn't give a new state, so there's less undercounting
                if (c==d) symmetry_factor *= 0.5;

                // Figure out which range of twoJp and twoJpp we will need
                int twoJp_min  = std::max(  std::min(std::abs(2*Jab - ji2),std::abs(2*Jab-jj2)), std::min(std::abs(2*Jcd - jl2),std::abs(2*Jcd-jk2)) );
                int twoJpp_min = std::max(  std::min(std::abs(2*Jcd - jj2),std::abs(2*Jcd-ji2)), std::min(std::abs(2*Jab - jk2),std::abs(2*Jab-jl2)) );
                int twoJp_max  = std::min(  2*Jab + std::max(ji2,jj2),  2*Jcd + std::max(jk2,jl2) );
                int twoJpp_max = std::min(  2*Jcd + std::max(ji2,jj2),  2*Jab + std::max(jk2,jl2) );

                if (twoJpp_max<twoJpp_min) continue;
                for (int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
                {

                    double xabicdl = (keep_abi and keep_cdl) ? X3.GetME_pn(Jab,Jcd,twoJp, a,b,i,c,d,l) : 0;
                    double xabicdk = (keep_abi and keep_cdk) ? X3.GetME_pn(Jab,Jcd,twoJp, a,b,i,c,d,k) : 0;
                    double xabjcdl = (keep_abj and keep_cdl) ? X3.GetME_pn(Jab,Jcd,twoJp, a,b,j,c,d,l) : 0;
                    double xabjcdk = (keep_abj and keep_cdk) ? X3.GetME_pn(Jab,Jcd,twoJp, a,b,j,c,d,k) : 0;

                  for (int twoJpp=twoJpp_min; twoJpp<=twoJpp_max; twoJpp+=2)
                  {
                    double Jp = 0.5 * twoJp;
                    double Jpp = 0.5 * twoJpp;
                    double hatfactor = (twoJp+1)*(twoJpp+1);

                    // I think having these in the inner loop may be disastrous for performance
                    double ninej1 = (keep_abi and keep_cdl and keep_cdj and keep_abk) ? Z.modelspace->GetNineJ(Jab,jk,Jpp, ji,J,jj, Jp,jl,Jcd) : 0;
                    double ninej2 = (keep_abj and keep_cdl and keep_cdi and keep_abk) ? Z.modelspace->GetNineJ(Jab,jk,Jpp, jj,J,ji, Jp,jl,Jcd) : 0; // permute i<->j
                    double ninej3 = (keep_abi and keep_cdk and keep_cdj and keep_abl) ? Z.modelspace->GetNineJ(Jab,jl,Jpp, ji,J,jj, Jp,jk,Jcd) : 0; // permute k<->l
                    double ninej4 = (keep_abj and keep_cdk and keep_cdi and keep_abl) ? Z.modelspace->GetNineJ(Jab,jl,Jpp, jj,J,ji, Jp,jk,Jcd) : 0; // permute i<->j and k<->l

                    // These phase factors account for the minus signs associated with the permutations
                    // so that all the permuted terms should just be added with their phase (no extra minus sign).
                    int phase1 = Z.modelspace->phase( (ji2+jk2+twoJp+twoJpp)/2 );
                    int phase2 = Z.modelspace->phase( (ji2+jk2+twoJp+twoJpp)/2 -J ); // from permuting i<->j
                    int phase3 = Z.modelspace->phase( (ji2+jk2+twoJp+twoJpp)/2 -J ); // from permuting k<->l
                    int phase4 = Z.modelspace->phase( (ji2+jk2+twoJp+twoJpp)/2 ); // from permuting i<->j and k<->l

                    double ycdjabk = (keep_cdj and keep_abk) ? Y3.GetME_pn(Jcd,Jab,twoJpp,c,d,j,a,b,k) : 0;
                    double ycdjabl = (keep_cdj and keep_abl) ? Y3.GetME_pn(Jcd,Jab,twoJpp,c,d,j,a,b,l) : 0;
                    double ycdiabk = (keep_cdi and keep_abk) ? Y3.GetME_pn(Jcd,Jab,twoJpp,c,d,i,a,b,k) : 0;
                    double ycdiabl = (keep_cdi and keep_abl) ? Y3.GetME_pn(Jcd,Jab,twoJpp,c,d,i,a,b,l) : 0;
                    zijkl += symmetry_factor * hatfactor * occupation_factor * (
                                      1*phase1 * ninej1 * xabicdl * ycdjabk
                                    + 1*phase2 * ninej2 * xabjcdl * ycdiabk // i<->j
                                    + 1*phase3 * ninej3 * xabicdk * ycdjabl // k<->l
                                    + 1*phase4 * ninej4 * xabjcdk * ycdiabl // i<->j and k<->l
                                   );

                  }// for twoJpp
                }// for twoJp
                
              }// for iket_cd
            }// for ch_cd
          }// for iket_ab
        }// for ch_ab
        // make it a normalized TBME
        zijkl /=  sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
        // the AddToTBME routine automatically takes care of the hermitian conjugate as well
        Z2.AddToTBME(ch,ch,ibra,iket,zijkl);
      }// for iket
    }// for ibra
  }// for ch
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}

*/



//*****************************************************************************************
//
// i|  j|  k|       Uncoupled expression:
//  |   |   *--[X]     Z_ijklmn = sum_a P(ij/k) (X_ka * Y_ijalmn) - P(lm/n) (Y_ijklma * X_an)
//  |   |   |a
//  *~~[Y]~~*       Coupled expression:
//  |   |   |          Z_{ijklmn}^{J1,J2,J} = sum_a P(ij/k)^{J1,J} (X_ka * Y_{ijalmn}^{J1,J2,J})
// l|  m|  n|                                      -P(lm/n)^{J2,J} (Y_{ijklma}^{J1,J2,J} * X_an)
//
//  Checked with UnitTest and passed
//
void comm133ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  auto& X1 = X.OneBody;
  auto& Y1 = Y.OneBody;

  int hermX = X.IsHermitian() ? 1 : -1;
  int hermY = Y.IsHermitian() ? 1 : -1;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  int norbs = Z.modelspace->GetNumberOrbits();
  double X3NORM = X3.Norm();
  double Y3NORM = Y3.Norm();
  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();
//  bool x3_allocated = X3.is_allocated;
//  bool y3_allocated = Y3.is_allocated;
//  if (X3NORM<1e-6 and Y3NORM<1e-6 ) return;
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    int twoJ = Tbc.twoJ;
    size_t nkets = Tbc.GetNumberKets();
    std::vector<size_t> kets_kept;
    std::map<size_t,size_t> kept_lookup;
    for (size_t iket=0; iket<nkets; iket++)
    {
      Ket3& ket = Tbc.GetKet(iket);
      double d_ei = std::abs(2*ket.op->n + ket.op->l - e_fermi.at(ket.op->tz2));
      double d_ej = std::abs(2*ket.oq->n + ket.oq->l - e_fermi.at(ket.oq->tz2));
      double d_ek = std::abs(2*ket.oR->n + ket.oR->l - e_fermi.at(ket.oR->tz2));
      double occnat_i = ket.op->occ_nat;
      double occnat_j = ket.oq->occ_nat;
      double occnat_k = ket.oR->occ_nat;
      if (  (d_ei+d_ej+d_ek ) > Z.modelspace->GetdE3max() ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      kets_kept.push_back( iket );
      kept_lookup[iket] = kets_kept.size()-1;
    }
    size_t nkets_kept = kets_kept.size();


    arma::mat X1MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat Y1MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat X3MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat Y3MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat Z3MAT( nkets_kept, nkets_kept, arma::fill::zeros);

    for (size_t index_bra=0; index_bra<nkets_kept; index_bra++)
    {
      size_t ibra = kets_kept[index_bra];
      Ket3& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      int Jij = bra.Jpq;
      double d_ei = std::abs(2*oi.n + oi.l - e_fermi.at(oi.tz2));
      double d_ej = std::abs(2*oj.n + oj.l - e_fermi.at(oj.tz2));
      double d_ek = std::abs(2*ok.n + ok.l - e_fermi.at(ok.tz2));
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;


      for (auto a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
        Orbit& oa = X.modelspace->GetOrbit(a);
        double d_ea = std::abs(2*oa.n+oa.l - e_fermi.at(oa.tz2));
        double occnat_a = oa.occ_nat;
        if (  (d_ea+d_ej+d_ek ) > Z.modelspace->GetdE3max() ) continue;
        if ( (occnat_a*(1-occnat_a) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        std::vector<size_t> ket_list;
        std::vector<double> recouple_list;
        
        size_t ch_check = Z3.GetKetIndex_withRecoupling( Jij, twoJ, a, j, k,  ket_list,  recouple_list );
        for (size_t ilist=0; ilist<ket_list.size(); ilist++)
        {
          auto iter_find = kept_lookup.find( ket_list[ilist] );
          if (iter_find == kept_lookup.end() ) continue;
          size_t index_ket = iter_find->second;
          double recouple = recouple_list[ilist];
          X1MAT(index_bra,index_ket) += X1(i,a) * recouple;
          Y1MAT(index_bra,index_ket) += Y1(i,a) * recouple;
        }
      }
      for (auto a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
      {

        Orbit& oa = Z.modelspace->GetOrbit(a);
        double d_ea = std::abs(2*oa.n+oa.l - e_fermi.at(oa.tz2));
        double occnat_a = oa.occ_nat;
        if (  (d_ei+d_ea+d_ek ) > Z.modelspace->GetdE3max() ) continue;
        if ( (occnat_i*(1-occnat_i) * occnat_a*(1-occnat_a) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        std::vector<size_t> ket_list;
        std::vector<double> recouple_list;
        size_t ch_check = Z3.GetKetIndex_withRecoupling( Jij, twoJ, i, a, k,  ket_list,  recouple_list );
        for (size_t ilist=0; ilist<ket_list.size(); ilist++)
        {
          auto iter_find = kept_lookup.find( ket_list[ilist] );
          if (iter_find == kept_lookup.end() ) continue;
          size_t index_ket = iter_find->second;
          double recouple = recouple_list[ilist];
          X1MAT(index_bra,index_ket) += X1(j,a) * recouple;
          Y1MAT(index_bra,index_ket) += Y1(j,a) * recouple;
        }
      }
      for (auto a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
      {
        Orbit& oa = Z.modelspace->GetOrbit(a);
        double d_ea = std::abs(2*oa.n+oa.l - e_fermi.at(oa.tz2));
        double occnat_a = oa.occ_nat;
        if (  (d_ei+d_ej+d_ea ) > Z.modelspace->GetdE3max() ) continue;
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_a*(1-occnat_a) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        std::vector<size_t> ket_list;
        std::vector<double> recouple_list;
        size_t ch_check = Z3.GetKetIndex_withRecoupling( Jij, twoJ, i, j, a,  ket_list,  recouple_list );
        for (size_t ilist=0; ilist<ket_list.size(); ilist++)
        {
          auto iter_find = kept_lookup.find( ket_list[ilist] );
          if (iter_find == kept_lookup.end() ) continue;
          size_t index_ket = iter_find->second;
          double recouple = recouple_list[ilist];
          X1MAT(index_bra,index_ket) += X1(k,a) * recouple;
          Y1MAT(index_bra,index_ket) += Y1(k,a) * recouple;
        }
      }
    }// for ibra

    // kept_lookup is a map   Full index => Kept index, so iter_bra.first gives the full index, and iter_bra.second is the
    // index for the 3-body state we keep in this commutator
    for ( auto& iter_bra : kept_lookup )
    {
      for ( auto& iter_ket : kept_lookup )
      {
//        if ( X3NORM > 1e-6)
        if ( x3_allocated )
           X3MAT( iter_bra.second, iter_ket.second) = X3.GetME_pn_ch(ch3,ch3, iter_bra.first, iter_ket.first );
//        if ( Y3NORM > 1e-6)
        if ( y3_allocated )
           Y3MAT( iter_bra.second, iter_ket.second) = Y3.GetME_pn_ch(ch3,ch3, iter_bra.first, iter_ket.first );
      }
    }


    // Do the matrix multiplication
//    Z3MAT = X1MAT*Y3MAT - Y1MAT*X3MAT +  hermX*hermY * ( X3MAT.t()*Y1MAT.t() - Y3MAT.t()*X1MAT.t() );
    Z3MAT = X1MAT*Y3MAT - Y1MAT*X3MAT;
    Z3MAT -=  hermX*hermY * Z3MAT.t();


    // unpack the result
    for ( auto& iter_bra : kept_lookup )
    {
      for ( auto& iter_ket : kept_lookup )
      {
        if ( iter_ket.first < iter_bra.first ) continue;
        Z3.AddToME_pn_ch(ch3,ch3, iter_bra.first,iter_ket.first,  Z3MAT(iter_bra.second,iter_ket.second) );
      }
    }
 
  }// for ch3
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}


// the old slow way
/*
void comm133ss( const Operator& X, const Operator& Y, Operator& Z )
{
//  double normY = Y.ThreeBodyNorm();
//  int e3maxcut = 96;
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  auto& X1 = X.OneBody;
  auto& Y1 = Y.OneBody;
//  std::cout << "Y1 =" << std::endl << Y1 << std::endl;
  int norbs = Z.modelspace->GetNumberOrbits();
  if (X3.Norm()<1e-6 and Y3.Norm()<1e-6 ) return;
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
//  #pragma omp parallel for schedule(dynamic,1)
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<=nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    int twoJ = Tbc.twoJ;
    size_t nkets = Tbc.GetNumberKets();
    for (size_t ibra=0; ibra<nkets; ibra++)
    {
      Ket3& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      int Jij = bra.Jpq;
      for (size_t iket=ibra; iket<nkets; iket++)
      {
        Ket3& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        int Jlm = ket.Jpq;
        
        double zsum =0;
        // First, connect on the bra side
        for (auto a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
        {
          zsum += X1(i,a) * Y3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
          zsum -= Y1(i,a) * X3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
        }
        for (auto a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
        {
          zsum += X1(j,a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
          zsum -= Y1(j,a) * X3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
        }
        for (auto a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
        {
          zsum += X1(k,a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
          zsum -= Y1(k,a) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
        }
        // Now connect on the ket side
        for (auto a : X.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
        {
          zsum -= X1(a,l) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
          zsum += Y1(a,l) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
        }
        for (auto a : X.OneBodyChannels.at({om.l,om.j2,om.tz2}) )
        {
          zsum -= X1(a,m) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
          zsum += Y1(a,m) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
        }
        for (auto a : X.OneBodyChannels.at({on.l,on.j2,on.tz2}) )
        {
          zsum -= X1(a,n) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
          zsum += Y1(a,n) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
        }
  
//        Z3.AddToME_pn(Jij, Jlm, twoJ, i,j,k,l,m,n, zsum );
        Z3.AddToME_pn_ch(ch3,ch3, ibra, iket, zsum );

      }// for iket
    }// for ibra
  }// for ch3
            
}
*/




//*****************************************************************************************
//
//  i|  j|  k|   Uncoupled expression:  
//   |~X~|   |      Z_ijklmn  =  P(ij/k)P(lm/n) sum_a  (X_ijla * Y_akmn - Y_ijla * X_akmn)
//   |   |a  |
//   |   |~Y~|   Coupled expression:
//  l|  m|  n|     Z_{ijklmn}^{J1,J2,J}  = Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9
//                 where each of those terms corresponds to a permutation from the uncoupled expression
//
// We cast this as mat-mult by flipping things around and recoupling:
// i|  j|  k|           l\ |i j|
//  |~X~|   |             \|~X~|
//  |   |a  |     -->          |a
//  |   |~Y~|                  |~Y~|\
// l|  m|  n|                 m|  n| \k
//  checked with UnitTest, and it looks good. 
//
void comm223ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  double t_internal = omp_get_wtime();
  int norbs = Z.modelspace->GetNumberOrbits();
  auto& Z3 = Z.ThreeBody;
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  int hX = X.IsHermitian() ? 1 : -1;
  int hY = Y.IsHermitian() ? 1 : -1;
  int hZ = Z.IsHermitian() ? 1 : -1;
  if ( (std::abs( X2.Norm() * Y2.Norm() ) < 1e-6 ) and not Z.modelspace->scalar3b_transform_first_pass) return;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();

  // create a function object for hashing 4 size_t (unsigned long int) to a single size_t for lookup in an unordered_map
  int max_twoJph = 6*Z.modelspace->GetEmax()+3;
  auto ch_pph_hash = [](int twoJ_ph, int parity, int twoTz_ph){ return 8*(twoJ_ph/2) + parity + (twoTz_ph+3);};
  std::vector<int> channels_pph(ch_pph_hash(max_twoJph,1,3)+1, -1); // map twoJ_ph, parity_ph, twoTz_ph  -> channel index
  std::vector< std::array<int,3>> channel_list_pph; // vector  channel_index -> twoJ-ph, parity_ph, twoTz_ph.

  auto hash_key_ijnJ = [](size_t i,size_t j, size_t n, size_t Jij){return ( i + (j<<12) + (n<<24) + (Jij<<36) );};

   //  0xF = 15 = 1111 (4bits), so 0xFFF is 12 bits of 1's. 0xFFFL makes it a long int
  auto unhash_key_ijnJ = [](size_t& i,size_t& j, size_t& n, size_t& Jij, size_t key){ i=( key     &0xFFFL);
                                                                                      j=((key>>12)&0xFFFL);
                                                                                      n=((key>>24)&0xFFFL);
                                                                                    Jij=((key>>36)&0xFFFL); }; 
  std::vector< std::unordered_map<size_t, size_t>> ket_lookup_pph; // in a given channel, map  i,j,n,Jij -> matrix index

  std::vector<float> Zbar;  // store the transformed Z in a 1D vector.
  std::vector<size_t> Zbar_start_pointers;  // ch -> index in Zbar where that channel starts

  // Different approach: store by ijklmn index, then for each ijklmn combo, loop over Jij,Jlm,twoJ
  // In this case, we need a map {i,j,k,l,m,n} => size_t, where the size_t points to the start locating in Zbar


//  std::cout << "Made it as far as line " << __LINE__ << std::endl;

  // we're looking at pph type states. We don't make a cut on the occupations in this channel
  // but if the first two orbits ij can't possibly make it past the occnat cut later on, we don't bother including them
  // To figure out if they'll make it, we want to know the biggest contribution that can possibly be made by a third orbit.
  double occnat_factor_max = 0;
  for ( auto i : Z.modelspace->all_orbits)
  {
     double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
     occnat_factor_max = std::max( occnat_factor_max, occnat_i*(1-occnat_i) );
  }

//  std::cout << "BeginLoop over ph channels" << std::endl;
  size_t nch_pph = 0;
  size_t Zbar_n_elements = 0; // total number of elements stored in Zbar
  // Generate all the pph type channels J,parity,Tz. Note that these will be
  // not the same as the standard 3-body channels since we aren't enforcing antisymmetry with the 3rd particle
  for ( auto& iter_obc : Z.modelspace->OneBodyChannels )
  {

    int parity_ph = iter_obc.first[0]%2;  // parity = l%2
    int twoJph = iter_obc.first[1];
    int twoTz_ph = iter_obc.first[2];
//     std::cout << "Made it as far as line " << __LINE__ << "   " << parity_ph << " " << twoJph << " " << twoTz_ph << std::endl;

    size_t ngood_ijn=0;
    std::unordered_map<size_t, size_t> good_kets_ijn;
    for ( size_t chij=0; chij<nch2; chij++)
    {
      TwoBodyChannel& tbc_ij = Z.modelspace->GetTwoBodyChannel(chij);
      int Jij = tbc_ij.J;
      int parity_ij = tbc_ij.parity;
      int Tz_ij = tbc_ij.Tz;
//     std::cout << "Made it as far as line " << __LINE__ << "   " << Jij << " " << parity_ij << " " << Tz_ij << std::endl;
      size_t nkets_ij = tbc_ij.GetNumberKets();
      int tz2n = 2*Tz_ij - twoTz_ph;
      if ( std::abs(tz2n)>1 ) continue;
      int j2n_min = std::abs(2*Jij-twoJph);
      int j2n_max = (2*Jij+twoJph);
      int parity_n = (parity_ij + parity_ph)%2;
      
      std::vector<size_t> good_n;
      for ( int j2n=j2n_min; j2n<=j2n_max; j2n+=2)
      {
        int l_n = ((j2n+1)/2)%2 == parity_n ? (j2n+1)/2  :  (j2n-1)/2;
        if (l_n > Z.modelspace->GetEmax() or l_n > Z.modelspace->GetLmax() ) continue;
        if ( Z.OneBodyChannels.find({l_n,j2n,tz2n}) == Z.OneBodyChannels.end()) continue;
        for ( size_t n : Z.OneBodyChannels.at({l_n,j2n,tz2n}) )
        {
          Orbit& on = Z.modelspace->GetOrbit(n);
          double occnat_n = on.occ_nat;
          double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
          if ( (occnat_n*(1-occnat_n) * occnat_factor_max * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
          if ( (d_en) > Z.modelspace->dE3max ) continue;
          good_n.push_back(n);
        }// for n
      }// for j2n
//      std::cout << " done looping j2n" << std::endl;

      std::vector<std::array<size_t,2>> good_ij;
      for ( size_t iket_ij=0; iket_ij<nkets_ij; iket_ij++)
      {
        // check that we pass our cuts on E3max, occupations, etc.
        Ket& ket_ij = tbc_ij.GetKet(iket_ij);
        size_t i = ket_ij.p;
        size_t j = ket_ij.q;
        double occnat_i = ket_ij.op->occ_nat;
        double occnat_j = ket_ij.oq->occ_nat;
        int ei = 2*ket_ij.op->n + ket_ij.op->l;
        int ej = 2*ket_ij.oq->n + ket_ij.oq->l;
        double d_ei = std::abs( ei - e_fermi[ket_ij.op->tz2]);
        double d_ej = std::abs( ej - e_fermi[ket_ij.oq->tz2]);
        // if i and j cant make it past the OccNat and dE3max cuts, don't bother including it
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( (d_ei+d_ej) > Z.modelspace->GetdE3max() ) continue;
        if ( (ei + ej) > Z.modelspace->GetE3max() ) continue;
        if ( perturbative_triples and ( (ket_ij.op->cvq==0 and ket_ij.oq->cvq!=0) or (ket_ij.op->cvq!=0 and ket_ij.oq->cvq==0)) ) continue;
        good_ij.push_back({i,j}); 
      }//for iket_ij

      for ( auto& ij : good_ij )
      {
        for ( size_t n : good_n )
        {
          if ( perturbative_triples ) // we just want <ppp|hhh> terms in the end, so for ijn we want pph or hhp.
          {
            Orbit& oi = Z.modelspace->GetOrbit(ij[0]);
            Orbit& oj = Z.modelspace->GetOrbit(ij[1]);
            Orbit& on = Z.modelspace->GetOrbit(n);
            if ( (oi.cvq==0 and oj.cvq==0 and on.cvq==0) or ( oi.cvq!=0 and oj.cvq!=0 and on.cvq!=0) ) continue;
          }
          size_t key = hash_key_ijnJ( ij[0],ij[1],n,size_t(Jij));
          good_kets_ijn[ key ] = ngood_ijn;
          ngood_ijn++;
        }
      }// for iket_ij
    }// for chij
   
//    std::cout << "Made it as far as line " << __LINE__ << std::endl;
    if ( (ngood_ijn < 1) ) continue;

    channels_pph[ ch_pph_hash(twoJph,parity_ph,twoTz_ph) ] = nch_pph;

    channel_list_pph.push_back({twoJph, parity_ph,twoTz_ph});
    ket_lookup_pph.push_back( good_kets_ijn );
    Zbar_start_pointers.push_back( Zbar_n_elements ) ;     
    Zbar_n_elements += ngood_ijn * (ngood_ijn+1) / 2 ;
    nch_pph++;
  }// for iter_obc

  std::cout << __func__ << "  allocating Zbar   " << Zbar_n_elements << " elements  ~ " << Zbar_n_elements * sizeof(float) /(1024.*1024*1024) << " GB" << std::endl;
  Zbar.resize( Zbar_n_elements, 0. );
  
  Z.profiler.timer["comm223_setup_loop"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();
//  std::cout << "done with setup" << std::endl;


 // Now we should fill those Zbar matrices.
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
 for (size_t ch_pph=0; ch_pph<nch_pph; ch_pph++)
 {

   auto& ch_info = channel_list_pph[ch_pph];
   int twoJph    = ch_info[0];
   int parity_ph = ch_info[1];
   int twoTz_ph  = ch_info[2];
   int l_a = ((twoJph+1)/2)%2 == parity_ph ? (twoJph+1)/2  :  (twoJph-1)/2;


//   auto& a_set = Z.modelspace->OneBodyChannels.at({l_a,twoJph,twoTz_ph}); // list of one-body orbits in this channel
   std::vector<size_t> a_list;
   for (auto a : Z.modelspace->OneBodyChannels.at({l_a,twoJph,twoTz_ph}) ) a_list.push_back(a);
   size_t number_a = a_list.size();

   auto& good_kets_ijn = ket_lookup_pph[ch_pph];
   size_t ngood_ijn = good_kets_ijn.size();
   arma::mat Xph( ngood_ijn, number_a );
   arma::mat Yph( number_a, ngood_ijn );
   for (auto& iter_ijn : good_kets_ijn )
   {
      size_t Iijn = iter_ijn.second;

      size_t key = iter_ijn.first;
      size_t i,j,n,Jij_tmp;
      unhash_key_ijnJ(i,j,n,Jij_tmp, key);
      int Jij = (int)Jij_tmp;

      for (size_t index_a=0; index_a<number_a; index_a++ )
      {
        size_t a = a_list[index_a];
        Orbit& oa = Z.modelspace->GetOrbit(a);
        Xph(Iijn, index_a) = X.TwoBody.GetTBME_J(Jij,i,j,n,a);
        Yph(index_a, Iijn) = Y.TwoBody.GetTBME_J(Jij,n,a,i,j);

      }// for index_a

   }//for iter_ijn
   arma::mat z_tmp = Xph*Yph;
   z_tmp -= hX*hY * z_tmp.t();
   // Fold the full (anti-)symmetric matrix into a more compact storage
   // we  use row-major ordering so that the inner loop runs over adjacent elements
   size_t start_ptr = Zbar_start_pointers[ch_pph];
   for (size_t II=0; II<ngood_ijn; II++)
   {
     size_t row_index = start_ptr + ( 2*ngood_ijn - II - 1) * II / 2 ;
     for (size_t JJ=II; JJ<ngood_ijn; JJ++)
     {
       Zbar[ row_index+JJ ] = z_tmp(II,JJ);
     }
   }
 }// for ch_pph

// std::cout << "Done filling matrices " << std::endl;


  Z.profiler.timer["comm223_fill_loop"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      if ( perturbative_triples and  not ( (oi.cvq + oj.cvq + ok.cvq)==0 or (oi.cvq>0 and oj.cvq>0 and ok.cvq>0)) ) continue;
      int J1 = bra.Jpq;

      // Set up the permutation stuff for ijk
      std::vector<std::array<size_t,3>> ijk = { {i,j,k}, {k,j,i}, {i,k,j} };
      std::vector<int> J1p_min = {J1,  std::max(std::abs(ok.j2-oj.j2),std::abs(twoJ-oi.j2) )/2,   std::max(std::abs(oi.j2-ok.j2), std::abs(twoJ-oj.j2) )/2 };
      std::vector<int> J1p_max = {J1,  std::min(ok.j2+oj.j2, twoJ+oi.j2)/2 , std::min(oi.j2+ok.j2, twoJ+oj.j2)/2 };
      std::vector<std::vector<double>> recouple_ijk = {{1},{},{} };

      for (int J1p=J1p_min[1]; J1p<=J1p_max[1]; J1p++)
           recouple_ijk[1].push_back( sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p) );

      for (int J1p=J1p_min[2]; J1p<=J1p_max[2]; J1p++)
           recouple_ijk[2].push_back( -Z.modelspace->phase((oj.j2+ok.j2)/2+J1+J1p)*sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p) );

      for (size_t iket=ibra; iket<nkets3; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( perturbative_triples and  not ( (ol.cvq + om.cvq + on.cvq)==0 or (ol.cvq>0 and om.cvq>0 and on.cvq>0)) ) continue;
        if ( perturbative_triples and (  (oi.cvq==0 and ol.cvq==0) or (oi.cvq!=0 and ol.cvq!=0) ) ) continue;
        int J2 = ket.Jpq;

        // Set up the permutation stuff for lmn
        std::vector<std::array<size_t,3>> lmn = { {l,m,n}, {n,m,l}, {l,n,m} };
        std::vector<int> J2p_min = {J2,  std::max(std::abs(on.j2-om.j2), std::abs(twoJ-ol.j2) )/2,   std::max( std::abs(ol.j2-on.j2), std::abs(twoJ-om.j2) )/2 };
        std::vector<int> J2p_max = {J2,  std::min(on.j2+om.j2, twoJ+ol.j2)/2 , std::min(ol.j2+on.j2,twoJ+om.j2)/2 };
        std::vector<std::vector<double>> recouple_lmn = {{1},{},{} };
            
        for (int J2p=J2p_min[1]; J2p<=J2p_max[1]; J2p++)
           recouple_lmn[1].push_back( sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p) );

        for (int J2p=J2p_min[2]; J2p<=J2p_max[2]; J2p++)
             recouple_lmn[2].push_back( -Z.modelspace->phase((om.j2+on.j2)/2+J2+J2p)*sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p) );


        double zijklmn = 0;
/// BEGIN THE SLOW BIT...

        for ( int perm_ijk=0; perm_ijk<3; perm_ijk++ )
        {
          size_t I1 = ijk[perm_ijk][0];
          size_t I2 = ijk[perm_ijk][1];
          size_t I3 = ijk[perm_ijk][2];
          Orbit& o1 = Z.modelspace->GetOrbit( I1 );
          Orbit& o2 = Z.modelspace->GetOrbit( I2 );
          Orbit& o3 = Z.modelspace->GetOrbit( I3 );

          double j3 = 0.5*o3.j2;
          for (int J1p=J1p_min[perm_ijk]; J1p<=J1p_max[perm_ijk]; J1p++)
          {
            if ((I1==I2) and J1p%2!=0 ) continue;
            double rec_ijk = recouple_ijk[perm_ijk].at(J1p-J1p_min[perm_ijk]);


            for ( int perm_lmn=0; perm_lmn<3; perm_lmn++ )
            {
              size_t I4 = lmn[perm_lmn][0];
              size_t I5 = lmn[perm_lmn][1];
              size_t I6 = lmn[perm_lmn][2];
              Orbit& o4 = Z.modelspace->GetOrbit( I4 );
              Orbit& o5 = Z.modelspace->GetOrbit( I5 );
              Orbit& o6 = Z.modelspace->GetOrbit( I6 );

              size_t key_126 = hash_key_ijnJ(std::min(I1,I2),std::max(I1,I2),I6,J1p);

              double j6 = 0.5*o6.j2;
              for (int J2p=J2p_min[perm_lmn]; J2p<=J2p_max[perm_lmn]; J2p++)
              {
                if ((I4==I5) and J2p%2!=0 ) continue;
                double rec_lmn = recouple_lmn[perm_lmn].at(J2p-J2p_min[perm_lmn]);

                int phase_12 = I1>I2 ? -Z.modelspace->phase( (o1.j2+o2.j2)/2 - J1p) : 1;
                int phase_45 = I4>I5 ? -Z.modelspace->phase( (o4.j2+o5.j2)/2 - J2p) : 1;
                double hat_factor = sqrt( (2*J1p+1)*(2*J2p+1) );
                double z_123456 = 0;

                size_t key_453 = hash_key_ijnJ(std::min(I4,I5),std::max(I4,I5),I3,J2p);

                int parity_ph = (o1.l+o2.l+o6.l)%2;
                int twoTz_ph = (o1.tz2+o2.tz2-o6.tz2);
                int twoJph_min = std::max( std::abs(2*J1p-o6.j2), std::abs(2*J2p-o3.j2) );
                int twoJph_max = std::min( 2*J1p+o6.j2 , 2*J2p+o3.j2 );

// This inner loop is slow
                for ( int twoJph=twoJph_min; twoJph<=twoJph_max; twoJph++)
                {
                   double zbar_126453=0;
                   double Jph = 0.5 * twoJph;


                   // This lookup is expensive. Replace map with vector makes things better
                   int ch_pph = channels_pph[ ch_pph_hash(twoJph,parity_ph,twoTz_ph)]; 
                   if ( ch_pph<0 ) continue;
 

                   size_t index_126 = ket_lookup_pph[ch_pph].at(key_126 );
                   size_t index_453 = ket_lookup_pph[ch_pph].at(key_453 );

                   size_t start_ptr = Zbar_start_pointers[ch_pph];
                   size_t ngood_ijn = ket_lookup_pph[ch_pph].size();
                   if ( index_126<= index_453 )
                   {
                     size_t zbar_index = start_ptr + ( 2*ngood_ijn - index_126 - 1) * index_126 / 2  + index_453 ;
                     zbar_126453 = Zbar[ zbar_index ]  * phase_12 * phase_45;
                   }
                   else
                   {
                     size_t zbar_index = start_ptr + ( 2*ngood_ijn - index_453 - 1) * index_453 / 2  + index_126 ;
                     zbar_126453 = hZ * Zbar[ zbar_index ]  * phase_12 * phase_45;
                   }

                   double sixj1 = Z.modelspace->GetSixJ(J1p,j3,Jtot, J2p,j6,Jph);
                   z_123456 +=   sixj1 * zbar_126453  ; 

                }

                zijklmn += hat_factor * rec_ijk * rec_lmn * z_123456;

              }//for J2p
            }// for perm_lmn
          }// for J1p
        }// for perm_ijk


        Z3.AddToME_pn_ch( ch3,ch3,ibra,iket, zijklmn );
      }// for iket
    }// for ibra
  }// for ch3
  Z.profiler.timer["comm223_compute_loop"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();
  if (Z.modelspace->scalar3b_transform_first_pass)
  Z.profiler.timer["comm223_firstpass"] += omp_get_wtime() - tstart;
  else
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}



void comm223ss_new( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
  double t_internal = omp_get_wtime();
  int norbs = Z.modelspace->GetNumberOrbits();
  auto& Z3 = Z.ThreeBody;
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  int hX = X.IsHermitian() ? 1 : -1;
  int hY = Y.IsHermitian() ? 1 : -1;
  int hZ = Z.IsHermitian() ? 1 : -1;
  if ( (std::abs( X2.Norm() * Y2.Norm() ) < 1e-6 ) and not Z.modelspace->scalar3b_transform_first_pass) return;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();

  // create a function object for hashing 4 size_t (unsigned long int) to a single size_t for lookup in an unordered_map
  int max_twoJph = 6*Z.modelspace->GetEmax()+3;
  auto ch_pph_hash = [](int twoJ_ph, int parity, int twoTz_ph){ return 8*(twoJ_ph/2) + parity + (twoTz_ph+3);};
  std::vector<int> channels_pph(ch_pph_hash(max_twoJph,1,3)+1, -1); // map twoJ_ph, parity_ph, twoTz_ph  -> channel index
  std::vector< std::array<int,3>> channel_list_pph; // vector  channel_index -> twoJ-ph, parity_ph, twoTz_ph.

  auto hash_key_ijnJ = [](size_t i,size_t j, size_t n, size_t Jij){return ( i + (j<<12) + (n<<24) + (Jij<<36) );};

   //  0xF = 15 = 1111 (4bits), so 0xFFF is 12 bits of 1's. 0xFFFL makes it a long int
  auto unhash_key_ijnJ = [](size_t& i,size_t& j, size_t& n, size_t& Jij, size_t key){ i=( key     &0xFFFL);
                                                                                      j=((key>>12)&0xFFFL);
                                                                                      n=((key>>24)&0xFFFL);
                                                                                    Jij=((key>>36)&0xFFFL); }; 
  std::vector< std::unordered_map<size_t, size_t>> ket_lookup_pph; // in a given channel, map  i,j,n,Jij -> matrix index

  std::vector<float> Zbar;  // store the transformed Z in a 1D vector.
  std::vector<size_t> Zbar_start_pointers;  // ch -> index in Zbar where that channel starts
  std::unordered_map<size_t,size_t> Zbar_start_pointers_ijklmn;  // ijklmn => index in Zbar where that channel starts

  // Different approach: store by ijklmn index, then for each ijklmn combo, loop over Jij,Jlm,twoJ
  // In this case, we need a map {i,j,k,l,m,n} => size_t, where the size_t points to the start locating in Zbar

  // This is the hash {i,j,k,l,m,n} => size_t, for use in Zbar_start_pointers, or maybe a different name?
  auto hash_key_ijklmn = [](size_t i, size_t j, size_t k, size_t l, size_t m, size_t n){ return i+(j<<10)+(k<<20)+(l<<30)+(m<<40)+(n<<50);};

  // we're looking at pph type states. We don't make a cut on the occupations in this channel
  // but if the first two orbits ij can't possibly make it past the occnat cut later on, we don't bother including them
  // To figure out if they'll make it, we want to know the biggest contribution that can possibly be made by a third orbit.
  double occnat_factor_max = 0;
  for ( auto i : Z.modelspace->all_orbits)
  {
     double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
     occnat_factor_max = std::max( occnat_factor_max, occnat_i*(1-occnat_i) );
  }

//  std::cout << "BeginLoop over ph channels" << std::endl;
  size_t nch_pph = 0;
  size_t Zbar_n_elements = 0; // total number of elements stored in Zbar
  // Generate all the pph type channels J,parity,Tz. Note that these will be
  // not the same as the standard 3-body channels since we aren't enforcing antisymmetry with the 3rd particle
  for ( auto& iter_obc : Z.modelspace->OneBodyChannels )
  {
    int parity_ph = iter_obc.first[0]%2;  // parity = l%2
    int twoJph = iter_obc.first[1];
    int twoTz_ph = iter_obc.first[2];

    size_t ngood_ijn=0;
    std::unordered_map<size_t, size_t> good_kets_ijn;
    for ( size_t chij=0; chij<nch2; chij++)
    {
      TwoBodyChannel& tbc_ij = Z.modelspace->GetTwoBodyChannel(chij);
      int Jij = tbc_ij.J;
      int parity_ij = tbc_ij.parity;
      int Tz_ij = tbc_ij.Tz;
      size_t nkets_ij = tbc_ij.GetNumberKets();
      int tz2n = 2*Tz_ij - twoTz_ph;
      if ( std::abs(tz2n)>1 ) continue;
      int j2n_min = std::abs(2*Jij-twoJph);
      int j2n_max = (2*Jij+twoJph);
      int parity_n = (parity_ij + parity_ph)%2;
      
      std::vector<size_t> good_n;
      for ( int j2n=j2n_min; j2n<=j2n_max; j2n+=2)
      {
        int l_n = ((j2n+1)/2)%2 == parity_n ? (j2n+1)/2  :  (j2n-1)/2;
        if (l_n > Z.modelspace->GetEmax() or l_n > Z.modelspace->GetLmax() ) continue;
        for ( size_t n : Z.OneBodyChannels.at({l_n,j2n,tz2n}) )
        {
          Orbit& on = Z.modelspace->GetOrbit(n);
          double occnat_n = on.occ_nat;
          double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
          if ( (occnat_n*(1-occnat_n) * occnat_factor_max * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
          if ( (d_en) > Z.modelspace->dE3max ) continue;
          good_n.push_back(n);
        }// for n
      }// for j2n

      std::vector<std::array<size_t,2>> good_ij;
      for ( size_t iket_ij=0; iket_ij<nkets_ij; iket_ij++)
      {
        // check that we pass our cuts on E3max, occupations, etc.
        Ket& ket_ij = tbc_ij.GetKet(iket_ij);
        size_t i = ket_ij.p;
        size_t j = ket_ij.q;
        double occnat_i = ket_ij.op->occ_nat;
        double occnat_j = ket_ij.oq->occ_nat;
        int ei = 2*ket_ij.op->n + ket_ij.op->l;
        int ej = 2*ket_ij.oq->n + ket_ij.oq->l;
        double d_ei = std::abs( ei - e_fermi[ket_ij.op->tz2]);
        double d_ej = std::abs( ej - e_fermi[ket_ij.oq->tz2]);
        // if i and j cant make it past the OccNat and dE3max cuts, don't bother including it
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( (d_ei+d_ej) > Z.modelspace->GetdE3max() ) continue;
        if ( (ei + ej) > Z.modelspace->GetE3max() ) continue;
        if ( perturbative_triples and ( (ket_ij.op->cvq==0 and ket_ij.oq->cvq!=0) or (ket_ij.op->cvq!=0 and ket_ij.oq->cvq==0)) ) continue;
        good_ij.push_back({i,j}); 
      }//for iket_ij

      for ( auto& ij : good_ij )
      {
        for ( size_t n : good_n )
        {
          if ( perturbative_triples ) // we just want <ppp|hhh> terms in the end, so for ijn we want pph or hhp.
          {
            Orbit& oi = Z.modelspace->GetOrbit(ij[0]);
            Orbit& oj = Z.modelspace->GetOrbit(ij[1]);
            Orbit& on = Z.modelspace->GetOrbit(n);
            if ( (oi.cvq==0 and oj.cvq==0 and on.cvq==0) or ( oi.cvq!=0 and oj.cvq!=0 and on.cvq!=0) ) continue;
          }
          size_t key = hash_key_ijnJ( ij[0],ij[1],n,size_t(Jij));
          good_kets_ijn[ key ] = ngood_ijn;
          ngood_ijn++;
        }
      }// for iket_ij
    }// for chij
    if ( (ngood_ijn < 1) ) continue;

    channels_pph[ ch_pph_hash(twoJph,parity_ph,twoTz_ph) ] = nch_pph;

    channel_list_pph.push_back({twoJph, parity_ph,twoTz_ph});
    ket_lookup_pph.push_back( good_kets_ijn );
    Zbar_start_pointers.push_back( Zbar_n_elements ) ;     
    Zbar_n_elements += ngood_ijn * (ngood_ijn+1) / 2 ;
    nch_pph++;
  }// for iter_obc

  // Layout of Zbar:  {ijklmn}, twoJ, Jij, Jlm
  //  looping over ijklmn, we keep i<=j,  l<=m,  i<=l, and if i==l then j<=m, and if i==l,j==m, then k<=n.
  // for for a given ijklmn,
  // twoJmax = min( ji+jj+jk,  jl+jm+jn)  and twoJmin = max(twoJmin_ijk, twoJmin_lmn)
  //    with twoJmin_ijk =  max ({abs(ji-jj)-k, jk-ji-jj, 1}) 
// 
//   then Jij_min = abs(twoJ-jk),  Jij_max = twoJ+jk,   and likewise for Jlm.
  int max_single_j2 = 2*Z.modelspace->GetEmax() + 1;

  size_t norbits = Z.modelspace->GetNumberOrbits();

  size_t ngood_ijklmn = 0;
  for ( size_t i=0; i<norbits; i++ )
  {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    int ei = 2*oi.n+oi.l;
    double d_ei = std::abs( ei - e_fermi[oi.tz2]);
    double occnat_i = oi.occ_nat;
    for ( size_t j=i; j<norbits; j++)
    {
      Orbit& oj = Z.modelspace->GetOrbit(j);
      int ej = 2*oj.n+oj.l;
      double d_ej = std::abs( ej - e_fermi[oj.tz2]);
      double occnat_j = oj.occ_nat;

      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
      if ( (d_ei+d_ej) > Z.modelspace->GetdE3max() ) continue;
      if ( (ei + ej) > Z.modelspace->GetE3max() ) continue;
      if ( perturbative_triples and ( (oi.cvq==0 and oj.cvq!=0) or (oi.cvq!=0 and oj.cvq==0)) ) continue;
      
      for (size_t k=0; k<norbits; k++)
      {
        Orbit& ok = Z.modelspace->GetOrbit(k);
        int ek = 2*ok.n+ok.l;
        double d_ek = std::abs( ek - e_fermi[ok.tz2]);
        double occnat_k = ok.occ_nat;

        if ( std::abs( oi.tz2+oj.tz2-ok.tz2 ) >1 ) continue;
        
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( (d_ei+d_ej+d_ek) > Z.modelspace->GetdE3max() ) continue;
        if ( (ei + ej + ek) > Z.modelspace->GetE3max() ) continue;
        if ( perturbative_triples and ( (oi.cvq==0 and oj.cvq==0 and ok.cvq==0) or ( oi.cvq!=0 and oj.cvq!=0 and ok.cvq!=0) ) ) continue;

        int twoJmin_ijk = std::max( { std::abs(oi.j2-oj.j2)-ok.j2, ok.j2-oi.j2-oj.j2, 1});
        int twoJmax_ijk = oi.j2+oj.j2+ok.j2;

        // now loop over the ket orbits lmn
        for (size_t l=i; l<norbits; l++)
        {
          Orbit& ol = Z.modelspace->GetOrbit(l);
          int el = 2*ol.n+ol.l;
          double d_el = std::abs( el - e_fermi[ol.tz2]);
          double occnat_l = ol.occ_nat;

          for (size_t m=( (i==l) ? j : l); m<norbits; m++)
          {
            Orbit& om = Z.modelspace->GetOrbit(m);
            int em = 2*om.n+om.l;
            double d_em = std::abs( em - e_fermi[om.tz2]);
            double occnat_m = om.occ_nat;

            if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
            if ( (d_el+d_em) > Z.modelspace->GetdE3max() ) continue;
            if ( (el + em) > Z.modelspace->GetE3max() ) continue;
            if ( perturbative_triples and ( (ol.cvq==0 and om.cvq!=0) or (ol.cvq!=0 and om.cvq==0)) ) continue;

            for (size_t n=( (i==l and j==m) ? k : 0 ); n<norbits; n++)
            {
              Orbit& on = Z.modelspace->GetOrbit(n);
              int en = 2*on.n+on.l;
              double d_en = std::abs( en - e_fermi[on.tz2]);
              double occnat_n = on.occ_nat;
              if ( (oi.l+oj.l+ok.l + ol.l+om.l+on.l)%2>0) continue;
              if ( (oi.tz2+oj.tz2-ok.tz2) != (ol.tz2+om.tz2-on.tz2) ) continue; // Note the minus signs! k and n are hole states

              if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
              if ( (d_el+d_em+d_en) > Z.modelspace->GetdE3max() ) continue;
              if ( (el + em + en) > Z.modelspace->GetE3max() ) continue;
              if ( perturbative_triples and ( (ol.cvq==0 and om.cvq==0 and on.cvq==0) or ( ol.cvq!=0 and om.cvq!=0 and on.cvq!=0) ) ) continue;
              
              int twoJmin_lmn = std::max( { std::abs(ol.j2-om.j2)-on.j2, on.j2-ol.j2-om.j2, 1} );
              int twoJmax_lmn = ol.j2+om.j2+on.j2;
              int twoJmin = std::max(twoJmin_ijk,twoJmin_lmn);
              int twoJmax = std::min({ twoJmax_ijk, twoJmax_lmn, max_single_j2});
              if (twoJmax>=twoJmin)
              {
                 size_t key = hash_key_ijklmn(i,j,k,l,m,n);
                 Zbar_start_pointers_ijklmn[ key ] = ngood_ijklmn;
//                 std::cout << " HASHING " << i << " " << j << " " << k << " " << l << " " << m << " " << n << "   to " << key << "   and pointing to " << ngood_ijklmn << std::endl;
              }
              for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
              {
                   int parity_ph = (oi.l + oj.l + ok.l)%2;
                   if ( (twoJ == max_single_j2) and (parity_ph != Z.modelspace->GetEmax()%2)) continue;
                   int twoTz_ph = (oi.tz2 + oj.tz2 - ok.tz2);
                   int ch_pph = channels_pph[ ch_pph_hash(twoJ,parity_ph,twoTz_ph)]; 
                int Jij_min = std::max( std::abs(oi.j2-oj.j2), std::abs( twoJ-ok.j2) )/2;
                int Jij_max = std::min( oi.j2+oj.j2, twoJ+ok.j2)/2;
                int Jlm_min = std::max( std::abs(ol.j2-om.j2), std::abs( twoJ-on.j2) )/2;
                int Jlm_max = std::min( ol.j2+om.j2, twoJ+on.j2)/2;
                for (int Jij=Jij_min; Jij<=Jij_max; Jij++)
                {
                  if ( i==j and Jij%2>0) continue;
                    size_t key_ijk = hash_key_ijnJ(std::min(i,j),std::max(i,j),k,Jij);
                     if ( ch_pph == -1 or ch_pph >= (int)ket_lookup_pph.size() )
                     {
                       std::cout << "TROUBLE line " << __LINE__ << " ijk = " << i << " " << j << " " << k << "   " << Jij
                                 << "   twoJ parity twoTz " << twoJ << " " << parity_ph << " " << twoTz_ph << "  ch = " << ch_pph << std::endl;
                     }
                     size_t index_ijk = ket_lookup_pph[ch_pph].at(key_ijk );


                  for (int Jlm=Jlm_min; Jlm<=Jlm_max; Jlm++)
                  {
                    if ( l==m and Jlm%2>0) continue;
                    if ( i==l and j==m and k==n and Jlm<Jij) continue;
                    ngood_ijklmn++;

                    size_t key_lmn = hash_key_ijnJ(std::min(l,m),std::max(l,m),n,Jlm);
                     size_t index_lmn = ket_lookup_pph[ch_pph].at(key_lmn );


                  }// for Jlm
                }// for Jij
              }// for twoJ

            }// for n
          }// for m
        }// for l
      }// for k
    }// for j
  }// for i
  std::cout << "Counted ngood_ijklmn = " << ngood_ijklmn << std::endl;



  Zbar.resize( Zbar_n_elements, 0. );
  std::cout << __func__ << "  allocating Zbar   " << Zbar_n_elements << " elements  ~ " << Zbar_n_elements * sizeof(float) /(1024.*1024*1024) << " GB" << std::endl;
  
  Z.profiler.timer["comm223_setup_loop"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();
//  std::cout << "done with setup" << std::endl;


 // Now we should fill those Zbar matrices.
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
 for (size_t ch_pph=0; ch_pph<nch_pph; ch_pph++)
 {

   auto& ch_info = channel_list_pph[ch_pph];
   int twoJph    = ch_info[0];
   int parity_ph = ch_info[1];
   int twoTz_ph  = ch_info[2];
   int l_a = ((twoJph+1)/2)%2 == parity_ph ? (twoJph+1)/2  :  (twoJph-1)/2;


//   auto& a_set = Z.modelspace->OneBodyChannels.at({l_a,twoJph,twoTz_ph}); // list of one-body orbits in this channel
   std::vector<size_t> a_list;
   for (auto a : Z.modelspace->OneBodyChannels.at({l_a,twoJph,twoTz_ph}) ) a_list.push_back(a);
   size_t number_a = a_list.size();

   auto& good_kets_ijn = ket_lookup_pph[ch_pph];
   size_t ngood_ijn = good_kets_ijn.size();
   arma::mat Xph( ngood_ijn, number_a );
   arma::mat Yph( number_a, ngood_ijn );
   for (auto& iter_ijn : good_kets_ijn )
   {
      size_t Iijn = iter_ijn.second;

      size_t key = iter_ijn.first;
      size_t i,j,n,Jij_tmp;
      unhash_key_ijnJ(i,j,n,Jij_tmp, key);
      int Jij = (int)Jij_tmp;

      for (size_t index_a=0; index_a<number_a; index_a++ )
      {
        size_t a = a_list[index_a];
        Orbit& oa = Z.modelspace->GetOrbit(a);
        Xph(Iijn, index_a) = X.TwoBody.GetTBME_J(Jij,i,j,n,a);
        Yph(index_a, Iijn) = Y.TwoBody.GetTBME_J(Jij,n,a,i,j);

      }// for index_a

   }//for iter_ijn

   // Do the matrix multiplication
   arma::mat z_tmp = Xph*Yph;
   z_tmp -= hX*hY * z_tmp.t();


   for (auto& iter_ijn : good_kets_ijn )
   {
      size_t Iijn = iter_ijn.second;

      size_t key_ijn = iter_ijn.first;
      size_t i,j,n,Jij_tmp;
      unhash_key_ijnJ(i,j,n,Jij_tmp, key_ijn);
      int Jij_this = (int)Jij_tmp;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& on = Z.modelspace->GetOrbit(n);
     
      int twoJmin_ijn = std::max( { std::abs(oi.j2-oj.j2)-on.j2, on.j2-oi.j2-oj.j2, 1});
      int twoJmax_ijn = oi.j2+oj.j2+on.j2;

 
      for (auto& iter_lmk : good_kets_ijn )
      {
         size_t Ilmk = iter_lmk.second;
//         if (Ilmk<Iijn) continue;
         size_t key_lmk = iter_lmk.first;
         size_t l,m,k,Jlm_tmp;
         unhash_key_ijnJ(l,m,k,Jlm_tmp, key_lmk);
         if ( (l<i)  or (l==i and m<j) or (l==i and m==j and k<n) ) continue;
         int Jlm_this = (int)Jlm_tmp;
         Orbit& ol = Z.modelspace->GetOrbit(l);
         Orbit& om = Z.modelspace->GetOrbit(m);
         Orbit& ok = Z.modelspace->GetOrbit(k);
         

         size_t key_ijnlmk = hash_key_ijklmn(i,j,n,l,m,k);
         if ( Zbar_start_pointers_ijklmn.find( key_ijnlmk) == Zbar_start_pointers_ijklmn.end() )
         {
           std::cout << "OOPS " << __LINE__ << " II: " << Iijn << " " << Ilmk << "  ijnlmk " << i << " " << j << " " << n << " " << l << " " << m << " " << k << " parity: " << (oi.l+oj.l+on.l)%2 << " vs " << (ol.l+om.l+ok.l)%2 << "  Tz " << (oi.tz2+oj.tz2-on.tz2) << " vs " << (ol.tz2+om.tz2-ok.tz2)  << std::endl;
         }
         size_t start_ptr = Zbar_start_pointers_ijklmn[ key_ijnlmk ];
         

         int twoJmin_lmk = std::max( { std::abs(ol.j2-om.j2)-on.j2, on.j2-ol.j2-om.j2, 1} );
         int twoJmax_lmk = ol.j2+om.j2+on.j2;
         int twoJmin = std::max(twoJmin_ijn,twoJmin_lmk);
         int twoJmax = std::min({ twoJmax_ijn, twoJmax_lmk, max_single_j2});

         bool found_it = false;
         size_t counter=0;

         // wind through the J values to find the one we want
         for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
         {
           int parity_ph = (oi.l + oj.l + on.l)%2;
           if ( (twoJ == max_single_j2) and (parity_ph != Z.modelspace->GetEmax()%2)) continue;
           int Jij_min = std::max( std::abs(oi.j2-oj.j2), std::abs( twoJ-on.j2) )/2;
           int Jij_max = std::min( oi.j2+oj.j2, twoJ+ok.j2)/2;
           int Jlm_min = std::max( std::abs(ol.j2-om.j2), std::abs( twoJ-ok.j2) )/2;
           int Jlm_max = std::min( ol.j2+om.j2, twoJ+on.j2)/2;
           for (int Jij=Jij_min; Jij<=Jij_max; Jij++)
           {
             if ( i==j and Jij%2>0) continue;

             for (int Jlm=Jlm_min; Jlm<=Jlm_max; Jlm++)
             {
               if ( l==m and Jlm%2>0) continue;
               if ( i==l and j==m and k==n and Jlm<Jij) continue;
               if ( (twoJ==twoJph) and (Jij==Jij_this) and (Jlm==Jlm_this) )
               {
                  if ( (start_ptr + counter) >= Zbar.size() )
                  {
                      std::cout << "TROUBLE " << __LINE__ << "  ijnlmk " << i << " " << j << " " << n << " " << l << " " << m << " " << k << std::endl;
                  }
                  Zbar[ start_ptr + counter ] = z_tmp(Iijn,Ilmk);
                  found_it = true;
               }
               counter++;
               if (found_it) break;

             }// for Jlm
               if (found_it) break;
           }// for Jij
               if (found_it) break;
         }// for twoJ


      }
   }


//   // Fold the full (anti-)symmetric matrix into a more compact storage
//   // we  use row-major ordering so that the inner loop runs over adjacent elements
//   size_t start_ptr = Zbar_start_pointers[ch_pph];
//   for (size_t II=0; II<ngood_ijn; II++)
//   {
//     size_t row_index = start_ptr + ( 2*ngood_ijn - II - 1) * II / 2 ;
//     for (size_t JJ=II; JJ<ngood_ijn; JJ++)
//     {
//       Zbar[ row_index+JJ ] = z_tmp(II,JJ);
//     }
//   }

 }// for ch_pph

// std::cout << "Done filling matrices " << std::endl;


  Z.profiler.timer["comm223_fill_loop"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      if ( perturbative_triples and  not ( (oi.cvq + oj.cvq + ok.cvq)==0 or (oi.cvq>0 and oj.cvq>0 and ok.cvq>0)) ) continue;
      int J1 = bra.Jpq;

      // Set up the permutation stuff for ijk
      std::vector<std::array<size_t,3>> ijk = { {i,j,k}, {k,j,i}, {i,k,j} };
      std::vector<int> J1p_min = {J1,  std::max(std::abs(ok.j2-oj.j2),std::abs(twoJ-oi.j2) )/2,   std::max(std::abs(oi.j2-ok.j2), std::abs(twoJ-oj.j2) )/2 };
      std::vector<int> J1p_max = {J1,  std::min(ok.j2+oj.j2, twoJ+oi.j2)/2 , std::min(oi.j2+ok.j2, twoJ+oj.j2)/2 };
      std::vector<std::vector<double>> recouple_ijk = {{1},{},{} };

      for (int J1p=J1p_min[1]; J1p<=J1p_max[1]; J1p++)
           recouple_ijk[1].push_back( sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p) );

      for (int J1p=J1p_min[2]; J1p<=J1p_max[2]; J1p++)
           recouple_ijk[2].push_back( -Z.modelspace->phase((oj.j2+ok.j2)/2+J1+J1p)*sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p) );

      for (size_t iket=ibra; iket<nkets3; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( perturbative_triples and  not ( (ol.cvq + om.cvq + on.cvq)==0 or (ol.cvq>0 and om.cvq>0 and on.cvq>0)) ) continue;
        if ( perturbative_triples and (  (oi.cvq==0 and ol.cvq==0) or (oi.cvq!=0 and ol.cvq!=0) ) ) continue;
        int J2 = ket.Jpq;

        // Set up the permutation stuff for lmn
        std::vector<std::array<size_t,3>> lmn = { {l,m,n}, {n,m,l}, {l,n,m} };
        std::vector<int> J2p_min = {J2,  std::max(std::abs(on.j2-om.j2), std::abs(twoJ-ol.j2) )/2,   std::max( std::abs(ol.j2-on.j2), std::abs(twoJ-om.j2) )/2 };
        std::vector<int> J2p_max = {J2,  std::min(on.j2+om.j2, twoJ+ol.j2)/2 , std::min(ol.j2+on.j2,twoJ+om.j2)/2 };
        std::vector<std::vector<double>> recouple_lmn = {{1},{},{} };
            
        for (int J2p=J2p_min[1]; J2p<=J2p_max[1]; J2p++)
           recouple_lmn[1].push_back( sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p) );

        for (int J2p=J2p_min[2]; J2p<=J2p_max[2]; J2p++)
             recouple_lmn[2].push_back( -Z.modelspace->phase((om.j2+on.j2)/2+J2+J2p)*sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p) );


        double zijklmn = 0;
/// BEGIN THE SLOW BIT...

        for ( int perm_ijk=0; perm_ijk<3; perm_ijk++ )
        {
          size_t I1 = ijk[perm_ijk][0];
          size_t I2 = ijk[perm_ijk][1];
          size_t I3 = ijk[perm_ijk][2];
          Orbit& o1 = Z.modelspace->GetOrbit( I1 );
          Orbit& o2 = Z.modelspace->GetOrbit( I2 );
          Orbit& o3 = Z.modelspace->GetOrbit( I3 );

          double j3 = 0.5*o3.j2;


            for ( int perm_lmn=0; perm_lmn<3; perm_lmn++ )
            {
              size_t I4 = lmn[perm_lmn][0];
              size_t I5 = lmn[perm_lmn][1];
              size_t I6 = lmn[perm_lmn][2];
              Orbit& o4 = Z.modelspace->GetOrbit( I4 );
              Orbit& o5 = Z.modelspace->GetOrbit( I5 );
              Orbit& o6 = Z.modelspace->GetOrbit( I6 );


              double j6 = 0.5*o6.j2;

              int parity_ph = (o1.l+o2.l+o6.l)%2;
              if ( std::abs( o1.tz2+o2.tz2-o6.tz2 ) >1 ) continue;
              if ( ( o1.tz2+o2.tz2-o6.tz2 ) !=  ( o4.tz2+o5.tz2-o3.tz2 ) ) continue;


              int twoJmin_126 = std::max( { std::abs(o1.j2-o2.j2)-o6.j2, o6.j2-o1.j2-o2.j2, 1});
              int twoJmax_126 = o1.j2+o2.j2+o6.j2;
              int twoJmin_453 = std::max( { std::abs(o4.j2-o5.j2)-o3.j2, on.j2-ol.j2-om.j2, 1} );
              int twoJmax_453 = o4.j2+o5.j2+o3.j2;
              int twoJph_min = std::max(twoJmin_126,twoJmin_453);
              int twoJph_max = std::min({ twoJmax_126, twoJmax_453, max_single_j2});

              size_t key_126453;
              int phase_126453 = 1;
              if ( (I4<I1)  or (I4==I1 and I5<I2) or (I4==I1 and I5==I2 and I3<I6) ) 
              {
                key_126453 = hash_key_ijklmn( std::min(I4,I5),std::max(I4,I5),I3 , std::min(I1,I2),std::max(I1,I2),I6 );
                phase_126453 *= hZ;
              }
              else
              {
                key_126453 = hash_key_ijklmn(std::min(I1,I2),std::max(I1,I2),I6,std::min(I4,I5),std::max(I4,I5),I3);
              }
              size_t start_ptr = Zbar_start_pointers_ijklmn.at( key_126453 );

                size_t counter = 0;
                for ( int twoJph=twoJph_min; twoJph<=twoJph_max; twoJph++)
                {
                   double zbar_126453=0;
                   double Jph = 0.5 * twoJph;
                   if ( (twoJph == max_single_j2) and (parity_ph != Z.modelspace->GetEmax()%2)) continue;


                   int J12_min = std::max( std::abs(o1.j2-o2.j2), std::abs( twoJph-o6.j2) )/2;
                   int J12_max = std::min( o1.j2+o2.j2, twoJph+o6.j2)/2;
                   int J45_min = std::max( std::abs(o4.j2-o5.j2), std::abs( twoJph-o3.j2) )/2;
                   int J45_max = std::min( o4.j2+o5.j2, twoJph+o3.j2)/2;

                   for (int J12=J12_min; J12<=J12_max; J12++)
                   {
                     if ( I1==I2 and J12%2>0) continue;
                     int phase_12 = I1>I2 ? -Z.modelspace->phase( (o1.j2+o2.j2)/2 - J12) : 1;

                     double rec_ijk = recouple_ijk[perm_ijk].at(J12-J1p_min[perm_ijk]);
                     for (int J45=J45_min; J45<=J45_max; J45++)
                     {
                       if ( I4==I5 and J45%2>0) continue;
//                       if ( I1==I4 and I2==I5 and I6==I3 and J45<J12) continue; // maybe we don't want this...
                     if (  AngMom::Triangle( J12,j3,Jtot) and AngMom::Triangle( J12, j6,Jph) and   
                           AngMom::Triangle( J45,j6,Jtot) and AngMom::Triangle( J45, j3,Jph) ) 
                     { 
                       double rec_lmn = recouple_lmn[perm_lmn].at(J45-J2p_min[perm_lmn]);
                       int phase_45 = I4>I5 ? -Z.modelspace->phase( (o4.j2+o5.j2)/2 - J45) : 1;

                       zbar_126453 = Zbar[ start_ptr + counter ] * phase_126453 * phase_12 * phase_45;

                       double hat_factor = sqrt( (2*J12+1.)*(2*J45+1) );
                       double sixj1 = Z.modelspace->GetSixJ(J12,j3,Jtot, J45,j6,Jph);
//                       z_123456 +=   sixj1 * zbar_126453  ; 
                       double z_123456 =   sixj1 * zbar_126453  ; 

                       zijklmn += hat_factor * rec_ijk * rec_lmn * z_123456;
                     }
                       counter++;

                     }// for J45
                   }// for J12
               }// for twoJph

            }// for perm_lmn
        }// for perm_ijk


//          for (int J1p=J1p_min[perm_ijk]; J1p<=J1p_max[perm_ijk]; J1p++)
//          {
//            if ((I1==I2) and J1p%2!=0 ) continue;
//            double rec_ijk = recouple_ijk[perm_ijk].at(J1p-J1p_min[perm_ijk]);
//
//
//              size_t key_126 = hash_key_ijnJ(std::min(I1,I2),std::max(I1,I2),I6,J1p);
//
//
//              for (int J2p=J2p_min[perm_lmn]; J2p<=J2p_max[perm_lmn]; J2p++)
//              {
//                if ((I4==I5) and J2p%2!=0 ) continue;
//                double rec_lmn = recouple_lmn[perm_lmn].at(J2p-J2p_min[perm_lmn]);
//
//                int phase_12 = I1>I2 ? -Z.modelspace->phase( (o1.j2+o2.j2)/2 - J1p) : 1;
//                int phase_45 = I4>I5 ? -Z.modelspace->phase( (o4.j2+o5.j2)/2 - J2p) : 1;
//                double hat_factor = sqrt( (2*J1p+1)*(2*J2p+1) );
//                double z_123456 = 0;
//
//                size_t key_453 = hash_key_ijnJ(std::min(I4,I5),std::max(I4,I5),I3,J2p);
//
//                int parity_ph = (o1.l+o2.l+o6.l)%2;
//                int twoTz_ph = (o1.tz2+o2.tz2-o6.tz2);
//                int twoJph_min = std::max( std::abs(2*J1p-o6.j2), std::abs(2*J2p-o3.j2) );
////                int twoJph_max = std::min( 2*J1p+o6.j2 , 2*J2p+o3.j2 );
//
//                int twoJph_max = std::min( { 2*J1p+o6.j2 , 2*J2p+o3.j2, max_single_j2 } );
//                if ( std::abs(twoTz_ph)>1 ) continue;
//
//// This inner loop is slow
////                for ( int twoJph=twoJph_min; twoJph<=twoJph_max; twoJph++)
////                {
////                   double zbar_126453=0;
////                   double Jph = 0.5 * twoJph;
//
//
//                   // This lookup is expensive. Replace map with vector makes things better
//                   int ch_pph = channels_pph[ ch_pph_hash(twoJph,parity_ph,twoTz_ph)]; 
//                   if ( ch_pph<0 ) continue;
// 
//
//                   size_t index_126 = ket_lookup_pph[ch_pph].at(key_126 );
//                   size_t index_453 = ket_lookup_pph[ch_pph].at(key_453 );
//
//                   size_t start_ptr = Zbar_start_pointers[ch_pph];
//                   size_t ngood_ijn = ket_lookup_pph[ch_pph].size();
//                   if ( index_126<= index_453 )
//                   {
//                     size_t zbar_index = start_ptr + ( 2*ngood_ijn - index_126 - 1) * index_126 / 2  + index_453 ;
//                     zbar_126453 = Zbar[ zbar_index ]  * phase_12 * phase_45;
//                   }
//                   else
//                   {
//                     size_t zbar_index = start_ptr + ( 2*ngood_ijn - index_453 - 1) * index_453 / 2  + index_126 ;
//                     zbar_126453 = hZ * Zbar[ zbar_index ]  * phase_12 * phase_45;
//                   }
//
//                   double sixj1 = Z.modelspace->GetSixJ(J1p,j3,Jtot, J2p,j6,Jph);
//                   z_123456 +=   sixj1 * zbar_126453  ; 
//
//                }
//
//                zijklmn += hat_factor * rec_ijk * rec_lmn * z_123456;
//
//              }//for J2p
//          }// for J1p
//            }// for perm_lmn
//        }// for perm_ijk


        Z3.AddToME_pn_ch( ch3,ch3,ibra,iket, zijklmn );
      }// for iket
    }// for ibra
  }// for ch3
  Z.profiler.timer["comm223_compute_loop"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();
  if (Z.modelspace->scalar3b_transform_first_pass)
  Z.profiler.timer["comm223_firstpass"] += omp_get_wtime() - tstart;
  else
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}







bool check_2b_channel_Tz_parity( const Operator& Op, Orbit& o1, Orbit&o2, Orbit& o3, Orbit& o4 ) 
{
  return (    ((o1.l+o2.l+o3.l+o4.l)%2==Op.parity)  and (  std::abs(o1.tz2+o2.tz2-o3.tz2-o4.tz2)==2*Op.rank_T)  );
}

/*
void comm223ss_debug( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
//  int e3maxcut = 999;
  int norbs = Z.modelspace->GetNumberOrbits();
  auto& Z3 = Z.ThreeBody;
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  if ( (std::abs( X2.Norm() * Y2.Norm() ) < 1e-6 ) and not Z.modelspace->scalar3b_transform_first_pass) return;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();
  // we loop over i, j<=i, k<=j    l<=i, m<=l, n<=m, and Jij, Jmn free unless  ijk == lmn, then we take Jmn <= Jij.

//  std::vector< std::map<std::array<size_t,2>,size_t>::iterator> channel_iterators;
  std::vector< std::pair<const std::array<size_t,2>,size_t>> channel_iterators;
  for (auto iter : Z.ThreeBody.ch_start ) channel_iterators.push_back(iter);
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  size_t n_iter = channel_iterators.size();
//  for (size_t ch3=0; ch3<nch3; ch3++ )
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t iiter=0; iiter<n_iter; iiter++ )
  {
    size_t ch_bra = channel_iterators[iiter].first[0];
    size_t ch_ket = channel_iterators[iiter].first[1];
//    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    auto& Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch_bra);
    auto& Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch_ket);
//    size_t nkets3 = Tbc.GetNumberKets();
    size_t nbras3 = Tbc_bra.GetNumberKets();
    size_t nkets3 = Tbc_ket.GetNumberKets();
//    int twoJ = Tbc.twoJ;
    int twoJ = Tbc_bra.twoJ;
    double Jtot = 0.5*twoJ;
    for (size_t ibra=0; ibra<nbras3; ibra++)
    {
//      Ket3& bra = Tbc.GetKet(ibra);
      Ket3& bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      int ei = 2*oi.n + oi.l;
      int ej = 2*oj.n + oj.l;
      int ek = 2*ok.n + ok.l;
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (std::abs(ei-e_fermi[oi.tz2]) + std::abs(ej-e_fermi[oj.tz2]) + std::abs(ek-e_fermi[ok.tz2])) > Z.modelspace->GetdE3max() ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      double ji = 0.5*oi.j2;
      double jj = 0.5*oj.j2;
      double jk = 0.5*ok.j2;
      int J1 = bra.Jpq;


          int Jx_min,Jx_max,Jy_min,Jy_max;

          std::vector<int> J_Pik;
          std::vector<double> recouple_Pik;
          Jx_min = std::max( std::abs(oi.j2-twoJ), std::abs(oj.j2-ok.j2) )/2;
          Jx_max = std::min( oi.j2+twoJ, oj.j2+ok.j2)/2;
          for (int Jx=Jx_min; Jx<=Jx_max; Jx++)
          {
             double recouple_coeff = (2*Jx+1) * Z.modelspace->GetSixJ(ji,jj,J1, jk,Jtot,Jx);
             if (std::abs( recouple_coeff ) <1e-7) continue;
             J_Pik.push_back( Jx );
             recouple_Pik.push_back(recouple_coeff);
          }


          std::vector<int> J_Pjk;
          std::vector<double> recouple_Pjk;
          Jx_min = std::max( std::abs(oj.j2-twoJ), std::abs(oi.j2-ok.j2) )/2;
          Jx_max = std::min( oj.j2+twoJ, oi.j2+ok.j2)/2;
          for (int Jx=Jx_min; Jx<=Jx_max; Jx++)
          {
             double recouple_coeff = (2*Jx+1) * Z.modelspace->GetSixJ(jj,ji,J1, jk,Jtot,Jx);
             if (std::abs( recouple_coeff ) <1e-7) continue;
             J_Pjk.push_back( Jx );
             recouple_Pjk.push_back(recouple_coeff);
          }


      size_t ket_min = (ch_bra==ch_ket) ? ibra : 0;
//      for (size_t iket=ibra; iket<nkets3; iket++)

      for (size_t iket=ket_min; iket<nkets3; iket++)
      {
//        Ket3& ket = Tbc.GetKet(iket);
        Ket3& ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        int el = 2*ol.n + ol.l;
        int em = 2*om.n + om.l;
        int en = 2*on.n + on.l;
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (std::abs(el-e_fermi[ol.tz2]) + std::abs(em-e_fermi[om.tz2]) + std::abs(en-e_fermi[on.tz2])) > Z.modelspace->GetdE3max() ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        double jl = 0.5*ol.j2;
        double jm = 0.5*om.j2;
        double jn = 0.5*on.j2;
        int J2 = ket.Jpq;

        double zijklmn = 0;




          std::vector<int> J_Pln;
          std::vector<double> recouple_Pln;
          Jx_min = std::max( std::abs(ol.j2-twoJ), std::abs(om.j2-on.j2) )/2;
          Jx_max = std::min( ol.j2+twoJ, om.j2+on.j2)/2;
          for (int Jx=Jx_min; Jx<=Jx_max; Jx++)
          {
             double recouple_coeff = (2*Jx+1) * Z.modelspace->GetSixJ(jl,jm,J2, jn,Jtot,Jx);
             if (std::abs( recouple_coeff ) <1e-7) continue;
             J_Pln.push_back( Jx );
             recouple_Pln.push_back(recouple_coeff);
          }


          std::vector<int> J_Pmn;
          std::vector<double> recouple_Pmn;
          Jx_min = std::max( std::abs(om.j2-twoJ), std::abs(ol.j2-on.j2) )/2;
          Jx_max = std::min( om.j2+twoJ, ol.j2+on.j2)/2;
          for (int Jx=Jx_min; Jx<=Jx_max; Jx++)
          {
             double recouple_coeff = (2*Jx+1) * Z.modelspace->GetSixJ(jm,jl,J2, jn,Jtot,Jx);
             if (std::abs( recouple_coeff ) <1e-7) continue;
             J_Pmn.push_back( Jx );
             recouple_Pmn.push_back(recouple_coeff);
          }


//         for (int a=0; a<norbs; a++) // TODO: this can maybe be made more efficient by first looping through ja and computing the 6js?
         for ( auto OBC_a : Z.OneBodyChannels ) 
         {
//          if ( OBC_a.second.size() < 1) continue;
          Orbit& oa1 = Z.modelspace->GetOrbit(  *(OBC_a.second.begin()) );
//          Orbit& oa = Z.modelspace->GetOrbit(a);

//          double ja = 0.5*oa.j2;
          double ja = 0.5*oa1.j2;


          // TODO: recast these 9 terms into a loop as done in the other []->3 commutators

          //Z5   remapped to ZZ1 
          if (   (check_2b_channel_Tz_parity(X,oi,oj,on,oa1) and check_2b_channel_Tz_parity(Y,ok,oa1,om,ol))
              or (check_2b_channel_Tz_parity(Y,oi,oj,on,oa1) and check_2b_channel_Tz_parity(X,ok,oa1,om,ol)) )
          {
            double sixj1 = Z.modelspace->GetSixJ(jn,ja,J1, jk,Jtot,J2);
            if (std::abs(sixj1)>1e-7)
            {
              for ( auto a : OBC_a.second )
              {
                double xijna = X2.GetTBME_J(J1,i,j,n,a);
                double yijna = Y2.GetTBME_J(J1,i,j,n,a);
                double xkalm = X2.GetTBME_J(J2,k,a,l,m);
                double ykalm = Y2.GetTBME_J(J2,k,a,l,m);
                zijklmn +=  sixj1 * (xijna*ykalm-yijna*xkalm);
              }
            }
          }



          //Z7  remapped to Pik ZZ1
          if (  ( (check_2b_channel_Tz_parity(X,ok,oj,on,oa1) and check_2b_channel_Tz_parity(Y,oa1,oi,om,ol))
              or (check_2b_channel_Tz_parity(Y,ok,oj,on,oa1) and check_2b_channel_Tz_parity(X,oa1,oi,om,ol)) )
              and ( (std::abs(oa1.j2-oi.j2)<=2*J2 ) and (oa1.j2+oi.j2>=2*J2)  )   )
          {
            Jx_min = std::max(std::abs(ok.j2-oj.j2), std::abs(oa1.j2-on.j2) )/2;
            Jx_max = std::min((ok.j2+oj.j2), (oa1.j2+on.j2) )/2;
            double phase = Z.modelspace->phase( (ol.j2+om.j2+oi.j2+oa1.j2)/2);
            for (auto a : OBC_a.second )
            {
             double xaiml = X2.GetTBME_J(J2,a,i,m,l);
             double yaiml = Y2.GetTBME_J(J2,a,i,m,l);
             for ( size_t i_Pik=0; i_Pik<J_Pik.size(); i_Pik++)
             {
               int Jx = J_Pik[i_Pik];
               double recouple_ik = recouple_Pik[i_Pik];
               double sixj2 = Z.modelspace->GetSixJ(jn,ja,Jx, ji,Jtot,J2);
               if (std::abs(sixj2)<1e-7) continue;
               double xkjna = X2.GetTBME_J(Jx,k,j,n,a);
               double ykjna = Y2.GetTBME_J(Jx,k,j,n,a);
               zijklmn += recouple_ik * phase * sixj2 * (xkjna*yaiml-ykjna*xaiml);
              }
            }
          }


          //Z9  remapped to Pjk ZZ1
          if ( (  (check_2b_channel_Tz_parity(X,oi,ok,on,oa1) and check_2b_channel_Tz_parity(Y,oa1,oj,om,ol))
               or (check_2b_channel_Tz_parity(Y,oi,ok,on,oa1) and check_2b_channel_Tz_parity(X,oa1,oj,om,ol)) ) 
              and ( (std::abs(oa1.j2-oj.j2)<=2*J2 ) and (oa1.j2+oj.j2>=2*J2)  )  )
          {
            Jx_min = std::max(std::abs(ok.j2-oi.j2), std::abs(oa1.j2-on.j2) )/2;
            Jx_max = std::min((ok.j2+oi.j2), (oa1.j2+on.j2) )/2;
            for (auto a : OBC_a.second )
            {
             double xajml = X2.GetTBME_J(J2,a,j,m,l);
             double yajml = Y2.GetTBME_J(J2,a,j,m,l);
             for ( size_t i_Pjk=0; i_Pjk<J_Pjk.size(); i_Pjk++)
             {
              int Jx = J_Pjk[i_Pjk];
              double recouple_jk = recouple_Pjk[i_Pjk];
              double phase = Z.modelspace->phase( (oa1.j2+ok.j2+ol.j2+om.j2)/2 -J1-Jx);
              double sixj2 = Z.modelspace->GetSixJ(jj,ja,J2, jn,Jtot,Jx);
              if (std::abs(sixj2)<1e-7) continue;
              double xikna = X2.GetTBME_J(Jx,i,k,n,a);
              double yikna = Y2.GetTBME_J(Jx,i,k,n,a);
              zijklmn +=  phase * recouple_jk * sixj2 * (xikna*yajml-yikna*xajml);
             }
            }
          }


          //Z4 remapped to Pmn ZZ1
          if (  ( (check_2b_channel_Tz_parity(X,oi,oj,om,oa1) and check_2b_channel_Tz_parity(Y,oa1,ok,ol,on))
               or (check_2b_channel_Tz_parity(Y,oi,oj,om,oa1) and check_2b_channel_Tz_parity(X,oa1,ok,ol,on)) )
              and ( (std::abs(om.j2-oa1.j2)<=2*J1 ) and (om.j2+oa1.j2>=2*J1)  ) )
          {
            Jx_min = std::max(std::abs(oa1.j2-ok.j2), std::abs(ol.j2-on.j2) )/2;
            Jx_max = std::min((oa1.j2+ok.j2), (ol.j2+on.j2) )/2;
            double phase = -Z.modelspace->phase( (om.j2-on.j2-oa1.j2-ok.j2)/2-J2);
            for (auto a : OBC_a.second )
            {
              double xijma = X2.GetTBME_J(J1,i,j,m,a);
              double yijma = Y2.GetTBME_J(J1,i,j,m,a);
              for ( size_t i_Pmn=0; i_Pmn<J_Pmn.size(); i_Pmn++)
              {
                int Jx = J_Pmn[i_Pmn];
                double recouple_mn = recouple_Pmn[i_Pmn];
                double sixj2 = Z.modelspace->GetSixJ(jm,ja,J1, jk,Jtot,Jx);

                if (std::abs(sixj2)<1e-7) continue;
                double xakln = X2.GetTBME_J(Jx,a,k,l,n);
                double yakln = Y2.GetTBME_J(Jx,a,k,l,n);
                zijklmn +=  phase * recouple_mn * sixj2 * (xijma*yakln-yijma*xakln);
              }
            }
          }



          //Z1   remapped to Pln ZZ1
          if (  (   (check_2b_channel_Tz_parity(X,oi,oj,ol,oa1) and check_2b_channel_Tz_parity(Y,oa1,ok,om,on))
                 or (check_2b_channel_Tz_parity(Y,oi,oj,ol,oa1) and check_2b_channel_Tz_parity(X,oa1,ok,om,on)) )
             and (std::abs(ol.j2-oa1.j2)<=2*J1) and ((ol.j2+oa1.j2)>=2*J1) )
          {
            Jx_min = std::max(std::abs(oa1.j2-ok.j2), std::abs(om.j2-on.j2) )/2;
            Jx_max = std::min((oa1.j2+ok.j2), (om.j2+on.j2) )/2;
            double phase = -Z.modelspace->phase( (om.j2-on.j2-oa1.j2-ok.j2)/2);
            for (auto a : OBC_a.second )
            {
              double xijla = X2.GetTBME_J(J1,i,j,l,a);
              double yijla = Y2.GetTBME_J(J1,i,j,l,a);
              for ( size_t i_Pln=0; i_Pln<J_Pln.size(); i_Pln++)
              {
                int Jx = J_Pln[i_Pln];
                double recouple_ln = recouple_Pln[i_Pln];
                double sixj2 = Z.modelspace->GetSixJ(jl,ja,J1, jk,Jtot,Jx);
                if (std::abs(sixj2)<1e-7) continue;
                double xakmn = X2.GetTBME_J(Jx,a,k,m,n);
                double yakmn = Y2.GetTBME_J(Jx,a,k,m,n);
                zijklmn +=  phase * recouple_ln * sixj2 * (xijla*yakmn-yijla*xakmn);
              }
            }
          }


          //Z2  remapped to   Pik Pln  ZZ1
          if (   (check_2b_channel_Tz_parity(X,ok,oj,ol,oa1) and check_2b_channel_Tz_parity(Y,oa1,oi,om,on))
              or (check_2b_channel_Tz_parity(Y,ok,oj,ol,oa1) and check_2b_channel_Tz_parity(X,oa1,oi,om,on)) )
          {
            Jx_min = std::max(std::abs(oa1.j2-ol.j2), std::abs(ok.j2-oj.j2) )/2;
            Jx_max = std::min((oa1.j2+ol.j2), (ok.j2+oj.j2) )/2;
            Jy_min = std::max(std::abs(oa1.j2-oi.j2), std::abs(om.j2-on.j2) )/2;
            Jy_max = std::min((oa1.j2+oi.j2), (om.j2+on.j2) )/2;
            double phase = Z.modelspace->phase( (om.j2+on.j2-oa1.j2-oi.j2)/2);
            for ( size_t i_Pik=0; i_Pik<J_Pik.size(); i_Pik++)
            {
              int Jx = J_Pik[i_Pik];
              if ( Jx<Jx_min or Jx>Jx_max) continue;
              double recouple_ik = recouple_Pik[i_Pik];
              for (auto a : OBC_a.second )
              {
                double xkjla = X2.GetTBME_J(Jx,k,j,l,a);
                double ykjla = Y2.GetTBME_J(Jx,k,j,l,a);
                for ( size_t i_Pln=0; i_Pln<J_Pln.size(); i_Pln++)
                {
                 int Jy = J_Pln[i_Pln];
                 double recouple_ln = recouple_Pln[i_Pln];
                 double sixj3 = Z.modelspace->GetSixJ(ji,ja,Jy, jl,Jtot,Jx);
                 if (std::abs(sixj3)<1e-7) continue;
                 double xaimn = X2.GetTBME_J(Jy,a,i,m,n);
                 double yaimn = Y2.GetTBME_J(Jy,a,i,m,n);
                 zijklmn +=  phase * recouple_ik  * recouple_ln *sixj3 * (xkjla*yaimn -ykjla*xaimn);
                }
              }
            }
          }


          //Z3  remapped to  Pjk Pln ZZ1
          if (   (check_2b_channel_Tz_parity(X,oi,ok,ol,oa1) and check_2b_channel_Tz_parity(Y,oa1,oj,om,on))
              or (check_2b_channel_Tz_parity(Y,oi,ok,ol,oa1) and check_2b_channel_Tz_parity(X,oa1,oj,om,on)) )
          {
            Jx_min = std::max(std::abs(oa1.j2-ol.j2), std::abs(oi.j2-ok.j2) )/2;
            Jx_max = std::min((oa1.j2+ol.j2), (oi.j2+ok.j2) )/2;
            Jy_min = std::max(std::abs(oa1.j2-oj.j2), std::abs(om.j2-on.j2) )/2;
            Jy_max = std::min((oa1.j2+oj.j2), (om.j2+on.j2) )/2;
            for ( size_t i_Pjk=0; i_Pjk<J_Pjk.size(); i_Pjk++)
            {
              int Jx = J_Pjk[i_Pjk];
              if ( Jx<Jx_min or Jx>Jx_max) continue;
              double recouple_jk = recouple_Pjk[i_Pjk];
              double phase =-Z.modelspace->phase( (ok.j2+om.j2+on.j2-oa1.j2)/2-J1-Jx);
               for (auto a : OBC_a.second )
               {
                 double xikla = X2.GetTBME_J(Jx,i,k,l,a);
                 double yikla = Y2.GetTBME_J(Jx,i,k,l,a);
                 for ( size_t i_Pln=0; i_Pln<J_Pln.size(); i_Pln++)
                 {
                  int Jy = J_Pln[i_Pln];
                  double recouple_ln = recouple_Pln[i_Pln];
                  double sixj3 = Z.modelspace->GetSixJ(jj,ja,Jy, jl,Jtot,Jx);
                  if (std::abs(sixj3)<1e-7) continue;
                  double xajmn = X2.GetTBME_J(Jy,a,j,m,n);
                  double yajmn = Y2.GetTBME_J(Jy,a,j,m,n);
                  zijklmn += phase * recouple_jk * recouple_ln *sixj3 * (xikla*yajmn-yikla*xajmn);
                 }
              }
            }
          }


          //Z6  remapped to Pik Pmn ZZ1
          if (   (check_2b_channel_Tz_parity(X,ok,oj,om,oa1) and check_2b_channel_Tz_parity(Y,oa1,oi,ol,on))
              or (check_2b_channel_Tz_parity(Y,ok,oj,om,oa1) and check_2b_channel_Tz_parity(X,oa1,oi,ol,on)) )
          {
            Jx_min = std::max(std::abs(oa1.j2-om.j2), std::abs(oj.j2-ok.j2) )/2;
            Jx_max = std::min((oa1.j2+om.j2), (oj.j2+ok.j2) )/2;
            Jy_min = std::max(std::abs(oa1.j2-oi.j2), std::abs(ol.j2-on.j2) )/2;
            Jy_max = std::min((oa1.j2+oi.j2), (ol.j2+on.j2) )/2;
            double phase =-Z.modelspace->phase( (on.j2-om.j2-oi.j2-oa1.j2)/2-J2);
            for ( size_t i_Pik=0; i_Pik<J_Pik.size(); i_Pik++)
            {
              int Jx = J_Pik[i_Pik];
              if ( Jx<Jx_min or Jx>Jx_max) continue;
              double recouple_ik = recouple_Pik[i_Pik];
              for (auto a : OBC_a.second )
              {
                double xkjma = X2.GetTBME_J(Jx,k,j,m,a);
                double ykjma = Y2.GetTBME_J(Jx,k,j,m,a);
                for ( size_t i_Pmn=0; i_Pmn<J_Pmn.size(); i_Pmn++)
                {
                  int Jy = J_Pmn[i_Pmn];
                  double recouple_mn = recouple_Pmn[i_Pmn];
                  double sixj3 = Z.modelspace->GetSixJ(ji,ja,Jy, jm,Jtot,Jx);
                  if (std::abs(sixj3)<1e-7) continue;
                  double xailn = X2.GetTBME_J(Jy,a,i,l,n);
                  double yailn = Y2.GetTBME_J(Jy,a,i,l,n);
                  zijklmn +=  phase * recouple_ik * recouple_mn * sixj3 * (xkjma*yailn-ykjma*xailn);
                }
             }
            }
          }


          //Z8  remapped to Pjk Pmn ZZ1
          if (   (check_2b_channel_Tz_parity(X,oi,ok,om,oa1) and check_2b_channel_Tz_parity(Y,oa1,oj,ol,on))
              or (check_2b_channel_Tz_parity(Y,oi,ok,om,oa1) and check_2b_channel_Tz_parity(X,oa1,oj,ol,on)) )
          {
            Jx_min = std::max(std::abs(oa1.j2-om.j2), std::abs(oi.j2-ok.j2) )/2;
            Jx_max = std::min((oa1.j2+om.j2), (oi.j2+ok.j2) )/2;
            Jy_min = std::max(std::abs(oa1.j2-oj.j2), std::abs(ol.j2-on.j2) )/2;
            Jy_max = std::min((oa1.j2+oj.j2), (ol.j2+on.j2) )/2;
            for ( size_t i_Pjk=0; i_Pjk<J_Pjk.size(); i_Pjk++)
            {
              int Jx = J_Pjk[i_Pjk];
              if ( Jx<Jx_min or Jx>Jx_max) continue;
              double recouple_jk = recouple_Pjk[i_Pjk];
              double phase =-Z.modelspace->phase( (ok.j2+om.j2+on.j2-oa1.j2)/2-J1-J2-Jx);
              for (auto a : OBC_a.second )
              {
                double xikma = X2.GetTBME_J(Jx,i,k,m,a);
                double yikma = Y2.GetTBME_J(Jx,i,k,m,a);
                for ( size_t i_Pmn=0; i_Pmn<J_Pmn.size(); i_Pmn++)
                {
                  int Jy = J_Pmn[i_Pmn];
                  double recouple_mn = recouple_Pmn[i_Pmn];
                  double sixj3 = Z.modelspace->GetSixJ(jj,ja,Jy, jm,Jtot,Jx);
                  double xajln = X2.GetTBME_J(Jy,a,j,l,n);
                  double yajln = Y2.GetTBME_J(Jy,a,j,l,n);
                  zijklmn +=  phase * recouple_jk * recouple_mn * sixj3 * (xikma*yajln-yikma*xajln);
                }
              }
            }
          }



         }// for OBC_a

         zijklmn *=  sqrt( (2*J1+1.)*(2*J2+1.));
         Z3.AddToME_pn_ch( ch_bra,ch_ket,ibra,iket, zijklmn );

    }// for iket
   }// for ibra
  }// for ch3
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
*/




/*
void comm223ss( const Operator& X, const Operator& Y, Operator& Z )
{
  double tstart = omp_get_wtime();
//  int e3maxcut = 999;
  int norbs = Z.modelspace->GetNumberOrbits();
  auto& Z3 = Z.ThreeBody;
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  if ( (std::abs( X2.Norm() * Y2.Norm() ) < 1e-6 ) and not Z.modelspace->scalar3b_transform_first_pass) return;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();
  // we loop over i, j<=i, k<=j    l<=i, m<=l, n<=m, and Jij, Jmn free unless  ijk == lmn, then we take Jmn <= Jij.

//  std::vector< std::map<std::array<size_t,2>,size_t>::iterator> channel_iterators;
  std::vector< std::pair<const std::array<size_t,2>,size_t>> channel_iterators;
  for (auto iter : Z.ThreeBody.ch_start ) channel_iterators.push_back(iter);
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  size_t n_iter = channel_iterators.size();
//  for (size_t ch3=0; ch3<nch3; ch3++ )
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t iiter=0; iiter<n_iter; iiter++ )
  {
    size_t ch_bra = channel_iterators[iiter].first[0];
    size_t ch_ket = channel_iterators[iiter].first[1];
//    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    auto& Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch_bra);
    auto& Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch_ket);
//    size_t nkets3 = Tbc.GetNumberKets();
    size_t nbras3 = Tbc_bra.GetNumberKets();
    size_t nkets3 = Tbc_ket.GetNumberKets();
//    int twoJ = Tbc.twoJ;
    int twoJ = Tbc_bra.twoJ;
    double Jtot = 0.5*twoJ;
    for (size_t ibra=0; ibra<nbras3; ibra++)
    {
//      Ket3& bra = Tbc.GetKet(ibra);
      Ket3& bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      int ei = 2*oi.n + oi.l;
      int ej = 2*oj.n + oj.l;
      int ek = 2*ok.n + ok.l;
      if ( (std::abs(ei-e_fermi[oi.tz2]) + std::abs(ej-e_fermi[oj.tz2]) + std::abs(ek-e_fermi[ok.tz2])) > Z.modelspace->GetdE3max() ) continue;
      double ji = 0.5*oi.j2;
      double jj = 0.5*oj.j2;
      double jk = 0.5*ok.j2;
      int J1 = bra.Jpq;
      size_t ket_min = (ch_bra==ch_ket) ? ibra : 0;
//      for (size_t iket=ibra; iket<nkets3; iket++)
      for (size_t iket=ket_min; iket<nkets3; iket++)
      {
//        Ket3& ket = Tbc.GetKet(iket);
        Ket3& ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        int el = 2*ol.n + ol.l;
        int em = 2*om.n + om.l;
        int en = 2*on.n + on.l;
        if ( (std::abs(el-e_fermi[ol.tz2]) + std::abs(em-e_fermi[om.tz2]) + std::abs(en-e_fermi[on.tz2])) > Z.modelspace->GetdE3max() ) continue;
        double jl = 0.5*ol.j2;
        double jm = 0.5*om.j2;
        double jn = 0.5*on.j2;
        int J2 = ket.Jpq;

        double zijklmn = 0;

         for (int a=0; a<norbs; a++) // TODO: this can maybe be made more efficient by first looping through ja and computing the 6js?
         {
          Orbit& oa = Z.modelspace->GetOrbit(a);

          double ja = 0.5*oa.j2;

          int Jx_min,Jx_max,Jy_min,Jy_max;


          //Z1
          if (  (   (check_2b_channel_Tz_parity(X,oi,oj,ol,oa) and check_2b_channel_Tz_parity(Y,oa,ok,om,on))
                 or (check_2b_channel_Tz_parity(Y,oi,oj,ol,oa) and check_2b_channel_Tz_parity(X,oa,ok,om,on)) )
             and (std::abs(ol.j2-oa.j2)<=2*J1) and ((ol.j2+oa.j2)>=2*J1) )
          {
            Jx_min = std::max(std::abs(oa.j2-ok.j2), std::abs(om.j2-on.j2) )/2;
            Jx_max = std::min((oa.j2+ok.j2), (om.j2+on.j2) )/2;
            double phase = -Z.modelspace->phase( (om.j2-on.j2-oa.j2-ok.j2)/2);
            double xijla = X2.GetTBME_J(J1,i,j,l,a);
            double yijla = Y2.GetTBME_J(J1,i,j,l,a);
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++) 
            {
              double hats  =  (2*Jx+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
              double sixj1 = Z.modelspace->GetSixJ(jl,jm,J2, jn,Jtot,Jx);
              double sixj2 = Z.modelspace->GetSixJ(jl,ja,J1, jk,Jtot,Jx);
              double xakmn = X2.GetTBME_J(Jx,a,k,m,n);
              double yakmn = Y2.GetTBME_J(Jx,a,k,m,n);
              zijklmn += hats * phase * sixj1*sixj2 * (xijla*yakmn-yijla*xakmn);
            }
          }


          //Z2
          if (   (check_2b_channel_Tz_parity(X,ok,oj,ol,oa) and check_2b_channel_Tz_parity(Y,oa,oi,om,on))
              or (check_2b_channel_Tz_parity(Y,ok,oj,ol,oa) and check_2b_channel_Tz_parity(X,oa,oi,om,on)) )
          {
            Jx_min = std::max(std::abs(oa.j2-ol.j2), std::abs(ok.j2-oj.j2) )/2;
            Jx_max = std::min((oa.j2+ol.j2), (ok.j2+oj.j2) )/2;
            Jy_min = std::max(std::abs(oa.j2-oi.j2), std::abs(om.j2-on.j2) )/2;
            Jy_max = std::min((oa.j2+oi.j2), (om.j2+on.j2) )/2;
            double phase = Z.modelspace->phase( (om.j2+on.j2-oa.j2-oi.j2)/2);
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
            {
              double sixj1 = Z.modelspace->GetSixJ(ji,jj,J1, jk,Jtot,Jx);
              double xkjla = X2.GetTBME_J(Jx,k,j,l,a);
              double ykjla = Y2.GetTBME_J(Jx,k,j,l,a);
              for ( int Jy=Jy_min; Jy<=Jy_max; Jy++)
              {
               double hats  =  (2*Jx+1) * (2*Jy+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
               double sixj2 = Z.modelspace->GetSixJ(jl,jm,J2, jn,Jtot,Jy);
               double sixj3 = Z.modelspace->GetSixJ(ji,ja,Jy, jl,Jtot,Jx);
               double xaimn = X2.GetTBME_J(Jy,a,i,m,n);
               double yaimn = Y2.GetTBME_J(Jy,a,i,m,n);
               zijklmn += hats * phase * sixj1*sixj2*sixj3 * (xkjla*yaimn -ykjla*xaimn);

              }
            }
          }


          //Z3
          if (   (check_2b_channel_Tz_parity(X,oi,ok,ol,oa) and check_2b_channel_Tz_parity(Y,oa,oj,om,on))
              or (check_2b_channel_Tz_parity(Y,oi,ok,ol,oa) and check_2b_channel_Tz_parity(X,oa,oj,om,on)) )
          {
            Jx_min = std::max(std::abs(oa.j2-ol.j2), std::abs(oi.j2-ok.j2) )/2;
            Jx_max = std::min((oa.j2+ol.j2), (oi.j2+ok.j2) )/2;
            Jy_min = std::max(std::abs(oa.j2-oj.j2), std::abs(om.j2-on.j2) )/2;
            Jy_max = std::min((oa.j2+oj.j2), (om.j2+on.j2) )/2;
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
            {
              double phase =-Z.modelspace->phase( (ok.j2+om.j2+on.j2-oa.j2)/2-J1-Jx);
              double sixj1 = Z.modelspace->GetSixJ(jj,ji,J1, jk,Jtot,Jx);
              double xikla = X2.GetTBME_J(Jx,i,k,l,a);
              double yikla = Y2.GetTBME_J(Jx,i,k,l,a);
              for ( int Jy=Jy_min; Jy<=Jy_max; Jy++)
              {
               double hats  = (2*Jx+1) * (2*Jy+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
               double sixj2 = Z.modelspace->GetSixJ(jl,jm,J2, jn,Jtot,Jy);
               double sixj3 = Z.modelspace->GetSixJ(jj,ja,Jy, jl,Jtot,Jx);
               double xajmn = X2.GetTBME_J(Jy,a,j,m,n);
               double yajmn = Y2.GetTBME_J(Jy,a,j,m,n);
               zijklmn += hats * phase * sixj1*sixj2*sixj3 * (xikla*yajmn-yikla*xajmn);
              }
            }
          }

          //Z4
          if (  ( (check_2b_channel_Tz_parity(X,oi,oj,om,oa) and check_2b_channel_Tz_parity(Y,oa,ok,ol,on))
               or (check_2b_channel_Tz_parity(Y,oi,oj,om,oa) and check_2b_channel_Tz_parity(X,oa,ok,ol,on)) )
              and ( (std::abs(om.j2-oa.j2)<=2*J1 ) and (om.j2+oa.j2>=2*J1)  ) )
          {
            Jx_min = std::max(std::abs(oa.j2-ok.j2), std::abs(ol.j2-on.j2) )/2;
            Jx_max = std::min((oa.j2+ok.j2), (ol.j2+on.j2) )/2;
            double phase = -Z.modelspace->phase( (om.j2-on.j2-oa.j2-ok.j2)/2-J2);
            double xijma = X2.GetTBME_J(J1,i,j,m,a);
            double yijma = Y2.GetTBME_J(J1,i,j,m,a);
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++) 
            {
              double hats  =  (2*Jx+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
              double sixj1 = Z.modelspace->GetSixJ(jm,jl,J2, jn,Jtot,Jx);
              double sixj2 = Z.modelspace->GetSixJ(jm,ja,J1, jk,Jtot,Jx);
              double xakln = X2.GetTBME_J(Jx,a,k,l,n);
              double yakln = Y2.GetTBME_J(Jx,a,k,l,n);
              zijklmn += hats * phase * sixj1*sixj2 * (xijma*yakln-yijma*xakln);
            }
          }

          //Z5
          if (   (check_2b_channel_Tz_parity(X,oi,oj,on,oa) and check_2b_channel_Tz_parity(Y,ok,oa,om,ol))
              or (check_2b_channel_Tz_parity(Y,oi,oj,on,oa) and check_2b_channel_Tz_parity(X,ok,oa,om,ol)) )
          {
            double hats  =  sqrt( (2*J1+1.)*(2*J2+1.) );
            double phase =  1;
            double sixj1 = Z.modelspace->GetSixJ(jn,ja,J1, jk,Jtot,J2);
            double xijna = X2.GetTBME_J(J1,i,j,n,a);
            double yijna = Y2.GetTBME_J(J1,i,j,n,a);
            double xkalm = X2.GetTBME_J(J2,k,a,l,m);
            double ykalm = Y2.GetTBME_J(J2,k,a,l,m);
            zijklmn += hats * phase * sixj1 * (xijna*ykalm-yijna*xkalm);
          }

          //Z6
          if (   (check_2b_channel_Tz_parity(X,ok,oj,om,oa) and check_2b_channel_Tz_parity(Y,oa,oi,ol,on))
              or (check_2b_channel_Tz_parity(Y,ok,oj,om,oa) and check_2b_channel_Tz_parity(X,oa,oi,ol,on)) )
          {
            Jx_min = std::max(std::abs(oa.j2-om.j2), std::abs(oj.j2-ok.j2) )/2;
            Jx_max = std::min((oa.j2+om.j2), (oj.j2+ok.j2) )/2;
            Jy_min = std::max(std::abs(oa.j2-oi.j2), std::abs(ol.j2-on.j2) )/2;
            Jy_max = std::min((oa.j2+oi.j2), (ol.j2+on.j2) )/2;
            double phase =-Z.modelspace->phase( (on.j2-om.j2-oi.j2-oa.j2)/2-J2);
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
            {
              double sixj1 = Z.modelspace->GetSixJ(ji,jj,J1, jk,Jtot,Jx);
              double xkjma = X2.GetTBME_J(Jx,k,j,m,a);
              double ykjma = Y2.GetTBME_J(Jx,k,j,m,a);
             for ( int Jy=Jy_min; Jy<=Jy_max; Jy++)
             {
              double hats  =  (2*Jx+1) * (2*Jy+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
              double sixj2 = Z.modelspace->GetSixJ(jm,jl,J2, jn,Jtot,Jy);
              double sixj3 = Z.modelspace->GetSixJ(ji,ja,Jy, jm,Jtot,Jx);
              double xailn = X2.GetTBME_J(Jy,a,i,l,n);
              double yailn = Y2.GetTBME_J(Jy,a,i,l,n);
              zijklmn += hats * phase * sixj1*sixj2*sixj3 * (xkjma*yailn-ykjma*xailn);
             }
            }
          }


          //Z7
          if (  ( (check_2b_channel_Tz_parity(X,ok,oj,on,oa) and check_2b_channel_Tz_parity(Y,oa,oi,om,ol))
              or (check_2b_channel_Tz_parity(Y,ok,oj,on,oa) and check_2b_channel_Tz_parity(X,oa,oi,om,ol)) )
              and ( (std::abs(oa.j2-oi.j2)<=2*J2 ) and (oa.j2+oi.j2>=2*J2)  )   )
          {
            Jx_min = std::max(std::abs(ok.j2-oj.j2), std::abs(oa.j2-on.j2) )/2;
            Jx_max = std::min((ok.j2+oj.j2), (oa.j2+on.j2) )/2;
            double phase = Z.modelspace->phase( (ol.j2+om.j2+oi.j2+oa.j2)/2);
            double xaiml = X2.GetTBME_J(J2,a,i,m,l);
            double yaiml = Y2.GetTBME_J(J2,a,i,m,l);
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++) 
            {
              double hats  =  (2*Jx+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
              double sixj1 = Z.modelspace->GetSixJ(ji,jj,J1, jk,Jtot,Jx);
              double sixj2 = Z.modelspace->GetSixJ(ji,ja,J2, jn,Jtot,Jx);
              double xkjna = X2.GetTBME_J(Jx,k,j,n,a);
              double ykjna = Y2.GetTBME_J(Jx,k,j,n,a);
              zijklmn += hats * phase * sixj1*sixj2 * (xkjna*yaiml-ykjna*xaiml);
            }
          }



          //Z8
          if (   (check_2b_channel_Tz_parity(X,oi,ok,om,oa) and check_2b_channel_Tz_parity(Y,oa,oj,ol,on))
              or (check_2b_channel_Tz_parity(Y,oi,ok,om,oa) and check_2b_channel_Tz_parity(X,oa,oj,ol,on)) )
          {
            Jx_min = std::max(std::abs(oa.j2-om.j2), std::abs(oi.j2-ok.j2) )/2;
            Jx_max = std::min((oa.j2+om.j2), (oi.j2+ok.j2) )/2;
            Jy_min = std::max(std::abs(oa.j2-oj.j2), std::abs(ol.j2-on.j2) )/2;
            Jy_max = std::min((oa.j2+oj.j2), (ol.j2+on.j2) )/2;
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
            {
              double phase =-Z.modelspace->phase( (ok.j2+om.j2+on.j2-oa.j2)/2-J1-J2-Jx);
              double sixj1 = Z.modelspace->GetSixJ(jj,ji,J1, jk,Jtot,Jx);
              double xikma = X2.GetTBME_J(Jx,i,k,m,a);
              double yikma = Y2.GetTBME_J(Jx,i,k,m,a);
             for ( int Jy=Jy_min; Jy<=Jy_max; Jy++)
             {
              double hats  =  (2*Jx+1) * (2*Jy+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
              double sixj2 = Z.modelspace->GetSixJ(jm,jl,J2, jn,Jtot,Jy);
              double sixj3 = Z.modelspace->GetSixJ(jj,ja,Jy, jm,Jtot,Jx);
              double xajln = X2.GetTBME_J(Jy,a,j,l,n);
              double yajln = Y2.GetTBME_J(Jy,a,j,l,n);
              zijklmn += hats * phase * sixj1*sixj2*sixj3 * (xikma*yajln-yikma*xajln);
             }
            }
          }


          //Z9
          if ( (  (check_2b_channel_Tz_parity(X,oi,ok,on,oa) and check_2b_channel_Tz_parity(Y,oa,oj,om,ol))
               or (check_2b_channel_Tz_parity(Y,oi,ok,on,oa) and check_2b_channel_Tz_parity(X,oa,oj,om,ol)) ) 
              and ( (std::abs(oa.j2-oj.j2)<=2*J2 ) and (oa.j2+oj.j2>=2*J2)  )  )
          {
            Jx_min = std::max(std::abs(ok.j2-oi.j2), std::abs(oa.j2-on.j2) )/2;
            Jx_max = std::min((ok.j2+oi.j2), (oa.j2+on.j2) )/2;
            double xajml = X2.GetTBME_J(J2,a,j,m,l);
            double yajml = Y2.GetTBME_J(J2,a,j,m,l);
            for ( int Jx=Jx_min; Jx<=Jx_max; Jx++) 
            {
              double phase = Z.modelspace->phase( (oa.j2+ok.j2+ol.j2+om.j2)/2 -J1-Jx);
              double hats  =  (2*Jx+1) * sqrt( (2*J1+1.)*(2*J2+1.) );
              double sixj1 = Z.modelspace->GetSixJ(jj,ji,J1, jk,Jtot,Jx);
              double sixj2 = Z.modelspace->GetSixJ(jj,ja,J2, jn,Jtot,Jx);
              double xikna = X2.GetTBME_J(Jx,i,k,n,a);
              double yikna = Y2.GetTBME_J(Jx,i,k,n,a);
              zijklmn += hats * phase * sixj1*sixj2 * (xikna*yajml-yikna*xajml);
            }
          }

         }// for a

//             Z3.AddToME_pn_ch( ch3,ch3,ibra,iket, zijklmn );
             Z3.AddToME_pn_ch( ch_bra,ch_ket,ibra,iket, zijklmn );

    }// for iket
   }// for ibra
  }// for ch3
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
*/



//*****************************************************************************************
//
//  |     |    |     Uncoupled expression:
// i|    j|   k|       Z_ijklmn = 1/2 sum_{ab} (n`an`b - nanb) P(ij/k) X_{ijab} Y_{abklmn} - P(lm/n) Y_{ijkabn} X_{ablm}
//  *~~X~~*    |                   
// a|    b|    | 
//  *~~~~[Y]~~~*     Coupled expression:
// l|    m|   n|     Z_{ijklmn}^{J1,J2,J} = 1/2 sum_{ab} (n`an`b - nanb) ( P(ij/k)^{J1,J} X_{ijab}^{J1} Y_{abklmn}^{J1,J2,J} 
//  |     |    |                                                         - P(lm/n)^{J2,J} Y_{ijkabn}^{J1,J2,J} X_{ablm}^{J2} )
//  |     |    |                                                                                                                   
//                                                                                          
//                                                           
//      Checked with UnitTest and passed                                                                      
//
void comm233_pp_hhss( const Operator& X, const Operator& Y, Operator& Z )
{
//  std::cout << "ENTER " <<__func__ << std::endl;
  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;

  int hermX = X.IsHermitian() ? 1 : -1;
  int hermY = Y.IsHermitian() ? 1 : -1;

  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  int norbs = Z.modelspace->GetNumberOrbits();
  double X3NORM = X3.Norm();
  double Y3NORM = Y3.Norm();
//  bool x3_allocated = X3.is_allocated;
//  bool y3_allocated = Y3.is_allocated;
  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();
//  std::cout << "Begin the ch3 loop " << std::endl;

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    size_t nkets = Tbc.GetNumberKets();
    std::vector<size_t> kets_kept;
    std::map<size_t,size_t> kept_lookup;
    for (size_t iket=0; iket<nkets; iket++)
    {
      Ket3& ket = Tbc.GetKet(iket);
      double d_ei = std::abs(2*ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
      double d_ej = std::abs(2*ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
      double d_ek = std::abs(2*ket.oR->n + ket.oR->l - e_fermi[ket.oR->tz2]);
      double occnat_i = ket.op->occ_nat;
      double occnat_j = ket.oq->occ_nat;
      double occnat_k = ket.oR->occ_nat;
      if (  (d_ei+d_ej+d_ek ) > Z.modelspace->GetdE3max() ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      kets_kept.push_back( iket );
      kept_lookup[iket] = kets_kept.size()-1;
    }
    size_t nkets_kept = kets_kept.size();


    arma::mat X2MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat Y2MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat X3MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat Y3MAT( nkets_kept, nkets_kept, arma::fill::zeros);
    arma::mat Z3MAT( nkets_kept, nkets_kept, arma::fill::zeros);

//    std::cout << "   fill the 2b matrices ch = " << ch3 << std::endl;

    for (size_t index_bra=0; index_bra<nkets_kept; index_bra++)
    {
      size_t ibra = kets_kept[index_bra];
      Ket3& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double d_ei = std::abs(2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs(2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs(2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      int Jij = bra.Jpq;
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;


      for (size_t ch_ab=0; ch_ab<nch2; ch_ab++)
      {
        auto& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
        size_t nkets_ab = tbc_ab.GetNumberKets();
        int Jab = tbc_ab.J;

//        std::cout << "   direct " << std::endl;
        if ( ((2*tbc_ab.Tz + ok.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + ok.l + Tbc.parity)%2==0) and Jab == Jij)
        {
          for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
          {
              auto& ket_ab = tbc_ab.GetKet(iket_ab);
              size_t a = ket_ab.p;
              size_t b = ket_ab.q;
              double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              double occnat_a = ket_ab.op->occ_nat;
              double occnat_b = ket_ab.oq->occ_nat;
              if ( (d_ea + d_eb + d_ek) > Z.modelspace->dE3max ) continue;
              if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
              double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
              if (std::abs(occfactor)<1e-6) continue;
              double symm_factor =  (a==b) ? 1.0 : 2.0;
              double xijab = X2.GetTBME_J(Jab,Jab,i,j,a,b);
              double yijab = Y2.GetTBME_J(Jab,Jab,i,j,a,b);


              std::vector<size_t> ket_list;
              std::vector<double> recouple_list;
              
              size_t ch_check = Z3.GetKetIndex_withRecoupling( Jab, twoJ, a, b, k,  ket_list,  recouple_list );
              for (size_t ilist=0; ilist<ket_list.size(); ilist++)
              {
                auto iter_find = kept_lookup.find( ket_list[ilist] );
                if (iter_find == kept_lookup.end() ) continue;
                size_t index_ket = iter_find->second;
                double recouple = recouple_list[ilist];
                X2MAT(index_bra,index_ket) += 0.5 * symm_factor * occfactor * recouple * xijab ;
                Y2MAT(index_bra,index_ket) += 0.5 * symm_factor * occfactor * recouple * yijab ;
              }// for ilist

          }// for iket_ab
        } // if J, Tz and parity check for term 1

//        std::cout << "   Pik " << std::endl;
        // Permute  Pik
        if ( ((2*tbc_ab.Tz + oi.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + oi.l + Tbc.parity)%2==0)
                   and (std::abs(2*Jab-oi.j2)<=twoJ) and ((2*Jab+oi.j2)>=twoJ) 
                   and (std::abs(2*Jab-oj.j2)<=ok.j2) and ((2*Jab+oj.j2)>=ok.j2) )
        {
           double sixj = Z.modelspace->GetSixJ(ji,jj,Jij,jk,Jtot,Jab);
           if (std::abs(sixj)>1e-6) 
           {
            double hats = sqrt( (2*Jij+1)*(2*Jab+1));
            for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
            {
              auto& ket_ab = tbc_ab.GetKet(iket_ab);
              size_t a = ket_ab.p;
              size_t b = ket_ab.q;
              double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              double occnat_a = ket_ab.op->occ_nat;
              double occnat_b = ket_ab.oq->occ_nat;
              if ( (d_ea + d_eb + d_ei) > Z.modelspace->dE3max ) continue;
              if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) < Z.modelspace->GetOccNat3Cut() ) continue;
              double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
              if (std::abs(occfactor)<1e-6) continue;
              double symm_factor =  (a==b) ? 1.0 : 2.0;
              double xkjab = X2.GetTBME_J(Jab,Jab,k,j,a,b);
              double ykjab = Y2.GetTBME_J(Jab,Jab,k,j,a,b);

              std::vector<size_t> ket_list;
              std::vector<double> recouple_list;
              
//              std::cout << " line " << __LINE__ << "  calling GetKetIndex_withRecoupling " << Jab << " " << twoJ << " " << a << " " << b << " " << i << std::endl;
              size_t ch_check = Z3.GetKetIndex_withRecoupling( Jab, twoJ, a, b, i,  ket_list,  recouple_list );
//              std::cout << "    size of lists " << ket_list.size() << " " << recouple_list.size() << std::endl;
              for (size_t ilist=0; ilist<ket_list.size(); ilist++)
              {
                auto iter_find = kept_lookup.find( ket_list[ilist] );
                if (iter_find == kept_lookup.end() ) continue;
                size_t index_ket = iter_find->second;
                double recouple = recouple_list[ilist];
                X2MAT(index_bra,index_ket) += 0.5 * hats * sixj * symm_factor * occfactor * recouple * xkjab ;
                Y2MAT(index_bra,index_ket) += 0.5 * hats * sixj * symm_factor * occfactor * recouple * ykjab ;
              }// for ilist
//              std::cout << "  ok. line " << __LINE__ << std::endl;

            }// for iket_ab
           }// if sixj nonzero
        } // if J, Tz and parity check for term 2


//        std::cout << "   Pjk " << std::endl;
        // Permute  Pjk
        if ( ((2*tbc_ab.Tz + oj.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + oj.l + Tbc.parity)%2==0) 
                   and (std::abs(2*Jab-oj.j2)<=twoJ) and ((2*Jab+oj.j2)>=twoJ) 
                   and (std::abs(2*Jab-oi.j2)<=ok.j2) and ((2*Jab+oi.j2)>=ok.j2) )
        {
           double sixj = Z.modelspace->GetSixJ(jj,ji,Jij,jk,Jtot,Jab);
           if (std::abs(sixj)>1e-6)
           {
             double hats = sqrt( (2*Jij+1)*(2*Jab+1));
             int phase = Z.modelspace->phase( (oj.j2+ok.j2)/2-Jij-Jab);
             for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
             {
               auto& ket_ab = tbc_ab.GetKet(iket_ab);
               size_t a = ket_ab.p;
               size_t b = ket_ab.q;
               double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
               double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
               double occnat_a = ket_ab.op->occ_nat;
               double occnat_b = ket_ab.oq->occ_nat;
               if ( (d_ea + d_eb + d_ej) > Z.modelspace->dE3max ) continue;
               if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) < Z.modelspace->GetOccNat3Cut() ) continue;
               double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
               if (std::abs(occfactor)<1e-6) continue;
               double symm_factor =  (a==b) ? 1.0 : 2.0;
               double xikab = X2.GetTBME_J(Jab,Jab,i,k,a,b);
               double yikab = Y2.GetTBME_J(Jab,Jab,i,k,a,b);

              std::vector<size_t> ket_list;
              std::vector<double> recouple_list;
              
              size_t ch_check = Z3.GetKetIndex_withRecoupling( Jab, twoJ, a, b, j,  ket_list,  recouple_list );
              for (size_t ilist=0; ilist<ket_list.size(); ilist++)
              {
                auto iter_find = kept_lookup.find( ket_list[ilist] );
                if (iter_find == kept_lookup.end() ) continue;
                size_t index_ket = iter_find->second;
                double recouple = recouple_list[ilist];
                X2MAT(index_bra,index_ket) -= 0.5 * phase * hats * sixj * symm_factor * occfactor * recouple * xikab ;
                Y2MAT(index_bra,index_ket) -= 0.5 * phase * hats * sixj * symm_factor * occfactor * recouple * yikab ;

              }// for ilist

             }// for iket_ab
          }// if sixj nonzero
        } // if  J, Tz and parity check for term 3

      }// for ch_ab

    }// for index_bra

//     std::cout << "Fill the 3b matrices " << std::endl;

     // Fill X3 and Y3
    // kept_lookup is a map   Full index => Kept index, so iter_bra.first gives the full index, and iter_bra.second is the
    // index for the 3-body state we keep in this commutator
    for ( auto& iter_bra : kept_lookup )
    {
      for ( auto& iter_ket : kept_lookup )
      {
        if ( x3_allocated )
           X3MAT( iter_bra.second, iter_ket.second) = X3.GetME_pn_ch(ch3,ch3, iter_bra.first, iter_ket.first );

        if ( y3_allocated )
           Y3MAT( iter_bra.second, iter_ket.second) = Y3.GetME_pn_ch(ch3,ch3, iter_bra.first, iter_ket.first );
      }
    }


    // Do the matrix multiplication
//    Z3MAT = X2MAT*Y3MAT - Y2MAT*X3MAT  +  hermX*hermY * ( X3MAT.t()*Y2MAT.t() - Y3MAT.t()*X2MAT.t() );
    Z3MAT = X2MAT*Y3MAT - Y2MAT*X3MAT ;
    Z3MAT -= hermX * hermY * Z3MAT.t(); 


//    std::cout << "Unpack..." << std::endl;

    // unpack the result
    for ( auto& iter_bra : kept_lookup )
    {
      for ( auto& iter_ket : kept_lookup )
      {
        if ( iter_ket.first > iter_bra.first ) continue;
        Z3.AddToME_pn_ch(ch3,ch3, iter_bra.first,iter_ket.first,  Z3MAT(iter_bra.second,iter_ket.second) );
      }
    }
 
  }// for ch3
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}


// old much slower way
/*

void comm233_pp_hhss_debug( const Operator& X, const Operator& Y, Operator& Z )
{

  double tstart = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch=0; ch<=nch3; ch++)
  {
    auto Tbc = Z.modelspace->GetThreeBodyChannel(ch);
    size_t nket3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nket3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      int J1 = bra.Jpq;
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5*oi.j2;
      double jj = 0.5*oj.j2;
      double jk = 0.5*ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      for (size_t iket=0; iket<=ibra; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        int J2 = ket.Jpq;
        int l = ket.p;
        int m = ket.q;
        int n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5*ol.j2;
        double jm = 0.5*om.j2;
        double jn = 0.5*on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;

        double z_ijklmn = 0;

        for (size_t ch_ab=0; ch_ab<nch2; ch_ab++)
        {
          auto& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
          size_t nkets_ab = tbc_ab.GetNumberKets();
          int Jab = tbc_ab.J;

          if ( ((2*tbc_ab.Tz + ok.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + ok.l + Tbc.parity)%2==0) and Jab == J1)
          {
            for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
            {
              auto& ket_ab = tbc_ab.GetKet(iket_ab);
              size_t a = ket_ab.p;
              size_t b = ket_ab.q;
              double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              if ( (d_ea + d_eb + d_ek) > Z.modelspace->dE3max ) continue;
              double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
              if (std::abs(occfactor)<1e-6) continue;
              double symm_factor =  (a==b) ? 1.0 : 2.0;
              double xijab = X2.GetTBME_J(Jab,Jab,i,j,a,b);
              double yijab = Y2.GetTBME_J(Jab,Jab,i,j,a,b);
              double xabklmn = X3.GetME_pn(J1,J2,twoJ,a,b,k,l,m,n);
              double yabklmn = Y3.GetME_pn(J1,J2,twoJ,a,b,k,l,m,n);
              z_ijklmn += 0.5 * symm_factor * occfactor * (xijab * yabklmn - yijab * xabklmn);
//              if (ch==0 and ibra==4 and iket==4)
//              {
//                 std::cout << "  iket_ab = " << iket_ab << "   adding to z   0.5 * " << symm_factor << " * " << occfactor << " * " << xijab << " " << yabklmn << "  =>  " << z_ijklmn << std::endl;;
//              }
            }// for iket_ab
          } // if J, Tz and parity check for term 1


          // Permute  Pik
          if ( ((2*tbc_ab.Tz + oi.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + oi.l + Tbc.parity)%2==0)
                     and (std::abs(2*Jab-oi.j2)<=twoJ) and ((2*Jab+oi.j2)>=twoJ) )
          {
             double sixj = Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,Jab);
             if (std::abs(sixj)>1e-6) 
             {
              double hats = sqrt( (2*J1+1)*(2*Jab+1));
              for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
              {
                auto& ket_ab = tbc_ab.GetKet(iket_ab);
                size_t a = ket_ab.p;
                size_t b = ket_ab.q;
                double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                if ( (d_ea + d_eb + d_ei) > Z.modelspace->dE3max ) continue;
                double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                if (std::abs(occfactor)<1e-6) continue;
                double symm_factor =  (a==b) ? 1.0 : 2.0;
                double xkjab = X2.GetTBME_J(Jab,Jab,k,j,a,b);
                double ykjab = Y2.GetTBME_J(Jab,Jab,k,j,a,b);
                double xabilmn = X3.GetME_pn(Jab,J2,twoJ,a,b,i,l,m,n);
                double yabilmn = Y3.GetME_pn(Jab,J2,twoJ,a,b,i,l,m,n);
                z_ijklmn += 0.5  * symm_factor * occfactor * hats * sixj * (xkjab * yabilmn - ykjab * xabilmn);
//              if (ch==0 and ibra==4 and iket==4)
//              {
//                 std::cout << "Pik  iket_ab = " << iket_ab << "   adding to z   0.5 * " << symm_factor << " * " << occfactor << " * " << hats
//                           << " * " << sixj << " * " << xkjab << " " << yabilmn << "  =>  " << z_ijklmn << std::endl;;
//              }
              }// for iket_ab
            }// if sixj nonzero
          } // if J, Tz and parity check for term 2

          // Permute  Pjk
          if ( ((2*tbc_ab.Tz + oj.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + oj.l + Tbc.parity)%2==0) 
                     and (std::abs(2*Jab-oj.j2)<=twoJ) and ((2*Jab+oj.j2)>=twoJ) )
          {
             double sixj = Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,Jab);
             if (std::abs(sixj)>1e-6)
             {
               double hats = sqrt( (2*J1+1)*(2*Jab+1));
               int phase = Z.modelspace->phase( (oj.j2+ok.j2)/2-J1-Jab);
               for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
               {
                 auto& ket_ab = tbc_ab.GetKet(iket_ab);
                 size_t a = ket_ab.p;
                 size_t b = ket_ab.q;
                 double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                 double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                 if ( (d_ea + d_eb + d_ej) > Z.modelspace->dE3max ) continue;
                 double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                 if (std::abs(occfactor)<1e-6) continue;
                 double symm_factor =  (a==b) ? 1.0 : 2.0;
                 double xikab = X2.GetTBME_J(Jab,Jab,i,k,a,b);
                 double yikab = Y2.GetTBME_J(Jab,Jab,i,k,a,b);
                 double xabjlmn = X3.GetME_pn(Jab,J2,twoJ,a,b,j,l,m,n);
                 double yabjlmn = Y3.GetME_pn(Jab,J2,twoJ,a,b,j,l,m,n);
                 z_ijklmn -= 0.5 * symm_factor  * occfactor * phase * hats * sixj * (xikab * yabjlmn - yikab * xabjlmn);
//              if (ch==0 and ibra==4 and iket==4)
//              {
//                 std::cout << "Pjk  iket_ab = " << iket_ab << "   adding to z   0.5 * " << symm_factor << " * " << occfactor << " * "
//                           << phase << " * " << hats << " * " << sixj << " * " << xikab << " " << yabjlmn << "  =>  " << z_ijklmn << std::endl;
//              }
               }// for iket_ab
            }// if sixj nonzero
          } // if  J, Tz and parity check for term 3
        
          // direct Y3 X2 term
          if ( ((2*tbc_ab.Tz + on.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + on.l + Tbc.parity)%2==0) and Jab == J2 )
          {
            for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
            {
              auto& ket_ab = tbc_ab.GetKet(iket_ab);
              size_t a = ket_ab.p;
              size_t b = ket_ab.q;
              double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
              if ( (d_ea + d_eb + d_en) > Z.modelspace->dE3max ) continue;
              double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
              if (std::abs(occfactor)<1e-6) continue;
              double symm_factor =  (a==b) ? 1.0 : 2.0;
              double xablm = X2.GetTBME_J(Jab,Jab,a,b,l,m);
              double yablm = Y2.GetTBME_J(Jab,Jab,a,b,l,m);
              double xijkabn = X3.GetME_pn(J1,J2,twoJ,i,j,k,a,b,n);
              double yijkabn = Y3.GetME_pn(J1,J2,twoJ,i,j,k,a,b,n);
              z_ijklmn -= 0.5 * symm_factor  * occfactor * (xablm * yijkabn - yablm * xijkabn);
            }// for iket_ab
          } // if  Tz and parity check for term 4


          // Permute  Pln
          if ( ((2*tbc_ab.Tz + ol.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + ol.l + Tbc.parity)%2==0) 
                     and (std::abs(2*Jab-ol.j2)<=twoJ) and ((2*Jab+ol.j2)>=twoJ) )
          {
             double sixj = Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,Jab);
             if (std::abs(sixj)>1e-6) 
             {
               double hats = sqrt( (2*J2+1)*(2*Jab+1));
               for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
               {
                 auto& ket_ab = tbc_ab.GetKet(iket_ab);
                 size_t a = ket_ab.p;
                 size_t b = ket_ab.q;
                 double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                 double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                 if ( (d_ea + d_eb + d_el) > Z.modelspace->dE3max ) continue;
                 double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                 if (std::abs(occfactor)<1e-6) continue;
                 double symm_factor =  (a==b) ? 1.0 : 2.0;
                 double xabnm = X2.GetTBME_J(Jab,Jab,a,b,n,m);
                 double yabnm = Y2.GetTBME_J(Jab,Jab,a,b,n,m);
                 double xijkabl = X3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,l);
                 double yijkabl = Y3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,l);
                 z_ijklmn -= 0.5 * symm_factor  * occfactor * hats * sixj * (xabnm * yijkabl - yabnm * xijkabl);
               }// for iket_ab
             }// if sixj nonzero
          } // if  Tz and parity check for term 5


          // Permute  Pmn
          if ( ((2*tbc_ab.Tz + om.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + om.l + Tbc.parity)%2==0) 
                     and (std::abs(2*Jab-om.j2)<=twoJ) and ((2*Jab+om.j2)>=twoJ) )
          {
              double sixj = Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,Jab);
             if (std::abs(sixj)>1e-6) 
             {
              double hats = sqrt( (2*J2+1)*(2*Jab+1));
              int phase = Z.modelspace->phase( (om.j2+on.j2)/2-J2-Jab);
              for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
              {
                auto& ket_ab = tbc_ab.GetKet(iket_ab);
                size_t a = ket_ab.p;
                size_t b = ket_ab.q;
                double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                if ( (d_ea + d_eb + d_em) > Z.modelspace->dE3max ) continue;
                double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                if (std::abs(occfactor)<1e-6) continue;
                double symm_factor =  (a==b) ? 1.0 : 2.0;
                double xabln = X2.GetTBME_J(Jab,Jab,a,b,l,n);
                double yabln = Y2.GetTBME_J(Jab,Jab,a,b,l,n);
                double xijkabm = X3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,m);
                double yijkabm = Y3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,m);
                z_ijklmn += 0.5 * symm_factor * phase  * occfactor * hats * sixj * (xabln * yijkabm - yabln * xijkabm);
              }// for iket_ab
            }// if sixj nonzero
          } // if  Tz and parity check for term 6

        }// for ch_ab

//        if (ch==0 and ibra<5 and iket<5)
//        {
//          std::cout << __func__ << "   ch = " << ch << "  ibra =  " <<  ibra << "  iket = " << iket 
//                    << "   ijklmn " << bra.p << " " << bra.q << " " << bra.r << " " << ket.p << " " << ket.q << " " << ket.r
//                    << "  Jij Jlm twoJ " << bra.Jpq << " " << ket.Jpq << "  " << twoJ 
//                    << "  zijklmn = " << z_ijklmn << std::endl;
//        }
        // no normalization to be done. we store un-normalized 3body matrix elements.
        if (std::abs(z_ijklmn)>1e-6)
        {
//            Z3.AddToME_pn( J1, J2, twoJ, i,j,k,l,m,n, z_ijklmn );
            Z3.AddToME_pn_ch( ch,ch,ibra,iket, z_ijklmn );
        }

      }// for iket
    }// for ibra
  }// for ch

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
*/


//*****************************************************************************************
//
// i|   j|        k|  Uncoupled expression:
//  |    |         |    C_ijklmn = sum_{ab} (n`anb -nan`b) P(ij/k)P(lm/n) Y_{ijalmb} X_{bkan}
//  |    |   *~~X~~|                
//  |    | a/ \b   |
//  |    |  \ /    |
//  |~~~[Y]~~*     |  Coupled expression:
//  |    |         |  C_{ijklmn}^{J1,J2,J} =  P^{J1,J}(ij/k)P^{J2,j}_{lm/n) sum_{ab} (n`anb -nan`b)  
// l|   m|        n|                          * sum_{J',J3} (2J'+1)(2J3+1) (-1)^{2J+k-n+J1-J2} 
//  |    |         |                          { J   J2  n  }                                                                       
//                                          * { J1  J'  a  } * Y_{ijalmb}^{J1,J2,J'} X_{bkan}^{J3}                                 
//                                            { k   b   J3 }
//                                              
//      Checked with UnitTest and passed                                                                      
//                                                                           
void comm233_phss( const Operator& X, const Operator& Y, Operator& Z )
{

  double tstart = omp_get_wtime();
  double t_internal = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();
  int hX = X.IsHermitian() ? 1 : -1;
  int hY = Y.IsHermitian() ? 1 : -1;

  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels(); // number of regular 2b channels
  size_t nch2_CC = Z.modelspace->TwoBodyChannels_CC.size(); // number of particle hole channels
  std::deque<std::vector<size_t>> ph_kets_cc(nch2_CC); // list of ph 2b kets for each CC channel
  std::deque<std::map<size_t,size_t>> ph_kets_cc_lookup(nch2_CC); // map i j1 to  matrix index
  std::deque<std::vector<double>> occfactors(nch2_CC);  // occupation factors for ab`


  auto hash_key_lmij = [](size_t l,size_t m, size_t Jlm, size_t i, size_t j, size_t Jij){return ( l + (m<<12) + (i<<24) + (j<<36) + (Jlm<<48) + (Jij<<56) );};
  auto unhash_key_lmij = [](size_t& l, size_t& m, size_t& Jlm, size_t& i,size_t& j, size_t& Jij, size_t key){ l=(key & 0xFFFL);m=((key>>12)&0xFFFL);i=((key>>24)&0xFFFL);j=((key>>36)&0xFFFL); Jlm=((key>>48)&0xFFL); Jij=((key>>56)&0xFFL); }; //  0xF = 15 = 1111 (4bits), so 0xFFF is 12 bits of 1's. 0xFFFL makes it a long
  std::deque<std::unordered_map<size_t,size_t>> ph_kets3_lookup(nch2_CC); // map  lm Jlm  (ij Jij)` to matrix index
//  std::deque<std::map<std::array<int,6>,size_t>> ph_kets3_lookup(nch2_CC); // map  lm Jlm  (ij Jij)` to matrix index

  std::deque<arma::mat> Z3_ph(nch2_CC); // matrices for ph transformed Z3. This will be the big one.

  double occnat_factor_max = 0;
  for ( auto i : Z.modelspace->all_orbits)
  {
     double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
     occnat_factor_max = std::max( occnat_factor_max, occnat_i*(1-occnat_i) );
  }


//  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch_cc=0; ch_cc<nch2_CC; ch_cc++)
  {
    TwoBodyChannel_CC& tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int Jph = tbc_CC.J;
    size_t nkets_CC = tbc_CC.GetNumberKets();
    // Figure out which 2b ph kets will contribute to ab`
    for ( size_t iket=0; iket<nkets_CC; iket++)
    {
      Ket& ket = tbc_CC.GetKet(iket);
      double occfactor = ket.op->occ - ket.oq->occ;
      if (std::abs(occfactor) > 1e-7 )
      {
        ph_kets_cc[ch_cc].push_back( iket );
        ph_kets_cc_lookup[ch_cc][ iket] =  ph_kets_cc[ch_cc].size()-1;
        occfactors[ch_cc].push_back( occfactor );
      }
    }// for iket

//    size_t nkets_ph = ph_kets_cc[ch_cc].size();


   auto& ph_kets3 = ph_kets3_lookup[ch_cc];
   size_t nkets3 = 0;
   // next, figure out how many states we'll need for the 3b matrices
   for (size_t ch_lm=0; ch_lm<nch2; ch_lm++)
   {
     TwoBodyChannel& tbc_lm = Z.modelspace->GetTwoBodyChannel(ch_lm);
     size_t nkets_lm = tbc_lm.GetNumberKets();
     for (size_t ch_ij=0; ch_ij<nch2; ch_ij++)
     {
       TwoBodyChannel& tbc_ij = Z.modelspace->GetTwoBodyChannel(ch_ij);
       if ( (tbc_lm.parity + tbc_ij.parity)%2 != tbc_CC.parity ) continue;
       if ( std::abs(tbc_lm.Tz - tbc_ij.Tz) != tbc_CC.Tz ) continue; // This abs bit may make things expensive...
       if ( (std::abs( tbc_lm.J - tbc_ij.J) > tbc_CC.J) or ( tbc_lm.J + tbc_ij.J) < tbc_CC.J) continue;

       size_t nkets_ij = tbc_ij.GetNumberKets();
       for ( size_t iket_lm=0; iket_lm<nkets_lm; iket_lm++)
       {
         Ket& ket_lm = tbc_lm.GetKet(iket_lm);
         size_t l = ket_lm.p;
         size_t m = ket_lm.q;
         int Jlm  = tbc_lm.J;
         Orbit& ol = Z.modelspace->GetOrbit(l);
         Orbit& om = Z.modelspace->GetOrbit(m);
         double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
         double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
         double occnat_l = ol.occ_nat;
         double occnat_m = om.occ_nat;
         if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
         if ( (d_el + d_em ) > Z.modelspace->dE3max ) continue;
 
         for (size_t iket_ij=0; iket_ij<nkets_ij; iket_ij++)
         {
           Ket& ket_ij = tbc_ij.GetKet(iket_ij);
           size_t i = ket_ij.p;
           size_t j = ket_ij.q;
           int Jij  = tbc_ij.J;
           // Only compute one ordering of   lmij  vs ijlm  since we can get the other ordering from symmetry
           if ( (std::min(i,j) > std::min(l,m)) or ((std::min(i,j)==std::min(l,m)) and (std::max(i,j) > std::max(l,m)))
               or ((std::min(i,j)==std::min(l,m)) and (std::max(i,j)==std::max(l,m)) and (Jij>Jlm)) ) continue;

           Orbit& oi = Z.modelspace->GetOrbit(i);
           Orbit& oj = Z.modelspace->GetOrbit(j);
           double occnat_i = oi.occ_nat;
           double occnat_j = oj.occ_nat;
           double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
           double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
           if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
           if ( (d_ei + d_ej ) > Z.modelspace->dE3max ) continue;

           size_t key = hash_key_lmij( l,m,(size_t)Jlm,i,j,(size_t)Jij);
           ph_kets3[ key ] = nkets3;
//           ph_kets3[{l,m,Jlm,i,j,Jij}] = nkets3;
           nkets3++;
         }
       }

     }// for ch_ij
   }// for ch_lm

   // so now we can allocate the ph 3-body matrices

   Z3_ph[ch_cc].zeros( 2*nkets_CC, nkets3 ); // Zbar_knlmij
  }// for ch_cc

  Z.profiler.timer["comm233_setup_lists"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // Compute the particle-hole recoupled matrix elements of Z
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch_cc=0; ch_cc<nch2_CC; ch_cc++)
  {
    TwoBodyChannel_CC& tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int Jph = tbc_CC.J;
    size_t nkets_CC = tbc_CC.GetNumberKets();
    size_t nkets_ph = ph_kets_cc[ch_cc].size();
    size_t nkets3 = Z3_ph[ch_cc].n_cols;

     arma::mat X2_ph( 2*nkets_CC, 2*nkets_ph, arma::fill::zeros ); // we need to store both permutations ij ji because the ph transformed 2bmes don't have as nice a symmetry
     arma::mat Y2_ph( 2*nkets_CC, 2*nkets_ph, arma::fill::zeros );

     arma::mat X3_ph( 2*nkets_ph,   nkets3, arma::fill::zeros ); // Xbar_ablmij
     arma::mat Y3_ph( 2*nkets_ph,   nkets3, arma::fill::zeros );

    // now we fill the matrices
    for (size_t ibra_ij=0; ibra_ij<nkets_CC; ibra_ij++)
    {
      Ket& bra_ij = tbc_CC.GetKet(ibra_ij);
      size_t i = bra_ij.p;
      size_t j = bra_ij.q;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      for ( size_t iket_ab=0; iket_ab<nkets_ph; iket_ab++)
      {
        Ket& ket_ab = tbc_CC.GetKet( ph_kets_cc[ch_cc].at(iket_ab));
        size_t a = ket_ab.p;
        size_t b = ket_ab.q;
        Orbit& oa = Z.modelspace->GetOrbit(a);
        Orbit& ob = Z.modelspace->GetOrbit(b);
        double ja = 0.5 * oa.j2;
        double jb = 0.5 * ob.j2;
//        double occ_factor = occfactors[ch_cc][iket_ab];

        double occ_factor = oa.occ-ob.occ;

        // XJph_ij`ab` = -sum J'  (2J'+1) { i j J_ph } XJ'_ibaj
        //                                { a b J'   }

        double xijab = 0;
        double xjiab = 0;
        double yijab = 0;
        double yjiab = 0;
        int Jprime_min = std::max( std::abs(oi.j2-ob.j2), std::abs(oa.j2-oj.j2) )/2;
        int Jprime_max = std::min( oi.j2+ob.j2, oa.j2+oj.j2 )/2;
        for (int Jprime=Jprime_min; Jprime<=Jprime_max; Jprime++)
        {
          double sixj_ij = Z.modelspace->GetSixJ(ji, jj, Jph, ja,jb, Jprime);
          xijab -= (2*Jprime+1) * sixj_ij * X2.GetTBME_J( Jprime, i,b,a,j );
          yijab -= (2*Jprime+1) * sixj_ij * Y2.GetTBME_J( Jprime, i,b,a,j );
        }
        Jprime_min = std::max( std::abs(oj.j2-ob.j2), std::abs(oa.j2-oi.j2) )/2;
        Jprime_max = std::min( oj.j2+ob.j2, oa.j2+oi.j2 )/2;
        for (int Jprime=Jprime_min; Jprime<=Jprime_max; Jprime++)
        {
          double sixj_ji = Z.modelspace->GetSixJ(jj, ji, Jph, ja,jb, Jprime);
          xjiab -= (2*Jprime+1) * sixj_ji * X2.GetTBME_J( Jprime, j,b,a,i );
          yjiab -= (2*Jprime+1) * sixj_ji * Y2.GetTBME_J( Jprime, j,b,a,i );
        }
        int phase_X = hX * Z.modelspace->phase( (oi.j2+oj.j2+oa.j2+ob.j2)/2);
        int phase_Y = hY * Z.modelspace->phase( (oi.j2+oj.j2+oa.j2+ob.j2)/2);

        X2_ph( ibra_ij,          iket_ab ) = xijab * occ_factor;
        Y2_ph( ibra_ij,          iket_ab ) = yijab * occ_factor;
        X2_ph( ibra_ij+nkets_CC, iket_ab ) = xjiab * occ_factor;
        Y2_ph( ibra_ij+nkets_CC, iket_ab ) = yjiab * occ_factor;
        X2_ph( ibra_ij,          iket_ab+nkets_ph) = phase_X * xjiab * (-occ_factor);
        Y2_ph( ibra_ij,          iket_ab+nkets_ph) = phase_Y * yjiab * (-occ_factor);
        X2_ph( ibra_ij+nkets_CC, iket_ab+nkets_ph) = phase_X * xijab * (-occ_factor);
        Y2_ph( ibra_ij+nkets_CC, iket_ab+nkets_ph) = phase_Y * yijab * (-occ_factor);

      }//for iket_ab
    }// for ibra_ij


    // now fill the 3-body matrices
    auto& ph_kets3 = ph_kets3_lookup[ch_cc];
    for (auto iter_lmij : ph_kets3 )
    {
      size_t index_lmij = iter_lmij.second;

      size_t key = iter_lmij.first;
      size_t l,m,i,j,Jlm_tmp,Jij_tmp;
      unhash_key_lmij(l,m,Jlm_tmp,i,j,Jij_tmp, key);
      int Jlm = (int)Jlm_tmp;
      int Jij = (int)Jij_tmp;


//      size_t index_ijlm = ph_kets3.at({i,j,Jij,l,m,Jlm});
      Orbit& ol = Z.modelspace->GetOrbit(l);
      Orbit& om = Z.modelspace->GetOrbit(m);
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      double occnat_l = ol.occ_nat;
      double occnat_m = om.occ_nat;
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
      double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);

//      for ( size_t iket_ab=0; iket_ab<nkets_ph; iket_ab++)
      for ( size_t iket_ab=0; iket_ab<2*nkets_ph; iket_ab++)
      {
        Ket& ket_ab = tbc_CC.GetKet( ph_kets_cc[ch_cc].at(iket_ab%nkets_ph));

        size_t a = (iket_ab < nkets_ph) ? ket_ab.p : ket_ab.q;
        size_t b = (iket_ab < nkets_ph) ? ket_ab.q : ket_ab.p;
        Orbit& oa = Z.modelspace->GetOrbit(a);
        Orbit& ob = Z.modelspace->GetOrbit(b);
        double ja = 0.5 * oa.j2;
        double jb = 0.5 * ob.j2;

        double occnat_a = oa.occ_nat;
        double occnat_b = ob.occ_nat;
        double d_ea = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
        double d_eb = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
        double occ_factor = oa.occ-ob.occ;
        if (std::abs(occ_factor)<1e-6) continue;

        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_a*(1-occnat_a) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_b*(1-occnat_b) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        if ( (d_ei + d_ej + d_ea) > Z.modelspace->dE3max ) continue;
        if ( (d_el + d_em + d_eb) > Z.modelspace->dE3max ) continue;


        double X3bar_ablmij = 0;
        double Y3bar_ablmij = 0;

        if ( (ob.l+ol.l+om.l+oa.l+oi.l+oj.l)%2>0 ) continue;
        if ( (ob.tz2+ol.tz2+om.tz2)!=(oa.tz2+oi.tz2+oj.tz2) ) continue;


        int twoJp_min = std::max( std::abs( oa.j2-2*Jij), std::abs(ob.j2-2*Jlm));
        int twoJp_max = std::min( ( oa.j2+2*Jij), (ob.j2+2*Jlm));
        for ( int twoJp=twoJp_min; twoJp<=twoJp_max; twoJp+=2)
        {
          int phase_y = Z.modelspace->phase( (oa.j2+twoJp)/2 );

          double xablmij = 0;
          double yablmij = 0;
          // get the X3 and Y3 elements in one go (to avoid doing the recoupling twice)
          if ( x3_allocated and y3_allocated )
          {
            auto xandy = Y3.GetME_pn_TwoOps( Jij,Jlm,twoJp, i,j,a, l,m,b, X3,Y3 );
            xablmij = xandy[0];
            yablmij = xandy[1];
          }
          else if (y3_allocated)
          {
             yablmij = Y3.GetME_pn(Jij, Jlm, twoJp,i,j,a,l,m,b);
          }
          else if (x3_allocated)
          {
             xablmij = X3.GetME_pn(Jij, Jlm, twoJp,i,j,a,l,m,b);
          }

          double sixjy = Z.modelspace->GetSixJ( Jij, Jlm, Jph, jb, ja, 0.5*twoJp);
          X3bar_ablmij -= phase_y * (twoJp+1) * sixjy * xablmij;
          Y3bar_ablmij -= phase_y * (twoJp+1) * sixjy * yablmij;
        }// for twoJp

        size_t index_ab =  iket_ab;
        size_t index_ba =  (iket_ab + nkets_ph) % (2*nkets_ph);
        X3_ph(index_ab, index_lmij) = X3bar_ablmij ;
        Y3_ph(index_ab, index_lmij) = Y3bar_ablmij ;

       }// for iket_ab
      }// for iter_lmij

      Z3_ph[ch_cc] = X2_ph * Y3_ph - Y2_ph * X3_ph;

  }// for ch_cc
//  std::cout << "DONE WITH THE CC BIT" << std::endl;

  Z.profiler.timer["comm233_fll_matrices"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // Now transform it back into standard coupling
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      int J1 = bra.Jpq;


      std::vector<std::array<size_t,3>> ijk = { {i,j,k}, {k,j,i}, {i,k,j} };

      std::vector<int> J1p_min = {J1,  std::max(std::abs(ok.j2-oj.j2),std::abs(twoJ-oi.j2) )/2,   std::max(std::abs(oi.j2-ok.j2), std::abs(twoJ-oj.j2) )/2 };
      std::vector<int> J1p_max = {J1,  std::min(ok.j2+oj.j2, twoJ+oi.j2)/2 , std::min(oi.j2+ok.j2, twoJ+oj.j2)/2 };
      std::vector<std::vector<double>> recouple_ijk = {{1},{},{} };
      for (int J1p=J1p_min[1]; J1p<=J1p_max[1]; J1p++)
           recouple_ijk[1].push_back( sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p) );

      for (int J1p=J1p_min[2]; J1p<=J1p_max[2]; J1p++)
           recouple_ijk[2].push_back( -Z.modelspace->phase((oj.j2+ok.j2)/2+J1+J1p)*sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p) );


      for (size_t iket=ibra; iket<nkets3; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        int J2 = ket.Jpq;


        std::vector<std::array<size_t,3>> lmn = { {l,m,n}, {n,m,l}, {l,n,m} };
        std::vector<int> J2p_min = {J2,  std::max(std::abs(on.j2-om.j2),std::abs(twoJ-ol.j2) )/2,   std::max(std::abs(ol.j2-on.j2), std::abs(twoJ-om.j2) )/2 };
        std::vector<int> J2p_max = {J2,  std::min(on.j2+om.j2, twoJ+ol.j2)/2 , std::min(ol.j2+on.j2, twoJ+om.j2)/2 };
        std::vector<std::vector<double>> recouple_lmn = {{1},{},{} };
            
        for (int J2p=J2p_min[1]; J2p<=J2p_max[1]; J2p++)
           recouple_lmn[1].push_back( sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p) );

        for (int J2p=J2p_min[2]; J2p<=J2p_max[2]; J2p++)
             recouple_lmn[2].push_back( -Z.modelspace->phase((om.j2+on.j2)/2+J2+J2p)*sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p) );


        double z_ijklmn = 0;


        for ( int perm_ijk=0; perm_ijk<3; perm_ijk++ )
        {
          size_t I1 = ijk[perm_ijk][0];
          size_t I2 = ijk[perm_ijk][1];
          size_t I3 = ijk[perm_ijk][2];
          Orbit& o1 = Z.modelspace->GetOrbit( I1 );
          Orbit& o2 = Z.modelspace->GetOrbit( I2 );
          Orbit& o3 = Z.modelspace->GetOrbit( I3 );
          double occnat_1 = o1.occ_nat;
          double occnat_2 = o2.occ_nat;
          double d_e1 = std::abs( 2*o1.n + o1.l - e_fermi[o1.tz2]);
          double d_e2 = std::abs( 2*o2.n + o2.l - e_fermi[o2.tz2]);
          double j3 = 0.5*o3.j2;
          for (int J1p=J1p_min[perm_ijk]; J1p<=J1p_max[perm_ijk]; J1p++)
          {
            if ((I1==I2) and J1p%2!=0 ) continue;
            double rec_ijk = recouple_ijk[perm_ijk].at(J1p-J1p_min[perm_ijk]);
            for ( int perm_lmn=0; perm_lmn<3; perm_lmn++ )
            {
              size_t I4 = lmn[perm_lmn][0];
              size_t I5 = lmn[perm_lmn][1];
              size_t I6 = lmn[perm_lmn][2];
              Orbit& o4 = Z.modelspace->GetOrbit( I4 );
              Orbit& o5 = Z.modelspace->GetOrbit( I5 );
              Orbit& o6 = Z.modelspace->GetOrbit( I6 );
              double occnat_4 = o4.occ_nat;
              double occnat_5 = o5.occ_nat;
              double d_e4 = std::abs( 2*o4.n + o4.l - e_fermi[o4.tz2]);
              double d_e5 = std::abs( 2*o5.n + o5.l - e_fermi[o5.tz2]);
              double j6 = 0.5*o6.j2;
              for (int J2p=J2p_min[perm_lmn]; J2p<=J2p_max[perm_lmn]; J2p++)
              {
                if ((I4==I5) and J2p%2!=0 ) continue;
                double rec_lmn = recouple_lmn[perm_lmn].at(J2p-J2p_min[perm_lmn]);

                double z_123456 = 0;
                int Jph_min = std::max( std::abs(2*J1p-2*J2p), std::abs(o3.j2-o6.j2) )/2;
                int Jph_max = std::min( 2*J1p+2*J2p , o3.j2+o6.j2 )/2;
                int Tz_ph = std::abs(o3.tz2 - o6.tz2)/2;
                int parity_ph = (o3.l + o6.l)%2;

                for ( int Jph=Jph_min; Jph<=Jph_max; Jph++)
                {
                   size_t ch_ph = Z.modelspace->GetTwoBodyChannelIndex(Jph,parity_ph,Tz_ph);
                   TwoBodyChannel_CC& tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch_ph);
                   size_t nkets_CC = tbc_CC.GetNumberKets();


                   // we only compute one ordering of the right side of Zbar, but we can get the other ordering by symmetry
                   // so we check here which ordering we need, and pick up the corresponding phase factor if needed.
                   size_t index_36 = tbc_CC.GetLocalIndex(I3,I6);
                   size_t key_4512 = hash_key_lmij(std::min(I4,I5),std::max(I4,I5),(size_t)J2p,std::min(I1,I2),std::max(I1,I2),(size_t)J1p);
                   int phase_45 = (I4>I5) ? -Z.modelspace->phase( (o4.j2+o5.j2)/2 - J2p) : 1;
                   int phase_12 = (I1>I2) ? -Z.modelspace->phase( (o1.j2+o2.j2)/2 - J1p) : 1;
                   int flip_phase = 1;

                   // if any of these are true, we're looking for the not-stored ordering
                   if ( (std::min(I1,I2) > std::min(I4,I5)) or ((std::min(I1,I2)==std::min(I4,I5)) and (std::max(I1,I2) > std::max(I4,I5)))
                        or ((std::min(I1,I2)==std::min(I4,I5)) and (std::max(I1,I2)==std::max(I4,I5)) and (J1p>J2p)) )
                   {
                     index_36 = tbc_CC.GetLocalIndex(I6,I3);
                     key_4512 = hash_key_lmij(std::min(I1,I2),std::max(I1,I2),(size_t)J1p,std::min(I4,I5),std::max(I4,I5),(size_t)J2p);
                     flip_phase = Z.modelspace->phase( (o3.j2-o6.j2)/2 );
                   }
                   size_t index_4512 = ph_kets3_lookup[ch_ph].at(key_4512);

                   // < 36` Jph | Z | 45 J2p  (12 J1p)` ; Jph >
                   double z_364512 =  flip_phase * phase_45*phase_12 * Z3_ph[ch_ph](index_36,index_4512);

                   double sixj_ph = Z.modelspace->GetSixJ(J1p,j3,Jtot,  j6,J2p,Jph);
                   int phase_ph = -Z.modelspace->phase( (twoJ +o3.j2)/2);
                   z_123456 -=  phase_ph * (2*Jph+1) * sixj_ph * z_364512;
                 }// for Jph

                // include the recoupling coefficients needed for the permutations
                z_ijklmn += rec_ijk * rec_lmn * z_123456;

              }// for J2p
            }// for perm_lmn
          }// for J1p
        }// for perm_ijk


        Z3.AddToME_pn_ch( ch3, ch3, ibra,iket, z_ijklmn);
      }// for iket
    }// for ibra
  }//for ch3
  Z.profiler.timer["comm233_pph_recouple"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}





void comm233_phss_debug( const Operator& X, const Operator& Y, Operator& Z )
{

  double tstart = omp_get_wtime();
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();


  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      int J1 = bra.Jpq;
      for (size_t iket=ibra; iket<nkets3; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        int J2 = ket.Jpq;

        double z_ijklmn = 0;

        for ( auto& a : Z.modelspace->all_orbits )
        {
         Orbit& oa = Z.modelspace->GetOrbit(a);
         double ja = 0.5 * oa.j2;
         double occnat_a = oa.occ_nat;
         double d_ea = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
         bool keep_ija = ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_a*(1-occnat_a) ) >= Z.modelspace->GetOccNat3Cut() ) ;
         bool keep_kja = ( (occnat_k*(1-occnat_k) * occnat_j*(1-occnat_j) * occnat_a*(1-occnat_a) ) >= Z.modelspace->GetOccNat3Cut() ) ;
         bool keep_ika = ( (occnat_i*(1-occnat_i) * occnat_k*(1-occnat_k) * occnat_a*(1-occnat_a) ) >= Z.modelspace->GetOccNat3Cut() ) ;
         keep_ija = keep_ija and (d_ei+d_ej+d_ea <= Z.modelspace->dE3max );
         keep_kja = keep_kja and (d_ek+d_ej+d_ea <= Z.modelspace->dE3max );
         keep_ika = keep_ika and (d_ei+d_ek+d_ea <= Z.modelspace->dE3max );
         if (not ( keep_ija or keep_kja or keep_ika) ) continue;

         for ( auto& b : Z.modelspace->all_orbits )
         {
//           if (b<a) continue;
           Orbit& ob = Z.modelspace->GetOrbit(b);
           double jb = 0.5 * ob.j2;
           double occnat_b = ob.occ_nat;
           double d_eb = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
           double occ_factor = oa.occ-ob.occ;
           if (std::abs(occ_factor)<1e-6) continue;

           bool keep_lmb = ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_b*(1-occnat_b) ) >= Z.modelspace->GetOccNat3Cut() ) ;
           bool keep_lnb = ( (occnat_l*(1-occnat_l) * occnat_n*(1-occnat_n) * occnat_b*(1-occnat_b) ) >= Z.modelspace->GetOccNat3Cut() ) ;
           bool keep_nmb = ( (occnat_n*(1-occnat_n) * occnat_m*(1-occnat_m) * occnat_b*(1-occnat_b) ) >= Z.modelspace->GetOccNat3Cut() ) ;
           keep_lmb = keep_lmb and (d_el+d_em+d_eb <= Z.modelspace->dE3max );
           keep_lnb = keep_lnb and (d_el+d_en+d_eb <= Z.modelspace->dE3max );
           keep_nmb = keep_nmb and (d_en+d_em+d_eb <= Z.modelspace->dE3max );
           if (not ( keep_lmb or keep_lnb or keep_nmb) ) continue;


           // Direct term i.e. Z1
           if ( (ob.l+ok.l+oa.l+on.l)%2==0   and (  (ob.tz2+ok.tz2)==(oa.tz2+on.tz2) ) and keep_ija and keep_lmb )
           {
            int Jx_min = std::max( std::abs(ob.j2-ok.j2),std::abs(oa.j2-on.j2) )/2;
            int Jx_max = std::min( ob.j2+ok.j2 ,  oa.j2+on.j2 )/2;
            int twoJJ_min = std::max( std::abs(oa.j2-2*J1) , std::abs(ob.j2-2*J2) );
            int twoJJ_max = std::min( oa.j2+2*J1 ,  ob.j2+2*J2 );
            int phase = -Z.modelspace->phase((ok.j2+on.j2)/2+J1+J2);
            for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
            {
              double xijalmb = X3.GetME_pn( J1,J2,twoJJ, i,j,a,l,m,b);
              double yijalmb = Y3.GetME_pn( J1,J2,twoJJ, i,j,a,l,m,b);
              double JJtot = 0.5 * twoJJ;
              for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
              {
                double xbkan = X2.GetTBME_J(Jx, b,k,a,n);
                double ybkan = Y2.GetTBME_J(Jx, b,k,a,n);
                double hats = (twoJJ+1)*(2*Jx+1);
                double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2,  jk,J1,Jtot, Jx,ja,jn);
                z_ijklmn -= occ_factor * hats * phase * ninej * (xbkan * yijalmb - ybkan*xijalmb);
              }// for Jx
            }// for twoJ
           }// Z1 block


           // Z2  ~ Pik Z1     Xbian Ykjalmb
           if ( (ob.l+oi.l+oa.l+on.l)%2==0   and (  (ob.tz2+oi.tz2)==(oa.tz2+on.tz2) ) and keep_kja and keep_lmb)
           {
            int J1p_min = std::abs(ok.j2-oj.j2)/2;
            int J1p_max = ( ok.j2+oj.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-oi.j2),std::abs(oa.j2-on.j2) )/2;
            int Jx_max = std::min( ob.j2+oi.j2 ,  oa.j2+on.j2 )/2;
            for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
            {
              double sixj = Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p);
              if (std::abs(sixj)<1e-6) continue;
              int twoJJ_min = std::max( std::abs(ob.j2-2*J2),std::abs(oa.j2-2*J1p) );
              int twoJJ_max = std::min( ob.j2+2*J2 ,  oa.j2+2*J1p );
              int phase = Z.modelspace->phase((oi.j2+on.j2)/2+J1p+J2);
              for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
              {
                double xkjalmb = X3.GetME_pn( J1p,J2,twoJJ, k,j,a,l,m,b);
                double ykjalmb = Y3.GetME_pn( J1p,J2,twoJJ, k,j,a,l,m,b);
                double JJtot = 0.5 * twoJJ;
                for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                {
                  double xbian = X2.GetTBME_J(Jx, b,i,a,n);
                  double ybian = Y2.GetTBME_J(Jx, b,i,a,n);
                  double hats = (twoJJ+1)*(2*Jx+1) * sqrt( (2*J1+1)*(2*J1p+1) );
                  double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2,  ji,J1p,Jtot, Jx,ja,jn);
                  z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbian*ykjalmb - ybian*xkjalmb);
                }// for Jx
              }// for twoJ
            }// for J1p
           }// Z2 block



           // Z3  ~ Pjk Z1      Xbjan Yikalmb
           if ( (ob.l+oj.l+oa.l+on.l)%2==0   and (  (ob.tz2+oj.tz2)==(oa.tz2+on.tz2) ) and keep_ika and keep_lmb)
           {
            int J1p_min = std::abs(oi.j2-ok.j2)/2;
            int J1p_max = ( oi.j2+ok.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-oj.j2),std::abs(oa.j2-on.j2) )/2;
            int Jx_max = std::min( ob.j2+oj.j2 ,  oa.j2+on.j2 )/2;
            for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
            {
              double sixj = Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p);
              if (std::abs(sixj)<1e-6) continue;
              int twoJJ_min = std::max( std::abs(ob.j2-2*J2),std::abs(oa.j2-2*J1p) );
              int twoJJ_max = std::min( ob.j2+2*J2 ,  oa.j2+2*J1p );
              int phase = Z.modelspace->phase((ok.j2+on.j2)/2+J1+J2);
              for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
              {
                double xikalmb = X3.GetME_pn( J1p,J2,twoJJ, i,k,a,l,m,b);
                double yikalmb = Y3.GetME_pn( J1p,J2,twoJJ, i,k,a,l,m,b);
                double JJtot = 0.5 * twoJJ;
                for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                {
                  double xbjan = X2.GetTBME_J(Jx, b,j,a,n);
                  double ybjan = Y2.GetTBME_J(Jx, b,j,a,n);
                  double hats = (twoJJ+1)*(2*Jx+1) * sqrt( (2*J1+1)*(2*J1p+1));
                  double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2,  jj,J1p,Jtot, Jx,ja,jn);
                  z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbjan*yikalmb - ybjan*xikalmb);
                }// for Jx
              }// for twoJ
            }// for J1p
           }// Z3 block




           // Z4  ~ Pln Z1    Xbkal Y ijanmb
           if ( (ob.l+ok.l+oa.l+ol.l)%2==0   and (  (ob.tz2+ok.tz2)==(oa.tz2+ol.tz2) ) and keep_ija and keep_nmb)
           {
            int J2p_min = std::abs(on.j2-om.j2)/2;
            int J2p_max = ( on.j2+om.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-ok.j2),std::abs(oa.j2-ol.j2) )/2;
            int Jx_max = std::min( ob.j2+ok.j2 ,  oa.j2+ol.j2 )/2;
            for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
            {
              double sixj = Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p);
              if (std::abs(sixj)<1e-6) continue;
              int twoJJ_min = std::max( std::abs(ob.j2-2*J2p),std::abs(oa.j2-2*J1) );
              int twoJJ_max = std::min( ob.j2+2*J2p ,  oa.j2+2*J1 );
              int phase = Z.modelspace->phase((ok.j2+ol.j2)/2+J1+J2p);
              for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
              {
                double xijanmb = X3.GetME_pn( J1,J2p,twoJJ, i,j,a,n,m,b);
                double yijanmb = Y3.GetME_pn( J1,J2p,twoJJ, i,j,a,n,m,b);
                double JJtot = 0.5 * twoJJ;
                for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                {
                  double xbkal = X2.GetTBME_J(Jx, b,k,a,l);
                  double ybkal = Y2.GetTBME_J(Jx, b,k,a,l);
                  double hats = (twoJJ+1)*(2*Jx+1)*sqrt((2*J2+1)*(2*J2p+1));
                  double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2p,  jk,J1,Jtot, Jx,ja,jl);
                  z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbkal*yijanmb - ybkal*xijanmb);
                }// for Jx
              }// for twoJ
            }// for J1p
           }// Z4 block



           // Z5  ~ Pmn Z1    Xbkam Yijalnb
           if ( (ob.l+ok.l+oa.l+om.l)%2==0   and (  (ob.tz2+ok.tz2)==(oa.tz2+om.tz2) ) and keep_ija and keep_lnb)
           {
            int J2p_min = std::abs(ol.j2-on.j2)/2;
            int J2p_max = ( ol.j2+on.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-ok.j2),std::abs(oa.j2-om.j2) )/2;
            int Jx_max = std::min( ob.j2+ok.j2 ,  oa.j2+om.j2 )/2;
            for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
            {
              double sixj = Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p);
              if (std::abs(sixj)<1e-6) continue;
              int twoJJ_min = std::max( std::abs(ob.j2-2*J2p),std::abs(oa.j2-2*J1) );
              int twoJJ_max = std::min( ob.j2+2*J2p ,  oa.j2+2*J1 );
              int phase = Z.modelspace->phase((ok.j2+on.j2)/2+J1+J2);
              for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
              {
                double xijalnb = X3.GetME_pn( J1,J2p,twoJJ, i,j,a,l,n,b);
                double yijalnb = Y3.GetME_pn( J1,J2p,twoJJ, i,j,a,l,n,b);
                double JJtot = 0.5 * twoJJ;
                for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                {
                  double xbkam = X2.GetTBME_J(Jx, b,k,a,m);
                  double ybkam = Y2.GetTBME_J(Jx, b,k,a,m);
                  double hats = (twoJJ+1)*(2*Jx+1)*sqrt((2*J2+1)*(2*J2p+1));
                  double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2p,  jk,J1,Jtot, Jx,ja,jm);
                  z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbkam*yijalnb - ybkam*xijalnb);
                }// for Jx
              }// for twoJ
            }// for J1p
           }// Z5 block



           // Z6  ~ Pik Pln Z1     Xbial Ykjanmb
           if ( (ob.l+oi.l+oa.l+ol.l)%2==0   and (  (ob.tz2+oi.tz2)==(oa.tz2+ol.tz2) ) and keep_kja and keep_nmb)
           {
            int J1p_min = std::abs(ok.j2-oj.j2)/2;
            int J1p_max = ( ok.j2+oj.j2 )/2;
            int J2p_min = std::abs(on.j2-om.j2)/2;
            int J2p_max = ( on.j2+om.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-oi.j2),std::abs(oa.j2-ol.j2) )/2;
            int Jx_max = std::min( ob.j2+oi.j2 ,  oa.j2+ol.j2 )/2;
            for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
            {
              double sixj1 = Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p);
              if (std::abs(sixj1)<1e-6) continue;
              for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
              {
                double sixj2 = Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p);
                if (std::abs(sixj2)<1e-6) continue;
                int twoJJ_min = std::max( std::abs(ob.j2-2*J2p),std::abs(oa.j2-2*J1p) );
                int twoJJ_max = std::min( ob.j2+2*J2p ,  oa.j2+2*J1p );
                int phase = -Z.modelspace->phase((oi.j2+ol.j2)/2+J1p+J2p);
                for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
                {
                  double xkjanmb = X3.GetME_pn( J1p,J2p,twoJJ, k,j,a,n,m,b);
                  double ykjanmb = Y3.GetME_pn( J1p,J2p,twoJJ, k,j,a,n,m,b);
                  double JJtot = 0.5 * twoJJ;
                  for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                  {
                    double xbial = X2.GetTBME_J(Jx, b,i,a,l);
                    double ybial = Y2.GetTBME_J(Jx, b,i,a,l);
                    double hats = (twoJJ+1)*(2*Jx+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1));
                    double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2p,  ji,J1p,Jtot, Jx,ja,jl);
                    z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbial*ykjanmb - ybial*xkjanmb);
                  }// for Jx
                }// for twoJ
              }// for J2p
            }// for J1p
           }// Z6 block


           // Z7  ~ Pjk Pln Z1     Xbjal Yikanmb
           if ( (ob.l+oj.l+oa.l+ol.l)%2==0   and (  (ob.tz2+oj.tz2)==(oa.tz2+ol.tz2) )  and keep_ika and keep_nmb)
           {
            int J1p_min = std::abs(oi.j2-ok.j2)/2;
            int J1p_max = ( oi.j2+ok.j2 )/2;
            int J2p_min = std::abs(on.j2-om.j2)/2;
            int J2p_max = ( on.j2+om.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-oj.j2),std::abs(oa.j2-ol.j2) )/2;
            int Jx_max = std::min( ob.j2+oj.j2 ,  oa.j2+ol.j2 )/2;
            for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
            {
              double sixj1 = Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p);
              if (std::abs(sixj1)<1e-6) continue;
              for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
              {
                double sixj2 = Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p);
                if (std::abs(sixj2)<1e-6) continue;
                int twoJJ_min = std::max( std::abs(ob.j2-2*J2p),std::abs(oa.j2-2*J1p) );
                int twoJJ_max = std::min( ob.j2+2*J2p ,  oa.j2+2*J1p );
                int phase = -Z.modelspace->phase((ok.j2+ol.j2)/2+J1+J2p);
                for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
                {
                  double xikanmb = X3.GetME_pn( J1p,J2p,twoJJ, i,k,a,n,m,b);
                  double yikanmb = Y3.GetME_pn( J1p,J2p,twoJJ, i,k,a,n,m,b);
                  double JJtot = 0.5 * twoJJ;
                  for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                  {
                    double xbjal = X2.GetTBME_J(Jx, b,j,a,l);
                    double ybjal = Y2.GetTBME_J(Jx, b,j,a,l);
                    double hats = (twoJJ+1)*(2*Jx+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1));
                    double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2p,  jj,J1p,Jtot, Jx,ja,jl);
                    z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbjal*yikanmb - ybjal*xikanmb);
                  }// for Jx
                }// for twoJ
              }// for J2p
            }// for J1p
           }// Z7 block



           // Z8  ~ Pik mln Z1     Xbiam Ykjalnb
           if ( (ob.l+oi.l+oa.l+om.l)%2==0   and (  (ob.tz2+oi.tz2)==(oa.tz2+om.tz2) )  and keep_kja and keep_lnb)
           {
            int J1p_min = std::abs(ok.j2-oj.j2)/2;
            int J1p_max = ( ok.j2+oj.j2 )/2;
            int J2p_min = std::abs(ol.j2-on.j2)/2;
            int J2p_max = ( ol.j2+on.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-oi.j2),std::abs(oa.j2-om.j2) )/2;
            int Jx_max = std::min( ob.j2+oi.j2 ,  oa.j2+om.j2 )/2;
            for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
            {
              double sixj1 = Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p);
              if (std::abs(sixj1)<1e-6) continue;
              for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
              {
                double sixj2 = Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p);
                if (std::abs(sixj2)<1e-6) continue;
                int twoJJ_min = std::max( std::abs(ob.j2-2*J2p),std::abs(oa.j2-2*J1p) );
                int twoJJ_max = std::min( ob.j2+2*J2p ,  oa.j2+2*J1p );
                int phase = -Z.modelspace->phase((oi.j2+on.j2)/2+J1p+J2);
                for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
                {
                  double xkjalnb = X3.GetME_pn( J1p,J2p,twoJJ, k,j,a,l,n,b);
                  double ykjalnb = Y3.GetME_pn( J1p,J2p,twoJJ, k,j,a,l,n,b);
                  double JJtot = 0.5 * twoJJ;
                  for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                  {
                    double xbiam = X2.GetTBME_J(Jx, b,i,a,m);
                    double ybiam = Y2.GetTBME_J(Jx, b,i,a,m);
                    double hats = (twoJJ+1)*(2*Jx+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1));
                    double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2p,  ji,J1p,Jtot, Jx,ja,jm);
                    z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbiam*ykjalnb - ybiam*xkjalnb);
                  }// for Jx
                }// for twoJ
              }// for J2p
            }// for J1p
           }// Z8 block


           // Z9  ~ Pjk mln Z1     Xbjam Yikalnb
           if ( (ob.l+oj.l+oa.l+om.l)%2==0   and (  (ob.tz2+oj.tz2)==(oa.tz2+om.tz2) )  and keep_ika and keep_lnb)
           {
            int J1p_min = std::abs(oi.j2-ok.j2)/2;
            int J1p_max = ( oi.j2+ok.j2 )/2;
            int J2p_min = std::abs(ol.j2-on.j2)/2;
            int J2p_max = ( ol.j2+on.j2 )/2;
            int Jx_min = std::max( std::abs(ob.j2-oj.j2),std::abs(oa.j2-om.j2) )/2;
            int Jx_max = std::min( ob.j2+oj.j2 ,  oa.j2+om.j2 )/2;
            for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
            {
              double sixj1 = Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p);
              if (std::abs(sixj1)<1e-6) continue;
              for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
              {
                double sixj2 = Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p);
                if (std::abs(sixj2)<1e-6) continue;
                int twoJJ_min = std::max( std::abs(ob.j2-2*J2p),std::abs(oa.j2-2*J1p) );
                int twoJJ_max = std::min( ob.j2+2*J2p ,  oa.j2+2*J1p );
                int phase = -Z.modelspace->phase((ok.j2+on.j2)/2+J1+J2);
                for (int twoJJ=twoJJ_min; twoJJ<=twoJJ_max; twoJJ+=2)
                {
                  double xikalnb = X3.GetME_pn( J1p,J2p,twoJJ, i,k,a,l,n,b);
                  double yikalnb = Y3.GetME_pn( J1p,J2p,twoJJ, i,k,a,l,n,b);
                  double JJtot = 0.5 * twoJJ;
                  for ( int Jx=Jx_min; Jx<=Jx_max; Jx++)
                  {
                    double xbjam = X2.GetTBME_J(Jx, b,j,a,m);
                    double ybjam = Y2.GetTBME_J(Jx, b,j,a,m);
                    double hats = (twoJJ+1)*(2*Jx+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1));
                    double ninej = Z.modelspace->GetNineJ(jb,JJtot,J2p,  jj,J1p,Jtot, Jx,ja,jm);
                    z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbjam*yikalnb - ybjam*xikalnb);
                  }// for Jx
                }// for twoJ
              }// for J2p
            }// for J1p
           }// Z9 block

         }// for b
        }// for a
        Z3.AddToME_pn_ch( ch3, ch3, ibra,iket, z_ijklmn);
      }// for iket
    }// for ibra
  }//for ch3

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}



//*****************************************************************************************
//
//  |    |    |     Uncoupled expression:
// i|   j|   k|       Z_ijklmn = 1/6 sum_{abc} (nanbnc + n`an`bn`c) ( X_{ijkabc} Y_{abclmn} - Y_{ijkabc} X_{abclmn}
//  *~~~[X]~~~*                                                                                                                     
// a|   b|   c| 
//  *~~~[Y]~~~*     Coupled expression:
// l|   m|   n|     Z_{ijklmn}^{J1,J2,J} = 1/6 sum_{abc} sum_{J'} (nanbnc + n`an`bn`c) ( X_{ijkabc}^{J1,J',J} Y_{abclmn}^{J',J2,J} - X<->Y )
//  |    |    |                                                                                                                       
//                                                           
//  Checked with UnitTest and passed
//
void comm333_ppp_hhhss( const Operator& X, const Operator& Y, Operator& Z ) 
{

  double tstart = omp_get_wtime();
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  int hX = X.IsHermitian() ? 1 : -1;
  int hY = Y.IsHermitian() ? 1 : -1;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();


  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();
  if ( not ( x3_allocated and y3_allocated) ) return;

  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
  #pragma omp parallel for schedule(dynamic,1)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    std::vector<size_t> abc_kets;
    std::vector<double> abc_factors;

    for (size_t iket_abc=0; iket_abc<nkets3; iket_abc++)
    {
       auto& ket_abc = Tbc.GetKet(iket_abc);
       double na = ket_abc.op->occ;
       double nb = ket_abc.oq->occ;
       double nc = ket_abc.oR->occ;
       double occ_factor = na*nb*nc - (1-na)*(1-nb)*(1-nc);
       if (std::abs(occ_factor)<1e-6) continue;
       double d_ea = std::abs( 2*ket_abc.op->n + ket_abc.op->l - e_fermi[ket_abc.op->tz2]);
       double d_eb = std::abs( 2*ket_abc.oq->n + ket_abc.oq->l - e_fermi[ket_abc.oq->tz2]);
       double d_ec = std::abs( 2*ket_abc.oR->n + ket_abc.oR->l - e_fermi[ket_abc.oR->tz2]);
       double occnat_a = ket_abc.op->occ_nat;
       double occnat_b = ket_abc.oq->occ_nat;
       double occnat_c = ket_abc.oR->occ_nat;
       if ( (d_ea + d_eb + d_ec) > Z.modelspace->GetdE3max() ) continue;
       if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
       double symm_factor = 1;
       if ( (ket_abc.p == ket_abc.q ) and (ket_abc.p == ket_abc.r) )
          symm_factor = 1./6;
       else if ( (ket_abc.p == ket_abc.q) or ( ket_abc.q == ket_abc.r) )
          symm_factor = 3./6;
       abc_kets.push_back(iket_abc);
//       abc_factors.push_back( sqrt( occ_factor * symm_factor) );
       abc_factors.push_back( occ_factor * symm_factor );
    }// for iket_abc

    std::vector<size_t> bras_kept;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      Ket3& bra = Tbc.GetKet(ibra);
      double d_ei = std::abs( 2*bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
      double d_ej = std::abs( 2*bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
      double d_ek = std::abs( 2*bra.oR->n + bra.oR->l - e_fermi[bra.oR->tz2]);
      double occnat_i = bra.op->occ_nat;
      double occnat_j = bra.oq->occ_nat;
      double occnat_k = bra.oR->occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->GetdE3max() ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      bras_kept.push_back( ibra );
    }

    size_t dim_abc = abc_kets.size();
    arma::mat XMAT(nkets3,dim_abc,arma::fill::zeros);
    arma::mat YMAT(nkets3,dim_abc,arma::fill::zeros);
    for (size_t index_abc=0; index_abc<dim_abc; index_abc++)
    {
      size_t iket_abc = abc_kets[index_abc];
      double factor = abc_factors[index_abc];
//      for (size_t ibra=0; ibra<nkets3; ibra++)
      for (size_t ibra : bras_kept )
      {
         double xijkabc = X3.GetME_pn_ch( ch3,ch3, ibra, iket_abc);
         double yijkabc = Y3.GetME_pn_ch( ch3,ch3, ibra, iket_abc);
         XMAT(ibra,index_abc) = factor * xijkabc ;
//         YMAT(ibra,index_abc) = factor * yijkabc ;
         YMAT(ibra,index_abc) =  yijkabc ;
      }// for ibra
    }// for index_abc


    // Now do the mat mult
    arma::mat ZMAT = hY * XMAT * YMAT.t() - hX*(YMAT * XMAT.t());

    // Store the results
//    for (size_t ibra=0; ibra<nkets3; ibra++)
    for (size_t ibra : bras_kept )
    {
//      for (size_t iket=ibra; iket<nkets3; iket++)
      for (size_t iket : bras_kept )
      {
        if ( iket < ibra ) continue;
        Z3.AddToME_pn_ch( ch3, ch3, ibra,iket, ZMAT(ibra,iket) );
        
      }// for iket
    }// for ibra
  }//for ch3

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}


//*****************************************************************************************
//
//  |i  |j      k/     Uncoupled expression:
//  |   |       /      Z_ijklmn =-1/2 P(ij/k)P(lm/n) sum_{abc} (nanbn`c + na`n`bnc) ( X_{ijcabn} Y_{abklmc} - Y_{ijcabn} X_{abklmc} )
//  *~~[X]~*   /                   
//  |   |  |\ /        Coupled expression:
// a|  b| c| /         Z_{ijklmn}^{J1,J2,J} = +1/2 P^{J1,J}(ij/k)P^{J2,j}_{lm/n) sum_{abc} (nanbn`c+n`an`bnc) sum_{J',J",J3} (2J'+1)(2J"+1)
//  |   |  |/ \                                  { k   J1  J  }
//  *~~[Y]~*   \                               * { J3  J'  n  } * ( X_{ijcabn}^{J1,J3,J'} Y_{abklmc}^{J3,J2,J"} - X<->Y )     
//  |   |       \                                { J"  c   J2 }                                                                            
//  |l  |m       \ n
//                    
//                    
//  Checked with UnitTest and passed
// TODO: Go back and implement the dE3max cut on the internal terms as well
// 
void comm333_pph_hhpss( const Operator& X, const Operator& Y, Operator& Z ) 
{

  double tstart = omp_get_wtime();
  double t_internal = omp_get_wtime();

  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  int hX = X.IsHermitian() ? 1 : -1;
  int hY = Y.IsHermitian() ? 1 : -1;

  bool x3_allocated = X3.IsAllocated();
  bool y3_allocated = Y3.IsAllocated();
  if ( not ( x3_allocated and y3_allocated) ) return;

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

  std::map<std::array<int,3>,size_t> channels_pph; // map twoJ_ph, parity_ph, twoTz_ph  -> channel index
  std::vector< std::array<int,3>> channel_list_pph; // vector  channel_index -> twoJ-ph, parity_ph, twoTz_ph.

//  std::vector< std::map<std::array<int,4>, size_t>> ket_lookup_pph; // in a given channel, map  i,j,k,Jij -> matrix index

  auto hash_key_ijnJ = [](size_t i,size_t j, size_t n, size_t Jij){return ( i + (j<<12) + (n<<24) + (Jij<<36) );};
  auto unhash_key_ijnJ = [](size_t& i,size_t& j, size_t& n, size_t& Jij, size_t key){ i=(key & 0xFFFL);j=((key>>12)&0xFFFL);n=((key>>24)&0xFFFL);Jij=((key>>36)&0xFFFL); }; //  0xF = 15 = 1111 (4bits), so 0xFFF is 12 bits of 1's. 0xFFFL makes it a long
  std::vector< std::unordered_map<size_t, size_t>> ket_lookup_pph; // in a given channel, map  i,j,n,Jij -> matrix index


  std::vector<std::vector<std::array<int,4>>> ket_lookup_abc;
  std::vector<std::vector<double>> occs_abc_list;

  std::vector< arma::mat > Zbar;

  // we're looking at pph type states. We don't make a cut on the occupations in this channel
  // but if the first two orbits ij can't possibly make it past the occnat cut later on, we don't bother including them
  // To figure out if they'll make it, we want to know the biggest contribution that can possibly be made by a third orbit.
  double occnat_factor_max = 0;
  for ( auto i : Z.modelspace->all_orbits)
  {
     double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
     occnat_factor_max = std::max( occnat_factor_max, occnat_i*(1-occnat_i) );
  }

  // Generate all the pph type channels J,parity,Tz. Note that these will be
  // not the same as the standard 3-body channels since we aren't enforcing antisymmetry with the 3rd particle
  int twoJph_min = 1;
  int twoJph_max = 6*Z.modelspace->GetEmax() + 3;
  int nch_pph = 0;
  for ( int twoJph=twoJph_min; twoJph<=twoJph_max; twoJph+=2)
  {
    for (int parity_ph=0; parity_ph<=1; parity_ph++)
    {
      for (int twoTz_ph=-3; twoTz_ph<=3; twoTz_ph+=2)
      {
//        std::vector<std::array<int,4>> good_kets_ijn;
        size_t ngood_ijn=0;

        std::unordered_map<size_t, size_t> good_kets_ijn;
//        std::map<std::array<int,4>, size_t> good_kets_ijn;

        std::vector<std::array<int,4>> good_kets_abc;
//        std::vector<size_t> kets_abc;
        std::vector<double> occs_abc;
        for ( size_t chij=0; chij<nch2; chij++)
        {
          TwoBodyChannel& tbc_ij = Z.modelspace->GetTwoBodyChannel(chij);
          int Jij = tbc_ij.J;
          int parity_ij = tbc_ij.parity;
          int Tz_ij = tbc_ij.Tz;
//          for (size_t n : Z.modelspace->all_orbits )
          for (int n : Z.modelspace->all_orbits )
          {
            Orbit& on = Z.modelspace->GetOrbit(n);
            if ( (on.l + parity_ij)%2 != parity_ph ) continue;
            if ( (2*Tz_ij - on.tz2) != twoTz_ph ) continue;  // note the minus sign because n is a hole
            if ( (std::abs(2*Jij-on.j2)> twoJph)  or  ((2*Jij+on.j2)<twoJph) ) continue;

            size_t nkets_ij = tbc_ij.GetNumberKets();
            for ( size_t iket_ij=0; iket_ij<nkets_ij; iket_ij++)
            {

              // check that we pass our cuts on E3max, occupations, etc.
              Ket& ket_ij = tbc_ij.GetKet(iket_ij);
//              size_t i = ket_ij.p;
//              size_t j = ket_ij.q;
              int i = ket_ij.p;
              int j = ket_ij.q;
              double occnat_i = ket_ij.op->occ_nat;
              double occnat_j = ket_ij.oq->occ_nat;
              double d_ei = std::abs( 2*ket_ij.op->n + ket_ij.op->l - e_fermi[ket_ij.op->tz2]);
              double d_ej = std::abs( 2*ket_ij.oq->n + ket_ij.oq->l - e_fermi[ket_ij.oq->tz2]);
              // if i and j cant make it past the OccNat and dE3max cuts, don't bother including it
              if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_factor_max ) < Z.modelspace->GetOccNat3Cut() ) continue;
              if ( (d_ei+d_ej) > Z.modelspace->dE3max ) continue;
              
              size_t key = hash_key_ijnJ( i,j,n,size_t(Jij));
              good_kets_ijn[ key ] = ngood_ijn;
//              good_kets_ijn[ {i,j,n,Jij} ] = ngood_ijn;
              ngood_ijn++;
              double occ_factor = ket_ij.op->occ * ket_ij.oq->occ * (1-on.occ) + (1-ket_ij.op->occ)*(1-ket_ij.oq->occ)*on.occ;

              if (i==j) occ_factor *=0.5; // because we only sum b<a
              if (std::abs(occ_factor)>1e-8 )
              {
                good_kets_abc.push_back( {i,j,n,Jij} );
                occs_abc.push_back( occ_factor );
              }
            }
          }
          
        }
        size_t ngood_abc = good_kets_abc.size();
        if ( (ngood_ijn < 1) or (ngood_abc<1)) continue;

        channels_pph[{twoJph, parity_ph,twoTz_ph}] = nch_pph;
        channel_list_pph.push_back({twoJph, parity_ph,twoTz_ph});
        ket_lookup_pph.push_back( good_kets_ijn );
        ket_lookup_abc.push_back( good_kets_abc );
        occs_abc_list.push_back( occs_abc );
        Zbar.push_back( arma::mat( ngood_ijn, ngood_ijn, arma::fill::zeros )  );
        nch_pph++;
      }// for twoTz_ph
    }// for parity_ph
  }// for twoJph
//  std::cout << "Nominally " << nch_pph << "  pph channels " << std::endl;


  Z.profiler.timer["comm333_pph_setup_lists"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

 // Now that we've made our lists and allocated the matrices, we fill the matrices inside the parallel block.
 #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
 for (size_t ch_pph=0; ch_pph<nch_pph; ch_pph++)
 {
   int twoJph =  channel_list_pph[ch_pph][0];
   double Jph = 0.5 * twoJph;
   auto& good_kets_ijn = ket_lookup_pph[ch_pph];
   auto& good_kets_abc = ket_lookup_abc[ch_pph];
   auto& occs_abc = occs_abc_list[ch_pph];
   size_t ngood_ijn = good_kets_ijn.size();
   size_t ngood_abc = good_kets_abc.size();
   arma::mat Xbar_ijnabc( ngood_ijn, ngood_abc, arma::fill::zeros );
   arma::mat Ybar_abcijn( ngood_abc, ngood_ijn, arma::fill::zeros );
   
   // now we fill Xbar an Ybar. This is where the heavy lifting is.
   for (auto& iter_ijn : good_kets_ijn )
   {
     size_t Iijn = iter_ijn.second;

     size_t key = iter_ijn.first;
     size_t i,j,n,Jij_tmp;
     unhash_key_ijnJ(i,j,n,Jij_tmp, key);
     int Jij = (int)Jij_tmp;
//     size_t i = iter_ijn.first[0];
//     size_t j = iter_ijn.first[1];
//     size_t n = iter_ijn.first[2];
//     int  Jij = iter_ijn.first[3];

     Orbit& oi = Z.modelspace->GetOrbit(i);
     Orbit& oj = Z.modelspace->GetOrbit(j);
     Orbit& on = Z.modelspace->GetOrbit(n);

     double occnat_i = oi.occ_nat;
     double occnat_j = oj.occ_nat;
     double occnat_n = on.occ_nat;
     double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
     double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
     double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
     double jn = 0.5*on.j2;
     for (size_t Iabc=0; Iabc<ngood_abc; Iabc++)
     {
       auto& abc_info = good_kets_abc[Iabc];
       size_t a = abc_info[0];
       size_t b = abc_info[1];
       size_t c = abc_info[2];
       double occ_factor = occs_abc[Iabc];
       Orbit& oa = Z.modelspace->GetOrbit(a);
       Orbit& ob = Z.modelspace->GetOrbit(b);
       Orbit& oc = Z.modelspace->GetOrbit(c);
       double jc = 0.5*oc.j2;
       int Jab  = abc_info[3];

       double occnat_a = oa.occ_nat;
       double occnat_b = ob.occ_nat;
       double occnat_c = oc.occ_nat;
       double d_ea = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
       double d_eb = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
       double d_ec = std::abs( 2*oc.n + oc.l - e_fermi[oc.tz2]);

       if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
       if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_c*(1-occnat_c) ) < Z.modelspace->GetOccNat3Cut() ) continue;
       if ( (d_ea+d_eb+d_en) > Z.modelspace->dE3max ) continue;
       if ( (d_ei+d_ej+d_ec) > Z.modelspace->dE3max ) continue;


       double xbar = 0;
       double ybar = 0;
       int twoJprime_min = std::max( std::abs(2*Jij - oc.j2), std::abs(2*Jab - on.j2) );
       int twoJprime_max = std::min( (2*Jij + oc.j2), (2*Jab + on.j2) );
       for (int twoJprime=twoJprime_min; twoJprime<=twoJprime_max; twoJprime+=2)
       {
         double sixj_ijn = Z.modelspace->GetSixJ( Jij, jn, Jph,  Jab, jc, 0.5*twoJprime);
        
         // TODO This comes up often enough that we should just make ThreeBodyMEpn do this.
         std::vector<double> xandy = Y3.GetME_pn_TwoOps( Jij,Jab,twoJprime, i,j,c, a,b,n, X3,Y3 );
         double xijcabn = xandy[0];
         double yijcabn = xandy[1];

         // this slightly verbose block is equivalent to
         // xbar += (twoJprime+1) * sixj_ijn * X3.GetME_pn( Jij, Jab, twoJprime, i,j,c, a,b,n);
         // ybar += (twoJprime+1) * sixj_ijn * Y3.GetME_pn( Jij, Jab, twoJprime, i,j,c, a,b,n) * occ_factor;
         // but it avoids doing the recoupling twice (if i,j,c and a,b,n aren't in the storage order)

//         std::vector<size_t> ibra_list, iket_list;
//         std::vector<double> recouple_bra_list, recouple_ket_list;
//         size_t ch_check_bra = Y.ThreeBody.GetKetIndex_withRecoupling( Jij, twoJprime, i, j, c, ibra_list, recouple_bra_list) ;
//         size_t ch_check_ket = Y.ThreeBody.GetKetIndex_withRecoupling( Jab, twoJprime, a, b, n, iket_list, recouple_ket_list) ;
//         double xijcabn = 0;
//         double yijcabn = 0;
//         for ( size_t ind_bra=0; ind_bra<ibra_list.size(); ind_bra++)
//         {
//           double rec_bra = recouple_bra_list[ind_bra];
//           size_t Ibra = ibra_list[ind_bra];
//           for ( size_t ind_ket=0; ind_ket<iket_list.size(); ind_ket++)
//           {
//             double rec_ket = recouple_ket_list[ind_ket];
//             size_t Iket = iket_list[ind_ket];
//              xijcabn += rec_bra * rec_ket * X3.GetME_pn_ch(ch_check_bra,ch_check_ket,Ibra,Iket);
//              yijcabn += rec_bra * rec_ket * Y3.GetME_pn_ch(ch_check_bra,ch_check_ket,Ibra,Iket);
//           }
//         }
         xbar += (twoJprime+1) * sixj_ijn * xijcabn;
         ybar += (twoJprime+1) * sixj_ijn * yijcabn ;
 
       }// for twoJprime
       
       Xbar_ijnabc( Iijn, Iabc ) = xbar;
       Ybar_abcijn( Iabc, Iijn ) = hY * ybar * occ_factor;
     }// for Iabc
   }// for iter_ijn

   auto& zbar = Zbar[ch_pph];
   zbar =  Xbar_ijnabc * Ybar_abcijn;
   zbar -= hX *hY * zbar.t();

  }// for ch_pph , end of parallel block

  Z.profiler.timer["comm333_pph_fill_matrices"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // Now that we've computed Zbar (i.e. the pph transformed commutator), we need to
  // transform that back to Z, and do all the permutations to ensure antisymmetry
  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      int J1 = bra.Jpq;

      // Here are the permutations for ijk
      std::vector<std::array<size_t,3>> ijk = { {i,j,k}, {k,j,i}, {i,k,j} };

      std::vector<int> J1p_min = {J1,  std::max(std::abs(ok.j2-oj.j2),std::abs(twoJ-oi.j2) )/2,   std::max(std::abs(oi.j2-ok.j2), std::abs(twoJ-oj.j2) )/2 };
      std::vector<int> J1p_max = {J1,  std::min(ok.j2+oj.j2, twoJ+oi.j2)/2 , std::min(oi.j2+ok.j2, twoJ+oj.j2)/2 };
//      std::vector<int> J1p_min = {J1,  std::abs(ok.j2-oj.j2)/2,   std::abs(oi.j2-ok.j2)/2 };
//      std::vector<int> J1p_max = {J1,  (ok.j2+oj.j2)/2 , (oi.j2+ok.j2)/2 };
      std::vector<std::vector<double>> recouple_ijk = {{1},{},{} };
      for (int J1p=J1p_min[1]; J1p<=J1p_max[1]; J1p++)
           recouple_ijk[1].push_back( sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p) );

      for (int J1p=J1p_min[2]; J1p<=J1p_max[2]; J1p++)
           recouple_ijk[2].push_back( -Z.modelspace->phase((oj.j2+ok.j2)/2+J1+J1p)*sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p) );


      for (size_t iket=0; iket<=ibra; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        int J2 = ket.Jpq;

        // permutations for lmn
        std::vector<std::array<size_t,3>> lmn = { {l,m,n}, {n,m,l}, {l,n,m} };

        std::vector<int> J2p_min = {J2,  std::max(std::abs(on.j2-om.j2),std::abs(twoJ-ol.j2) )/2,   std::max(std::abs(ol.j2-on.j2), std::abs(twoJ-om.j2) )/2 };
        std::vector<int> J2p_max = {J2,  std::min(on.j2+om.j2, twoJ+ol.j2)/2 , std::min(ol.j2+on.j2, twoJ+om.j2)/2 };
//        std::vector<int> J2p_min = {J2,  std::abs(on.j2-om.j2)/2,   std::abs(ol.j2-on.j2)/2 };
//        std::vector<int> J2p_max = {J2,  (on.j2+om.j2)/2 , (ol.j2+on.j2)/2 };
        std::vector<std::vector<double>> recouple_lmn = {{1},{},{} };
            
        for (int J2p=J2p_min[1]; J2p<=J2p_max[1]; J2p++)
           recouple_lmn[1].push_back( sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p) );

        for (int J2p=J2p_min[2]; J2p<=J2p_max[2]; J2p++)
             recouple_lmn[2].push_back( -Z.modelspace->phase((om.j2+on.j2)/2+J2+J2p)*sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p) );


        double z_ijklmn = 0;

        for ( int perm_ijk=0; perm_ijk<3; perm_ijk++ )
        {
          size_t I1 = ijk[perm_ijk][0];
          size_t I2 = ijk[perm_ijk][1];
          size_t I3 = ijk[perm_ijk][2];
          Orbit& o1 = Z.modelspace->GetOrbit( I1 );
          Orbit& o2 = Z.modelspace->GetOrbit( I2 );
          Orbit& o3 = Z.modelspace->GetOrbit( I3 );
          double j3 = 0.5*o3.j2;
          for (int J1p=J1p_min[perm_ijk]; J1p<=J1p_max[perm_ijk]; J1p++)
          {
            if ( (I1==I2) and J1p%2>0) continue;
            double rec_ijk = recouple_ijk[perm_ijk].at(J1p-J1p_min[perm_ijk]);
            for ( int perm_lmn=0; perm_lmn<3; perm_lmn++ )
            {
              size_t I4 = lmn[perm_lmn][0];
              size_t I5 = lmn[perm_lmn][1];
              size_t I6 = lmn[perm_lmn][2];
              Orbit& o4 = Z.modelspace->GetOrbit( I4 );
              Orbit& o5 = Z.modelspace->GetOrbit( I5 );
              Orbit& o6 = Z.modelspace->GetOrbit( I6 );
              double j6 = 0.5*o6.j2;

              int twoTz_ph = (o1.tz2 + o2.tz2 - o6.tz2);
              int parity_ph = (o1.l + o2.l + o6.l)%2;

              for (int J2p=J2p_min[perm_lmn]; J2p<=J2p_max[perm_lmn]; J2p++)
              {
                double rec_lmn = recouple_lmn[perm_lmn].at(J2p-J2p_min[perm_lmn]);
                if ( (I4==I5) and J2p%2>0) continue;

                size_t key_126 = hash_key_ijnJ(std::min(I1,I2),std::max(I1,I2),I6,J1p);
                size_t key_453 = hash_key_ijnJ(std::min(I4,I5),std::max(I4,I5),I3,J2p);

                int twoJph_min = std::max( std::abs(2*J1p - o6.j2), std::abs(2*J2p-o3.j2) );
                int twoJph_max = std::min( (2*J1p + o6.j2), (2*J2p+o3.j2) );
                if ( twoJph_min > twoJph_max) continue;

//                std::vector<double> Zbar_126453 ( (twoJph_max-twoJph_min)/2+1, 0. );
                int phase_12 = I1>I2 ? -Z.modelspace->phase( (o1.j2+o2.j2)/2 - J1p) : 1;
                int phase_45 = I4>I5 ? -Z.modelspace->phase( (o4.j2+o5.j2)/2 - J2p) : 1;

                double z_123456 =0;
                for (int twoJph = twoJph_min; twoJph<=twoJph_max; twoJph+=2)
                {
                  auto iter_ch = channels_pph.find({ twoJph, parity_ph, twoTz_ph});
                  if ( iter_ch == channels_pph.end() ) continue;
                  size_t ch_pph = iter_ch->second;

                   size_t index_126 = ket_lookup_pph[ch_pph].at(key_126 );
                   size_t index_453 = ket_lookup_pph[ch_pph].at(key_453 );
//                  size_t index_126 = ket_lookup_pph[ch_pph].at({ std::min(I1,I2),std::max(I1,I2),I6,J1p} );
//                  size_t index_453 = ket_lookup_pph[ch_pph].at({ std::min(I4,I5),std::max(I4,I5),I3,J2p} );
                  double zbar =  phase_12 * phase_45 * Zbar[ch_pph](index_126,index_453);

                  double Jph = 0.5 * twoJph;
                  double sixj_ph = Z.modelspace->GetSixJ( J1p, j3, Jtot,  J2p, j6, Jph );
                  z_123456 += (twoJph+1) * sixj_ph  * zbar ;
                }// for Jph 
                z_ijklmn +=  rec_ijk * rec_lmn * z_123456;

                 }// for J2p
               }// for perm_lmn
             }// for J1p
           }// for perm_ijk

        Z3.AddToME_pn_ch( ch3, ch3, ibra, iket, z_ijklmn);
      }// for iket
    }// for ibra
  }//for ch3

  Z.profiler.timer["comm333_pph_recouple"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}

/*
//
// The old slow way, just to make sure we didn't screw anything up by going faster
//
void comm333_pph_hhpss_debug( const Operator& X, const Operator& Y, Operator& Z ) 
{

  double tstart = omp_get_wtime();

  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      int J1 = bra.Jpq;



      std::vector<std::array<size_t,3>> ijk = { {i,j,k}, {k,j,i}, {i,k,j} };
      std::vector<int> J1p_min = {J1,  std::abs(ok.j2-oj.j2)/2,   std::abs(oi.j2-ok.j2)/2 };
      std::vector<int> J1p_max = {J1,  (ok.j2+oj.j2)/2 , (oi.j2+ok.j2)/2 };
      std::vector<std::vector<double>> recouple_ijk = {{1},{},{} };
      for (int J1p=J1p_min[1]; J1p<=J1p_max[1]; J1p++)
           recouple_ijk[1].push_back( sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p) );

      for (int J1p=J1p_min[2]; J1p<=J1p_max[2]; J1p++)
           recouple_ijk[2].push_back( -Z.modelspace->phase((oj.j2+ok.j2)/2+J1+J1p)*sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p) );


//      for (size_t iket=ibra; iket<nkets3; iket++)
      for (size_t iket=0; iket<=ibra; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        int J2 = ket.Jpq;


        std::vector<std::array<size_t,3>> lmn = { {l,m,n}, {n,m,l}, {l,n,m} };
        std::vector<int> J2p_min = {J2,  std::abs(on.j2-om.j2)/2,   std::abs(ol.j2-on.j2)/2 };
        std::vector<int> J2p_max = {J2,  (on.j2+om.j2)/2 , (ol.j2+on.j2)/2 };
        std::vector<std::vector<double>> recouple_lmn = {{1},{},{} };
            
        for (int J2p=J2p_min[1]; J2p<=J2p_max[1]; J2p++)
           recouple_lmn[1].push_back( sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p) );

        for (int J2p=J2p_min[2]; J2p<=J2p_max[2]; J2p++)
             recouple_lmn[2].push_back( -Z.modelspace->phase((om.j2+on.j2)/2+J2+J2p)*sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p) );


        double z_ijklmn = 0;


        for ( int perm_ijk=0; perm_ijk<3; perm_ijk++ )
        {
          size_t I1 = ijk[perm_ijk][0];
          size_t I2 = ijk[perm_ijk][1];
          size_t I3 = ijk[perm_ijk][2];
          Orbit& o1 = Z.modelspace->GetOrbit( I1 );
          Orbit& o2 = Z.modelspace->GetOrbit( I2 );
          Orbit& o3 = Z.modelspace->GetOrbit( I3 );
          double occnat_1 = o1.occ_nat;
          double occnat_2 = o2.occ_nat;
          double occnat_3 = o3.occ_nat;
          double d_e1 = std::abs( 2*o1.n + o1.l - e_fermi[o1.tz2]);
          double d_e2 = std::abs( 2*o2.n + o2.l - e_fermi[o2.tz2]);
          double d_e3 = std::abs( 2*o3.n + o3.l - e_fermi[o3.tz2]);
          double j3 = 0.5*o3.j2;
          for (int J1p=J1p_min[perm_ijk]; J1p<=J1p_max[perm_ijk]; J1p++)
          {
            double rec_ijk = recouple_ijk[perm_ijk].at(J1p-J1p_min[perm_ijk]);
            for ( int perm_lmn=0; perm_lmn<3; perm_lmn++ )
            {
              size_t I4 = lmn[perm_lmn][0];
              size_t I5 = lmn[perm_lmn][1];
              size_t I6 = lmn[perm_lmn][2];
              Orbit& o4 = Z.modelspace->GetOrbit( I4 );
              Orbit& o5 = Z.modelspace->GetOrbit( I5 );
              Orbit& o6 = Z.modelspace->GetOrbit( I6 );
              double occnat_4 = o4.occ_nat;
              double occnat_5 = o5.occ_nat;
              double occnat_6 = o6.occ_nat;
              double d_e4 = std::abs( 2*o4.n + o4.l - e_fermi[o4.tz2]);
              double d_e5 = std::abs( 2*o5.n + o5.l - e_fermi[o5.tz2]);
              double d_e6 = std::abs( 2*o6.n + o6.l - e_fermi[o6.tz2]);
              double j6 = 0.5*o6.j2;
              for (int J2p=J2p_min[perm_lmn]; J2p<=J2p_max[perm_lmn]; J2p++)
              {
                double rec_lmn = recouple_lmn[perm_lmn].at(J2p-J2p_min[perm_lmn]);

                for (size_t ch2=0; ch2<nch2; ch2++)
                {
                  auto& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch2);
                  if ( std::abs(Tbc.twoTz-2*tbc_ab.Tz)==5) continue; // TODO there are probably other checks at the channel level...
                  size_t nkets_ab = tbc_ab.GetNumberKets();
                  int Jab = tbc_ab.J;

                  for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
                  {
                    Ket& ket_ab = tbc_ab.GetKet(iket_ab);
                    size_t a = ket_ab.p;
                    size_t b = ket_ab.q;
                    Orbit& oa = Z.modelspace->GetOrbit(a);
                    Orbit& ob = Z.modelspace->GetOrbit(b);
                    if (std::abs(oa.occ * ob.occ)<1e-6 and std::abs( (1-oa.occ)*(1-ob.occ))<1e-6) continue;

                    double occnat_a = oa.occ_nat;
                    double occnat_b = ob.occ_nat;
                    double d_ea = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
                    double d_eb = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
                    bool keep_ab3 = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_3*(1-occnat_3) ) >= Z.modelspace->GetOccNat3Cut(); 
                    bool keep_ab6 = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_6*(1-occnat_6) ) >= Z.modelspace->GetOccNat3Cut(); 
                    keep_ab3 = keep_ab3 and (d_ea+d_eb+d_e3 <= Z.modelspace->dE3max );
                    keep_ab6 = keep_ab6 and (d_ea+d_eb+d_e6 <= Z.modelspace->dE3max );

                    if ( not ( keep_ab3 or keep_ab6) ) continue;
                    
                    for (auto c : Z.modelspace->all_orbits)
                    {
                      Orbit& oc = Z.modelspace->GetOrbit(c);
                      double occ_factor = oa.occ * ob.occ * (1-oc.occ) + (1-oa.occ)*(1-ob.occ)*oc.occ;
                      if (std::abs(occ_factor)<1e-6) continue;
                      if (a==b) occ_factor *=0.5; // because we only sum b<a
                      double occnat_c = oc.occ_nat;
                      double d_ec = std::abs( 2*oc.n + oc.l - e_fermi[oc.tz2]);
                      double jc = 0.5 * oc.j2;
                      bool keep_12c = (occnat_1*(1-occnat_1) * occnat_2*(1-occnat_2) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
                      bool keep_45c = (occnat_4*(1-occnat_4) * occnat_5*(1-occnat_5) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
                      keep_12c = keep_12c and (d_e1+d_e2+d_ec <= Z.modelspace->dE3max );
                      keep_45c = keep_45c and (d_e4+d_e5+d_ec <= Z.modelspace->dE3max );


                      if ( ((oa.l+ob.l+o3.l+o4.l+o5.l+oc.l)%2==0) and ((o1.l+o2.l+oc.l+oa.l+ob.l+o6.l)%2==0)
                       and  ( (oa.tz2+ob.tz2+o3.tz2)==(o4.tz2+o5.tz2+oc.tz2) ) and ( (o1.tz2+o2.tz2+oc.tz2)==(oa.tz2+ob.tz2+o6.tz2))
                        and keep_ab3 and keep_45c and keep_12c and keep_ab6)
//                        )//and keep_ab3 and keep_45c and keep_12c and keep_ab6)
                      {

                        int twoJx_min = std::max( std::abs(2*Jab - o3.j2), std::abs(2*J2p - oc.j2) );
                        int twoJx_max = std::min( 2*Jab+o3.j2 , 2*J2p + oc.j2 );
                        int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - o6.j2) );
                        int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + o6.j2 );
                        if (twoJx_min <= twoJx_max and twoJy_min<=twoJy_max) 
                        {
                          std::vector<double> xab345c( (twoJx_max-twoJx_min)/2+1, 0);
                          std::vector<double> yab345c( (twoJx_max-twoJx_min)/2+1, 0);
                          std::vector<double> x12cab6( (twoJy_max-twoJy_min)/2+1, 0);
                          std::vector<double> y12cab6( (twoJy_max-twoJy_min)/2+1, 0);
                          for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                          {
                            size_t iJx = (twoJx-twoJx_min)/2;
                            xab345c[iJx] = X3.GetME_pn(Jab,J2p,twoJx, a,b,I3,I4,I5,c);
                            yab345c[iJx] = Y3.GetME_pn(Jab,J2p,twoJx, a,b,I3,I4,I5,c);
                          }
                          for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                          {
                             size_t iJy = (twoJy-twoJy_min)/2;
                             x12cab6[iJy] = X3.GetME_pn(J1p,Jab,twoJy, I1,I2,c,a,b,I6);
                             y12cab6[iJy] = Y3.GetME_pn(J1p,Jab,twoJy, I1,I2,c,a,b,I6);
                          }
                          for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                          {
                            double JJx = 0.5 * twoJx;
                            size_t iJx = (twoJx-twoJx_min)/2;
                            for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                            {
                               double JJy = 0.5 * twoJy;
                               size_t iJy = (twoJy-twoJy_min)/2;
                               double hats = (twoJx+1)*(twoJy+1);
                               double ninej = Z.modelspace->GetNineJ( j3,Jab,JJx, J1p,JJy,jc, Jtot,j6,J2p);
                                z_ijklmn +=  rec_ijk * rec_lmn * occ_factor * hats * ninej * ( xab345c[iJx]*y12cab6[iJy] - yab345c[iJx]*x12cab6[iJy] );
                            }// for twoJy
                          }// for twoJx
                        }
                      }// Z1 block

                     }// for c
                   }// for iket_ab
                 }// for ch2


                    }// for J2p
                  }// for perm_lmn
                }// for J1p
              }// for perm_ijk

        Z3.AddToME_pn_ch( ch3, ch3, ibra, iket, z_ijklmn);
      }// for iket
    }// for ibra
  }//for ch3

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
*/




//void comm333_pph_hhpss( const Operator& X, const Operator& Y, Operator& Z ) 
void comm333_pph_hhpss_debug( const Operator& X, const Operator& Y, Operator& Z ) 
{

  double tstart = omp_get_wtime();

  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

  size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
  for (size_t ch3=0; ch3<nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    size_t nkets3 = Tbc.GetNumberKets();
    int twoJ = Tbc.twoJ;
    double Jtot = 0.5 * twoJ;
    for (size_t ibra=0; ibra<nkets3; ibra++)
    {
      auto& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;
      if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
      if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
      int J1 = bra.Jpq;
//      for (size_t iket=ibra; iket<nkets3; iket++)
      for (size_t iket=0; iket<=ibra; iket++)
      {
        auto& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
        if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        int J2 = ket.Jpq;

//              if ( not ((i==2 and j==4 and k==5 and l==4 and m==4 and n==5) 
//              if ( not ((i==4 and j==4 and k==5 and l==2 and m==4 and n==5) 
//               or (i==5 and j==4 and k==2 and l==5 and m==4 and n==4)) ) continue;
        double z_ijklmn = 0;

        for (size_t ch2=0; ch2<nch2; ch2++)
        {
          auto& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch2);
          if ( std::abs(Tbc.twoTz-2*tbc_ab.Tz)==5) continue; // TODO there are probably other checks at the channel level...
          size_t nkets_ab = tbc_ab.GetNumberKets();
          int Jab = tbc_ab.J;

          for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
          {
            Ket& ket_ab = tbc_ab.GetKet(iket_ab);
            size_t a = ket_ab.p;
            size_t b = ket_ab.q;
            Orbit& oa = Z.modelspace->GetOrbit(a);
            Orbit& ob = Z.modelspace->GetOrbit(b);
            if (std::abs(oa.occ * ob.occ)<1e-6 and std::abs( (1-oa.occ)*(1-ob.occ))<1e-6) continue;

            double occnat_a = oa.occ_nat;
            double occnat_b = ob.occ_nat;
            double d_ea = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
            double d_eb = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
            bool keep_abi = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_i*(1-occnat_i) ) >= Z.modelspace->GetOccNat3Cut(); 
            bool keep_abj = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_j*(1-occnat_j) ) >= Z.modelspace->GetOccNat3Cut(); 
            bool keep_abk = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_k*(1-occnat_k) ) >= Z.modelspace->GetOccNat3Cut(); 
            bool keep_abl = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_l*(1-occnat_l) ) >= Z.modelspace->GetOccNat3Cut(); 
            bool keep_abm = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_m*(1-occnat_m) ) >= Z.modelspace->GetOccNat3Cut(); 
            bool keep_abn = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_n*(1-occnat_n) ) >= Z.modelspace->GetOccNat3Cut(); 
            keep_abi = keep_abi and (d_ea+d_eb+d_ei <= Z.modelspace->dE3max );
            keep_abj = keep_abj and (d_ea+d_eb+d_ej <= Z.modelspace->dE3max );
            keep_abk = keep_abk and (d_ea+d_eb+d_ek <= Z.modelspace->dE3max );
            keep_abl = keep_abl and (d_ea+d_eb+d_el <= Z.modelspace->dE3max );
            keep_abm = keep_abm and (d_ea+d_eb+d_em <= Z.modelspace->dE3max );
            keep_abn = keep_abn and (d_ea+d_eb+d_en <= Z.modelspace->dE3max );
            
            for (auto c : Z.modelspace->all_orbits)
            {
              Orbit& oc = Z.modelspace->GetOrbit(c);
              double occ_factor = oa.occ * ob.occ * (1-oc.occ) + (1-oa.occ)*(1-ob.occ)*oc.occ;
              if (std::abs(occ_factor)<1e-6) continue;
              if (a==b) occ_factor *=0.5; // because we only sum b<a
              double occnat_c = oc.occ_nat;
              double d_ec = std::abs( 2*oc.n + oc.l - e_fermi[oc.tz2]);
              double jc = 0.5 * oc.j2;
              bool keep_ijc = (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
              bool keep_kjc = (occnat_k*(1-occnat_k) * occnat_j*(1-occnat_j) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
              bool keep_ikc = (occnat_i*(1-occnat_i) * occnat_k*(1-occnat_k) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
              bool keep_lmc = (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
              bool keep_nmc = (occnat_n*(1-occnat_n) * occnat_m*(1-occnat_m) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
              bool keep_lnc = (occnat_l*(1-occnat_l) * occnat_n*(1-occnat_n) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut(); 
              keep_ijc = keep_ijc and (d_ei+d_ej+d_ec <= Z.modelspace->dE3max );
              keep_kjc = keep_kjc and (d_ek+d_ej+d_ec <= Z.modelspace->dE3max );
              keep_ikc = keep_ikc and (d_ei+d_ek+d_ec <= Z.modelspace->dE3max );
              keep_lmc = keep_lmc and (d_el+d_em+d_ec <= Z.modelspace->dE3max );
              keep_nmc = keep_nmc and (d_en+d_em+d_ec <= Z.modelspace->dE3max );
              keep_lnc = keep_lnc and (d_el+d_en+d_ec <= Z.modelspace->dE3max );


              // Direct Z1 term
              if ( ((oa.l+ob.l+ok.l+ol.l+om.l+oc.l)%2==0) and ((oi.l+oj.l+oc.l+oa.l+ob.l+on.l)%2==0)
               and  ( (oa.tz2+ob.tz2+ok.tz2)==(ol.tz2+om.tz2+oc.tz2) ) and ( (oi.tz2+oj.tz2+oc.tz2)==(oa.tz2+ob.tz2+on.tz2))
                and keep_abk and keep_lmc and keep_ijc and keep_abn)
              {
                int twoJx_min = std::max( std::abs(2*Jab - ok.j2), std::abs(2*J2 - oc.j2) );
                int twoJx_max = std::min( 2*Jab+ok.j2 , 2*J2 + oc.j2 );
                int twoJy_min = std::max( std::abs(2*J1 - oc.j2), std::abs(2*Jab - on.j2) );
                int twoJy_max = std::min( 2*J1+oc.j2 , 2*Jab + on.j2 );
                if (twoJx_min <= twoJx_max and twoJy_min<=twoJy_max) 
                {
                  std::vector<double> xabklmc( (twoJx_max-twoJx_min)/2+1, 0);
                  std::vector<double> yabklmc( (twoJx_max-twoJx_min)/2+1, 0);
                  std::vector<double> xijcabn( (twoJy_max-twoJy_min)/2+1, 0);
                  std::vector<double> yijcabn( (twoJy_max-twoJy_min)/2+1, 0);
                  for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                  {
  //                  double JJx = 0.5 * twoJx;
                    size_t iJx = (twoJx-twoJx_min)/2;
                    xabklmc[iJx] = X3.GetME_pn(Jab,J2,twoJx, a,b,k,l,m,c);
                    yabklmc[iJx] = Y3.GetME_pn(Jab,J2,twoJx, a,b,k,l,m,c);
  //                  double xabklmc = X3.GetME_pn(Jab,J2,twoJx, a,b,k,l,m,c);
  //                  double yabklmc = Y3.GetME_pn(Jab,J2,twoJx, a,b,k,l,m,c);
                  }
                  for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                  {
                     size_t iJy = (twoJy-twoJy_min)/2;
                     xijcabn[iJy] = X3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,n);
                     yijcabn[iJy] = Y3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,n);
                  }
                  for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                  {
                    double JJx = 0.5 * twoJx;
                    size_t iJx = (twoJx-twoJx_min)/2;
                    for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                    {
                       double JJy = 0.5 * twoJy;
                       size_t iJy = (twoJy-twoJy_min)/2;
                       double hats = (twoJx+1)*(twoJy+1);
                       double ninej = Z.modelspace->GetNineJ( jk,Jab,JJx, J1,JJy,jc, Jtot,jn,J2);
  //                     double xijcabn = X3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,n);
  //                     double yijcabn = Y3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,n);
  //                     z_ijklmn +=  occ_factor * hats * ninej * ( xabklmc*yijcabn - yabklmc*xijcabn );
                       z_ijklmn +=  occ_factor * hats * ninej * ( xabklmc[iJx]*yijcabn[iJy] - yabklmc[iJx]*xijcabn[iJy] );
                    }// for twoJy
                  }// for twoJx
                }
              }// Z1 block
//              std::cout << "  after Z1  " << z_ijklmn << std::endl;


              // Z2 term  Pik    Xabilmc Ykjcabn
              if ( ((oa.l+ob.l+oi.l+ol.l+om.l+oc.l)%2==0) and ((ok.l+oj.l+oc.l+oa.l+ob.l+on.l)%2==0)
               and  ( (oa.tz2+ob.tz2+oi.tz2)==(ol.tz2+om.tz2+oc.tz2) ) and ( (ok.tz2+oj.tz2+oc.tz2)==(oa.tz2+ob.tz2+on.tz2)) 
                and keep_abi and keep_lmc and keep_kjc and keep_abn )
              {
                int twoJx_min = std::max( std::abs(2*Jab - oi.j2), std::abs(2*J2 - oc.j2) );
                int twoJx_max = std::min( 2*Jab+oi.j2 , 2*J2 + oc.j2 );
                int J1p_min = std::abs( ok.j2 - oj.j2 )/2;
                int J1p_max = ( ok.j2 + oj.j2 )/2;
                for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                {
                  double JJx = 0.5 * twoJx;
                  double xabilmc = X3.GetME_pn(Jab,J2,twoJx, a,b,i,l,m,c);
                  double yabilmc = Y3.GetME_pn(Jab,J2,twoJx, a,b,i,l,m,c);

                  for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
                  {
                    int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - on.j2) );
                    int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + on.j2 );
                    double sixj = Z.modelspace->GetSixJ(ji,jj,J1, jk,Jtot,J1p);
                    if (std::abs(sixj)<1e-6) continue;
                    int phase = -1;
                    for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                    {
                       double JJy = 0.5 * twoJy;
                       double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J1+1)*(2*J1p+1));
                       double ninej = Z.modelspace->GetNineJ( ji,Jab,JJx, J1p,JJy,jc, Jtot,jn,J2);
                       double xkjcabn = X3.GetME_pn(J1p,Jab,twoJy, k,j,c,a,b,n);
                       double ykjcabn = Y3.GetME_pn(J1p,Jab,twoJy, k,j,c,a,b,n);
                       z_ijklmn -=  occ_factor * phase * hats * sixj * ninej * ( xabilmc*ykjcabn - yabilmc*xkjcabn );
                    }// for twoJy
                  }// for J1p
                }// for twoJx
              }// Z2 block
//              std::cout << "  after Z2  " << z_ijklmn << std::endl;


              // Z3 term  Pjk    Xabjlmc Yikcabn
              if ( ((oa.l+ob.l+oj.l+ol.l+om.l+oc.l)%2==0) and ((oi.l+ok.l+oc.l+oa.l+ob.l+on.l)%2==0)
               and  ( (oa.tz2+ob.tz2+oj.tz2)==(ol.tz2+om.tz2+oc.tz2) ) and ( (oi.tz2+ok.tz2+oc.tz2)==(oa.tz2+ob.tz2+on.tz2)) 
                and keep_abj and keep_lmc and keep_ikc and keep_abn )
              {
                int twoJx_min = std::max( std::abs(2*Jab - oj.j2), std::abs(2*J2 - oc.j2) );
                int twoJx_max = std::min( 2*Jab+oj.j2 , 2*J2 + oc.j2 );
                int J1p_min = std::abs( oi.j2 - ok.j2 )/2;
                int J1p_max = ( oi.j2 + ok.j2 )/2;
                for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                {
                  double JJx = 0.5 * twoJx;
                  double xabjlmc = X3.GetME_pn(Jab,J2,twoJx, a,b,j,l,m,c);
                  double yabjlmc = Y3.GetME_pn(Jab,J2,twoJx, a,b,j,l,m,c);

                  for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
                  {
                    int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - on.j2) );
                    int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + on.j2 );
                    double sixj = Z.modelspace->GetSixJ(jj,ji,J1, jk,Jtot,J1p);
                    if (std::abs(sixj)<1e-6) continue;
                    int phase = Z.modelspace->phase( (oj.j2+ok.j2)/2+J1+J1p);
                    for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                    {
                       double JJy = 0.5 * twoJy;
                       double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J1+1)*(2*J1p+1));
                       double ninej = Z.modelspace->GetNineJ( jj,Jab,JJx, J1p,JJy,jc, Jtot,jn,J2);
                       double xikcabn = X3.GetME_pn(J1p,Jab,twoJy, i,k,c,a,b,n);
                       double yikcabn = Y3.GetME_pn(J1p,Jab,twoJy, i,k,c,a,b,n);
                       z_ijklmn -=  occ_factor * phase * hats * sixj * ninej * ( xabjlmc*yikcabn - yabjlmc*xikcabn );
                    }// for twoJy
                  }// for J1p
                }// for twoJx
              }// Z3 block
//              std::cout << "  after Z3  " << z_ijklmn << std::endl;


              // Z4 term  Pln    Xabknmc Yijcabl
              if ( ((oa.l+ob.l+ok.l+on.l+om.l+oc.l)%2==0) and ((oi.l+oj.l+oc.l+oa.l+ob.l+ol.l)%2==0)
               and  ( (oa.tz2+ob.tz2+ok.tz2)==(on.tz2+om.tz2+oc.tz2) ) and ( (oi.tz2+oj.tz2+oc.tz2)==(oa.tz2+ob.tz2+ol.tz2)) 
                and keep_abk and keep_nmc and keep_ijc and keep_abl )
              {
                int twoJy_min = std::max( std::abs(2*J1 - oc.j2), std::abs(2*Jab - ol.j2) );
                int twoJy_max = std::min( 2*J1+oc.j2 , 2*Jab + ol.j2 );
                int J2p_min = std::abs( on.j2 - om.j2 )/2;
                int J2p_max = ( on.j2 + om.j2 )/2;
                for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                {
                  double xijcabl = X3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,l);
                  double yijcabl = Y3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,l);
                  double JJy = 0.5 * twoJy;
                  for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
                  {
                    double sixj = Z.modelspace->GetSixJ(jl,jm,J2, jn,Jtot,J2p);
                    int phase = -1;
                    if (std::abs(sixj)<1e-6) continue;
                    int twoJx_min = std::max( std::abs(2*Jab - ok.j2), std::abs(2*J2p - oc.j2) );
                    int twoJx_max = std::min( 2*Jab+ok.j2 , 2*J2p + oc.j2 );
                    for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                    {
                      double JJx = 0.5 * twoJx;
                      double xabknmc = X3.GetME_pn(Jab,J2p,twoJx, a,b,k,n,m,c);
                      double yabknmc = Y3.GetME_pn(Jab,J2p,twoJx, a,b,k,n,m,c);
                      double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J2+1)*(2*J2p+1));
                      double ninej = Z.modelspace->GetNineJ( jk,Jab,JJx, J1,JJy,jc, Jtot,jl,J2p);
                      z_ijklmn -=   occ_factor * phase * hats * sixj * ninej * ( xabknmc*yijcabl - yabknmc*xijcabl );
                    }// for twoJx
                  }// for J2p
                }// for twoJy
              }// Z4 block
//              std::cout << "  after Z4  " << z_ijklmn << std::endl;


              // Z5 term  Pmn    Xabklnc Yijcabm
              if ( ((oa.l+ob.l+ok.l+ol.l+on.l+oc.l)%2==0) and ((oi.l+oj.l+oc.l+oa.l+ob.l+om.l)%2==0)
               and  ( (oa.tz2+ob.tz2+ok.tz2)==(ol.tz2+on.tz2+oc.tz2) ) and ( (oi.tz2+oj.tz2+oc.tz2)==(oa.tz2+ob.tz2+om.tz2)) 
                and keep_abk and keep_lnc and keep_ijc and keep_abm )
              {
                int twoJy_min = std::max( std::abs(2*J1 - oc.j2), std::abs(2*Jab - om.j2) );
                int twoJy_max = std::min( 2*J1+oc.j2 , 2*Jab + om.j2 );
                int J2p_min = std::abs( ol.j2 - on.j2 )/2;
                int J2p_max = ( ol.j2 + on.j2 )/2;
                for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                {
                  double xijcabm = X3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,m);
                  double yijcabm = Y3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,m);
                  double JJy = 0.5 * twoJy;
                  for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
                  {
                    double sixj = Z.modelspace->GetSixJ(jm,jl,J2, jn,Jtot,J2p);
                    int phase = Z.modelspace->phase( (om.j2+on.j2)/2 + J2 + J2p);
                    if (std::abs(sixj)<1e-6) continue;
                    int twoJx_min = std::max( std::abs(2*Jab - ok.j2), std::abs(2*J2p - oc.j2) );
                    int twoJx_max = std::min( 2*Jab+ok.j2 , 2*J2p + oc.j2 );
                    for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                    {
                      double JJx = 0.5 * twoJx;
                      double xabklnc = X3.GetME_pn(Jab,J2p,twoJx, a,b,k,l,n,c);
                      double yabklnc = Y3.GetME_pn(Jab,J2p,twoJx, a,b,k,l,n,c);
                      double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J2+1)*(2*J2p+1));
                      double ninej = Z.modelspace->GetNineJ( jk,Jab,JJx, J1,JJy,jc, Jtot,jm,J2p);
                      z_ijklmn -=  occ_factor * phase * hats * sixj * ninej * ( xabklnc*yijcabm - yabklnc*xijcabm );
                    }// for twoJx
                  }// for J2p
                }// for twoJy
              }// Z5 block
//              std::cout << "  after Z5  " << z_ijklmn << std::endl;


              // Z6 term  Pik Pln    Xabinmc Ykjcabl
              if ( ((oa.l+ob.l+oi.l+on.l+om.l+oc.l)%2==0) and ((ok.l+oj.l+oc.l+oa.l+ob.l+ol.l)%2==0)
               and  ( (oa.tz2+ob.tz2+oi.tz2)==(on.tz2+om.tz2+oc.tz2) ) and ( (ok.tz2+oj.tz2+oc.tz2)==(oa.tz2+ob.tz2+ol.tz2)) 
                and keep_abi and keep_nmc and keep_kjc and keep_abl )
              {
                int J1p_min = std::abs( ok.j2 - oj.j2 )/2;
                int J1p_max = ( ok.j2 + oj.j2 )/2;
                int J2p_min = std::abs( on.j2 - om.j2 )/2;
                int J2p_max = ( on.j2 + om.j2 )/2;
                for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
                {
                  double sixj2 = Z.modelspace->GetSixJ(jl,jm,J2, jn,Jtot,J2p);
                  if (std::abs(sixj2)<1e-6) continue;
                  int twoJx_min = std::max( std::abs(2*Jab - oi.j2), std::abs(2*J2p - oc.j2) );
                  int twoJx_max = std::min( 2*Jab+oi.j2 , 2*J2p + oc.j2 );
                  for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                  {
                    double JJx = 0.5 * twoJx;
                    double xabinmc = X3.GetME_pn(Jab,J2p,twoJx, a,b,i,n,m,c);
                    double yabinmc = Y3.GetME_pn(Jab,J2p,twoJx, a,b,i,n,m,c);

                    for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
                    {
                      int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - ol.j2) );
                      int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + ol.j2 );
                      double sixj1 = Z.modelspace->GetSixJ(ji,jj,J1, jk,Jtot,J1p);
                      if (std::abs(sixj1)<1e-6) continue;
                      int phase = 1;
                      for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                      {
//                         if (twoJy==5) continue;
                         double JJy = 0.5 * twoJy;
                         double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1) );
                         double ninej = Z.modelspace->GetNineJ( ji,Jab,JJx, J1p,JJy,jc, Jtot,jl,J2p);
                         double xkjcabl = X3.GetME_pn(J1p,Jab,twoJy, k,j,c,a,b,l);
                         double ykjcabl = Y3.GetME_pn(J1p,Jab,twoJy, k,j,c,a,b,l);
                         z_ijklmn +=  occ_factor * phase * hats * sixj1 * sixj2 * ninej * ( xabinmc*ykjcabl - yabinmc*xkjcabl );
                      }// for twoJy
                    }// for J1p
                  }// for twoJx
                }// for J2p
              }// Z6 block
//              std::cout << "  after Z6  " << z_ijklmn << std::endl;


              // Z7 term  Pjk Pln    Xabjnmc Yikcabl
              if ( ((oa.l+ob.l+oj.l+on.l+om.l+oc.l)%2==0) and ((oi.l+ok.l+oc.l+oa.l+ob.l+ol.l)%2==0)
               and  ( (oa.tz2+ob.tz2+oj.tz2)==(on.tz2+om.tz2+oc.tz2) ) and ( (oi.tz2+ok.tz2+oc.tz2)==(oa.tz2+ob.tz2+ol.tz2)) 
                and keep_abj and keep_nmc and keep_ikc and keep_abl )
              {
                int J1p_min = std::abs( oi.j2 - ok.j2 )/2;
                int J1p_max = ( oi.j2 + ok.j2 )/2;
                int J2p_min = std::abs( on.j2 - om.j2 )/2;
                int J2p_max = ( on.j2 + om.j2 )/2;
                for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
                {
                  double sixj2 = Z.modelspace->GetSixJ(jl,jm,J2, jn,Jtot,J2p);
                  if (std::abs(sixj2)<1e-6) continue;
                  int twoJx_min = std::max( std::abs(2*Jab - oj.j2), std::abs(2*J2p - oc.j2) );
                  int twoJx_max = std::min( 2*Jab+oj.j2 , 2*J2p + oc.j2 );
                  for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                  {
                    double JJx = 0.5 * twoJx;
                    double xabjnmc = X3.GetME_pn(Jab,J2p,twoJx, a,b,j,n,m,c);
                    double yabjnmc = Y3.GetME_pn(Jab,J2p,twoJx, a,b,j,n,m,c);

                    for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
                    {
                      int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - ol.j2) );
                      int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + ol.j2 );
                      double sixj1 = Z.modelspace->GetSixJ(jj,ji,J1, jk,Jtot,J1p);
                      if (std::abs(sixj1)<1e-6) continue;
                      int phase = -Z.modelspace->phase( (oj.j2+ok.j2)/2+J1+J1p);
                      for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                      {
                         double JJy = 0.5 * twoJy;
                         double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1) );
                         double ninej = Z.modelspace->GetNineJ( jj,Jab,JJx, J1p,JJy,jc, Jtot,jl,J2p);
                         double xikcabl = X3.GetME_pn(J1p,Jab,twoJy, i,k,c,a,b,l);
                         double yikcabl = Y3.GetME_pn(J1p,Jab,twoJy, i,k,c,a,b,l);
                         z_ijklmn +=  occ_factor * phase * hats * sixj1 * sixj2 * ninej * ( xabjnmc*yikcabl - yabjnmc*xikcabl );
                      }// for twoJy
                    }// for J1p
                  }// for twoJx
                }// for J2p
              }// Z7 block
//              std::cout << "  after Z7  " << z_ijklmn << std::endl;





              // Z8 term  Pik Pmn    Xabilnc Ykjcabm
              if ( ((oa.l+ob.l+oi.l+ol.l+on.l+oc.l)%2==0) and ((ok.l+oj.l+oc.l+oa.l+ob.l+om.l)%2==0)
               and  ( (oa.tz2+ob.tz2+oi.tz2)==(ol.tz2+on.tz2+oc.tz2) ) and ( (ok.tz2+oj.tz2+oc.tz2)==(oa.tz2+ob.tz2+om.tz2)) 
                and keep_abi and keep_lnc and keep_kjc and keep_abm )
              {
                int J1p_min = std::abs( ok.j2 - oj.j2 )/2;
                int J1p_max = ( ok.j2 + oj.j2 )/2;
                int J2p_min = std::abs( ol.j2 - on.j2 )/2;
                int J2p_max = ( ol.j2 + on.j2 )/2;
                for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
                {
                  double sixj2 = Z.modelspace->GetSixJ(jm,jl,J2, jn,Jtot,J2p);
                  if (std::abs(sixj2)<1e-6) continue;
                  int twoJx_min = std::max( std::abs(2*Jab - oi.j2), std::abs(2*J2p - oc.j2) );
                  int twoJx_max = std::min( 2*Jab+oj.j2 , 2*J2p + oc.j2 );
                  for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                  {
                    double JJx = 0.5 * twoJx;
                    double xabilnc = X3.GetME_pn(Jab,J2p,twoJx, a,b,i,l,n,c);
                    double yabilnc = Y3.GetME_pn(Jab,J2p,twoJx, a,b,i,l,n,c);

                    for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
                    {
                      int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - om.j2) );
                      int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + om.j2 );
                      double sixj1 = Z.modelspace->GetSixJ(ji,jj,J1, jk,Jtot,J1p);
                      if (std::abs(sixj1)<1e-6) continue;
                      int phase = -Z.modelspace->phase( (om.j2+on.j2)/2+J2+J2p);
                      for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                      {
                         double JJy = 0.5 * twoJy;
                         double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1) );
                         double ninej = Z.modelspace->GetNineJ( ji,Jab,JJx, J1p,JJy,jc, Jtot,jm,J2p);
                         double xkjcabm = X3.GetME_pn(J1p,Jab,twoJy, k,j,c,a,b,m);
                         double ykjcabm = Y3.GetME_pn(J1p,Jab,twoJy, k,j,c,a,b,m);
                         z_ijklmn +=  occ_factor * phase * hats * sixj1 * sixj2 * ninej * ( xabilnc*ykjcabm - yabilnc*xkjcabm );
                      }// for twoJy
                    }// for J1p
                  }// for twoJx
                }// for J2p
              }// Z8 block
//              std::cout << "  after Z8  " << z_ijklmn << std::endl;





              // Z9 term  Pjk Pmn    Xabjlnc Yikcabm
              if ( ((oa.l+ob.l+oj.l+ol.l+on.l+oc.l)%2==0) and ((oi.l+ok.l+oc.l+oa.l+ob.l+om.l)%2==0)
               and  ( (oa.tz2+ob.tz2+oj.tz2)==(ol.tz2+on.tz2+oc.tz2) ) and ( (oi.tz2+ok.tz2+oc.tz2)==(oa.tz2+ob.tz2+om.tz2)) 
                and keep_abj and keep_lnc and keep_ikc and keep_abm )
              {
                int J1p_min = std::abs( oi.j2 - ok.j2 )/2;
                int J1p_max = ( oi.j2 + ok.j2 )/2;
                int J2p_min = std::abs( ol.j2 - on.j2 )/2;
                int J2p_max = ( ol.j2 + on.j2 )/2;
                for (int J2p=J2p_min; J2p<=J2p_max; J2p++)
                {
                  double sixj2 = Z.modelspace->GetSixJ(jm,jl,J2, jn,Jtot,J2p);
                  if (std::abs(sixj2)<1e-6) continue;
                  int twoJx_min = std::max( std::abs(2*Jab - oj.j2), std::abs(2*J2p - oc.j2) );
                  int twoJx_max = std::min( 2*Jab+oj.j2 , 2*J2p + oc.j2 );
                  for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                  {
                    double JJx = 0.5 * twoJx;
                    double xabjlnc = X3.GetME_pn(Jab,J2p,twoJx, a,b,j,l,n,c);
                    double yabjlnc = Y3.GetME_pn(Jab,J2p,twoJx, a,b,j,l,n,c);

                    for (int J1p=J1p_min; J1p<=J1p_max; J1p++)
                    {
                      int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - om.j2) );
                      int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + om.j2 );
                      double sixj1 = Z.modelspace->GetSixJ(jj,ji,J1, jk,Jtot,J1p);
                      if (std::abs(sixj1)<1e-6) continue;
                      int phase = Z.modelspace->phase( (oj.j2+ok.j2+om.j2+on.j2)/2+J1+J1p+J2+J2p);
                      for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                      {
                         double JJy = 0.5 * twoJy;
                         double hats = (twoJx+1)*(twoJy+1) * sqrt( (2*J1+1)*(2*J1p+1)*(2*J2+1)*(2*J2p+1) );
                         double ninej = Z.modelspace->GetNineJ( jj,Jab,JJx, J1p,JJy,jc, Jtot,jm,J2p);
                         double xikcabm = X3.GetME_pn(J1p,Jab,twoJy, i,k,c,a,b,m);
                         double yikcabm = Y3.GetME_pn(J1p,Jab,twoJy, i,k,c,a,b,m);
                         z_ijklmn +=  occ_factor * phase * hats * sixj1 * sixj2 * ninej * ( xabjlnc*yikcabm - yabjlnc*xikcabm );
                      }// for twoJy
                    }// for J1p
                  }// for twoJx
                }// for J2p
              }// Z9 block
//              std::cout << "  after Z9  " << z_ijklmn << std::endl;



            }// for c
          }// for iket_ab
        }// for ch2
//        std::cout << "  " << i << " " << j << " " << k << " " << l << " " << m << " " << n << "  J1 J2 two J " << J1 << " " << J2 << " " << twoJ << "   Z = " << z_ijklmn << std::endl;
        Z3.AddToME_pn_ch( ch3, ch3, ibra, iket, z_ijklmn);
      }// for iket
    }// for ibra
  }//for ch3

  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////   BEGIN SCALAR-TENSOR COMMUTATORS      //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// This is no different from the scalar-scalar version
void comm111st( const Operator & X, const Operator& Y, Operator& Z)
{
   double tstart = omp_get_wtime();
   comm111ss(X,Y,Z);
   X.profiler.timer["comm111st"] += omp_get_wtime() - tstart;
}


//*****************************************************************************************
//                                       |
//      i |              i |             |
//        |    ___.Y       |__X__        |
//        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)  
//                                       |      
//---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
// X is scalar one-body, Y is tensor two-body
// There must be a better way to do this looping. 
//
//void Operator::comm121st( Operator& Y, Operator& Z) 
void comm121st( const Operator& X, const Operator& Y, Operator& Z) 
{

   double tstart = omp_get_wtime();
//   int norbits = Z.modelspace->GetNumberOrbits();
   int norbits = Z.modelspace->all_orbits.size();
   int Lambda = Z.GetJRank();
   Z.modelspace->PreCalculateSixJ();
   std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.rank_J))
   for (int indexi=0;indexi<norbits;++indexi)
//   for (int i=0;i<norbits;++i)
   {
//      auto i = Z.modelspace->all_orbits[indexi];
      auto i = allorb_vec[indexi];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double ji = 0.5*oi.j2;
      for (auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          Orbit &oj = Z.modelspace->GetOrbit(j);
          double jj = 0.5*oj.j2;
          if (j<i) continue; // only calculate upper triangle
          double& Zij = Z.OneBody(i,j);
          for (auto a : Z.modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = Z.modelspace->GetOrbit(a);
             double ja = 0.5*oa.j2;
             for (auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
             {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ);
                  int J1min = std::abs(ji-ja);
                  int J1max = ji + ja;
                  for (int J1=J1min; J1<=J1max; ++J1)
                  {
                    int phasefactor = Z.modelspace->phase(jj+ja+J1+Lambda);
                    int J2min = std::max(std::abs(Lambda - J1),std::abs(int(ja-jj)));
                    int J2max = std::min(Lambda + J1,int(ja+jj));
                    for (int J2=J2min; J2<=J2max; ++J2)
                    {
                      if ( ! ( J2>=std::abs(ja-jj) and J2<=ja+jj )) continue;
                      double prefactor = nanb*phasefactor * sqrt((2*J1+1)*(2*J2+1)) * Z.modelspace->GetSixJ(J1,J2,Lambda,jj,ji,ja);
                      Zij +=  prefactor * ( X.OneBody(a,b) * Y.TwoBody.GetTBME_J(J1,J2,b,i,a,j) - X.OneBody(b,a) * Y.TwoBody.GetTBME_J(J1,J2,a,i,b,j ));
                    }
                  }
             }
             // Now, X is scalar two-body and Y is tensor one-body
             for (auto& b : Y.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
             {

                Orbit &ob = Z.modelspace->GetOrbit(b);
                double jb = 0.5*ob.j2;
                if (std::abs(ob.occ-1) < ModelSpace::OCC_CUT) continue;
                double nanb = oa.occ * (1-ob.occ);
                int J1min = std::max(std::abs(ji-jb),std::abs(jj-ja));
                int J1max = std::min(ji+jb,jj+ja);
                double zij = 0;
                for (int J1=J1min; J1<=J1max; ++J1)
                {
                  zij -= Z.modelspace->phase(ji+jb+J1) * (2*J1+1) * Z.modelspace->GetSixJ(ja,jb,Lambda,ji,jj,J1) * X.TwoBody.GetTBME_J(J1,J1,b,i,a,j);
                }

                J1min = std::max(std::abs(ji-ja),std::abs(jj-jb));
                J1max = std::min(ji+ja,jj+jb);
                for (int J1=J1min; J1<=J1max; ++J1)
                {
                  zij += Z.modelspace->phase(ji+jb+J1) * (2*J1+1) * Z.modelspace->GetSixJ(jb,ja,Lambda,ji,jj,J1) * X.TwoBody.GetTBME_J(J1,J1,a,i,b,j) ;
                }

                Zij += nanb * Y.OneBody(a,b) * zij;
             }

             
          }
      }
   }
   
   X.profiler.timer["comm121st"] += omp_get_wtime() - tstart;
}




//*****************************************************************************************
//
//    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
//    |     |___.Y          |__X___|         
//    |     |         _     |      |          
//    |_____|               |      |_____.Y        
//    |  X  |               |      |            
//
// -- AGREES WITH NATHAN'S RESULTS
// Right now, this is the slowest one...
// Agrees with previous code in the scalar-scalar limit
//void Operator::comm122st( Operator& Y, Operator& Z ) 
//void Operator::comm122st( const Operator& X, const Operator& Y ) 
void comm122st( const Operator& X, const Operator& Y , Operator& Z) 
{
   double tstart = omp_get_wtime();
   int Lambda = Z.rank_J;

    std::vector< int > bra_channels;
    std::vector< int > ket_channels;
    for ( auto& itmat : Z.TwoBody.MatEl )
    {
      bra_channels.push_back( itmat.first[0] );
      ket_channels.push_back( itmat.first[1] );
    }
    int nmat = bra_channels.size();
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.rank_J))
    for (int ii=0; ii<nmat; ++ii)
    {
     int ch_bra = bra_channels[ii];
     int ch_ket = ket_channels[ii];

      TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      double hatfactor = sqrt((2*J1+1)*(2*J2+1));
      arma::mat& Z2 = Z.TwoBody.GetMatrix(ch_bra,ch_ket);

      for (int ibra = 0;ibra<nbras; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         double ji = oi.j2/2.0;
         double jj = oj.j2/2.0;
         for (int iket=0;iket<nkets; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit& ok = Z.modelspace->GetOrbit(k);
            Orbit& ol = Z.modelspace->GetOrbit(l);
            double jk = ok.j2/2.0;
            double jl = ol.j2/2.0;

            double cijkl = 0;
            double c1 = 0;
            double c2 = 0;
            double c3 = 0;
            double c4 = 0;

            for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
              c1 += X.OneBody(i,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,a,j,k,l);
            }
            if (i==j)
            {
              c2 = c1; // there should be a phase here, but if the ket exists, it'd better be +1.
            }
            else
            {
              for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
              {
                 c2 += X.OneBody(j,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,a,k,l);
              }
            }
            for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
               c3 += X.OneBody(a,k) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,a,l);
            }
            if (k==l)
            {
              c4 = c3 ;
            }
            else
            {
              for ( int a : X.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
              {
                 c4 += X.OneBody(a,l) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,k,a);
              }
            }

            cijkl = c1 + c2 - c3 - c4;


            c1=0;
            c2=0;
            c3=0;
            c4=0;
            int phase1 = Z.modelspace->phase(ji+jj+J2+Lambda);
            int phase2 = Z.modelspace->phase(J1-J2+Lambda);
            int phase3 = Z.modelspace->phase(J1-J2+Lambda);
            int phase4 = Z.modelspace->phase(jk+jl-J1+Lambda);


            for ( int a : Y.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
               double ja = Z.modelspace->GetOrbit(a).j2*0.5;
               c1 -=   Z.modelspace->GetSixJ(J2,J1,Lambda,ji,ja,jj) * Y.OneBody(i,a) * X.TwoBody.GetTBME(ch_ket,ch_ket,a,j,k,l) ;
            }
            if (i==j)
            {
              c2 = -c1;
            }
            else
            {
              for ( int a : Y.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
              {
                 double ja = Z.modelspace->GetOrbit(a).j2*0.5;
                 c2 +=  Z.modelspace->GetSixJ(J2,J1,Lambda,jj,ja,ji) * Y.OneBody(j,a) * X.TwoBody.GetTBME(ch_ket,ch_ket,a,i,k,l);
              }
            }
            for ( int a : Y.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
               double ja = Z.modelspace->GetOrbit(a).j2*0.5;
               c3 -=  Z.modelspace->GetSixJ(J1,J2,Lambda,jk,ja,jl) * Y.OneBody(a,k) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,l,a) ;
            }
            if (k==l)
            {
              c4 = -c3;
            }
            else
            {
              for ( int a : Y.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
              {
                 double ja = Z.modelspace->GetOrbit(a).j2*0.5;
                 c4 +=  Z.modelspace->GetSixJ(J1,J2,Lambda,jl,ja,jk) * Y.OneBody(a,l) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,k,a) ;
              }
            }
            cijkl += hatfactor*(phase1*c1+phase2*c2+phase3*c3+phase4*c4);


            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : PhysConst::SQRT2;
            Z2(ibra,iket) += cijkl /norm;
         }
      }
   }
   X.profiler.timer["comm122st"] += omp_get_wtime() - tstart;
}





// Since comm222_pp_hh and comm211 both require the construction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
// X is a scalar, Y is a tensor
//void Operator::comm222_pp_hh_221st( Operator& Y, Operator& Z )  
//void Operator::comm222_pp_hh_221st( const Operator& X, const Operator& Y )  
void comm222_pp_hh_221st( const Operator& X, const Operator& Y, Operator& Z )  
{

   double tstart = omp_get_wtime();
   int Lambda = Z.GetJRank();

   TwoBodyME Mpp = Y.TwoBody;
   TwoBodyME Mhh = Y.TwoBody;
   TwoBodyME Mff = Y.TwoBody;

   std::vector<int> vch_bra;
   std::vector<int> vch_ket;
   std::vector<const arma::mat*> vmtx;
   for ( auto& itmat : Y.TwoBody.MatEl )
   {
     vch_bra.push_back(itmat.first[0]);
     vch_ket.push_back(itmat.first[1]);
     vmtx.push_back(&(itmat.second));
   }
   size_t nchan = vch_bra.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (size_t i=0;i<nchan; ++i)
   {
    int ch_bra = vch_bra[i];
    int ch_ket = vch_ket[i];

    TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);

    auto& LHS1 = X.TwoBody.GetMatrix(ch_bra,ch_bra);
    auto& LHS2 = X.TwoBody.GetMatrix(ch_ket,ch_ket);

    auto& RHS  =  *vmtx[i];

    arma::mat& Matrixpp =  Mpp.GetMatrix(ch_bra,ch_ket);
    arma::mat& Matrixhh =  Mhh.GetMatrix(ch_bra,ch_ket);
   
    const arma::uvec& bras_pp = tbc_bra.GetKetIndex_pp();
    const arma::uvec& bras_hh = tbc_bra.GetKetIndex_hh();
    const arma::uvec& bras_ph = tbc_bra.GetKetIndex_ph();
    const arma::uvec& kets_pp = tbc_ket.GetKetIndex_pp();
    const arma::uvec& kets_hh = tbc_ket.GetKetIndex_hh();
    const arma::uvec& kets_ph = tbc_ket.GetKetIndex_ph();

    // the complicated-looking construct after the % signs just multiply the matrix elements by the proper occupation numbers (nanb, etc.)

    arma::mat MLeft  = join_horiz( LHS1.cols(bras_hh) , -RHS.cols(kets_hh) );
    arma::mat MRight = join_vert( RHS.rows(bras_hh)  % tbc_bra.Ket_occ_hh.cols( arma::uvec(RHS.n_cols,arma::fill::zeros ) ),
                                 LHS2.rows(kets_hh)  % tbc_ket.Ket_occ_hh.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) ));

    Matrixhh = MLeft * MRight;


    MLeft  = join_horiz( LHS1.cols(join_vert(bras_pp,join_vert(bras_hh,bras_ph))), -RHS.cols(join_vert(kets_pp,join_vert(kets_hh,kets_ph))) );
    MRight = join_vert( join_vert(     RHS.rows(bras_pp), 
                          join_vert( RHS.rows(bras_hh)  % tbc_bra.Ket_unocc_hh.cols( arma::uvec(RHS.n_cols,arma::fill::zeros) )  ,
                                     RHS.rows(bras_ph)  % tbc_bra.Ket_unocc_ph.cols( arma::uvec(RHS.n_cols,arma::fill::zeros) ) )),
                      join_vert(     LHS2.rows(kets_pp),
                          join_vert( LHS2.rows(kets_hh) % tbc_ket.Ket_unocc_hh.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) ),
                                     LHS2.rows(kets_ph) % tbc_ket.Ket_unocc_ph.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) )))
                     );

    Matrixpp = MLeft * MRight;
                                

    if (Z.GetParticleRank()>1)
    {
      Z.TwoBody.GetMatrix(ch_bra,ch_ket) += Matrixpp - Matrixhh;
    }

   }// for itmat

      // The one body part takes some additional work

//   int norbits = Z.modelspace->GetNumberOrbits();
   int norbits = Z.modelspace->all_orbits.size();
   std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(),Z.modelspace->all_orbits.end());
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.rank_J))
   for (int indexi=0;indexi<norbits;++indexi)
//   for (int i=0;i<norbits;++i)
   {
//      auto i = Z.modelspace->all_orbits[indexi];
      auto i = allorb_vec[indexi];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double ji = oi.j2/2.0;
      for (auto j : Z.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         if (j<i) continue;
         Orbit &oj = Z.modelspace->GetOrbit(j);
         double jj = oj.j2/2.0;
         double cijJ = 0;
         // Sum c over holes and include the nbar_a * nbar_b terms
           for (auto& c : Z.modelspace->holes)
           {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2/2.0;
              int j1min = std::abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = std::max( int(std::abs(jc-jj)), std::abs(Lambda-J1) );
               int j2max = std::min( int(jc+jj), J1+Lambda ); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = Z.modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * Z.modelspace->phase(jj + jc + J1 + Lambda)
                        * ( oc.occ * Mpp.GetTBME_J(J1,J2,c,i,c,j) + (1-oc.occ) * Mhh.GetTBME_J(J1,J2,c,i,c,j));
               }
              }
           // Sum c over particles and include the n_a * n_b terms
           }
           for (auto& c : Z.modelspace->particles)
           {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2/2.0;
              int j1min = std::abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = std::max( int(std::abs(jc-jj)), std::abs(Lambda-J1) );
               int j2max = std::min( int(jc+jj), J1+Lambda ); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = Z.modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * Z.modelspace->phase(jj + jc + J1 + Lambda) * Mhh.GetTBME_J(J1,J2,c,i,c,j);
               }
              }
           }
//         #pragma omp critical
         Z.OneBody(i,j) += cijJ ;
      } // for j
    } // for i
    X.profiler.timer["comm222_pp_hh_221st"] += omp_get_wtime() - tstart;
}






//**************************************************************************
//
//  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
//                        { k l J'}
// TENSOR VARIETY
/// The scalar Pandya transformation is defined as
/// \f[
///  \bar{X}^{J}_{i\bar{j}k\bar{l}} = - \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  X^{J}_{ilkj}
/// \f]
/// where the overbar indicates time-reversed orbits.
/// This function is designed for use with comm222_phss() and so it takes in
/// two arrays of matrices, one for hp terms and one for ph terms.
//void Operator::DoTensorPandyaTransformation(vector<arma::mat>& TwoBody_CC_hp, vector<arma::mat>& TwoBody_CC_ph)
//void Operator::DoTensorPandyaTransformation(map<array<int,2>,arma::mat>& TwoBody_CC_hp, map<array<int,2>,arma::mat>& TwoBody_CC_ph) const
void DoTensorPandyaTransformation( const Operator& Z, std::map<std::array<index_t,2>,arma::mat>& TwoBody_CC_ph)
{
   int Lambda = Z.rank_J;
   // loop over cross-coupled channels
   index_t nch = Z.modelspace->SortedTwoBodyChannels_CC.size();

   // Allocate map for matrices -- this needs to be serial.
   for ( index_t ch_bra_cc : Z.modelspace->SortedTwoBodyChannels_CC )
   {
      TwoBodyChannel& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
      index_t nph_bras = bras_ph.n_rows;
      for ( index_t ch_ket_cc : Z.modelspace->SortedTwoBodyChannels_CC )
      {
        TwoBodyChannel& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
        index_t nKets_cc = tbc_ket_cc.GetNumberKets();

         TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}] =  arma::mat(2*nph_bras,   nKets_cc, arma::fill::zeros);
      }
   }

   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Lambda))
   for (index_t ich=0;ich<nch;++ich)
   {
      index_t ch_bra_cc = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      int Jbra_cc = tbc_bra_cc.J;
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
//      arma::uvec& bras_ph = tbc_bra_cc.GetKetIndex_ph();
      index_t nph_bras = bras_ph.size();

      for ( index_t ch_ket_cc : Z.modelspace->SortedTwoBodyChannels_CC )
      {
        TwoBodyChannel& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
        int Jket_cc = tbc_ket_cc.J;
        if ( (Jbra_cc+Jket_cc < Z.GetJRank()) or std::abs(Jbra_cc-Jket_cc)>Z.GetJRank() ) continue;
        if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + Z.GetParity())%2>0 ) continue;

        index_t nKets_cc = tbc_ket_cc.GetNumberKets();

//        arma::mat& MatCC_hp = TwoBody_CC_hp[{ch_bra_cc,ch_ket_cc}];
        arma::mat& MatCC_ph = TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}];
        // loop over ph bras <ad| in this channel
        for (index_t ibra=0; ibra<nph_bras; ++ibra)
        {
           Ket & bra_cc = tbc_bra_cc.GetKet( bras_ph[ibra] );
           index_t a = bra_cc.p;
           index_t b = bra_cc.q;
           Orbit & oa = Z.modelspace->GetOrbit(a);
           Orbit & ob = Z.modelspace->GetOrbit(b);
           double ja = oa.j2*0.5;
           double jb = ob.j2*0.5;

           // loop over kets |bc> in this channel
           index_t iket_max =  nKets_cc ;
           for (index_t iket_cc=0; iket_cc<iket_max; ++iket_cc)
           {
              Ket & ket_cc = tbc_ket_cc.GetKet(iket_cc%nKets_cc);
              index_t c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
              index_t d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
              Orbit & oc = Z.modelspace->GetOrbit(c);
              Orbit & od = Z.modelspace->GetOrbit(d);
              double jc = oc.j2*0.5;
              double jd = od.j2*0.5;


              int j1min = std::abs(ja-jd);
              int j1max = ja+jd;
              double sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = std::max( int(std::abs(jc-jb)), std::abs(J1-Lambda) );
                int j2max = std::min(int(jc+jb),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
//                  double ninej = Z.modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
                    double ninej = 0;
                    if ( Lambda==0 )
                    {
                      ninej = AngMom::phase( jb + jd + J1 + Jbra_cc ) * Z.modelspace->GetSixJ( ja, jb, Jbra_cc, jc, jd, J1 ) / sqrt((2*J2+1)*(2*Jbra_cc+1) );
                    }
                    else 
                    {
                      ninej = Z.modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
                    }


                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = Z.TwoBody.GetTBME_J(J1,J2,a,d,c,b);
                  sm -= hatfactor * Z.modelspace->phase(jb+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
              MatCC_ph(ibra,iket_cc) = sm;

              // Exchange (a <-> b) to account for the (n_a - n_b) term
              // Get Tz,parity and range of J for <bd || ca > coupling
              j1min = std::abs(jb-jd);
              j1max = jb+jd;
              sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = std::max( int(std::abs(jc-ja)), std::abs(J1-Lambda) );
                int j2max = std::min(int(jc+ja),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
//                  double ninej = Z.modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  double ninej = 0;
                  if ( Lambda==0 )
                  {
                    ninej = AngMom::phase( ja + jd + J1 + Jbra_cc ) * Z.modelspace->GetSixJ( jb, ja, Jbra_cc, jc, jd, J1 ) / sqrt((2*J2+1)*(2*Jbra_cc+1) );
                  }
                  else 
                  {
                    ninej = Z.modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  }


                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = Z.TwoBody.GetTBME_J(J1,J2,b,d,c,a);
                  sm -= hatfactor * Z.modelspace->phase(ja+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
              MatCC_ph(ibra+nph_bras,iket_cc) = sm;

           }
        }
    }
   }
}



// This happens inside an OMP loop, and so everything here needs to be thread safe
//
//void Operator::DoTensorPandyaTransformation_SingleChannel( arma::mat& TwoBody_CC_ph, int ch_bra_cc, int ch_ket_cc) const
//void Operator::DoTensorPandyaTransformation_SingleChannel( arma::mat& MatCC_ph, int ch_bra_cc, int ch_ket_cc) const
void DoTensorPandyaTransformation_SingleChannel( const Operator& Z, arma::mat& MatCC_ph, int ch_bra_cc, int ch_ket_cc)
{
   int Lambda = Z.rank_J;

   TwoBodyChannel& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
   arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
   int nph_bras = bras_ph.n_rows;

   TwoBodyChannel& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
   int nKets_cc = tbc_ket_cc.GetNumberKets();

   // The Pandya-transformed (formerly cross-coupled) particle-hole type matrix elements
   // (this is the output of this method)
   MatCC_ph =  arma::mat(2*nph_bras,   nKets_cc, arma::fill::zeros);

   int Jbra_cc = tbc_bra_cc.J;
   int Jket_cc = tbc_ket_cc.J;
   if ( (Jbra_cc+Jket_cc < Z.GetJRank() ) or std::abs(Jbra_cc-Jket_cc)>Z.GetJRank() ) return;
   if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + Z.GetParity())%2>0 ) return;


   // loop over ph bras <ad| in this channel
   for (int ibra=0; ibra<nph_bras; ++ibra)
   {
      Ket & bra_cc = tbc_bra_cc.GetKet( bras_ph[ibra] );
      int a = bra_cc.p;
      int b = bra_cc.q;
      Orbit & oa = Z.modelspace->GetOrbit(a);
      Orbit & ob = Z.modelspace->GetOrbit(b);
      double ja = oa.j2*0.5;
      double jb = ob.j2*0.5;

      // loop over kets |bc> in this channel
      int iket_max =  nKets_cc ;
      for (int iket_cc=0; iket_cc<iket_max; ++iket_cc)
      {
         Ket & ket_cc = tbc_ket_cc.GetKet(iket_cc%nKets_cc);
         int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
         int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
         Orbit & oc = Z.modelspace->GetOrbit(c);
         Orbit & od = Z.modelspace->GetOrbit(d);
         double jc = oc.j2*0.5;
         double jd = od.j2*0.5;


         int j1min = std::abs(ja-jd);
         int j1max = ja+jd;
         double sm = 0;
         for (int J1=j1min; J1<=j1max; ++J1)
         {
           int j2min = std::max( int(std::abs(jc-jb)), std::abs(J1-Lambda) );
           int j2max = std::min( int(jc+jb), J1+Lambda) ;
           for (int J2=j2min; J2<=j2max; ++J2)
           {
//             double ninej = Z.modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
             double ninej = 0;
             if ( Lambda==0 )
             {
               ninej = AngMom::phase( jb + jd + J1 + Jbra_cc ) * Z.modelspace->GetSixJ( ja, jb, Jbra_cc, jc, jd, J1 ) / sqrt((2*J2+1)*(2*Jbra_cc+1) );
             }
             else 
             {
               ninej = Z.modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
             }

             if (std::abs(ninej) < 1e-8) continue;
             double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
             double tbme = Z.TwoBody.GetTBME_J(J1,J2,a,d,c,b);
             sm -= hatfactor * Z.modelspace->phase(jb+jd+Jket_cc+J2) * ninej * tbme ;
           }
         }
         MatCC_ph(ibra,iket_cc) = sm;

         // Exchange (a <-> b) to account for the (n_a - n_b) term
         if (a==b)
         {
           MatCC_ph(ibra+nph_bras,iket_cc) = sm;
         }
         else
         {

           // Get Tz,parity and range of J for <bd || ca > coupling
           j1min = std::abs(jb-jd);
           j1max = jb+jd;
           sm = 0;
           for (int J1=j1min; J1<=j1max; ++J1)
           {
             int j2min = std::max( int(std::abs(jc-ja)), std::abs(J1-Lambda) );
             int j2max = std::min( int(jc+ja), J1+Lambda );
             for (int J2=j2min; J2<=j2max; ++J2)
             {
//               double ninej = Z.modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
               double ninej = 0;
               if ( Lambda==0 )
               {
                 ninej = AngMom::phase( ja + jd + J1 + Jbra_cc ) * Z.modelspace->GetSixJ( jb, ja, Jbra_cc, jc, jd, J1 ) / sqrt((2*J2+1)*(2*Jbra_cc+1) );
               }
               else 
               {
                 ninej = Z.modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
               }


               if (std::abs(ninej) < 1e-8) continue;
               double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
               double tbme = Z.TwoBody.GetTBME_J(J1,J2,b,d,c,a);
               sm -= hatfactor * Z.modelspace->phase(ja+jd+Jket_cc+J2) * ninej * tbme ;
             }
           }
           MatCC_ph(ibra+nph_bras,iket_cc) = sm;
         }
      }
   }
}





//// Take Pandya-transformed matrix Zbar for a single channel, invert the Pandya transformation and add the result to the current operator.
////void Operator::AddInverseTensorPandyaTransformation_SingleChannel(arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc)
////void AddInverseTensorPandyaTransformation_SingleChannel(Operator& Z, arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc)
//void AddInverseTensorPandyaTransformation_SingleChannel(Operator& Z, arma::mat& Zbar, size_t ch_bra_cc, size_t ch_ket_cc)
//{
//    // Do the inverse Pandya transform
//   if (ch_bra_cc > ch_ket_cc)  // hopefully this won't happen
//   {
//      std::cout << "WARNING: Called Operator::AddInverseTensorPandyaTransformation_SingleChannel with ch_bra_cc > ch_ket_cc : " << ch_bra_cc << " > " << ch_ket_cc << std::endl;
//      std::cout << " Skipping this channel." << std::endl;
//      return;
//   }
////   int n_channels_kept = 0;
//   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup(Z.GetJRank(), Z.GetTRank(),Z.GetParity())[{ch_bra_cc,ch_ket_cc}];
//   int Lambda = Z.rank_J;
//   int hZ = Z.IsHermitian() ? 1 : -1;
//   TwoBodyChannel_CC& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
//   TwoBodyChannel_CC& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
//   int J3 = tbc_bra_cc.J;
//   int J4 = tbc_ket_cc.J;
//   int Tz_bra_cc = tbc_bra_cc.Tz;
//   int Tz_ket_cc = tbc_ket_cc.Tz;
//   int parity_bra_cc = tbc_bra_cc.parity;
//   int parity_ket_cc = tbc_ket_cc.parity;
//   int nkets_cc = tbc_ket_cc.GetNumberKets();
//
////   int nchannels = pandya_lookup[0].size();
//   int nchannels = pandya_lookup.size();
//   for (int ich=0;ich<nchannels;++ich)
//   {
//      int ch_bra = pandya_lookup[0][ich];
//      int ch_ket = pandya_lookup[1][ich];
//      TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
//      TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
//      int J1 = tbc_bra.J;
//      int J2 = tbc_ket.J;
//      int nBras = tbc_bra.GetNumberKets();
//      int nKets = tbc_ket.GetNumberKets();
//      arma::mat& Zijkl = Z.TwoBody.GetMatrix(ch_bra,ch_ket);
////      bool inner_loop = false;
//
//      double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
//
//      for (int ibra=0; ibra<nBras; ++ibra)
//      {
//         Ket & bra = tbc_bra.GetKet(ibra);
//         int i = bra.p;
//         int j = bra.q;
//         Orbit & oi = Z.modelspace->GetOrbit(i);
//         Orbit & oj = Z.modelspace->GetOrbit(j);
//         double ji = 0.5*oi.j2;
//         double jj = 0.5*oj.j2;
//         int ketmin = ch_bra==ch_ket ? ibra : 0;
//         for (int iket=ketmin; iket<nKets; ++iket)
//         {
//            Ket & ket = tbc_ket.GetKet(iket);
//            int k = ket.p;
//            int l = ket.q;
//            Orbit & ok = Z.modelspace->GetOrbit(k);
//            Orbit & ol = Z.modelspace->GetOrbit(l);
//            double jk = 0.5*ok.j2;
//            double jl = 0.5*ol.j2;
//
//            double commij = 0;
//            double commji = 0;
//
//            // Transform Z_ilkj
//            int j3min = std::abs(int(ji-jl));
//            int j3max = ji+jl;
//            int j4min = std::abs(int(jk-jj));
//            int j4max = int(jk+jj);
//
//            if (   (oi.l+ol.l)%2==parity_bra_cc             and (ok.l+oj.l)%2==parity_ket_cc
////                      and std::abs(oi.tz2+ol.tz2)==2*Tz_bra_cc   and std::abs(ok.tz2+oj.tz2)==2*Tz_ket_cc
//                      and std::abs(oi.tz2-ol.tz2)==2*Tz_bra_cc   and std::abs(ok.tz2-oj.tz2)==2*Tz_ket_cc
//                      and j3min<=J3 and J3<=j3max           and j4min<=J4 and J4<=j4max )
//            {
//               int indx_il = tbc_bra_cc.GetLocalIndex( std::min(i,l), std::max(i,l) );
//               int indx_kj = tbc_ket_cc.GetLocalIndex( std::min(j,k), std::max(j,k) );
//
////               double ninej = Z.modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
//               double ninej = 0;
//               if ( Lambda==0 )
//               {
//                 ninej = AngMom::phase( jj + jl + J1 + J3 ) * Z.modelspace->GetSixJ( ji, jj, J1, jk, jl, J3 ) / sqrt((2*J2+1)*(2*J4+1) );
//               }
//               else 
//               {
//                 ninej = Z.modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
//               }
//               double tbme = 0;
//
//               if (i<=l) tbme = Zbar( indx_il ,indx_kj+(k>j?nkets_cc:0) );
//               else      tbme = Zbar( indx_il ,indx_kj+(k>j?0:nkets_cc) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
//               commij += hatfactor * Z.modelspace->phase(jj+jl+J2+J4) * ninej * tbme ;
////               inner_loop = true;
//            }
//
//            if (  (ch_bra_cc != ch_ket_cc)  and ((oi.l+ol.l)%2==parity_ket_cc)           and ((ok.l+oj.l)%2==parity_bra_cc)
////                                            and (std::abs(oi.tz2+ol.tz2)==2*Tz_ket_cc)   and (std::abs(ok.tz2+oj.tz2)==2*Tz_bra_cc)
//                                            and (std::abs(oi.tz2-ol.tz2)==2*Tz_ket_cc)   and (std::abs(ok.tz2-oj.tz2)==2*Tz_bra_cc)
//                                            and (j3min<=J4) and (J4<=j3max)  and  (j4min<=J3) and (J3<=j4max )    )
//              {
//                 int indx_kj = tbc_bra_cc.GetLocalIndex( std::min(j,k), std::max(j,k));
//                 int indx_il = tbc_ket_cc.GetLocalIndex( std::min(i,l), std::max(i,l));
//
////                 double ninej = Z.modelspace->GetNineJ(ji,jl,J4,jj,jk,J3,J1,J2,Lambda);
//                 double ninej = 0;
//                 if ( Lambda==0 )
//                 {
//                   ninej = AngMom::phase( jj + jl + J1 + J4 ) * Z.modelspace->GetSixJ( ji, jj, J1, jk, jl, J4 ) / sqrt((2*J2+1)*(2*J4+1) );
//                 }
//                 else 
//                 {
//                   ninej = Z.modelspace->GetNineJ(ji,jl,J4,jj,jk,J3,J1,J2,Lambda);
//                 }
//
//                 double tbme = 0;
//
//                 if(k<=j) tbme = Zbar(indx_kj, indx_il+(i>l?nkets_cc:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
//                 else     tbme = Zbar(indx_kj, indx_il+(i>l?0:nkets_cc) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
//
//                 commij += hatfactor * Z.modelspace->phase(jj+jl+J2+J3) * ninej * tbme ;
////               inner_loop = true;
//              }
//
//
//            if (i==j)
//            {
//              commji = commij;
//            }
//            else
//            {
//              j3min = std::abs(int(jj-jl));
//              j3max = jj+jl;
//              j4min = std::abs(int(jk-ji));
//              j4max = int(jk+ji);
//              if (   ((oj.l+ol.l)%2==parity_bra_cc) and ((ok.l+oi.l)%2==parity_ket_cc)
////                      and (std::abs(oj.tz2+ol.tz2)==2*Tz_bra_cc) and (std::abs(ok.tz2+oi.tz2)==2*Tz_ket_cc)
//                      and (std::abs(oj.tz2-ol.tz2)==2*Tz_bra_cc) and (std::abs(ok.tz2-oi.tz2)==2*Tz_ket_cc)
//                      and (J3>=j3min) and (J3<=j3max) and (J4>=j4min) and (J4<=j4max)   )
//              {
//                // Transform Z_jlki
//                int indx_jl = tbc_bra_cc.GetLocalIndex( std::min(j,l), std::max(j,l) );
//                int indx_ki = tbc_ket_cc.GetLocalIndex( std::min(k,i), std::max(k,i) );
//
////                double ninej = Z.modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
//                double ninej = 0;
//                if ( Lambda==0 )
//                {
//                  ninej = AngMom::phase( ji + jl + J1 + J3 ) * Z.modelspace->GetSixJ( jj, ji, J1, jk, jl, J3 ) / sqrt((2*J2+1)*(2*J4+1) );
//                }
//                else 
//                {
//                  ninej = Z.modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
//                }
//
//
//                double tbme = 0;
//                if (j<=l) tbme = Zbar( indx_jl ,indx_ki+(k>i?nkets_cc:0) );
//                else      tbme = Zbar( indx_jl ,indx_ki+(k>i?0:nkets_cc) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
//                commji += hatfactor * Z.modelspace->phase(ji+jl+J2+J4) * ninej * tbme ;
////               inner_loop = true;
//              }
//              if ( (ch_bra_cc!=ch_ket_cc) and ((oj.l+ol.l)%2==parity_ket_cc) and ((ok.l+oi.l)%2==parity_bra_cc)
////                                       and (std::abs(oj.tz2+ol.tz2)==2*Tz_ket_cc) and (std::abs(ok.tz2+oi.tz2)==2*Tz_bra_cc)
//                                       and (std::abs(oj.tz2-ol.tz2)==2*Tz_ket_cc) and (std::abs(ok.tz2-oi.tz2)==2*Tz_bra_cc)
//                                       and (J4>=j3min) and (J4<=j3max) and (J3>=j4min) and (J3<=j4max)   )
//              {
//                int indx_jl = tbc_ket_cc.GetLocalIndex( std::min(j,l), std::max(j,l));
//                int indx_ki = tbc_bra_cc.GetLocalIndex( std::min(k,i), std::max(k,i));
////                double ninej = Z.modelspace->GetNineJ(jj,jl,J4,ji,jk,J3,J1,J2,Lambda);
//                double ninej = 0;
//                if ( Lambda==0 )
//                {
//                  ninej = AngMom::phase( ji + jl + J1 + J4 ) * Z.modelspace->GetSixJ( jj, ji, J1, jk, jl, J4 ) / sqrt((2*J2+1)*(2*J4+1) );
//                }
//                else 
//                {
//                  ninej = Z.modelspace->GetNineJ(jj,jl,J4,ji,jk,J3,J1,J2,Lambda);
//                }
//
//                double tbme = 0;
//
//                if(k<=i) tbme = Zbar(indx_ki, indx_jl+(j>l?nkets_cc:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
//                else     tbme = Zbar(indx_ki, indx_jl+(j>l?0:nkets_cc) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
//                commji += hatfactor * Z.modelspace->phase(ji+jl+J2+J3) * ninej * tbme ;
////               inner_loop = true;
//              }
//            }
//
//
//            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : PhysConst::SQRT2;
//            #pragma omp atomic
//            Zijkl(ibra,iket) +=  (commij - Z.modelspace->phase(ji+jj-J1)*commji) / norm;
//            if (ch_bra==ch_ket) 
//            {
//               #pragma omp atomic write
//               Zijkl(iket,ibra) = hZ * Zijkl(ibra,iket);
//            }
//         }
//      }
////      if (inner_loop)  // if we made it to the inner loop, count it.
////      {
////        n_channels_kept++;
////      }
//   }
//}





//void Operator::AddInverseTensorPandyaTransformation( const std::map<std::array<index_t,2>,arma::mat>&  Zbar )
void AddInverseTensorPandyaTransformation( Operator& Z, const std::map<std::array<index_t,2>,arma::mat>&  Zbar )
{
    // Do the inverse Pandya transform
   int Lambda = Z.rank_J;
//   std::vector<std::map<std::array<int,2>,arma::mat>::iterator> iteratorlist;
   std::vector<std::map<std::array<size_t,2>,arma::mat>::iterator> iteratorlist;
//   for (std::map<std::array<int,2>,arma::mat>::iterator iter= Z.TwoBody.MatEl.begin(); iter!= Z.TwoBody.MatEl.end(); ++iter) iteratorlist.push_back(iter);
   for (auto iter= Z.TwoBody.MatEl.begin(); iter!= Z.TwoBody.MatEl.end(); ++iter) iteratorlist.push_back(iter);
   int niter = iteratorlist.size();
   int hZ = Z.IsHermitian() ? 1 : -1;

    // Only go parallel if we've previously calculated the SixJs/NineJs. Otherwise, it's not thread safe.
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Lambda))
   for (int i=0; i<niter; ++i)
   {
      const auto iter = iteratorlist[i];
      int ch_bra = iter->first[0];
      int ch_ket = iter->first[1];
      arma::mat& Zijkl = iter->second;
      const TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      const TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      index_t nBras = tbc_bra.GetNumberKets();
      index_t nKets = tbc_ket.GetNumberKets();

      for (index_t ibra=0; ibra<nBras; ++ibra)
      {
         const Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         const Orbit & oi = Z.modelspace->GetOrbit(i);
         const Orbit & oj = Z.modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;
         index_t ketmin = ch_bra==ch_ket ? ibra : 0;
         for (index_t iket=ketmin; iket<nKets; ++iket)
         {
            const Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            const Orbit & ok = Z.modelspace->GetOrbit(k);
            const Orbit & ol = Z.modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;

            double commij = 0;
            double commji = 0;

            // Transform Z_ilkj
            int parity_bra_cc = (oi.l+ol.l)%2;
            int parity_ket_cc = (ok.l+oj.l)%2;
//            int Tz_bra_cc = std::abs(oi.tz2+ol.tz2)/2;
//            int Tz_ket_cc = std::abs(ok.tz2+oj.tz2)/2;
            int Tz_bra_cc = std::abs(oi.tz2-ol.tz2)/2;
            int Tz_ket_cc = std::abs(ok.tz2-oj.tz2)/2;
            int j3min = std::abs(int(ji-jl));
            int j3max = ji+jl;
            for (int J3=j3min; J3<=j3max; ++J3)
            {
              index_t ch_bra_cc = Z.modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
              const TwoBodyChannel_CC& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
              index_t nbras = tbc_bra_cc.GetNumberKets();
              index_t indx_il = tbc_bra_cc.GetLocalIndex( std::min(i,l), std::max(i,l) );
              int j4min = std::max( std::abs(int(jk-jj)), std::abs(J3-Lambda) );
              int j4max = std::min( int(jk+jj), J3+Lambda );
              for (int J4=j4min; J4<=j4max; ++J4)
              {

                 index_t ch_ket_cc = Z.modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                 const TwoBodyChannel_CC& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                 index_t nkets = tbc_ket_cc.GetNumberKets();
                 index_t indx_kj = tbc_ket_cc.GetLocalIndex( std::min(j,k), std::max(j,k) );

//                  double ninej = Z.modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
                  double ninej=0;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
                  double tbme = 0;


                  if (Lambda==0)
                  {
                    ninej =  AngMom::phase( jj+jl+J1+J3) * Z.modelspace->GetSixJ(ji,jj,J1, jk,jl, J3) / sqrt( (2*J2+1)*(2*J4+1))  ;
                  }
                  else
                  {
                     ninej = Z.modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
                  }
                  if (std::abs(ninej) < 1e-8) continue;
                  index_t ch_lo = std::min(ch_bra_cc,ch_ket_cc);
                  index_t ch_hi = std::max(ch_bra_cc,ch_ket_cc);
                  auto zbar_iter = Zbar.find({ch_lo,ch_hi});
                  if (zbar_iter == Zbar.end()) continue;
                  const auto& Zmat = zbar_iter->second;

                  if (ch_bra_cc <= ch_ket_cc)
                  {
                    if (i<=l) tbme = Zmat( indx_il ,indx_kj+(k>j?nkets:0) );
                    else      tbme = Zmat( indx_il ,indx_kj+(k>j?0:nkets) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
                  }
                  else
                  {
                      if(k<=j) tbme = Zmat(indx_kj, indx_il+(i>l?nbras:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                      else     tbme = Zmat(indx_kj, indx_il+(i>l?0:nbras) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                  }

/*
*/
                  commij += hatfactor * Z.modelspace->phase(jj+jl+J2+J4) * ninej * tbme ;
              }
            }

            if (i==j)
            {
              commji = commij;
            }
            else
            {
              // Transform Z_jlki

              parity_bra_cc = (oj.l+ol.l)%2;
              parity_ket_cc = (ok.l+oi.l)%2;
              Tz_bra_cc = std::abs(oj.tz2-ol.tz2)/2;
              Tz_ket_cc = std::abs(ok.tz2-oi.tz2)/2;
              j3min = std::abs(int(jj-jl));
              j3max = jj+jl;
  
              for (int J3=j3min; J3<=j3max; ++J3)
              {

                int ch_bra_cc = Z.modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
                const TwoBodyChannel_CC& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
                int nbras = tbc_bra_cc.GetNumberKets();
                int indx_jl = tbc_bra_cc.GetLocalIndex( std::min(j,l), std::max(j,l) );
                int j4min = std::max( std::abs(int(jk-ji)), std::abs(J3-Lambda) );
                int j4max = std::min( int(jk+ji), J3+Lambda );
                for (int J4=j4min; J4<=j4max; ++J4)
                {

                   int ch_ket_cc = Z.modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                   const TwoBodyChannel_CC& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                   int nkets = tbc_ket_cc.GetNumberKets();
                   int indx_ki = tbc_ket_cc.GetLocalIndex( std::min(k,i), std::max(k,i) );

//                    double ninej = Z.modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);

                   double ninej=0;
                   double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
                   double tbme = 0;

                   if (Lambda==0)
                   {
                     ninej =  AngMom::phase( ji+jl+J1+J3) * Z.modelspace->GetSixJ(jj,ji,J1, jk,jl, J3) / sqrt( (2*J2+1)*(2*J4+1))  ;
                   }
                   else
                   {
                      ninej = Z.modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
                   }


                    if (std::abs(ninej) < 1e-8) continue;
  
                    index_t ch_lo = std::min(ch_bra_cc,ch_ket_cc);
                    index_t ch_hi = std::max(ch_bra_cc,ch_ket_cc);
                    auto zbar_iter = Zbar.find({ch_lo,ch_hi});
                    if (zbar_iter == Zbar.end()) continue;
                    const auto& Zmat = zbar_iter->second;

                    if (ch_bra_cc <= ch_ket_cc)
                    {
                      if (j<=l) tbme = Zmat( indx_jl ,indx_ki+(k>i?nkets:0) );
                      else      tbme = Zmat( indx_jl ,indx_ki+(k>i?0:nkets) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
                    }
                    else
                    {
                        if(k<=i) tbme = Zmat(indx_ki, indx_jl+(j>l?nbras:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                        else     tbme = Zmat(indx_ki, indx_jl+(j>l?0:nbras) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                    }
/*
*/  
  
                      commji += hatfactor * Z.modelspace->phase(ji+jl+J2+J4) * ninej * tbme ;
                }
              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : PhysConst::SQRT2;
            Zijkl(ibra,iket) +=  (commij - Z.modelspace->phase(ji+jj-J1)*commji) / norm;
            if (ch_bra==ch_ket) Zijkl(iket,ibra) = hZ * Zijkl(ibra,iket);
         }
      }
   }
}




//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.     
//                                             
//   |          |      |          |           
//   |     __Y__|      |     __X__|            
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |            
//   |____\/    |      |____\/    |            
//   |  X       |      |  Y       |            
//           
//            
// -- This appears to agree with Nathan's results
//
/// Calculates the part of \f$ [X_{(2)},\mathbb{Y}^{\Lambda}_{(2)}]_{ijkl} \f$ which involves ph intermediate states, here indicated by \f$ \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} \f$
/// \f[
/// \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} = \sum_{abJ_3J_4}(n_a-n_b) \hat{J_1}\hat{J_2}\hat{J_3}\hat{J_4}
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_l  &  J_3 \\
///  j_j  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
/// \left( \bar{X}^{J3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
///  -(-1)^{j_i+j_j-J_1}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_l  &  J_3 \\
///  j_i  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
/// \left( \bar{X}^{J_3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
/// \right]
/// \f]
/// This is implemented by defining an intermediate matrix
/// \f[
/// \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{i\bar{l}k\bar{j}} \equiv \sum_{ab}(n_a\bar{n}_b)
/// \left[ \left( \bar{X}^{J3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
/// -\left( \bar{X}^{J_3}_{i\bar{l}b\bar{a}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{b\bar{a}k\bar{j}} - 
///    \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}b\bar{a}}\bar{X}^{J_4}_{b\bar{a}k\bar{j}} \right)\right]
/// \f]
/// The Pandya-transformed matrix elements are obtained with DoTensorPandyaTransformation().
/// The matrices \f$ \bar{X}^{J_3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} \f$
/// and \f$ \bar{\mathbb{Y}}^{J_4J_3\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_3}_{a\bar{b}k\bar{j}} \f$
/// are related by a Hermitian conjugation, which saves two matrix multiplications, provided we
/// take into account the phase \f$ (-1)^{J_3-J_4} \f$ from conjugating the spherical tensor.
/// The commutator is then given by
/// \f[
/// \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} = \sum_{J_3J_4} \hat{J_1}\hat{J_2}\hat{J_3}\hat{J_4}
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_l  &  J_3 \\
///  j_j  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
///  \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{i\bar{l}k\bar{j}}
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_l  &  J_3 \\
///  j_i  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
///  \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{j\bar{l}k\bar{i}}
///  \right]
///  \f]
///
//void Operator::comm222_phst( Operator& Y, Operator& Z ) 
//void Operator::comm222_phst( const Operator& X, const Operator& Y ) 
void comm222_phst( const Operator& X, const Operator& Y, Operator& Z ) 
{

   int hX = X.IsHermitian() ? 1 : -1;
   int hY = Y.IsHermitian() ? 1 : -1;

   double t_start = omp_get_wtime();
   // We reuse Xt_bar multiple times, so it makes sense to calculate them once and store them in a deque.
   std::deque<arma::mat> Xt_bar_ph = InitializePandya( Z, Z.nChannels, "transpose"); // We re-use the scalar part multiple times, so there's a significant speed gain for saving it
   std::map<std::array<index_t,2>,arma::mat> Y_bar_ph;
   DoPandyaTransformation(X, Xt_bar_ph, "transpose" );
   X.profiler.timer["DoTensorPandyaTransformationX"] += omp_get_wtime() - t_start;


   t_start = omp_get_wtime();
   // Construct the intermediate matrix Z_bar
   // First, we initialize the map Z_bar with empty matrices
   // to avoid problems in the parallel loop -- (do we even want a parallel loop here?)
   std::map<std::array<index_t,2>,arma::mat> Z_bar;

   t_start = omp_get_wtime();
   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup(Z.GetJRank(), Z.GetTRank(), Z.GetParity() );
   X.profiler.timer["PandyaLookup"] += omp_get_wtime() - t_start;

   std::vector<index_t> ybras;
   std::vector<index_t> ykets;
//   for (auto ich_bra : Z.modelspace->SortedTwoBodyChannels_CC)
   for (auto ich_bra : pandya_lookup )
   {
     auto tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ich_bra);
//     int n_rows = Z.modelspace->GetTwoBodyChannel_CC(ich_bra).GetNumberKets();
     int n_rows = tbc_bra_cc.GetNumberKets();
     if (n_rows<1) continue;
//     for (auto ich_ket : Z.modelspace->SortedTwoBodyChannels_CC)
     for (auto ich_ket : pandya_lookup)
     {
       if (ich_bra>ich_ket) continue;
//       if (pandya_lookup.at({(int)ich_bra,(int)ich_ket})[0].size()<1) continue;
//       if (pandya_lookup.at({ich_bra,ich_ket})[0].size()<1) continue;
         auto tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ich_ket);
         int n_cols = 2*tbc_ket_cc.GetNumberKets();
//         int n_cols = 2*Z.modelspace->GetTwoBodyChannel_CC(ich_ket).GetNumberKets();
         if (n_cols<1) continue;
         if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + Z.parity)%2>0 ) continue;
         if ( (tbc_bra_cc.J + tbc_ket_cc.J < Z.GetJRank() ) or ( std::abs(tbc_bra_cc.J - tbc_ket_cc.J) > Z.GetJRank() ) ) continue;
         // Important to remember. For CC channels, Tz is the magnitude of the difference of the isospin | tz1 - tz2|
         // For RankT=0, we can have <pn|pn>, <pp|nn>, <pp|pp>, <nn|nn>, (Tz_bra,Tz_ket) =>  (1,1) , (0,0)
         // For RankT=1, we can have <pn|pp>, <pn|nn>  (Tz_bra,Tz_ket) => (0,1) , (1,0)
         // For RankT=2, we can have <pn|pn>   (Tz_bra,Tz_ket) => (1,1)
         if ( not ( (tbc_bra_cc.Tz+tbc_ket_cc.Tz==Z.GetTRank()) or (std::abs(tbc_bra_cc.Tz-tbc_ket_cc.Tz)==Z.GetTRank()) ) ) continue;

         ybras.push_back(ich_bra);
         ykets.push_back(ich_ket);
         Z_bar[{ich_bra,ich_ket}] = arma::mat(n_rows,n_cols);
//         Z_bar[{ich_bra,ich_ket}] = arma::mat();
     }
   }
   int counter = ybras.size();


   X.profiler.timer["Allocate Z_bar_tensor"] += omp_get_wtime() - t_start;




   t_start = omp_get_wtime();


// BEGIN OLD WAY
   if ( Z.GetJRank()>0 )
   {
//      std::cout << "  in  " << __func__ << "  doing it the old way. Counter = " << counter << std::endl;
      #ifndef OPENBLAS_NOUSEOMP
      #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.GetJRank()))
      #endif
      for(int i=0;i<counter;++i)
      {
//         std::cout << "      i = " << i << std::endl;
         index_t ch_bra_cc = ybras[i];
         index_t ch_ket_cc = ykets[i];
   ////      const auto plookup = pandya_lookup.find({(int)ch_bra_cc,(int)ch_ket_cc});
   //      const auto plookup = pandya_lookup.find({ch_bra_cc,ch_ket_cc});
   //      if ( plookup == pandya_lookup.end() or plookup->second[0].size()<1 )
   //      {
   //       continue;
   //      }
   
         const auto& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
         const auto& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
         int Jbra = tbc_bra_cc.J;
         int Jket = tbc_ket_cc.J;
   
         arma::mat YJ1J2;
         arma::mat YJ2J1;
         const auto& XJ1 = Xt_bar_ph[ch_bra_cc];
         const auto& XJ2 = Xt_bar_ph[ch_ket_cc];
   
         arma::uvec kets_ph = arma::join_cols( tbc_ket_cc.GetKetIndex_hh(), tbc_ket_cc.GetKetIndex_ph() );
         arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
   
         DoTensorPandyaTransformation_SingleChannel(Y, YJ1J2, ch_bra_cc, ch_ket_cc);
         if (ch_bra_cc==ch_ket_cc)
         {
            YJ2J1 = YJ1J2;
         }
         else
         {
            DoTensorPandyaTransformation_SingleChannel(Y, YJ2J1, ch_ket_cc, ch_bra_cc);
         }
   
         int flipphaseY = hY * Z.modelspace->phase( Jbra - Jket ) ;
         // construct a matrix of phases (-1)^{k+j+p+h} used below to generate X_phkj for k>j
         arma::mat PhaseMatXJ2( tbc_ket_cc.GetNumberKets(), kets_ph.size(), arma::fill::ones) ;
         arma::mat PhaseMatYJ1J2( bras_ph.size(), tbc_ket_cc.GetNumberKets(), arma::fill::ones) ;
         for ( index_t iket=0;iket<(index_t)tbc_ket_cc.GetNumberKets();iket++)
         {
           const Ket& ket = tbc_ket_cc.GetKet(iket);
           if ( Z.modelspace->phase((ket.op->j2+ket.oq->j2)/2)<0)
           {
              PhaseMatXJ2.row(iket) *=-1;
              PhaseMatYJ1J2.col(iket) *=-1;
           }
         }
         for (index_t iph=0;iph<kets_ph.size();iph++)
         {
           const Ket& ket_ph = tbc_ket_cc.GetKet( kets_ph[iph] );
           if ( Z.modelspace->phase((ket_ph.op->j2+ket_ph.oq->j2)/2)<0)   PhaseMatXJ2.col(iph) *=-1;
         }
         for (index_t iph=0;iph<bras_ph.size();iph++)
         {
           const Ket& bra_ph = tbc_bra_cc.GetKet( bras_ph[iph] );
           if ( Z.modelspace->phase((bra_ph.op->j2+bra_ph.oq->j2)/2)<0)   PhaseMatYJ1J2.row(iph) *=-1;
         }
         PhaseMatYJ1J2 *= flipphaseY;
   
   
   
   //                J2                       J1         J2                       J2          J2
   //             k<=j     k>=j                hp  -ph    hp   ph                 k<=j       k<=j
   //      J1   [       |       ]       J1   [           |          ]         [hp        |ph        ]
   //     i<=j  [  Zbar | Zbar  ]  =   i<=j  [   Xbar    | -Ybar    ]   * J1  [   Ybar   |   Ybar'  ]
   //           [       |       ]            [           |          ]         [ph        |hp        ]      where Ybar'_phkj = Ybar_hpkj * (-1)^{p+h+k+j}*(-1)^{J1-J2}*hY
   //                                                                         [----------|----------]       and
   //                                                                     J2  [hp        |ph        ]            Xbar'_phkj = Xbar_hpkj * (-1)^{p+h+k+j}*hX
   //                                                                         [   Xbar   |   Xbar'  ]
   //                                                                         [-ph       |-hp       ]
   //    
   //
         int halfncx2 = XJ2.n_cols/2;
         int halfnry12 = YJ1J2.n_rows/2;
   
         arma::mat Mleft = join_horiz( XJ1,  -flipphaseY * YJ2J1.t() );
         arma::mat Mright = join_vert( join_horiz( YJ1J2 ,  join_vert( YJ1J2.tail_rows(halfnry12)%PhaseMatYJ1J2 ,
                                                                       YJ1J2.head_rows(halfnry12)%PhaseMatYJ1J2  )   ), 
                                     hX*join_vert( XJ2,    join_horiz(   XJ2.tail_cols(halfncx2)%PhaseMatXJ2 ,
                                                                         XJ2.head_cols(halfncx2)%PhaseMatXJ2     )   ).t() );
   
         auto& Zmat = Z_bar.at({ch_bra_cc,ch_ket_cc});
   
         Zmat = Mleft * Mright;
//         std::cout << "   ... line " << __LINE__ << "  " << ch_bra_cc << " " << ch_ket_cc << "  |Z| = " << arma::norm(Zmat,"fro") << std::endl;
   
      }


   }
   else  // faster, more memory hungry way
   {
//      std::cout << "  in  " << __func__ << "  doing it the new way" << std::endl;

      std::deque<arma::mat> YJ1J2_list(counter);
      std::deque<arma::mat> YJ2J1_list(counter);
   
   //   #ifndef OPENBLAS_NOUSEOMP
      #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.GetJRank()))
   //   #endif
      for(int i=0;i<counter;++i)
      {
         index_t ch_bra_cc = ybras[i];
         index_t ch_ket_cc = ykets[i];
   
         const auto& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
         const auto& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
//         int Jbra = tbc_bra_cc.J;
//         int Jket = tbc_ket_cc.J;
   
   //      arma::mat YJ1J2;
   //      arma::mat YJ2J1;
         arma::mat& YJ1J2 = YJ1J2_list[i];
         arma::mat& YJ2J1 = YJ2J1_list[i];
//         const auto& XJ1 = Xt_bar_ph[ch_bra_cc];
//         const auto& XJ2 = Xt_bar_ph[ch_ket_cc];
   
         arma::uvec kets_ph = arma::join_cols( tbc_ket_cc.GetKetIndex_hh(), tbc_ket_cc.GetKetIndex_ph() );
         arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
   
         DoTensorPandyaTransformation_SingleChannel(Y, YJ1J2, ch_bra_cc, ch_ket_cc);
         if (ch_bra_cc==ch_ket_cc)
         {
   //         YJ2J1 = YJ1J2;
             // Dont do nothing..
         }
         else
         {
            DoTensorPandyaTransformation_SingleChannel(Y, YJ2J1, ch_ket_cc, ch_bra_cc);
         }
      }
   
   
      X.profiler.timer["DoTensorPandyaTransformationY"] += omp_get_wtime() - t_start;
   
   
      t_start = omp_get_wtime();
   
      #ifndef OPENBLAS_NOUSEOMP
      #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.GetJRank()))
      #endif
      for(int i=0;i<counter;++i)
      {
         index_t ch_bra_cc = ybras[i];
         index_t ch_ket_cc = ykets[i];
   
         const auto& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
         const auto& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
         int Jbra = tbc_bra_cc.J;
         int Jket = tbc_ket_cc.J;
   
   //      arma::mat YJ1J2;
   //      arma::mat YJ2J1;
         arma::mat& YJ1J2 = YJ1J2_list[i];
         arma::mat& YJ2J1 = ( ch_bra_cc==ch_ket_cc) ? YJ1J2_list[i]  : YJ2J1_list[i];
   
   //      arma::mat& YJ2J1 = YJ2J1_list[i];
   
         const auto& XJ1 = Xt_bar_ph[ch_bra_cc];
         const auto& XJ2 = Xt_bar_ph[ch_ket_cc];
   
         arma::uvec kets_ph = arma::join_cols( tbc_ket_cc.GetKetIndex_hh(), tbc_ket_cc.GetKetIndex_ph() );
         arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
   
   
         int flipphaseY = hY * Z.modelspace->phase( Jbra - Jket ) ;
         // construct a matrix of phases (-1)^{k+j+p+h} used below to generate X_phkj for k>j
         arma::mat PhaseMatXJ2( tbc_ket_cc.GetNumberKets(), kets_ph.size(), arma::fill::ones) ;
         arma::mat PhaseMatYJ1J2( bras_ph.size(), tbc_ket_cc.GetNumberKets(), arma::fill::ones) ;
         for ( index_t iket=0;iket<(index_t)tbc_ket_cc.GetNumberKets();iket++)
         {
           const Ket& ket = tbc_ket_cc.GetKet(iket);
           if ( Z.modelspace->phase((ket.op->j2+ket.oq->j2)/2)<0)
           {
              PhaseMatXJ2.row(iket) *=-1;
              PhaseMatYJ1J2.col(iket) *=-1;
           }
         }
         for (index_t iph=0;iph<kets_ph.size();iph++)
         {
           const Ket& ket_ph = tbc_ket_cc.GetKet( kets_ph[iph] );
           if ( Z.modelspace->phase((ket_ph.op->j2+ket_ph.oq->j2)/2)<0)   PhaseMatXJ2.col(iph) *=-1;
         }
         for (index_t iph=0;iph<bras_ph.size();iph++)
         {
           const Ket& bra_ph = tbc_bra_cc.GetKet( bras_ph[iph] );
           if ( Z.modelspace->phase((bra_ph.op->j2+bra_ph.oq->j2)/2)<0)   PhaseMatYJ1J2.row(iph) *=-1;
         }
         PhaseMatYJ1J2 *= flipphaseY;
   
   
   
   //                J2                       J1         J2                       J2          J2
   //             k<=j     k>=j                hp  -ph    hp   ph                 k<=j       k<=j
   //      J1   [       |       ]       J1   [           |          ]         [hp        |ph        ]
   //     i<=j  [  Zbar | Zbar  ]  =   i<=j  [   Xbar    | -Ybar    ]   * J1  [   Ybar   |   Ybar'  ]
   //           [       |       ]            [           |          ]         [ph        |hp        ]      where Ybar'_phkj = Ybar_hpkj * (-1)^{p+h+k+j}*(-1)^{J1-J2}*hY
   //                                                                         [----------|----------]       and
   //                                                                     J2  [hp        |ph        ]            Xbar'_phkj = Xbar_hpkj * (-1)^{p+h+k+j}*hX
   //                                                                         [   Xbar   |   Xbar'  ]
   //                                                                         [-ph       |-hp       ]
   //    
   //
         int halfncx2 = XJ2.n_cols/2;
         int halfnry12 = YJ1J2.n_rows/2;
   
         arma::mat Mleft = join_horiz( XJ1,  -flipphaseY * YJ2J1.t() );
         arma::mat Mright = join_vert( join_horiz( YJ1J2 ,  join_vert( YJ1J2.tail_rows(halfnry12)%PhaseMatYJ1J2 ,
                                                                       YJ1J2.head_rows(halfnry12)%PhaseMatYJ1J2  )   ), 
                                     hX*join_vert( XJ2,    join_horiz(   XJ2.tail_cols(halfncx2)%PhaseMatXJ2 ,
                                                                         XJ2.head_cols(halfncx2)%PhaseMatXJ2     )   ).t() );
   
   
         Z_bar.at({ch_bra_cc,ch_ket_cc}) = Mleft * Mright;
   
      }
      X.profiler.timer["Build Z_bar_tensor"] += omp_get_wtime() - t_start;



  } //else J=0


   t_start = omp_get_wtime();
   AddInverseTensorPandyaTransformation(Z, Z_bar);

   X.profiler.timer["InverseTensorPandyaTransformation"] += omp_get_wtime() - t_start;

   Z.modelspace->tensor_transform_first_pass.at( Z.GetJRank() ) = false;

}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Begin Scalar-Dagger Commutator Terms////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///  In the below formulas, the number in parentheses indicates the number of external legs of a diagram, rather
///  than the usual particle rank. In this notation, a one-body operator has 2 external legs and so would be indicated
///  with a (2). The (1) then means a single creation/annihilation operator, and (3) means two creators and one annihilator.





//*****************************************************************************************
// [X^(2), Y^(1)]^(1)
//
//      |i                 i|   
//      |                   |   
//     (X)         ---      |  (Y)
//       a\                 |  /a
//         (Y)             (X)/
//           
//
// Sum_ia X_ia Y_a^{lambda}
// Adapted from
//    void Operator::comm111ss( const Operator & X, const Operator& Y) 
// This could be modified to be more efficent since we only need one column from Y to go to one column in Z.
// But this is so far from being the bottleneck that it makes more sense to be clear.
void comm211sd( const Operator& X, const Operator& Y, Operator& Z )
{
   Z.OneBody += X.OneBody*Y.OneBody;
}



//*****************************************************************************************
// [X^(2),Y^(3)]^(1)
//         |i             |i
//         |              |
//  (X)_ b |         a _(Y)
//    \_\  |   ---    /_/
//    a  (Y)         (X) b
//
// Sum_ab Sum_J (na n`b - n`a nb) X_ab Y_bia  hat(J)^2/hat(lambda)^2
// or
// Sum_ab Sum_J (2J+1)/(2lambda+1)  (na n`b) ( X_ab Y_bia - X_ba Y_aib )
//
// Adapted from
//    void Operator::comm121ss( const Operator& X, const Operator& Y) 
void comm231sd( const Operator& X, const Operator& Y, Operator& Z) 
{
   double t_start = omp_get_wtime();
   if (Y.legs<3) return;
//   index_t norbits = Z.modelspace->GetNumberOrbits();
   index_t Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
//   for (index_t i=0;i<norbits;++i)
//   #pragma omp parallel for 
   for (auto i : Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
   {
          for (auto& a : Z.modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = Z.modelspace->GetOrbit(a);
//             for (index_t b=0; b<norbits; ++b)
             for (auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2} ))
             {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ); // despite what the name suggests, nanb is n_a * (1-n_b).
                if (std::abs(nanb)<ModelSpace::OCC_CUT) continue;
                int Jmin = std::abs(oa.j2 - oQ.j2)/2;
                int Jmax = (oa.j2+oQ.j2)/2;
                double Ymon_bia = 0;
                double Ymon_aib = 0;
                for (int J=Jmin;J<=Jmax;J++)
                {
                  Ymon_bia += (2*J+1.0) / (oQ.j2+1.0) * Y.ThreeLeg.GetME_J(J,b,i,a);
                  Ymon_aib += (2*J+1.0) / (oQ.j2+1.0) * Y.ThreeLeg.GetME_J(J,a,i,b);
                }
//                Z.OneBody(i,Q) += nanb * X.OneBody(a,b) * Ymon_bia - X.OneBody(b,a) * Ymon_aib;
                Z.OneBody(i,0) += nanb * X.OneBody(a,b) * Ymon_bia - X.OneBody(b,a) * Ymon_aib;
//                  Z.OneBody(i,Q) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,Q) ;  // Is this still the right way to do this? Do we need to worry about normalization? (It looks ok).
//                  Z.OneBody(i,Q) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,Q) ;  // GetTBMEmonopole returns unnormalized TBME summed over J times (2J+1)/((2ji+1)*(2jj+1))
                
             }
          }
   }
   Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}





//*****************************************************************************************
// [X^(2),Y^(3)]^(3)  and [X^(4),Y^(1)]^(3), combined
//
//     |   |  /       |  | /             \ /         \ |  (Y)
//    (X)  | /   __   | (Y)              (X)    __    \| / 
//      \  |/         | /        and     / \          (X)
//       (Y)         (X)                /  (Y)         | 
//
//  ZJ_ijkQ = sum_a  X_ia YJ_ajkQ  + X_ja YJ_iakQ + X_ak YJ_ijaQ + XJ_ijka Y_aQ
//
//  I adapted this, but then gave up and brute forced it because the intermediate matrix
//  stuff was too obfuscated, and I'm the one who wrote it... -SRS 14/04/2018
// Adapted from
//  void Operator::comm122ss( const Operator& X, const Operator& Y ) 
void comm413_233sd( const Operator& X, const Operator& Y, Operator& Z) 
{
   double t_start = omp_get_wtime();
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;

   index_t Qorbit = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Qorbit);
//   int norb = Z.modelspace->GetNumberOrbits();

   int n_nonzero = Z.modelspace->SortedTwoBodyChannels.size();
//   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         double norm_ij = (i==j) ? PhysConst::INVSQRT2 : 1.0;
//         for ( int k=0; k<norb; ++k )
         for ( auto k : Z.modelspace->all_orbits )
         {
           Orbit& ok = Z.modelspace->GetOrbit(k);
//           if (not tbc.CheckChannel_ket(&ok,&oQ)) continue;
//           if ( ( (ok.tz2+oQ.tz2)!=2*tbc.Tz ) or ( (ok.l+oQ.l)%2!=tbc.parity ) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue;
           if ( ( (ok.tz2-oQ.tz2)!=2*tbc.Tz ) or ( (ok.l+oQ.l)%2!=tbc.parity ) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue;

//           double norm_kQ = (k==Qorbit) ? PhysConst::INVSQRT2 : 1.0;
           double norm_kQ = 1.0;
           double cijk = 0;

//           std::cout << "INNER LOOPs ijk = " << i << " " << j << " " << k << " Z has " << Z.GetNumberLegs() << " legs " << " and there are " << Z.ThreeLeg.MatEl.size() << " 3 leg channels " << std::endl;
           // The first three loops are the [2,3]->3 bit
           for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
           {
//             cijk += X1(i,a) * Y.TwoBody.GetTBME_norm(ch,ch,a,j,k,Qorbit);
//             cijk += X1(i,a) * Y.TwoBody.GetTBME(ch,ch,a,j,k,Qorbit);
             cijk += X1(i,a) * Y.ThreeLeg.GetME(ch,a,j,k);
           }
//           std::cout << "check 1. size of X.OneBodyChannels is " << X.OneBodyChannels.Size() << std::endl;
           for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
           {
//             cijk += X1(j,a) * Y.TwoBody.GetTBME_norm(ch,ch,i,a,k,Qorbit);
//             cijk += X1(j,a) * Y.TwoBody.GetTBME(ch,ch,i,a,k,Qorbit);
             cijk += X1(j,a) * Y.ThreeLeg.GetME(ch,i,a,k);
           }
//           std::cout << "check 2" << std::endl;
           for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
           {
//             cijk += X1(a,k) * Y.TwoBody.GetTBME_norm(ch,ch,i,j,a,Qorbit); // should this have a minus sign?
//             cijk -= X1(a,k) * Y.TwoBody.GetTBME_norm(ch,ch,i,j,a,Qorbit); 
//             cijk -= X1(a,k) * Y.TwoBody.GetTBME(ch,ch,i,j,a,Qorbit); 
             cijk -=  Y.ThreeLeg.GetME(ch,i,j,a) * X1(a,k); 
           }

//           std::cout << " DOING 413 bit "  << std::endl;
           // This here loop is the [4,1]->3 bit
           for ( int a : Y.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
           {
//             cijk += Y1(a,Qorbit) * X.TwoBody.GetTBME_norm(ch,ch,i,j,k,a);   // This determines the normalization for the (adagger adagger a) term.
//             cijk +=  X.TwoBody.GetTBME(ch,ch,i,j,k,a) * Y1(a,Qorbit) ;   // This determines the normalization for the (adagger adagger a) term.
             cijk +=  X.TwoBody.GetTBME(ch,ch,i,j,k,a) * Y1(a,0) ;   // This determines the normalization for the (adagger adagger a) term.
//             if (i==0 and j==8 and k==0)
//             {
//                std::cout << "   " << __func__ << "  a = " << a << "  Xijka, Ya = " << X.TwoBody.GetTBME(ch,ch,i,j,k,a) << " , " << Y1(a,Qorbit) << "  cijk = " << cijk << std::endl;
//             }

           }


           cijk *= (norm_ij * norm_kQ);  // We normalize like a scalar TBME so that the setter/getters make sense. We'll fix the wrong |kQ> normalization when writing to file.

           Z.ThreeLeg.AddToME(ch,i,j,k, cijk);   // AddToME  assumes a normalized matrix element.
         }
      }
   }
   Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}








void ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, ThreeLegME& Mpp, ThreeLegME& Mhh)
{
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

      auto& LHS = X.TwoBody.GetMatrix(ch,ch);
      auto& RHS = Y.ThreeLeg.GetMatrix(ch);

      auto& Matrixpp = Mpp.GetMatrix(ch);
      auto& Matrixhh = Mhh.GetMatrix(ch);

      auto& kets_pp = tbc.GetKetIndex_pp();
      auto& kets_hh = tbc.GetKetIndex_hh();
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb    = tbc.Ket_occ_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;


   } //for ch

}


/// Since comm443_pp_hhsd() and comm431sd() both require the construction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
/// We do the same thing in the standard scalar-scalar commutators.
void comm433_pp_hh_431sd( const Operator& X, const Operator& Y, Operator& Z )  
{

   double t = omp_get_wtime();

//   static TwoBodyME Mpp = Z.TwoBody;
//   static TwoBodyME Mhh = Z.TwoBody;
   ThreeLegME Mpp = Z.ThreeLeg;
   ThreeLegME Mhh = Z.ThreeLeg;
   Mpp *= 0;
   Mhh *= 0;

   ConstructDaggerMpp_Mhh( X, Y, Z, Mpp, Mhh);

   Z.ThreeLeg += Mpp;
   Z.ThreeLeg -= Mhh;

   Z.profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;

   // Now, the one body part
   t = omp_get_wtime();
   int Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
//   auto ilist = Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2});
   std::vector<index_t> ilist(Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}).begin(), Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}).end());
   int ni = ilist.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int i_ind=0; i_ind<ni; ++i_ind)
   {
      int i = ilist[i_ind];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double cijJ = 0;
      for (int ch=0;ch<Z.nChannels;++ch)
      {
         TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
         double Jfactor = (2*tbc.J+1.0);
         // Sum c over holes and include the nbar_a * nbar_b terms
         for (auto& c : Z.modelspace->holes)
         {
            Orbit& oc = Z.modelspace->GetOrbit(c);
//            cijJ += Jfactor * oc.occ * Mpp.GetTBME(ch,c,i,c,Q);     // We use the GetTBME, which returns an unnormalized matrix element, as required.
//            cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,Q); 
            cijJ += Jfactor * oc.occ * Mpp.GetME(ch,c,i,c);     // We use the GetTBME, which returns an unnormalized matrix element, as required.
            cijJ += Jfactor * (1-oc.occ) * Mhh.GetME(ch,c,i,c); 
         }
         // Sum c over particles and include the n_a * n_b terms
         for (auto& c : Z.modelspace->particles)
         {
//            cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,Q);
            cijJ += Jfactor * Mhh.GetME(ch,c,i,c);
         }
      }
      Z.OneBody(i,0) += cijJ /(oi.j2+1.0);   // The factor of 1/2 in the formula is absorbed by the fact that the mat-mult only sums a<=b.
//      Z.OneBody(i,Q) += cijJ /(oi.j2+1.0);   // The factor of 1/2 in the formula is absorbed by the fact that the mat-mult only sums a<=b.
   } // for i
   Z.profiler.timer["pphh 413sd"] += omp_get_wtime() - t;
}




//*****************************************************************************************
// [X^(4),Y^(3)]^(3)]  ph piece, the slow way
//                                             
//   |           |       |           |           
//   |      _(Y)_|       |_(X)_      |            
//   |     /\            |    /\     |
//   |    (  )      __   |   (  )    |            
//   |_(X)_\/            |    \/_(Y)_|            
//   |                   |            
//
//  Straightfoward and very slow implementation, only used for unit testing the
//  faster mat-mult implementation.  The two implementations agree as of Nov 2, 2018 - SRS 
//  Benchmarked against the m-scheme version in UnitTest.cc, and it agrees (although
//  some changes were required compared to the impoementation tested in Nov 2018).
//  So this again agrees with the more efficient way, and with the m-scheme verion as of Sep 30 2019 - SRS
void comm433sd_ph_dumbway( const Operator& X, const Operator& Y, Operator& Z)
{

   double t_start = omp_get_wtime();
//   int norb = Z.modelspace->GetNumberOrbits();
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   index_t Q = Y.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
   double jQ = 0.5*oQ.j2;
   auto kets_ph = Z.modelspace->KetIndex_ph;
   for ( auto khh : Z.modelspace->KetIndex_hh ) // ph kets needs to include hh list in case there are fractionally occupied states
   {
     kets_ph.push_back(khh);
   }

   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      size_t nkets = tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets; ibra++)
      {
        auto& bra = tbc.GetKet(ibra);
        std::vector<index_t> i_cases = { bra.p, bra.q };
        std::vector<index_t> j_cases = { bra.q, bra.p };
        std::vector<int> ijsign_cases = { +1, bra.Phase(J)};  // This accounts for the 1 - (-1)^(ji +jj-J) Pij  factor.
        double norm_ij = (bra.p==bra.q) ? PhysConst::INVSQRT2 : 1.0;

          for ( auto k : Z.modelspace->all_orbits )
          {
            Orbit& ok = Z.modelspace->GetOrbit(k);
//            if ( not tbc.CheckChannel_ket(&ok,&oQ) ) continue;
//            if ( ( (ok.tz2+oQ.tz2)!=2*tbc.Tz) or ( (ok.l+oQ.l)%2!=tbc.parity) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue; // if |kQ> doesn't live in this channel, move along.
            if ( ( (ok.tz2-oQ.tz2)!=2*tbc.Tz) or ( (ok.l+oQ.l)%2!=tbc.parity) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue; // if |kQ> doesn't live in this channel, move along.
            double zijk = 0.;
//            double norm_kQ = (k==Q) ? PhysConst::INVSQRT2 : 1.0;
            double norm_kQ =  1.0;
            double jk = 0.5* ok.j2;


        for (int ijcase=0; ijcase<=1; ijcase++)
        {
//         if (ijcase==1) continue;
          index_t i = i_cases[ijcase];
          index_t j = j_cases[ijcase];
//          std::cout << "~~~~ i,j = " << i << " " << j << std::endl;
          int ijsign = ijsign_cases[ijcase];
          Orbit& oi = Z.modelspace->GetOrbit(i);
          Orbit& oj = Z.modelspace->GetOrbit(j);
          double ji = 0.5* oi.j2;
          double jj = 0.5* oj.j2;
//          for ( size_t k=0; k<norb; k++)



            int Jprime_min = std::max(  std::abs(ji-jQ), std::abs(jj-jk) );
            int Jprime_max = std::min(  ji+jQ, jj+jk  );
            for (int Jprime=Jprime_min; Jprime<=Jprime_max; ++Jprime)
            {
              double sixjprime =  Z.modelspace->GetSixJ(ji,jj,J,jk,jQ,Jprime);
              if (std::abs(sixjprime) < 1e-8) continue;

              double XYprod = 0;

              for (auto& iketab : kets_ph )
              {
                auto& ketab = Z.modelspace->GetKet(iketab);
                  if (  (ketab.op->j2+ketab.oq->j2)<2*Jprime  or  std::abs(ketab.op->j2-ketab.oq->j2)>2*Jprime ) continue;
                std::vector<index_t> ab_cases = { ketab.p, ketab.q };
                for (int abcase = 0; abcase<=1; abcase++)
                {
                  index_t a  = ab_cases[abcase];
                  index_t b  = ab_cases[1-abcase];

                  Orbit& oa = Z.modelspace->GetOrbit(a);
                  Orbit& ob = Z.modelspace->GetOrbit(b);
                  double ja = 0.5*oa.j2;
                  double jb = 0.5*ob.j2;
                  double nanb = oa.occ - ob.occ;

//                  if ( not (a==0 and b==10) or (a==10 and b==0) ) continue;
//                  if ( not ((a==1 and b==8) or (a==8 and b==1)) ) continue;

                  // we can probably also check triangle for (a,b,Jprime)


                  int JA_min = std::max(  std::abs(ja-jj), std::abs(jk-jb)  );
                  int JA_max = std::min( ja+jj,  jk+jb );

                  int JB_min = std::max(  std::abs(ja-jQ), std::abs(ji-jb)  );
                  int JB_max = std::min( ja+jQ,  ji+jb );
                  double matelX = 0;
                  double matelY = 0;
                  for (int JA = JA_min; JA<=JA_max; ++JA)
                  {
                    double sixjA = Z.modelspace->GetSixJ( ja, jb, Jprime, jk, jj, JA );
                    matelX -= (2*JA+1) * sixjA * X.TwoBody.GetTBME_J(JA,a,j,k,b);   // GetTBME_J returns an un-normalized matrix element
                  }
                  for (int JB = JB_min; JB<=JB_max; ++JB)
                  {
                    double sixjB = Z.modelspace->GetSixJ( ja, jb, Jprime, ji, jQ, JB );
//                    matelY -= (2*JB+1) * sixjB * Y.TwoBody.GetTBME_J(JB,i,b,a,Q);   // GetTBME_J returns an un-normalized matrix element
                    matelY -= (2*JB+1) * sixjB * Y.ThreeLeg.GetME_J(JB,i,b,a);   // GetTBME_J returns an un-normalized matrix element
                  }
                  XYprod += nanb * matelX * matelY;
//                  if ( ( (i==0 and j==1) or (i==1 and j==0)) and k==1  and std::abs(nanb)>0)
//                  if ( ( (i==0 and j==1) or (i==1 and j==0)) and k==1 and ((a==0 and b==10) or (b==0 and a==10)) )
//                  if ( ( (i==0 and j==1) or (i==1 and j==0)) and k==1 and ((a==1 and b==8) or (b==1 and a==8)) )
//                  {
//                    std::cout << "   J = " << J << "  a,b = " << a << " " << b << "  Jprime = " << Jprime << "  nanb " << nanb 
//                              << "  Xbar_abkj Ybar_iQab " << matelX << " " << matelY << " => XYprod = " << XYprod << std::endl;
//                  }
                }
              } // loop over ab cases
               zijk -= ijsign * (2*Jprime+1) * sixjprime * XYprod;
//                  if (( (i==0 and j==1) or (i==1 and j==0)) and k==1)
//                  {
//                    std::cout << "i,j = " << i << " " << j << " ijsign(J=" << J << ") = " << ijsign << "  2Jprime+1 " << 2*Jprime+1 << "  sixjprime " << sixjprime << "  XYprod " << XYprod << "  => zijk " << zijk
//                               <<  std::endl;
//                  }
            }  // loop over ab kets
 
//                  if (( (i==0 and j==1) or (i==1 and j==0)) and k==1)
//                  {
//                    std::cout << " zijk = " << zijk <<  " norm_ij norm_kQ = " << norm_ij << " " << norm_kQ << "    before setting, the ME is " << Z.ThreeLeg.GetME_J(J,i,j,k) << std::endl << std::endl;
//                  }

        } // loop over ij cases

//            Z.TwoBody.AddToTBME_J(J,i,j,k,Q,  zijk * norm_ij*norm_kQ);   // AddToTBME_J assumes a normalized matrix element.
            Z.ThreeLeg.AddToME_J(J,bra.p,bra.q,k,  zijk * norm_ij*norm_kQ);   // AddToTBME_J assumes a normalized matrix element.
          }  // loop over k 
      } // for ibra
   } // for ich
   Z.profiler.timer[__func__] += omp_get_wtime() - t_start;

}




//*****************************************************************************************
// [X^(4),Y^(3)]^(3)]  ph piece
//                                             
//   |           |       |           |           
//   |      _(Y)_|       |_(X)_      |            
//   |     /\            |    /\     |
//   |    (  )      __   |   (  )    |            
//   |_(X)_\/            |    \/_(Y)_|            
//   |                   |            
//
// For formula and details of implementation, see Operator::comm222_phss. The only change
// made here to accommodate a dagger operator is to replace XY - YX  with -YX. This
// ensures that we don't include contributions from the Qspace orbit to X.
// Adapted from
//    void Operator::comm222_phss( const Operator& X, const Operator& Y ) 


void comm433sd_ph( const Operator& X, const Operator& Y, Operator& Z)  
{

   double t_start = omp_get_wtime();

   int hx = X.IsHermitian() ? 1 : -1; // assuming X is either hermitian or antihermitian. I'm sure this will come back to bite me some day.

   // Construct the intermediate matrix Z_bar and fill it with zeros.
//   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup( Z.GetJRank(), Z.GetTRank(), Z.GetParity() );
   size_t nch = Z.modelspace->SortedTwoBodyChannels_CC.size();
   size_t norb = Z.modelspace->GetNumberOrbits();
   std::deque<arma::mat> Z_bar (Z.nChannels );
   std::vector<bool> lookup_empty(Z.nChannels,true);
   for (size_t ich=0;ich<nch;++ich)
   {
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      index_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros( norb, 2*nKets_cc );  // Z_iQ`kj`   we want both orderings of k and j, but Q is fixed.
//      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
       lookup_empty[ich] = false;
   }


   // loop over channels for Z_bar
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   #endif
   for (size_t ich=0; ich<nch; ++ich )
   {
      if (lookup_empty.at(ich)) continue;
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC.at(ich);
      const TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      size_t nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

      arma::mat Y_bar_ph;
      arma::mat X_bar_ph;

      DoPandyaTransformation_SingleChannel_Dagger(Y,Y_bar_ph,ch);      // Generate  YbarJ'_iQ`ab` * (na-nb)      for a<=b and a>b
      DoPandyaTransformation_SingleChannel(X,X_bar_ph,ch,"normal");  // Generate XbarJ'_ab`kj`                  for a<=b and a>b
      auto& Zbar_ch = Z_bar.at(ch); // Marginally useful aliasing to avoid lookups...

     // Zbar will be
     // ZbarJ'_iQ`kj` = Sum_ab   (na-nb) <iQ`|YbarJ'|ab`> <ab`|XbarJ'|kj`>
     // So then we can perform an inverse Pandya transform (which is really just another Pandya transform) to get
     // ZJ_ijkQ = - Sum_J' (2J'+1) { ji  jj  J  }  ZbarJ'_iQ`kj`
     //                            { jk  jQ  J' }


     // The shape of Xt_bar_ph is  ( rows, columns) = ( 2*nph_kets,  nKets_cc)  =>    [     Xbar    ]  |  ab`, element of nph_kets, which has ph` and hp`
     //                                                                               [             ]  v
     //                                                                                    -> 
     //                                                                                    kj` is an element of kets_cc
     //
     // The shape of Y_bar_ph is   (rows, columns) = (norb, 2*nph_kets)         =>    [     Ybar    ]  | iQ`  runs over all orbits i, with Q fixed.
     //                                                                               [             ]  v
     //                                                                                    -> 
     //                                                                                    ab` is an element of nph_kets, here we have ph` and hp`
     //
     // So we should multiply Ybar * Xbar to sum over ab` and get out Zbar      =>      [     Zbar    ]  |  iQ`
     //                                                                                 [             ]  v
     //                                                                                      -> 
     //                                                                                      kj`
     //
     // Really, we need  <iQ`|Zbar|kj`> and <iQ`|Zbar|jk`> 



      if (Y_bar_ph.size()<1 or X_bar_ph.size()<1)   // for an armadillo matrix, .size()  gives the total number of elements, i.e. n_rows * n_cols
      {
//        Zbar_ch = arma::zeros( norb, 2*nKets_cc );  // This seems unnecessary since we initialized things earlier... try getting rid of it and see if things break.
        continue;
      }

      // get the phases for taking the transpose
      arma::mat PhaseMat(nKets_cc, nKets_cc, arma::fill::ones );
      for (index_t iket=0;iket<nKets_cc;iket++)
      {
         const Ket& ket = tbc_cc.GetKet(iket);
         if ( Z.modelspace->phase( (ket.op->j2 + ket.oq->j2)/2 ) > 0) continue;
         PhaseMat.col( iket ) *= -1;
         PhaseMat.row( iket ) *= -1;
      }
      arma::uvec phkets = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
      auto PhaseMatX = PhaseMat.rows(phkets) * hx;

//                                           [      |     ]
//     create full Y matrix from the half:   [  Yhp | Y'ph]   where the prime indicates multiplication by (-1)^(i+j+k+l) h_x
//                                           [      |     ]   Flipping hp <-> ph and multiplying by the phase is equivalent to
//                                           [  Yph | Y'hp]   having kets |kj> with k>j.
      Zbar_ch =  Y_bar_ph * join_horiz(X_bar_ph, join_vert(   X_bar_ph.tail_rows(nph_kets)%PhaseMatX,
                                                              X_bar_ph.head_rows(nph_kets)%PhaseMatX) );

   }

   Z.profiler.timer["Build Z_bar_Dagger"] += omp_get_wtime() - t_start;
//   std::cout << "Done building Z_bar_dagger" << std::endl;

   // Perform inverse Pandya transform on Z_bar to get Z
   t_start = omp_get_wtime();
   AddInversePandyaTransformation_Dagger(Z_bar, Z);

//   Z.modelspace->scalar_transform_first_pass = false;
   Z.profiler.timer["InversePandyaTransformation_Dagger"] += omp_get_wtime() - t_start;

   Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}






//**************************************************************************
// DAGGER VARIETY
//
//  X^J_iQ`ab` = - sum_J' { i Q J } (2J'+1) X^J'_ibaQ
//                        { a b J'}
/// A modification of the Pandya transform for the case of a dagger operator.
/// We embed the dagger operator as a scalar with a ficticious fourth index Q,
/// which is outside the Hilbert space. So our relevant transformation is
/// \f[
///  \bar{X}^{J}_{i\bar{Q}a\bar{b}} = - \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i  &  j_Q  &  J \\
///  j_a  &  j_b  &  J' \\
///  \end{array} \right\}
///  X^{J}_{ibaQ}
/// \f]
/// where the overbar indicates time-reversed orbits.
/// Since we use this for the commutator expression, the k and j orbits will be a and b,
/// i.e. the orbits we sum over with a factor (na-nb). We include the (na-nb) factor
/// in the output matrix.
void DoPandyaTransformation_SingleChannel_Dagger(const Operator& Z, arma::mat& TwoBody_CC_ph, int ch_cc)
{
   TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
//   int nKets_cc = tbc_cc.GetNumberKets();
   size_t norb = Z.modelspace->GetNumberOrbits();
   arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
   int nph_kets = kets_ph.n_rows;
   int J_cc = tbc_cc.J;

//   TwoBody_CC_ph.zeros( nKets_cc, 2*nph_kets);
   TwoBody_CC_ph.zeros( norb, 2*nph_kets);   // factor 2 because we want both orderings ab and ba.

   index_t Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
   double jQ = oQ.j2*0.5;

   // loop over cross-coupled ph kets |ab> in this channel
   // (this is the side that gets summed over in the matrix multiplication)
   for (int iket=0; iket<nph_kets; ++iket)
   {
     Ket & ket_cc = tbc_cc.GetKet( kets_ph[iket] );
     int a = ket_cc.p;
     int b = ket_cc.q;
     Orbit & oa = Z.modelspace->GetOrbit(a);
     Orbit & ob = Z.modelspace->GetOrbit(b);
     double ja = oa.j2*0.5;
     double jb = ob.j2*0.5;

     // Here we make some lists so that we don't need to repeat code when switching (a,b) -> (b,a).
     // Instead we can just iterate over the two cases.
     std::vector<int> ab     = { a, b };
     std::vector<double>jab  = { ja, jb };
     std::vector<double>nanb = { (oa.occ-ob.occ),  (ob.occ-oa.occ) };

     // we loop over both orderings, a<b and a>b. Here, this is done by exchanging a<->b and taking a minus sign due to the (na-nb) factor.
     for ( int ab_case=0; ab_case<=1; ab_case++)
     {
       size_t ab1 = ab[ab_case];
       size_t ab2 = ab[1-ab_case];
       double jab1  = jab[ab_case];
       double jab2  = jab[1-ab_case];
       size_t indx_ab = iket + ab_case*nph_kets;

       // loop over orbits i, and check if <iQ| is in the desired channel.
       for (size_t i=0; i<norb; i++)
       {
          Orbit& oi = Z.modelspace->GetOrbit(i);
          if (not tbc_cc.CheckChannel_ket( &oi, &oQ ) ) continue;  //  <iQ|  isn't in this channel, so move along. (Note, for this check, the ordering of i,Q doesn't matter).
          size_t indx_iQ = i;   // since Q is fixed, we can label the bra <iQ| by the index i.
          double ji = oi.j2*0.5;

          int Jmin = std::max( std::abs(jab1-jQ), std::abs(ji-jab2) );
          int Jmax = std::min( jab1+jQ, ji+jab2 );
          double Xbar = 0;
          for (int J_std=Jmin; J_std<=Jmax; ++J_std)
          {
             double sixj = Z.modelspace->GetSixJ(ji,jQ,J_cc,jab1,jab2,J_std);
             if (std::abs(sixj) < 1e-8) continue;
//             double tbme = Z.TwoBody.GetTBME_J(J_std,i,ab2,ab1,Q);
             double tbme = Z.ThreeLeg.GetME_J(J_std,i,ab2,ab1);
             Xbar -= (2*J_std+1) * sixj * tbme ;
          }
          TwoBody_CC_ph(indx_iQ, indx_ab) = Xbar * nanb[ab_case] ;
       }
     }
   }
}



// Same idea as the scalar commutator version, with enough modifications to make it not wrong (I hope).
// One big difference is that the Pandya-transformed operator Zbar has its bra indices numbered
// by orbit index rather than ket index, because it's   <iQ|Zbar|kj>  and Q is fixed.
// Otherwise, it should look pretty similar.
void AddInversePandyaTransformation_Dagger( const std::deque<arma::mat>& Zbar, Operator& Z )
{
    // Do the inverse Pandya transform
   int n_nonzeroChannels = Z.modelspace->SortedTwoBodyChannels.size();
   size_t norb = Z.modelspace->GetNumberOrbits();
   size_t Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
   double jQ = 0.5 * oQ.j2;

    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
//   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      size_t ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      size_t nKets = tbc.GetNumberKets();

      // the bra is the <ij| part of Z_ijkQ
      for (size_t ibra=0; ibra<nKets; ++ibra)
      {
        Ket & bra = tbc.GetKet(ibra);

        // Two-element lists for use below, because we need to antisymmetrize with 1-Pij, including a phase factor (-1)^(ji+jj-J)
        // The code would look like the exact same thing copy-pasted twice, with i and j exchanged in the second copy.
        std::vector<size_t> ij_switcheroo = {bra.p , bra.q};
        std::vector<int> phaseij = { +1, bra.Phase(J) }; 

        double norm_ij = (bra.p==bra.q) ? PhysConst::INVSQRT2 : 1.0;

        for (size_t k=0; k<norb; k++)
        {
          Orbit & ok = Z.modelspace->GetOrbit(k);
//          if ( not tbc.CheckChannel_ket(&ok, &oQ) ) continue; // if |kQ> doesn't live in this channel, move along.
//          if ( ( (ok.tz2+oQ.tz2)!=2*tbc.Tz) or ( (ok.l+oQ.l)%2!=tbc.parity) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue; // if |kQ> doesn't live in this channel, move along.
          if ( ( (ok.tz2-oQ.tz2)!=2*tbc.Tz) or ( (ok.l+oQ.l)%2!=tbc.parity) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue; // if |kQ> doesn't live in this channel, move along.

          double jk = ok.j2/2.;
//          double norm_kQ = (k==Q) ? PhysConst::INVSQRT2  : 1.0;
          double norm_kQ = 1.0;

          double commijk = 0;  // collect the sum in commijk 

          for (int ij_case=0; ij_case<=1; ij_case++) // loop over the two cases corresponding to 1-Pij (with the appropriate phase on the exchange term)
          {
            size_t i = ij_switcheroo[ij_case];
            size_t j = ij_switcheroo[1-ij_case];
            int Pij = phaseij[ij_case];

            Orbit & oi = Z.modelspace->GetOrbit(i);
            Orbit & oj = Z.modelspace->GetOrbit(j);
            double ji = oi.j2/2.;
            double jj = oj.j2/2.;

            // now loop over the Jprime entering in the Pandya transformation
            int parity_cc = (oi.l+oQ.l)%2;
//            int Tz_cc = std::abs(oi.tz2+oQ.tz2)/2;
            int Tz_cc = std::abs(oi.tz2-oQ.tz2)/2;
            int Jmin = std::max(std::abs(int(ji-jQ)),std::abs(int(jk-jj)));
            int Jmax =  std::min( ji+jQ, jk+jj );
            for (int Jprime=Jmin; Jprime<=Jmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(ji,jj,J,jk,jQ,Jprime);
              if (std::abs(sixj)<1e-8) continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
              TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              size_t indx_iQ = i;
              size_t indx_kj = tbc_cc.GetLocalIndex(std::min(j,k),std::max(j,k)) +(k>j?nkets_cc:0);
              double me1 = Zbar.at(ch_cc)(indx_iQ,indx_kj);
              commijk -= (2*Jprime+1.) * sixj * me1  * Pij;
            }

          } // for ij_case
          commijk *= norm_ij * norm_kQ;
//            Z.TwoBody.AddToTBME_J(J,i,j,k,Q,  commijk );   // AddToTBME_J assumes a normalized matrix element.
//            Z.ThreeLeg.AddToME_J(J,i,j,k,  commijk );   // AddToTBME_J assumes a normalized matrix element.
          Z.ThreeLeg.AddToME_J(J,bra.p,bra.q,k,  commijk );   // AddToTBME_J assumes a normalized matrix element.
        } // for k
      } // for ibra, which is <ij|
   } // for ich, which is J etc
}





} // namespace Commutator
