/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                  ____                                         ///
///        _________________           _____________/   /\               _________________        ///
///       /____/_____/_____/|         /____/_____/ /___/  \             /____/_____/_____/|       ///
///      /____/_____/__G_ /||        /____/_____/|/   /\  /\           /____/_____/____ /||       ///
///     /____/_____/__+__/|||       /____/_____/|/ G /  \/  \         /____/_____/_____/|||       ///
///    |     |     |     ||||      |     |     |/___/   /\  /\       |     |     |     ||||       ///
///    |  I  |  M  |     ||/|      |  I  |  M  /   /\  /  \/  \      |  I  |  M  |     ||/|       ///
///    |_____|_____|_____|/||      |_____|____/ + /  \/   /\  /      |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |   /___/   /\  /  \/       |     |     |     ||||       ///
///    |  S  |  R  |     ||/|      |  S  |   \   \  /  \/   /        |  S  |  R  |  G  ||/|       ///
///    |_____|_____|_____|/||      |_____|____\ __\/   /\  /         |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |     \   \  /  \/          |     |     |     ||||       ///
///    |     |  +  |     ||/       |     |  +  |\ __\/   /           |     |  +  |  +  ||/        ///
///    |_____|_____|_____|/        |_____|_____|/\   \  /            |_____|_____|_____|/         ///
///                                               \___\/                                          ///
///                                                                                               ///
///           imsrg++ : Interface for performing standard IMSRG calculations.                     ///
///                     Usage is imsrg++  option1=value1 option2=value2 ...                       ///
///                     To get a list of options, type imsrg++ help                               ///
///                                                                                               ///
///                                                      - Ragnar Stroberg 2016                   ///
///                                                                                               ///
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
//    imsrg++.cc, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"
#include "PhysicalConstants.hh"
#include "version.hh"

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
#ifdef BUILDVERSION
  std::cout << "######  imsrg++ build version: " << BUILDVERSION << std::endl;
#endif

  Parameters parameters(argc,argv);
  if (parameters.help_mode) return 0;

  std::string inputtbme = parameters.s("2bme");
  std::string input3bme = parameters.s("3bme");
  std::string input3bme_type = parameters.s("3bme_type");
  std::string no2b_precision = parameters.s("no2b_precision");
  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");
  std::string custom_valence_space = parameters.s("custom_valence_space");
  std::string basis = parameters.s("basis");
  std::string method = parameters.s("method");
  std::string flowfile = parameters.s("flowfile");
  std::string intfile = parameters.s("intfile");
  std::string core_generator = parameters.s("core_generator");
  std::string valence_generator = parameters.s("valence_generator");
  std::string fmt2 = parameters.s("fmt2");
  std::string fmt3 = parameters.s("fmt3");
  std::string input_op_fmt = parameters.s("input_op_fmt");
  std::string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
  std::string LECs = parameters.s("LECs");
  std::string scratch = parameters.s("scratch");
  std::string valence_file_format = parameters.s("valence_file_format");
  std::string occ_file = parameters.s("occ_file");
  std::string physical_system = parameters.s("physical_system");
  bool use_brueckner_bch = parameters.s("use_brueckner_bch") == "true";
  bool nucleon_mass_correction = parameters.s("nucleon_mass_correction") == "true";
  bool relativistic_correction = parameters.s("relativistic_correction") == "true";
  bool IMSRG3 = parameters.s("IMSRG3") == "true";
  bool imsrg3_n7 = parameters.s("imsrg3_n7") == "true";
  bool imsrg3_at_end = parameters.s("imsrg3_at_end") == "true";
  bool write_omega = parameters.s("write_omega") == "true";
  bool freeze_occupations = parameters.s("freeze_occupations")=="true";
  bool discard_no2b_from_3n = parameters.s("discard_no2b_from_3n")=="true";
  bool hunter_gatherer = parameters.s("hunter_gatherer") == "true";
  bool goose_tank = parameters.s("goose_tank") == "true";
  bool discard_residual_input3N = parameters.s("discard_residual_input3N")=="true";
  bool use_NAT_occupations = (parameters.s("use_NAT_occupations")=="true") ? true : false;
  bool store_3bme_pn = (parameters.s("store_3bme_pn")=="true");
  bool only_2b_eta = (parameters.s("only_2b_eta")=="true");
  bool only_2b_omega = (parameters.s("only_2b_omega")=="true");
  bool perturbative_triples = (parameters.s("perturbative_triples")=="true");
  bool brueckner_restart = false;

  int eMax = parameters.i("emax");
  int lmax = parameters.i("lmax"); // so far I only use this with atomic systems.
  int E3max = parameters.i("e3max");
  int lmax3 = parameters.i("lmax3");
  int targetMass = parameters.i("A");
  int nsteps = parameters.i("nsteps");
  int file2e1max = parameters.i("file2e1max");
  int file2e2max = parameters.i("file2e2max");
  int file2lmax = parameters.i("file2lmax");
  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");
  int atomicZ = parameters.i("atomicZ");
  int emax_unocc = parameters.i("emax_unocc");
  int eMax_imsrg = parameters.i("emax_imsrg");
  int e2Max_imsrg = parameters.i("e2max_imsrg");
  int e3Max_imsrg = parameters.i("e3max_imsrg");
  if (e2Max_imsrg==-1 and eMax_imsrg != -1) e2Max_imsrg = 2*eMax_imsrg;
  if (e3Max_imsrg==-1 and eMax_imsrg != -1) e3Max_imsrg = std::min(E3max, 3*eMax_imsrg);

  double hw = parameters.d("hw");
  double smax = parameters.d("smax");
  double ode_tolerance = parameters.d("ode_tolerance");
  double dsmax = parameters.d("dsmax");
  double ds_0 = parameters.d("ds_0");
  double domega = parameters.d("domega");
  double omega_norm_max = parameters.d("omega_norm_max");
  double denominator_delta = parameters.d("denominator_delta");
  double BetaCM = parameters.d("BetaCM");
  double hwBetaCM = parameters.d("hwBetaCM");
  double eta_criterion = parameters.d("eta_criterion");
  double hw_trap = parameters.d("hw_trap");
  double dE3max = parameters.d("dE3max");
  double OccNat3Cut = parameters.d("OccNat3Cut");

  std::vector<std::string> opnames = parameters.v("Operators");
  std::vector<std::string> opsfromfile = parameters.v("OperatorsFromFile");
  std::vector<std::string> opnamesPT1 = parameters.v("OperatorsPT1");
  std::vector<std::string> opnamesRPA = parameters.v("OperatorsRPA");
  std::vector<std::string> opnamesTDA = parameters.v("OperatorsTDA");

  std::vector<Operator> ops;
  std::vector<std::string> spwf = parameters.v("SPWF");


  // test 2bme file
  if (inputtbme != "none" and fmt2.find("oakridge")==std::string::npos and fmt2 != "schematic" )
  {
    if( not std::ifstream(inputtbme).good() )
    {
      std::cout << "trouble reading " << inputtbme << "  fmt2 = " << fmt2 << "   exiting. " << std::endl;
      return 1;
    }
  }
  // test 3bme file
  if (input3bme != "none")
  {
    if( not std::ifstream(input3bme).good() )
    {
      std::cout << "trouble reading " << input3bme << " exiting. " << std::endl;
      return 1;
    }
  }

  // If we're reading in other operators, make sure those are ok too
  for (auto& tag : opsfromfile)
  {
     std::istringstream ss(tag);
     std::string opname,qnumbers,f2name,f3name="";
  
     getline(ss,opname,'^');
     getline(ss,qnumbers,'^');
     getline(ss,f2name,'^');

     if( not std::ifstream(f2name).good() )
     {
       std::cout << "trouble reading " << f2name << " exiting. " << std::endl;
       return 1;
     }

     if ( not ss.eof() ) // is there a 3-body file too?
     {
       getline(ss,f3name,'^');
       if( not std::ifstream(f3name).good() )
       {
         std::cout << "trouble reading " << f3name << " exiting. " << std::endl;
         return 1;
       }
     }
  }



  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  rw.Set3NFormat( fmt3 );



  // deal with some short-hand method names
  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=50000;
    method = "magnus";
  }
  if (method.find("brueckner") != std::string::npos)
  {
    if (method=="brueckner2") brueckner_restart=true;
    if (method=="brueckner1step")
    {
       nsteps = 1;
       core_generator = valence_generator;
    }
    use_brueckner_bch = true;
    omega_norm_max=500;
    method = "magnus";
  }



  // Test whether the scratch directory exists and we can write to it.
  // This is necessary because otherwise you get garbage for transformed operators and it's
  // not obvious what went wrong.
  if ( (method == "magnus") and  ( (opnames.size() + opsfromfile.size()) > 0 )  )
  {
    if ( scratch=="/dev/null" or scratch=="/dev/null/")
    {
      std::cout << "ERROR!!! using Magnus with scratch = " << scratch << " but you're also trying to transform some operators. Dying now. " << std::endl;
      exit(EXIT_FAILURE);
    }
    else if ( scratch != "" )
    {
      std::string testfilename = scratch + "/_this_is_a_test_delete_me_ID" + std::to_string(getpid());
      std::ofstream testout(testfilename);
      testout << "PASSED" << std::endl;
      testout.close();

      // now read it back.
      std::ifstream testin(testfilename);
      std::string checkpassed;
      testin >> checkpassed;
      if ( (checkpassed != "PASSED") or ( not testout.good() ) or ( not testin.good() ) )
      {
        std::cout << "ERROR in " << __FILE__ <<  " failed test write to scratch directory " << scratch << " that's bad. Dying now." << std::endl;
        exit(EXIT_FAILURE);
      }

      if ( remove(testfilename.c_str()) != 0 )
      {
        std::cout << "Error when attempting to delete " << testfilename << std::endl;
      }

    }
  }

//  if ( method=="magnus" and  scratch != "" and scratch!= "/dev/null" and scratch != "/dev/null/")
//  {
//    std::string testfilename = scratch + "/_this_is_a_test_delete_me";
//    std::ofstream testout(testfilename);
//    testout << "PASSED" << std::endl;
//    testout.close();
//    std::remove( testfilename.c_str() );
//    if ( not testout.good() )
//    {
//      std::cout << "WARNING in " << __FILE__ <<  " failed test write to scratch directory " << scratch;
//      if (opnames.size()>0 )
//      {
//      std::cout << "   dying now. " << std::endl;
//      exit(EXIT_FAILURE);
//      }
//      else std::cout << std::endl;
//    }
//  }
//  if ( (method=="magnus") and (scratch=="/dev/null" or scratch=="/dev/null/") )
//  {
//    if ( opnames.size() > 0 )
//    {
//      std::cout << "WARNING!!! using Magnus with scratch = " << scratch << " but you're also trying to transform some operators: ";
//      for (auto opn : opnames ) std::cout << opn << " ";
//      std::cout << "   dying now." << std::endl;
//      exit(EXIT_FAILURE);
//    }
//  }


//  ModelSpace modelspace;

  if (custom_valence_space!="") // if a custom space is defined, the input valence_space is just used as a name
  {
    if (valence_space=="") // if no name is given, then just name it "custom"
    {
      parameters.string_par["valence_space"] = "custom";
      flowfile = parameters.DefaultFlowFile();
      intfile = parameters.DefaultIntFile();
    }
    valence_space = custom_valence_space;
  }


  ModelSpace modelspace = ( reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space) );

//  std::cout << __LINE__ << "  constructed modelspace " << std::endl;
  modelspace.SetE3max(E3max);
  modelspace.SetLmax(lmax);
//  std::cout << __LINE__ << "  done setting E3max and lmax " << std::endl;


  if (emax_unocc>0)
  {
    modelspace.SetEmaxUnocc(emax_unocc);
  }

  if (physical_system == "atomic")
  {
    modelspace.InitSingleSpecies(eMax, reference, valence_space);
  }

//bhu
  std::string jobname = intfile.substr( intfile.find_last_of("/\\")+1 ) + '_';
  int nOmega=0;  //bhu
  std::ostringstream tempname;
  tempname << rw.GetScratchDir().c_str() << "/OMEGA_" << jobname << "Number";
  std::ifstream Omegafile(tempname.str());
  bool Omega_input = false;
if ( Omegafile.good() )
{
  Omega_input = true;
  tempname.str("");
  tempname << rw.GetScratchDir().c_str() << "/OMEGA_" << jobname << "occ";
  occ_file = tempname.str();
  freeze_occupations = true;
  Omegafile>>nOmega; Omegafile.ignore(128,'\n');
  std::cout<< "Read number of " << nOmega <<" OMEGA from " << jobname << std::endl;
  Omegafile.close();
}

  if (occ_file != "none" and occ_file != "" )
  {
    modelspace.Init_occ_from_file(eMax,valence_space,occ_file);
  }


  if (nsteps < 0) // default to 1 step for single ref, 2 steps for valence decoupling
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);



// For both dagger operators and single particle wave functions, it's convenient to
// just get every orbit in the valence space. So if SPWF="valence" ,  we append all valence orbits
  if ( std::find( spwf.begin(), spwf.end(), "valence" ) != spwf.end() )
  {
    // this erase/remove idiom is needed because remove just shuffles things around rather than actually removing it.
    spwf.erase( std::remove( spwf.begin(), spwf.end(), "valence" ), std::end(spwf) );
    for ( auto v : modelspace.valence )
    {
      spwf.push_back( modelspace.Index2String(v) );
    }
  }

  if ( std::find( opnames.begin(), opnames.end(), "rhop_all") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "rhop_all"), std::end(opnames) );
    for ( double r=0.0; r<=10.0; r+=0.2 )
    {
       std::ostringstream opn;
       opn << "rhop_" << r;
       opnames.push_back( opn.str() );
    }
  }

  if ( std::find( opnames.begin(), opnames.end(), "rhon_all") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "rhon_all"), std::end(opnames) );
    for ( double r=0.0; r<=10.0; r+=0.2 )
    {
       std::ostringstream opn;
       opn << "rhon_" << r;
       opnames.push_back( opn.str() );
    }
  }

  if ( std::find( opnames.begin(), opnames.end(), "DaggerHF_valence") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "DaggerHF_valence"), std::end(opnames) );
    for ( auto v : modelspace.valence )
    {
      opnames.push_back( "DaggerHF_"+modelspace.Index2String(v) );
    }
    std::cout << "I found DaggerHF_valence, so I'm changing the opnames list to :" << std::endl;
    for ( auto opn : opnames ) std::cout << opn << " ,  ";
    std::cout << std::endl;
  }

  if ( std::find( opnames.begin(), opnames.end(), "DaggerAlln_valence") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "DaggerAlln_valence"), std::end(opnames) );
    for ( auto v : modelspace.valence )
    {
      opnames.push_back( "DaggerAlln_"+modelspace.Index2String(v) );
    }
    std::cout << "I found DaggerAlln_valence, so I'm changing the opnames list to :" << std::endl;
    for ( auto opn : opnames ) std::cout << opn << " ,  ";
    std::cout << std::endl;
  }

//bhu
  int particle_rank = 2;
if ( Omega_input )
{
  particle_rank = 2;
}
else
{
  //  std::cout << "Making the Hamiltonian..." << std::endl;
//bhu  int particle_rank = input3bme=="none" ? 2 : 3;
  particle_rank = input3bme=="none" ? 2 : 3;
}
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();


//  if ( goose_tank )
//  {
  Commutator::SetUseGooseTank(goose_tank);
//  }

if ( ! Omega_input )
{

  std::cout << "Reading interactions..." << std::endl;


  if (inputtbme != "none")
  {
    if (fmt2 == "me2j")
      rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
    else if (fmt2 == "navratil" or fmt2 == "Navratil")
      rw.ReadBareTBME_Navratil(inputtbme, Hbare);
    else if (fmt2 == "oslo" )
      rw.ReadTBME_Oslo(inputtbme, Hbare);
    else if (fmt2.find("oakridge") != std::string::npos )
    { // input format should be: singleparticle.dat,vnn.dat
      size_t comma_pos = inputtbme.find_first_of(",");
      if ( fmt2.find("bin") != std::string::npos )
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "binary");
      else
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "ascii");
    }
    else if (fmt2 == "takayuki" )
      rw.ReadTwoBody_Takayuki( inputtbme, Hbare);
    else if (fmt2 == "nushellx" )
      rw.ReadNuShellX_int( Hbare, inputtbme );
//bhu      rw.ReadTokyo( inputtbme, Hbare, "tokyo");
    else if (fmt2 == "schematic" )
    {
      std::cout << "using schematic potential " << inputtbme << std::endl;
      if ( inputtbme == "Minnesota") Hbare += imsrg_util::MinnesotaPotential( modelspace );
    }

    std::cout << "done reading 2N" << std::endl;
  }

  // Read in the 3-body file
  if (Hbare.particle_rank >=3)
  {
    if(input3bme_type == "full"){
      rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
      std::cout << "done reading 3N" << std::endl;
    }
    if(input3bme_type == "no2b"){

      Hbare.ThreeBody.SetMode("no2b");
      if (no2b_precision == "half")  Hbare.ThreeBody.SetMode("no2bhalf");

      Hbare.ThreeBody.ReadFile( {input3bme}, {file3e1max, file3e2max, file3e3max, file3e1max} );
      rw.File3N = input3bme;
      std::cout << "done reading 3N" << std::endl;

    }
  }

  if (store_3bme_pn)
  {
    Hbare.ThreeBody.TransformToPN();
  }




  if (inputtbme == "none" and physical_system == "atomic")
  {

    using PhysConst::M_ELECTRON;
    using PhysConst::M_NUCLEON;
    int Z = (atomicZ>=0) ?  atomicZ : modelspace.GetTargetZ() ;
    Hbare -= Z*imsrg_util::VCentralCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;
    Hbare += imsrg_util::VCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;  // convert oscillator length from fm with nucleon mass to nm with electon mass (in eV).
    Hbare += imsrg_util::KineticEnergy_Op(modelspace); // Don't need to rescale this, because it's related to the oscillator frequency, which we input.
    Hbare /= PhysConst::HARTREE; // Convert to Hartree
  }

  if (fmt2 != "nushellx" and physical_system != "atomic" and hw_trap < 0)  // Don't need to add kinetic energy if we read a shell model interaction
  {
    Hbare += imsrg_util::Trel_Op(modelspace);
    if (Hbare.OneBody.has_nan())
    {
       std::cout << "  Looks like the Trel op is hosed from the get go. Dying." << std::endl;
       std::exit(EXIT_FAILURE);
    }
  }

  // Add an external harmonic trap
  if ( hw_trap > 0 )
  {
    Hbare += 0.5 * (PhysConst::M_NUCLEON * hw_trap * hw_trap)/(PhysConst::HBARC*PhysConst::HBARC) * imsrg_util::RSquaredOp(modelspace); 
    Hbare += imsrg_util::KineticEnergy_Op(modelspace); // use lab-frame kinetic energy
  }

  // correction to kinetic energy because M_proton != M_neutron
  if ( nucleon_mass_correction)
  {
    Hbare += imsrg_util::Trel_Masscorrection_Op(modelspace);
  }

  if ( relativistic_correction)
  {
    Hbare += imsrg_util::KineticEnergy_RelativisticCorr(modelspace);
  }




  // Add a Lawson center of mass term. If hwBetaCM is specified, use that frequency, otherwise use the basis frequency
  if (std::abs(BetaCM)>1e-3)
  {
    if (hwBetaCM < 0) hwBetaCM = modelspace.GetHbarOmega();
    std::ostringstream hcm_opname;
    hcm_opname << "HCM_" << hwBetaCM;
    Hbare += BetaCM * imsrg_util::OperatorFromString( modelspace, hcm_opname.str());
  }

}



  std::cout << "Creating HF" << std::endl;
  HFMBPT hf(Hbare); // HFMBPT inherits from HartreeFock, so this works for HF and NAT bases.

  if (not freeze_occupations )  hf.UnFreezeOccupations();
  if ( discard_no2b_from_3n) hf.DiscardNO2Bfrom3N();
  std::cout << "Solving" << std::endl;

  if (basis!="oscillator" and (! Omega_input) )
  {
    hf.Solve();
  }

  // decide what to keep after normal ordering
  int hno_particle_rank = 2;
  if ((IMSRG3) and (Hbare.ThreeBodyNorm() > 1e-5))  hno_particle_rank = 3;
  if (discard_residual_input3N) hno_particle_rank = 2;

  Operator& HNO = Hbare; // The reference & means we overwrite Hbare and save some memory
  if (basis == "HF" and method !="HF")
  {
    HNO = hf.GetNormalOrderedH( hno_particle_rank );
    if ((IMSRG3 or perturbative_triples) and OccNat3Cut>0 ) hf.GetNaturalOrbitals();
  }
  else if (basis == "NAT") // we want to use the natural orbital basis
  {
    hf.UseNATOccupations( use_NAT_occupations );

//  GetNaturalOrbitals() calls GetDensityMatrix(), which computes the 1b density matrix up to MBPT2
//  using the NO2B Hamiltonian in the HF basis, obtained with GetNormalOrderedH().
//  Then it calls DiagonalizeRho() which diagonalizes the density matrix, yielding the natural orbital basis.
    hf.GetNaturalOrbitals();
    HNO = hf.GetNormalOrderedHNAT( hno_particle_rank );

//  SRS: I'm commenting this out because this is not reasonably-expected default behavior
//    // For now, even if we use the NAT occupations, we switch back to naive occupations after the normal ordering
//    // This should be investigated in more detail.
//    if (use_NAT_occupations)
//    {
//      hf.FillLowestOrbits();
//      std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
//      HNO = HNO.UndoNormalOrdering();
//      hf.UpdateReference();
//      modelspace.SetReference(modelspace.core); // change the reference
//      std::cout << "Doing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
//      HNO = HNO.DoNormalOrdering();
//    }

  }
  else if (basis == "oscillator")
  {
    HNO = Hbare.DoNormalOrdering();
  }

  if (perturbative_triples)
  {
    modelspace.SetdE3max(dE3max);
    modelspace.SetOccNat3Cut(OccNat3Cut);
    std::array<size_t,2> nstates = modelspace.CountThreeBodyStatesInsideCut();
    std::cout << "We will compute perturbative triples corrections" << std::endl;
    std::cout << "Truncations: dE3max = " << dE3max << "   OccNat3Cut = " << std::scientific << OccNat3Cut << "  ->  number of 3-body states kept:  " << nstates[0] << " out of " << nstates[1] << std::endl << std::fixed;
  }

  if (IMSRG3  )
  {
    modelspace.SetdE3max(dE3max);
    modelspace.SetOccNat3Cut(OccNat3Cut);
    std::array<size_t,2> nstates = modelspace.CountThreeBodyStatesInsideCut();
    std::cout << "You have chosen IMSRG3. good luck..." << std::endl;
    std::cout << "Truncations: dE3max = " << dE3max << "   OccNat3Cut = " << std::scientific << OccNat3Cut << "  ->  number of 3-body states kept:  " << nstates[0] << " out of " << nstates[1] << std::endl;

    if (hno_particle_rank<3 ) // if we're doing IMSRG3, we need a 3 body operator
    {
      Operator H3(modelspace,0,0,0,3);
      std::cout << "Constructed H3" << std::endl;
      H3.ZeroBody = HNO.ZeroBody;
      H3.OneBody = HNO.OneBody;
      H3.TwoBody = HNO.TwoBody;
      HNO = H3;
      std::cout << "Replacing HNO" << std::endl;
      std::cout << "Hbare Three Body Norm is " << Hbare.ThreeBodyNorm() << std::endl;
      HNO.ThreeBody.SwitchToPN_and_discard();
    }
  }



//  if ( spwf.size() > 0 )
//  {
    // If the length of spwf is zero, nothing happens
    imsrg_util::WriteSPWaveFunctions( spwf, hf, intfile);
//  }



  HNO -= BetaCM * 1.5*hwBetaCM; // This is just the zero-body piece. The other stuff was added earlier.
  std::cout << "Hbare 0b = " << std::setprecision(8) << HNO.ZeroBody << std::endl;

  if (method != "HF")
  {
    std::cout << "Perturbative estimates of gs energy:" << std::endl;
    double EMP2 = HNO.GetMP2_Energy();
    double EMP2_3B = HNO.GetMP2_3BEnergy();
    std::cout << "EMP2 = " << EMP2 << std::endl;
    std::cout << "EMP2_3B = " << EMP2_3B << std::endl;
    std::array<double,3> Emp_3 = HNO.GetMP3_Energy();
    double EMP3 = Emp_3[0]+Emp_3[1]+Emp_3[2];
    std::cout << "E3_pp = " << Emp_3[0] << "  E3_hh = " << Emp_3[1] << " E3_ph = " << Emp_3[2] << "   EMP3 = " << EMP3 << std::endl;
    std::cout << "To 3rd order, E = " << HNO.ZeroBody + EMP2 + EMP3 + EMP2_3B << std::endl;
  }



  std::cout << "done with pert stuff, method = " << method << std::endl;
  // Calculate all the desired operators. If we're using magnus, we'll do this after the flow is over
  if ( method != "magnus" )
  {

    std::cout << "MAKING THE OPERATORS LINE " << __LINE__ << std::endl;
    for (auto& opname : opnames)
    {
//        ops.emplace_back( imsrg_util::OperatorFromString(modelspace,opname) );
        Operator optmp = imsrg_util::OperatorFromString(modelspace,opname);

        if ((basis == "HF") and (opname.find("DaggerHF") == std::string::npos)  )
        {
          optmp = hf.TransformToHFBasis(optmp);
        }
        else if ((basis == "NAT") and (opname.find("DaggerHF") == std::string::npos)  )
        {
          optmp = hf.TransformHOToNATBasis(optmp);
        }
        optmp = optmp.DoNormalOrdering();
//bhu        rw.WriteTensorTokyo(intfile+"_"+opname+".snt",optmp);
        if (method == "HF" or method == "MP3")
        {
          std::cout << "HF expectation value  " << opname << "  " << optmp.ZeroBody << std::endl;
        }
        if (method == "MP3")
        {
          double dop = optmp.MP1_Eval( HNO );
          std::cout << "Operator 1st order correction  " << dop << "  ->  " << optmp.ZeroBody + dop << std::endl;
        }
        if ( opname == "Rp2" )
        {
          double Rp2 = optmp.ZeroBody;
          int Z = modelspace.GetTargetZ();
          int A = modelspace.GetTargetMass();
          std::cout << " HF point proton radius = " << sqrt( Rp2 ) << std::endl;
          std::cout << " HF charge radius = " << ( abs(Rp2)<1e-6 ? 0.0 : sqrt( Rp2 + r2p + r2n*(A-Z)/Z + DarwinFoldy) ) << std::endl;
        }
    }
    opnames.resize(0); // we already dealt with all the operators in opnames, so clear it out.

    // Calculate first order perturbative correction to some operators, if that's what we asked for.
    // Strictly speaking, it doesn't make much sense to do this and then proceed with the IMSRG calculation,
    // but I'm not here to tell people what to do...
    for (auto& opnamept1 : opnamesPT1 )
    {
      ops.emplace_back( imsrg_util::FirstOrderCorr_1b( imsrg_util::OperatorFromString(modelspace,opnamept1)   , HNO ) );
      opnames.push_back( opnamept1+"PT1" );
    }
    for (auto& opnametda : opnamesTDA )
    {  // passing the argument "TDA" just sets the phhp and hpph blocks to zero in the RPA calculation
      ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnametda)   , HNO, "TDA" ) );
      opnames.push_back( opnametda+"TDA" );
    }
    for (auto& opnamerpa : opnamesRPA )
    {
      ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnamerpa)   , HNO, "RPA" ) );
      opnames.push_back( opnamerpa+"RPA" );
    }
  
  
    // the format should look like OpName^j_t_p_r^/path/to/2bfile^/path/to/3bfile  if particle rank of Op is 2-body, then 3bfile is not needed.
    for (auto& tag : opsfromfile)
    {
      std::istringstream ss(tag);
      std::string opname,qnumbers,f2name,f3name="";
//      std::vector<int> qn(4);
      int j,t,p,r;
  
      getline(ss,opname,'^');
      getline(ss,qnumbers,'^');
      getline(ss,f2name,'^');
      if ( not ss.eof() )  getline(ss,f3name,'^');

      ss.str(qnumbers);
      ss.clear();
      std::string tmp;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> j;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> t;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> p;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> r;
      
      std::cout << "Parsed tag. opname = " << opname << "  " << j << " " << t << " " << p << " " << r << "   file2 = " << f2name   << "    file3 = " << f3name << std::endl;

      Operator op(modelspace,j,t,p,r);
      if (r>2) op.ThreeBody.Allocate();
  //    std::cout << "Reading operator " << opname << "  in " << input_op_fmt << "  format from files " << f2name << "  ,  " << f3name << std::endl;
  //    std::cout << "Operator has particle rank " << op.GetParticleRank() << std::endl;
      if ( input_op_fmt == "navratil" )
      {
        rw.Read2bCurrent_Navratil( f2name, op );
      }
      else if ( input_op_fmt == "miyagi" )
      {
        if (f2name != "")
        {   
            Operator optmp = rw.ReadOperator2b_Miyagi( f2name, modelspace );
            op.TwoBody = optmp.TwoBody;
        }
        if ( r>2 and f3name != "")  rw.Read_Darmstadt_3body( f3name, op,  file3e1max,file3e2max,file3e3max);
      }
      ops.push_back( op );
      opnames.push_back( opname );
    }



//  for (auto& op : ops)
   for (size_t i=0;i<ops.size();++i)
   {
//     std::cout << "Before transforming  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
      // We don't transform a DaggerHF, because we want the a^dagger to already refer to the HF basis.
     if ((basis == "HF") and (opnames[i].find("DaggerHF") == std::string::npos)  )
     {
       ops[i] = hf.TransformToHFBasis(ops[i]);
     }
     else if ((basis == "NAT") and (opnames[i].find("DaggerHF") == std::string::npos)  )
     {
       ops[i] = hf.TransformHOToNATBasis(ops[i]);
     }
//     std::cout << "After transforming  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
     ops[i] = ops[i].DoNormalOrdering();
//     std::cout << "Before normal ordering  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
     if (method == "HF" or method == "MP3")
     {
       std::cout << "HF expectation value  " << opnames[i] << "  " << ops[i].ZeroBody << std::endl;
     }
     if (method == "MP3")
     {
       double dop = ops[i].MP1_Eval( HNO );
       std::cout << "Operator 1st order correction  " << dop << "  ->  " << ops[i].ZeroBody + dop << std::endl;
     }
    if ( opnames[i] == "Rp2" )
    {
      double Rp2 = ops[i].ZeroBody;
      int Z = modelspace.GetTargetZ();
      int A = modelspace.GetTargetMass();
      std::cout << " HF point proton radius = " << sqrt( Rp2 ) << std::endl;
      std::cout << " HF charge radius = " << ( abs(Rp2)<1e-6 ? 0.0 : sqrt( Rp2 + r2p + r2n*(A-Z)/Z + DarwinFoldy) ) << std::endl;
    }
   }// for ops.size




  }// if method != "magnus"


  if (basis=="HF" or basis=="NAT")
  {
    std::cout << basis << " Single particle energies and wave functions:" << std::endl;
    hf.PrintSPEandWF();
    std::cout << std::endl;
    std::cout << basis << " diagonal matrix:" << std::endl;
    hf.PrintSPE();
    std::cout << std::endl;
  }

  if ( method == "HF" or method == "MP3")
  {
    HNO.PrintTimes();
    return 0;
  }


  if (method == "FCI")
  {
  std::cout << __func__ << "  line " << __LINE__ << std::endl;
   // we want the 1b piece to be diagonal in the vacuum NO representation
    HNO = HNO.UndoNormalOrdering();
    double previous_zero_body = HNO.ZeroBody;
    modelspace.SetReference("vacuum");
    HartreeFock hfvac(HNO);
    hfvac.Solve();

//    Operator Hvac = hfvac.GetNormalOrderedH();
    HNO = hfvac.GetNormalOrderedH();
    std::cout << "HNO had zero body = " << HNO.ZeroBody << "  and I add " << previous_zero_body << " to it. " << std::endl;
    HNO.ZeroBody += previous_zero_body;

    rw.WriteNuShellX_int(HNO,intfile+".int");
    rw.WriteNuShellX_sps(HNO,intfile+".sp");

    std::cout << "NO wrt vacuum. One Body term is hopfully still diagonal?" << std::endl << HNO.OneBody << std::endl;

    for (index_t i=0;i<ops.size();++i)
    {
      ops[i] = ops[i].UndoNormalOrdering();
      if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
      {
        rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
      }
      else
      {
        rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
        rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
      }
    }
    HNO.PrintTimes();
    return 0;
  }
  if (only_2b_omega)
  {
    std::cout << " Restricting the Magnus operator Omega to be 2b." << std::endl;
    Commutator::SetOnly2bOmega(only_2b_omega);
  }

//bhu
  if (write_omega and (! Omega_input))
  {
    std::ofstream fileocc;
    std::ostringstream occname;
    occname << rw.GetScratchDir().c_str() << "/OMEGA_" << jobname << "occ";
    fileocc.open( occname.str(), std::ofstream::out);
    int wint = 4; int wdouble = 20; int pdouble = 10;
    for (auto h : modelspace.holes)
    {
      Orbit& oi = modelspace.GetOrbit(h);
      if ( std::abs(oi.occ)>1e-6 )
      {
        fileocc << std::setw(wint) << oi.n << std::setw(wint) << oi.l << std::setw(wint) << oi.j2 << std::setw(wint) << oi.tz2 << std::setw(wdouble) << std::setiosflags(std::ios::fixed) << std::setprecision(pdouble) << std::scientific << oi.occ << std::endl;
      }
    }
    fileocc.close();
  }

  if( Omega_input )
  {
    std::string filehf = rw.GetScratchDir() + "/OMEGA_" + jobname + "HF.mat";
    bool filehfsucess = false;
    if (basis == "HF")
    {
       filehfsucess = hf.C.load(filehf);
    }
    else if (basis == "NAT")
    {
       filehfsucess = hf.C_HO2NAT.load(filehf);
    }
    if (filehfsucess == false)
    {
      std::cout<<"Couldn't load HF coefficient matrix: "<< filehf << std::endl;
      exit(EXIT_FAILURE);
    }
  }


  // We may want to use a smaller model space for the IMSRG evolution than we used for the HF step.
  // This is most effective when using natural orbitals.
//  ModelSpace modelspace_imsrg = ( reference=="default" ? ModelSpace(eMax_imsrg,e2Max_imsrg,e3Max_imsrg,valence_space) : ModelSpace(eMax_imsrg,e2Max_imsrg,e3Max_imsrg,reference,valence_space) );
  ModelSpace modelspace_imsrg = modelspace;
  if ( (eMax_imsrg != -1) or (e2Max_imsrg != -1) or (e3Max_imsrg) != -1)
  {
//     ModelSpace modelspace_imsrg = modelspace;
     std::cout << "Truncating modelspace for IMSRG calculation: emax e2max e3max  ->  " << eMax_imsrg << " " << e2Max_imsrg << " " << e3Max_imsrg << std::endl;
     modelspace_imsrg.SetEmax( eMax_imsrg);
     modelspace_imsrg.SetE2max( e2Max_imsrg);
     modelspace_imsrg.SetE3max( e3Max_imsrg);
     modelspace_imsrg.Init( eMax_imsrg, reference, valence_space);
   //  if (emax_unocc>0) modelspace_imsrg.SetEmaxUnocc(emax_unocc);
     if (physical_system == "atomic") modelspace_imsrg.InitSingleSpecies(eMax_imsrg, reference, valence_space);
     if (occ_file != "none" and occ_file != "" ) modelspace_imsrg.Init_occ_from_file(eMax_imsrg,valence_space,occ_file);
//     if (physical_system == "atomic") modelspace_imsrg.InitSingleSpecies(eMax_imsrg, eMax_imsrg, e3Max_imsrg, reference, valence_space);
//     if (occ_file != "none" and occ_file != "" ) modelspace_imsrg.Init_occ_from_file(eMax_imsrg,e2Max_imsrg,e3Max_imsrg,valence_space,occ_file);


     // If the occupations in modelspace were different from the naive filling, we want to keep those.
     std::map<index_t,double> hole_map;
     for ( auto& i_new : modelspace_imsrg.all_orbits )
     {
        Orbit& oi_new = modelspace_imsrg.GetOrbit(i_new);
        index_t i_old = modelspace.GetOrbitIndex( oi_new.n, oi_new.l, oi_new.j2, oi_new.tz2 );
        Orbit& oi_old = modelspace.GetOrbit(i_old);
        hole_map[i_new] = oi_old.occ;
     }
     modelspace_imsrg.SetReference( hole_map );


     HNO = HNO.Truncate(modelspace_imsrg);

//     modelspace = modelspace_imsrg;  // this could cause some confusion later on...
//    hf.PrintSPEandWF();
  }
  else
  {
    HNO.SetModelSpace(modelspace_imsrg);
  }



  IMSRGSolver imsrgsolver(HNO);
  imsrgsolver.SetHin(HNO); // necessary?
  imsrgsolver.SetReadWrite(rw);
  imsrgsolver.SetMethod(method);
  imsrgsolver.SetEtaCriterion(eta_criterion);
  imsrgsolver.GetGenerator().SetOnly2bEta(only_2b_eta);
  imsrgsolver.max_omega_written = 500;
  imsrgsolver.SetHunterGatherer( hunter_gatherer );
  imsrgsolver.SetPerturbativeTriples(perturbative_triples);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(dsmax);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

//  if (method == "NSmagnus") // "No split" magnus
//  {
//    omega_norm_max=50000;
//    method = "magnus";
//  }
//  if (method.find("brueckner") != std::string::npos)
//  {
//    if (method=="brueckner2") brueckner_restart=true;
//    if (method=="brueckner1step")
//    {
//       nsteps = 1;
//       core_generator = valence_generator;
//    }
//    use_brueckner_bch = true;
//    omega_norm_max=500;
//    method = "magnus";
//  }

  Commutator::SetUseBruecknerBCH(use_brueckner_bch);
  Commutator::SetUseIMSRG3(IMSRG3);
  Commutator::SetUseIMSRG3N7(imsrg3_n7);
  if (use_brueckner_bch)
  {
    std::cout << "Using Brueckner flavor of BCH" << std::endl;
  }
  if (IMSRG3)
  {
    std::cout << "Using IMSRG(3) commutators. This will probably be slow..." << std::endl;
  }
  if (imsrg3_n7)
  {
    std::cout << "  only including IMSRG3 commutator terms that scale up to n7" << std::endl;
  }


  // Here's where we need to have the operators
//  std::cout << "MADE IT TO LINE " << __LINE__ << std::endl;


  if (method == "flow" or method == "flow_RK4" )
  {
    for (auto& op : ops )  imsrgsolver.AddOperator( op );
  }

  imsrgsolver.SetGenerator(core_generator);
  if (core_generator.find("imaginary")!=std::string::npos or core_generator.find("wegner")!=std::string::npos )
  {
   if (ds_0>1e-2)
   {
     ds_0 = 1e-4;
     dsmax = 1e-2;
     imsrgsolver.SetDs(ds_0);
     imsrgsolver.SetDsmax(dsmax);
   }
  }

//bhu
if ( Omega_input )
{
  imsrgsolver.jobname = jobname;
  imsrgsolver.n_omega_written = nOmega;
}

  imsrgsolver.Solve();

  if (IMSRG3)
  {
    std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
  }
  if ( perturbative_triples and method=="magnus" )
  {
//    modelspace.SetdE3max(dE3max);
//    modelspace.SetOccNat3Cut(OccNat3Cut);
//    size_t nstates_kept = modelspace.CountThreeBodyStatesInsideCut();
//    std::array<size_t,2> nstates = modelspace.CountThreeBodyStatesInsideCut();
//    std::cout << "Truncations: dE3max = " << dE3max << "   OccNat3Cut = " << std::scientific << OccNat3Cut << "  ->  number of 3-body states kept:  " << nstates[0] << " out of " << nstates[1] << std::endl << std::fixed;
    double dE_triples = imsrgsolver.GetPerturbativeTriples();
    std::cout << "Perturbative triples:  " << std::setw(16) << std::setprecision(8) << dE_triples << " -> " << imsrgsolver.GetH_s().ZeroBody + dE_triples << std::endl;
  }



  if (brueckner_restart)
  {
     arma::mat newC = hf.C * arma::expmat( -imsrgsolver.GetOmega(0).OneBody  );
//     if (input3bme != "none") Hbare.SetParticleRank(3);
     HNO = hf.GetNormalOrderedH(newC);
     imsrgsolver.SetHin(HNO);
     imsrgsolver.s = 0;
     imsrgsolver.Solve();
  }

  if (nsteps > 1 and valence_space != reference) // two-step decoupling, do core first
  {
    if (method == "magnus") smax *= 2;

    imsrgsolver.SetGenerator(valence_generator);
    std::cout << "Setting generator to " << valence_generator << std::endl;
//    modelspace.ResetFirstPass();
    modelspace_imsrg.ResetFirstPass();
    if (valence_generator.find("imaginary")!=std::string::npos or valence_generator.find("wegner")!=std::string::npos)
    {
     if (ds_0>1e-2)
     {
       ds_0 = 1e-4;
       dsmax = 1e-2;
       imsrgsolver.SetDs(ds_0);
       imsrgsolver.SetDsmax(dsmax);
     }
    }
    imsrgsolver.SetSmax(smax);
    imsrgsolver.Solve();
  }


  if ( imsrg3_at_end )
  {
    if ( method.find("magnus") != std::string::npos )
    {
      std::cout << "Performing final BCH transformation at the IMSRG(3) level" << std::endl;

//      modelspace.SetdE3max(dE3max);
//      modelspace.SetOccNat3Cut(OccNat3Cut);
//      int new_E3max = std::min(modelspace.GetE3max(), int( std::ceil(3*std::max( modelspace.GetEFermi()[-1], modelspace.GetEFermi()[+1])+dE3max)));
      modelspace_imsrg.SetdE3max(dE3max);
      modelspace_imsrg.SetOccNat3Cut(OccNat3Cut);
      int new_E3max = std::min(modelspace_imsrg.GetE3max(), int( std::ceil(3*std::max( modelspace_imsrg.GetEFermi()[-1], modelspace_imsrg.GetEFermi()[+1])+dE3max)));
      std::cout << "Setting new E3max = " << new_E3max << std::endl;
//      modelspace.SetE3max(  new_E3max);
      modelspace_imsrg.SetE3max(  new_E3max);

//      Operator H3(modelspace,0,0,0,3);
      Operator H3(modelspace_imsrg,0,0,0,3);
      std::cout << "Constructed H3" << std::endl;
      H3.ZeroBody = HNO.ZeroBody;
      H3.OneBody = HNO.OneBody;
      H3.TwoBody = HNO.TwoBody;
      HNO = H3;
      std::cout << "Replacing HNO" << std::endl;
      std::cout << "Hbare Three Body Norm is " << Hbare.ThreeBodyNorm() << std::endl;
      HNO.ThreeBody.SwitchToPN_and_discard();


      Commutator::SetUseIMSRG3(true);
      Commutator::SetUseIMSRG3N7(imsrg3_n7);
//      Operator H_with_3 = imsrgsolver.Transform(  *(imsrgsolver.H_0) );
      Operator H_with_3 = imsrgsolver.Transform( HNO );

      // Now throw away the residual 3-body so we don't need to keep it after re-normal ordering
      H_with_3.SetNumberLegs(4);
      H_with_3.SetParticleRank(2);
      imsrgsolver.FlowingOps[0] = H_with_3;
    }
    else
    {
      std::cout << "selected imsrg3_at_end, but method != magnus, so I don't know what to do. Ignoring." << std::endl;
    }

  }


/*
  // Transform all the operators
  if (method == "magnus")
  {
    if (ops.size()>0) std::cout << "transforming operators" << std::endl;
    for (size_t i=0;i<ops.size();++i)
    {
      std::cout << opnames[i] << " " << std::endl;
      ops[i] = imsrgsolver.Transform(ops[i]);
      std::cout << " (" << ops[i].ZeroBody << " ) " << std::endl;
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");
    }
    std::cout << std::endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }
  if (method == "flow" or method == "flow_RK4" )
  {
    for (size_t i=0;i<ops.size();++i)
    {
      ops[i] = imsrgsolver.GetOperator(i+1);  // the zero-th operator is the Hamiltonian
    }
  }
*/


  // If we're doing targeted/ensemble normal ordering
  // we now re-normal order wrt to the core
  // and do any remaining flow.
//  ModelSpace ms2(modelspace);
  ModelSpace ms2(modelspace_imsrg);
  ms2.SetReference(ms2.core); // change the reference
  bool renormal_order = false;
//  if (modelspace.valence.size() > 0 )
  if (modelspace_imsrg.valence.size() > 0 )
//  if (modelspace.valence.size() > 0 or basis=="NAT")
  {
//    renormal_order = modelspace.holes.size() != modelspace.core.size();
    renormal_order = modelspace.holes.size() != modelspace_imsrg.core.size();
    if (not renormal_order)
    {
//      for (auto c : modelspace.core)
      for (auto c : modelspace_imsrg.core)
      {
//         if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (std::abs(1-modelspace.GetOrbit(c).occ)>1e-6))
         if ( (find( modelspace_imsrg.holes.begin(), modelspace_imsrg.holes.end(), c) == modelspace_imsrg.holes.end()) or (std::abs(1-modelspace_imsrg.GetOrbit(c).occ)>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
  if ( renormal_order )
  {

    HNO = imsrgsolver.GetH_s();

//    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
//    std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
    std::cout << "Undoing NO wrt A=" << modelspace_imsrg.GetAref() << " Z=" << modelspace_imsrg.GetZref() << std::endl;
    std::cout << "Before doing so, the spes are " << std::endl;
//    for ( auto i : modelspace.all_orbits ) std::cout << "  " << i << " : " << HNO.OneBody(i,i) << std::endl;
    for ( auto i : modelspace_imsrg.all_orbits ) std::cout << "  " << i << " : " << HNO.OneBody(i,i) << std::endl;
    if (IMSRG3)
    {
      std::cout << "Re-normal-ordering wrt the core. For now, we just throw away the 3N at this step." << std::endl;
      HNO.SetNumberLegs(4);
      HNO.SetParticleRank(2);
    }

    HNO = HNO.UndoNormalOrdering();

//bhu
    std::cout << std::endl;
    std::cout << "After undoing NO, the spes are " << std::endl;
    for ( auto i : modelspace_imsrg.all_orbits ) std::cout << "  " << i << " : " << HNO.OneBody(i,i) << std::endl;
    std::cout << std::endl;
    Operator HNO_1B = HNO;
    HNO_1B.EraseZeroBody();
    HNO_1B.EraseOneBody();
    HNO_1B = HNO_1B.DoNormalOrdering();
    std::cout << "One-body part only from two body is " << std::endl;
    for ( auto i : modelspace_imsrg.all_orbits ) std::cout << "  " << i << " : " << HNO_1B.OneBody(i,i) << std::endl;
    std::cout << std::endl;
//bhu

    HNO.SetModelSpace(ms2);
    std::cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << std::endl;
    HNO = HNO.DoNormalOrdering();

//bhu
    std::cout << "After Re-NO, the spes are " << std::endl;
    for ( auto i : ms2.all_orbits ) std::cout << "  " << i << " : " << HNO.OneBody(i,i) << std::endl;
    std::cout << std::endl;
//bhu

    imsrgsolver.FlowingOps[0] = HNO;

// More flowing is unnecessary, since things should stay decoupled.
//    imsrgsolver.SetHin(HNO);
//    imsrgsolver.SetEtaCriterion(1e-4);
//    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
//    std::cout << "Final transformation on the operators..." << std::endl;
//    int iop = 0;
//    for (auto& op : ops)
//    {
//      std::cout << opnames[iop++] << std::endl;
//      op = op.UndoNormalOrdering();
//      op.SetModelSpace(ms2);
//      op = op.DoNormalOrdering();
//      // transform using the remaining omegas
//      op = imsrgsolver.Transform_Partial(op,nOmega);
//    }
  }

//bhu
if( ! Omega_input )
{

  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
//  if (modelspace.valence.size() > 0)
  if (modelspace_imsrg.valence.size() > 0)
  {
    if (valence_file_format == "antoine") // this is still being tested...
    {
      rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
//      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace_imsrg.GetAref(),modelspace_imsrg.GetZref());
    }
    std::cout << "Writing files: " << intfile << std::endl;
    if (valence_file_format == "tokyo")
    {
     rw.WriteTokyo(imsrgsolver.GetH_s(),intfile+".snt", "");
    }
    else
    {
      rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
      rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
    }

//    if (method == "magnus" or method=="flow_RK4")
//    {
//       for (index_t i=0;i<ops.size();++i)
//       {
//          if ( ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1) and (ops[i].GetNumberLegs()%2==0) )
//          {
//            if (valence_file_format == "tokyo")
//            {
//              rw.WriteTokyo(ops[i],intfile+opnames[i]+".snt", "op");
//            }
//            else
//            {
//              rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
//            }
//          }
//          else if ( ops[i].GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
//          {
////            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
//            rw.WriteDaggerOperator( ops[i], intfile+opnames[i]+".dag",opnames[i]);
//          }
//          else
//          {
//            if (valence_file_format == "tokyo")
//            {
//              rw.WriteTensorTokyo(intfile+opnames[i]+"_2b.snt",ops[i]);
//            }
//            else
//            {
//              rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
//              rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
//            }
//          }
//       }
//    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    std::cout << "Core Energy = " << std::setprecision(6) << imsrgsolver.GetH_s().ZeroBody << std::endl;
    rw.WriteTokyoFull(imsrgsolver.GetH_s(),intfile+".snt");
    for (index_t i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      std::cout << opnames[i] << " = " << std::setprecision(6) << ops[i].ZeroBody << "  " << std::setprecision(6) << sqrt( ops[i].ZeroBody ) << std::endl;
      if ( opnames[i] == "Rp2" )
      {
//         int Z = modelspace.GetTargetZ();
//         int A = modelspace.GetTargetMass();
         int Z = modelspace_imsrg.GetTargetZ();
         int A = modelspace_imsrg.GetTargetMass();
         std::cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
         std::cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DarwinFoldy) << std::endl;
      }
      if ((op.GetJRank()>0) or (op.GetTRank()>0) or (opnames[i]=="ISM")) // if it's a tensor, you probably want the full operator
      {
        std::cout << "Writing operator to " << intfile+opnames[i]+".op" << std::endl;
        rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
      }
    }
  }

}



/////////////////////
/// Transform operators and write them



  if (method == "magnus")
  {
    if (opnames.size()>0) std::cout << "transforming operators" << std::endl;

    for (size_t i=0;i<opnames.size();++i)
    {
      std::cout << "i..." << i << std::endl;
      auto opname = opnames[i];
      std::cout << i << ": " << opname << " " << std::endl;

      Operator op = imsrg_util::OperatorFromString( modelspace, opname );

      if ( basis == "HF")
      {
        op = hf.TransformToHFBasis(op).DoNormalOrdering();
      }
      if ( basis == "NAT")
      {
        op = hf.TransformHOToNATBasis(op).DoNormalOrdering();
      }
      if ( basis == "oscillator")
      {
        op = op.DoNormalOrdering();
      }
      std::cout << "   HF: " << std::setprecision(6) << op.ZeroBody << std::endl;

      if ( (eMax_imsrg != -1) or (e2Max_imsrg != -1) or (e3Max_imsrg) != -1)
      {
//     ModelSpace modelspace_imsrg = modelspace;
        std::cout << "Truncating modelspace for IMSRG calculation: emax e2max e3max  ->  " << eMax_imsrg << " " << e2Max_imsrg << " " << e3Max_imsrg << std::endl;
        op = op.Truncate(modelspace_imsrg);
      }



      op = imsrgsolver.Transform(op);

      if (renormal_order) 
      {
        op = op.UndoNormalOrdering();
        op.SetModelSpace(ms2);
        op = op.DoNormalOrdering();
      }
//      std::cout << " (" << ops[i].ZeroBody << " ) " << std::endl;
      std::cout << "   IMSRG: " << std::setprecision(6) << op.ZeroBody << std::endl;
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");


    if (modelspace.valence.size() > 0)
    {

    std::cout << "      " << op.GetJRank() << " " << op.GetTRank() << " " << op.GetParity() << "   " << op.GetNumberLegs() << std::endl;
    if ( ((op.GetJRank()+op.GetTRank()+op.GetParity())<1) and (op.GetNumberLegs()%2==0) )
    {
       std::cout << "writing scalar files " << std::endl;
      if (valence_file_format == "tokyo")
      {
        rw.WriteTokyo(op,intfile+"_"+opnames[i]+".snt", "op");
      }
      else
      {
        rw.WriteNuShellX_op(op,intfile+opnames[i]+".int");
      }
    }
    else if ( op.GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
    {
//      rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
      rw.WriteDaggerOperator( op, intfile+opnames[i]+".dag",opnames[i]);
    }
    else
    {
       std::cout << "writing tensor files " << std::endl;
      if (valence_file_format == "tokyo")
      {
        rw.WriteTensorTokyo(intfile+"_"+opnames[i]+".snt",op);
      }
      else
      {
        rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",op,opnames[i]);
        rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",op,opnames[i]);
      }
    }

    }
    else // single ref. just print the zero body pieces out.
    {

      std::cout << opnames[i] << " = " << std::setprecision(6) << op.ZeroBody << "  " << std::setprecision(6) << sqrt( op.ZeroBody ) << std::endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         std::cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
         std::cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DarwinFoldy) << std::endl;
      }
      if ((op.GetJRank()>0) or (op.GetTRank()>0) or (opnames[i]=="ISM")) // if it's a tensor, you probably want the full operator
      {
        std::cout << "Writing operator to " << intfile+opnames[i]+".op" << std::endl;
        rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
      }
    }

    }// for opnames

  }// if method == "magnus"



  if (method == "flow" or method == "flow_RK4" )
  {
    for (size_t i=0;i<ops.size();++i)
    {
      ops[i] = imsrgsolver.GetOperator(i+1);  // the zero-th operator is the Hamiltonian
    }
  }











//  std::cout << "Made it here and write_omega is " << write_omega << std::endl;
  if ( write_omega and (! Omega_input) )
  {
    std::string scratch = rw.GetScratchDir();
    imsrgsolver.FlushOmegaToScratch();
    for (int i=0; i < imsrgsolver.GetNOmegaWritten() ; i++)
    {
       std::ostringstream inputfile,outputfile;
       inputfile << scratch << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
       outputfile << scratch << "/OMEGA_" << jobname << std::setw(3) << std::setfill('0') << i;
       rw.CopyFile( inputfile.str(), outputfile.str() );
    }
//    rw.WriteOmega(intfile,scratch, imsrgsolver.n_omega_written);
    std::string hffile = scratch + "/OMEGA_" + jobname + "HF.mat";
    bool filesucess = false;
    if (basis == "HF")
    {
       filesucess = hf.C.save(hffile);
    }
    else if (basis == "NAT")
    {
       filesucess = hf.C_HO2NAT.save(hffile);
    }
    if (filesucess == false)
    {
      std::cout<<"Couldn't save HF coefficient matrix."<<std::endl;
    }
    // std::cout << "writing Omega to " << intfile << "_omega.op" << std::endl;
    // rw.WriteOperatorHuman(imsrgsolver.Omega.back(),intfile+"_omega.op");

    std::ofstream filename;
    tempname.str("");
    tempname << rw.GetScratchDir().c_str() << "/OMEGA_" << jobname << "Number";
    filename.open( tempname.str(), std::ofstream::out);
    filename << imsrgsolver.GetNOmegaWritten() << std::endl;
    filename.close();

  }



  if (IMSRG3)
  {
    std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
  }
  Hbare.PrintTimes();

  return 0;
}

