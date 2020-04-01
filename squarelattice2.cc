
#include "itensor/all.h"
#include "tdvp_new.h"
#include <random>
#include <string>
#include <iostream>
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib>
using std::move;
using std::vector;

using namespace itensor;

int main()
{
    //system size
    int Nx = 4;
    int Ny = 4;
    int N = Nx*Ny;
    //total evolution time
    Real ttotal = 8;
    //timestep
    Real tstep = 0.02;
    Real tstepq;
    Real omega = 1.;
    Real V0 = 3.;
    double intpart;


    //number of total timesteps
    auto nt = int(ttotal / tstep + (1e-9 * (ttotal / tstep)));
    if (std::fabs(nt * tstep - ttotal) > 1E-9)
    {
        Error("Timestep not commensurate with total time");
    }

    //The part of the chain that we are going to calculate the correlators
    int C_Width = N;
    int C_Start = 1;

    //corelators measured every "skip" time steps
    int skip = 1;
    int data = int(nt/skip);   
    double Cor_z_t[data][C_Width];
    Print(data);

    //definitions of files that we save local expectation values
    std::string Cor_z_name,Fid_name;
    Cor_z_name = "C_z_" + std::to_string(N) + "11_.dat";
    Fid_name = "Fidelity_" + std::to_string(N) + "11_.dat";
    
    //This is just for streaming the data out as the simulation goes on
    ofstream outdata; 
    outdata.precision(16);



    //sites objects represent the Hilbert space,
    //a collection of "physical" indices
    auto sites = SpinHalf(N, {"ConserveQNs=", false});

    //Use AutoMPO to make Hamiltonian MPO
    auto ampo = AutoMPO(sites);
    for (int j = 1; j < N; ++j)
    {
        ampo += omega, "Sx", j;
        //NN
        if (j% Ny >0) ampo += V0, "projUp", j, "projUp", j + 1;
        if (j+Ny < N+1) ampo += V0, "projUp", j, "projUp", j + Ny;
        //NNN straight
        if (j% Ny >0 & j% Ny < Ny-1) ampo += V0/pow(2.,6.0), "projUp", j, "projUp", j + 2;
        if (j+2*Ny < N+1) ampo += V0/pow(2.,6.0), "projUp", j, "projUp", j + 2*Ny;
        //NNN diagonal
        if (j% Ny >0 & j+Ny < N+1) ampo += V0/pow(sqrt(2.),6.0), "projUp", j, "projUp", j + Ny+1;
        if (j% Ny != 1 & j+Ny < N+1)  ampo += V0/pow(sqrt(2.),6.0), "projUp", j, "projUp", j + Ny-1;
    }
    ampo += omega, "Sx", N;
    auto H = toMPO(ampo,{"Cutoff=",1E-10});
    //Print(maxLinkDim(H));
    //Create MPS
    //staggered starting state
    auto state = InitState(sites);
    for (int i = 1; i <= N; ++i)
    {
        double param = (i-1)/Ny+1;
        double fractpart = modf (param , &intpart);
        if (Ny % 2 == 0)
        {
        if (i % 2 == 1 & (int) intpart %2 == 1)  state.set(i, "Up");
        if (i % 2 == 1 & (int) intpart %2 == 0)  state.set(i, "Dn");
        if (i % 2 == 0 & (int) intpart %2 == 1)  state.set(i, "Dn");
        if (i % 2 == 0 & (int) intpart %2 == 0)  state.set(i, "Up");

        if (i % 2 == 1 & (int) intpart %2 == 1)  Print(1);
        if (i % 2 == 1 & (int) intpart %2 == 0)  Print(0);
        if (i % 2 == 0 & (int) intpart %2 == 1)  Print(0);
        if (i % 2 == 0 & (int) intpart %2 == 0)  Print(1);
        }
        

        if (Ny % 2 == 1)
        {
        if (i % 2 == 1 & (int) intpart %2 == 1)  state.set(i, "Up");
        if (i % 2 == 1 & (int) intpart %2 == 0)  state.set(i, "Up");
        if (i % 2 == 0 & (int) intpart %2 == 1)  state.set(i, "Dn");
        if (i % 2 == 0 & (int) intpart %2 == 0)  state.set(i, "Dn");

        if (i % 2 == 1 & (int) intpart %2 == 1)  Print(1);
        if (i % 2 == 1 & (int) intpart %2 == 0)  Print(1);
        if (i % 2 == 0 & (int) intpart %2 == 1)  Print(0);
        if (i % 2 == 0 & (int) intpart %2 == 0)  Print(0); 
        }
        

    }
    auto psi = MPS(state);
    auto psi0 = MPS(state);

//Perform 5 sweeps of DMRG
/*auto sweeps = Sweeps(40);
//Specify max number of states kept each sweep 
sweeps.maxdim() = 50, 50, 100, 100, 200;

//Run the DMRG algorithm
auto [energy,phi] = dmrg(H,psi,sweeps,"Quiet");

//Continue to analyze wavefunction afterward 
Print(inner(phi,H,phi)); //<psi|H|psi>*/

    //Cutoff of singular value decomposition and maximum bond dimension
    //Integration orders available 2,4
    auto sweeps = Sweeps(nt);
    auto args = Args("Cutoff", 1E-10, "MaxDim", 100, "IntegrationOrder", 2,"LocalDim", 2);
    //gauge fixing and normalization
    psi.position(1);
    psi.normalize();
    

     Print(real(innerC(psi,H,psi))); 
    for(auto sw : range1(sweeps.nsweep()))
        {
            TDVP2(psi,H, tstep,args);  

         //Print(maxLinkDim(psi));
         Print(sw);
      if(sw % skip == 0)
      {
      //loop over selected part of the chain
       for(int i = C_Start; i < (int) C_Start + C_Width; i++)
         {   
        //bring state to canonical form 
         psi.position(i);
        //local operator of interest 
         auto op_z_i = sites.op("projUp",i);
        //calculation of the correlation
         Cor_z_t[(int)(sw/skip-1)][i-1] = (dag(prime(psi(i), "Site"))*op_z_i*psi(i)).takeReal().real();
        //streaming the data to the folder
         outdata.open(Cor_z_name, std::ios::app);
         outdata << (double) Cor_z_t[(int)(sw/skip-1)][i-1] << endl;
         outdata.close();
         }  
		 psi.position(1);       
      }
        //fidelity
         outdata.open(Fid_name, std::ios::app);
         outdata << (double)pow(abs(innerC(psi,psi0)),2.0) << endl;
         outdata.close();
        }
    Print(cpu_mytime());
    return 0;
}
