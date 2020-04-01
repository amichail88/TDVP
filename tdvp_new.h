//The following ifndef/define/endif pattern is called a
//scope guard, and prevents the C++ compiler (actually, preprocessor)
//from including a header file more than once.
#ifndef __TDVPO2_H
#define __TDVPO2_H

#include "itensor/all.h"
namespace itensor
{

//two-site TDVP time-step algorithm
//options: 2nd and 4th order integrators

template <class LocalOpT>
Real TDVP2(MPS &psi,
           LocalOpT &H,
           Cplx t,
           Args const &args);

//one-site TDVP time-step algorithm
//options: 2nd and 4th order integrators

template <class LocalOpT>
Real TDVP1(MPS &psi,
           LocalOpT &H,
           Cplx t,
           Args const &args);

template <class LocalOpT>
Real TDVPh(MPS &psi,
           LocalOpT &H,
           Cplx t,
           Args const &args);

template <class LocalOpT>
Real TDVP2(MPS &psi,
           LocalOpT &H,
           Cplx t,
           Args const &args)
{
    const int N = length(psi);
    Real alpha = 1. / (2. - pow(2., 0.33333333333));
    Real beta = 1. - 2 * alpha;

    auto Heff = LocalMPO(H);

    auto sweeps = Sweeps(1);

    Cplx tstep = t;
    Cplx tstepq;

    if (args.getInt("IntegrationOrder") == 2)
    {
        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                //Effective Hamiltonian for sites b,b+1
                Heff.numCenter(2);
                Heff.position(b, psi);

                //create two site tensor
                auto phi = psi(b) * psi(b + 1);

                //Krylov evolution of two site tensor and effective Hamiltonian
                auto argsd = Args("DebugLevel", -1);
                applyExp(Heff, phi, -Cplx_i * tstep / 2., argsd);

                //splitting of the evolved two-site tensor to single sites using SVD with some truncation error and some maximum dimension
                auto u = commonInds(psi(b), phi);
                auto [A, D, B] = svd(phi, u, args);

                //the norm of the state is just the norm of the singular value matrix D,  norm(D)

                if (ha == 1) //sweeping right
                {
                    //moving gauge to the right
                    psi.ref(b) = A;
                    psi.ref(b + 1) = D * B;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    auto phi = psi(b + 1);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    if (b == N - 1)
                    {
                        tstepq = 0;
                    }
                    if (b < N - 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * tstepq / 2.);
                    psi.ref(b + 1) = phi;
                }
                else if (ha == 2) //sweeping left
                {
                    psi.ref(b) = D * A;
                    psi.ref(b + 1) = B;
                    Heff.numCenter(1);
                    Heff.position(b, psi);
                    auto phi = psi(b);
                    if (b == 1)
                    {
                        tstepq = 0;
                    }
                    if (b > 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * tstepq / 2.);
                    psi.ref(b) = phi;
                }

            } // for loop over b
        }
    }

    if (args.getInt("IntegrationOrder") == 4)
    {
        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                //Effective Hamiltonian for sites b,b+1
                Heff.numCenter(2);
                Heff.position(b, psi);

                //create two site tensor
                auto phi = psi(b) * psi(b + 1);

                //Krylov evolution of two site tensor and effective Hamiltonian
                auto argsd = Args("DebugLevel", -1);
                applyExp(Heff, phi, -Cplx_i * alpha * tstep / 2., argsd);

                //splitting of the evolved two-site tensor to single sites using SVD with some truncation error and some maximum dimension
                auto u = commonInds(psi(b), phi);
                auto [A, D, B] = svd(phi, u, args);

                //the norm of the state is just the norm of the singular value matrix D,  norm(D)

                if (ha == 1) //sweeping right
                {
                    //moving gauge to the right
                    psi.ref(b) = A;
                    psi.ref(b + 1) = D * B;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    auto phi = psi(b + 1);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    if (b == N - 1)
                    {
                        tstepq = 0;
                    }
                    if (b < N - 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * alpha * tstepq / 2.);
                    psi.ref(b + 1) = phi;
                }
                else if (ha == 2) //sweeping left
                {
                    psi.ref(b) = D * A;
                    psi.ref(b + 1) = B;
                    Heff.numCenter(1);
                    Heff.position(b, psi);
                    auto phi = psi(b);
                    if (b == 1)
                    {
                        tstepq = 0;
                    }
                    if (b > 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * alpha * tstepq / 2.);
                    psi.ref(b) = phi;
                }

            } // for loop over b
        }

        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                //Effective Hamiltonian for sites b,b+1
                Heff.numCenter(2);
                Heff.position(b, psi);

                //create two site tensor
                auto phi = psi(b) * psi(b + 1);

                //Krylov evolution of two site tensor and effective Hamiltonian
                auto argsd = Args("DebugLevel", -1);
                applyExp(Heff, phi, -Cplx_i * beta * tstep / 2., argsd);

                //splitting of the evolved two-site tensor to single sites using SVD with some truncation error and some maximum dimension
                auto u = commonInds(psi(b), phi);
                auto [A, D, B] = svd(phi, u, args);

                //the norm of the state is just the norm of the singular value matrix D,  norm(D)

                if (ha == 1) //sweeping right
                {
                    //moving gauge to the right
                    psi.ref(b) = A;
                    psi.ref(b + 1) = D * B;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    auto phi = psi(b + 1);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    if (b == N - 1)
                    {
                        tstepq = 0;
                    }
                    if (b < N - 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * beta * tstepq / 2.);
                    psi.ref(b + 1) = phi;
                }
                else if (ha == 2) //sweeping left
                {
                    psi.ref(b) = D * A;
                    psi.ref(b + 1) = B;
                    Heff.numCenter(1);
                    Heff.position(b, psi);
                    auto phi = psi(b);
                    if (b == 1)
                    {
                        tstepq = 0;
                    }
                    if (b > 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * beta * tstepq / 2.);
                    psi.ref(b) = phi;
                }

            } // for loop over b
        }

        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                //Effective Hamiltonian for sites b,b+1
                Heff.numCenter(2);
                Heff.position(b, psi);

                //create two site tensor
                auto phi = psi(b) * psi(b + 1);

                //Krylov evolution of two site tensor and effective Hamiltonian
                auto argsd = Args("DebugLevel", -1);
                applyExp(Heff, phi, -Cplx_i * alpha * tstep / 2., argsd);

                //splitting of the evolved two-site tensor to single sites using SVD with some truncation error and some maximum dimension
                auto u = commonInds(psi(b), phi);
                auto [A, D, B] = svd(phi, u, args);

                //the norm of the state is just the norm of the singular value matrix D,  norm(D)

                if (ha == 1) //sweeping right
                {
                    //moving gauge to the right
                    psi.ref(b) = A;
                    psi.ref(b + 1) = D * B;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    auto phi = psi(b + 1);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    if (b == N - 1)
                    {
                        tstepq = 0;
                    }
                    if (b < N - 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * alpha * tstepq / 2.);
                    psi.ref(b + 1) = phi;
                }
                else if (ha == 2) //sweeping left
                {
                    psi.ref(b) = D * A;
                    psi.ref(b + 1) = B;
                    Heff.numCenter(1);
                    Heff.position(b, psi);
                    auto phi = psi(b);
                    if (b == 1)
                    {
                        tstepq = 0;
                    }
                    if (b > 1)
                    {
                        tstepq = tstep;
                    }
                    applyExp(Heff, phi, +Cplx_i * alpha * tstepq / 2.);
                    psi.ref(b) = phi;
                }

            } // for loop over b
        }
    }
    return 0;
}

template <class LocalOpT>
Real TDVP1(MPS &psi,
           LocalOpT &H,
           Cplx t,
           Args const &args)
{
    const int N = length(psi);
    Real alpha = 1. / (2. - pow(2., 0.33333333333));
    Real beta = 1. - 2 * alpha;

    auto Heff = LocalMPO(H);

    auto sweeps = Sweeps(1);

    Cplx tstep = t;
    Cplx tstepq;
    if (args.getInt("IntegrationOrder") == 2)
    {
        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                if (ha == 1) //sweeping right
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b, psi);

                    //create two site tensor
                    auto phi = psi(b);

                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b == 1)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b > 1)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * tstepq / 2., argsd);
                   
                    //moving gauge to the right
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(0);
                    Heff.position(b + 1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * tstepq / 2.);
                    psi.ref(b + 1) *= BondTensor;
                }
                else if (ha == 2) //sweeping left
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    //create two site tensor
                    auto phi = psi(b + 1);
                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b + 1 == N)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b + 1 < N)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * tstepq / 2., argsd);

                    //moving gauge to the left
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b + 1) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b
                    Heff.numCenter(0);
                    Heff.position(b+1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * tstep / 2.);
                    psi.ref(b) *= BondTensor;
                }

            } // for loop over b
        }
    }

    if (args.getInt("IntegrationOrder") == 4)
    {

        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                if (ha == 1) //sweeping right
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b, psi);

                    //create two site tensor
                    auto phi = psi(b);

                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b == 1)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b > 1)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * alpha * tstepq / 2., argsd);
                    //moving gauge to the right
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(0);
                    Heff.position(b + 1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * alpha * tstepq / 2.);
                    psi.ref(b + 1) *= BondTensor;
                }
                else if (ha == 2) //sweeping left
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    //create two site tensor
                    auto phi = psi(b + 1);
                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b + 1 == N)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b + 1 < N)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * alpha * tstepq / 2., argsd);

                    //moving gauge to the left
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b + 1) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b
                    Heff.numCenter(0);
                    Heff.position(b+1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * alpha * tstep / 2.);
                    psi.ref(b) *= BondTensor;
                }

            } // for loop over b
        }

        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                if (ha == 1) //sweeping right
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b, psi);

                    //create two site tensor
                    auto phi = psi(b);

                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b == 1)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b > 1)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * beta * tstepq / 2., argsd);
                    //moving gauge to the right
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(0);
                    Heff.position(b + 1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * beta * tstepq / 2.);
                    psi.ref(b + 1) *= BondTensor;
                }
                else if (ha == 2) //sweeping left
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    //create two site tensor
                    auto phi = psi(b + 1);
                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b + 1 == N)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b + 1 < N)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * beta * tstepq / 2., argsd);

                    //moving gauge to the left
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b + 1) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b
                    Heff.numCenter(0);
                    Heff.position(b+1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * beta * tstep / 2.);
                    psi.ref(b) *= BondTensor;
                }

            } // for loop over b
        }

        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

                if (ha == 1) //sweeping right
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b, psi);

                    //create two site tensor
                    auto phi = psi(b);

                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b == 1)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b > 1)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * alpha * tstepq / 2., argsd);
                    //moving gauge to the right
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b+1
                    Heff.numCenter(0);
                    Heff.position(b + 1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * alpha * tstepq / 2.);
                    psi.ref(b + 1) *= BondTensor;
                }
                else if (ha == 2) //sweeping left
                {
                    // guarantee correct mixed canonical gauge around site b
                    Heff.numCenter(1);
                    Heff.position(b + 1, psi);
                    //create two site tensor
                    auto phi = psi(b + 1);
                    //Krylov evolution of two site tensor and effective Hamiltonian
                    auto argsd = Args("DebugLevel", -1);
                    if (b + 1 == N)
                    {
                        tstepq = 2 * tstep;
                    };
                    if (b + 1 < N)
                    {
                        tstepq = tstep;
                    };
                    applyExp(Heff, phi, -Cplx_i * alpha * tstepq / 2., argsd);

                    //moving gauge to the left
                    auto u = commonInds(psi(b), psi(b + 1));
                    auto [A, D, B] = svd(phi, u, args);
                    psi.ref(b + 1) = B;
                    auto BondTensor = D * A;

                    //Effective Hamiltonian for site b
                    Heff.numCenter(0);
                    Heff.position(b+1, psi);

                    //Backward evolution of site b+1
                    //The ifs guarantee that we dont evolve backwards at the boundary
                    tstepq = tstep;
                    applyExp(Heff, BondTensor, +Cplx_i * alpha * tstep / 2.);
                    psi.ref(b) *= BondTensor;
                }

            } // for loop over b
        }
    }

    return 0;
}


//works for spin-1/2 remember to define afterwards
template <class LocalOpT>
Real TDVPh(MPS &psi,
           LocalOpT &H,
           Cplx t,
           Args const &args)
{
    const int N = length(psi);
    Real alpha = 1. / (2. - pow(2., 0.33333333333));
    Real beta = 1. - 2 * alpha;
    //Print(1);
    auto Heff = LocalMPO(H);

    auto sweeps = Sweeps(1);

    Cplx tstep = t;
    Cplx tstepq;
    int Maxdim =args.getInt("MaxDim"); 
    int localdim = args.getInt("LocalDim");
    int localdimmax [N-1];
    int bondtest [N-1]; 

    for (int i = 0; i < N-1; i++)
    {     
        //min(i+1, N-1 - i) min distance from left-right boundary  
        localdimmax[i] = (int) pow(localdim, (int) fmin(i+1, N-1 - i));
        bondtest[i] = 0;
        //Print(localdimmax[i]);
    }
    
    //Print(1);
    if (args.getInt("IntegrationOrder") == 2)
    {
        for (auto sw : range1(sweeps.nsweep()))
        {
            //Loop over bonds
            for (int b = 1, ha = 1; ha != 3; sweepnext(b, ha, N))
            {
                //Print(localdimmax[b-1]);
                //printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);
                auto u = commonInds(psi(b), psi(b+1));
                //Print(b);
                if ((dim(u) < localdimmax[b-1] && dim(u) < Maxdim) ||(bondtest[b-1] == 0 && ha == 2)/*|| (b == 1 || b == N -1)*/)
                {
                    //Print(b);
                    //Print(localdimmax[b]);
                    //Effective Hamiltonian for sites b,b+1
                    Heff.numCenter(2);
                    Heff.position(b, psi);

                    //create two site tensor
                    auto phi = psi(b) * psi(b + 1);

                    //Krylov evolution of two site tensor and effective Hamiltonian
                    applyExp(Heff, phi, -Cplx_i * tstep / 2.);

                    //splitting of the evolved two-site tensor to single sites using SVD with some truncation error and some maximum dimension
                    auto u = commonInds(psi(b), phi);
                    auto [A, D, B] = svd(phi, u, args);

                    //the norm of the state is just the norm of the singular value matrix D,  norm(D)

                    if (ha == 1) //sweeping right
                    {
                        //moving gauge to the right
                        psi.ref(b) = A;
                        psi.ref(b + 1) = D * B;

                        //Effective Hamiltonian for site b+1
                        Heff.numCenter(1);
                        Heff.position(b + 1, psi);
                        auto phi = psi(b + 1);

                        //Backward evolution of site b+1
                        //The ifs guarantee that we dont evolve backwards at the boundary
                        if (b == N - 1)
                        {
                            tstepq = 0;
                        }
                        if (b < N - 1)
                        {
                            tstepq = tstep;
                        }
                        applyExp(Heff, phi, +Cplx_i * tstepq / 2.);
                        psi.ref(b + 1) = phi;
                    }
                    else if (ha == 2) //sweeping left
                    {
                        psi.ref(b) = D * A;
                        psi.ref(b + 1) = B;
                        Heff.numCenter(1);
                        Heff.position(b, psi);
                        auto phi = psi(b);
                        if (b == 1)
                        {
                            tstepq = 0;
                        }
                        if (b > 1)
                        {
                            tstepq = tstep;
                        }
                    applyExp(Heff, phi, +Cplx_i * tstepq / 2.);
                    psi.ref(b) = phi;
                    }

                }
                if ((dim(u) >=  localdimmax[b-1] || dim(u) >=  Maxdim) /*&& (b >1 && b< N -1)*/)
                {
                    
                    if (ha == 1) //sweeping right
                    {
                        //Print(b);
                        bondtest[b-1] = 1;
                        // guarantee correct mixed canonical gauge around site b
                        Heff.numCenter(1);
                        Heff.position(b, psi);
                        
                        //create two site tensor
                        auto phi = psi(b);

                        //Krylov evolution of two site tensor and effective Hamiltonian
                        if (b == 1)
                        {
                            tstepq = 2 * tstep;
                        };
                        if (b > 1)
                        {
                            tstepq = tstep;
                        };
                        
                        applyExp(Heff, phi, -Cplx_i * tstepq / 2.);
                        
                        //moving gauge to the right
                        auto u = commonInds(psi(b), psi(b + 1));
                        auto [A, D, B] = svd(phi, u, {"Truncate",false});
                        psi.ref(b) = B;
                        auto BondTensor = D * A;

                        //Effective Hamiltonian for site b+1
                        Heff.numCenter(0);
                        Heff.position(b + 1, psi);
                        //Print(2);
                        //Backward evolution of site b+1
                        tstepq = tstep;
                        applyExp(Heff, BondTensor, +Cplx_i * tstepq / 2.);
                        psi.ref(b + 1) *= BondTensor;
                    }
                    if (ha == 2 && bondtest[b-1]== 1) //sweeping left
                    {
                        //Print(-b);
                        // guarantee correct mixed canonical gauge around site b
                        Heff.numCenter(1);
                        Heff.position(b + 1, psi);
                        //create two site tensor
                        auto phi = psi(b + 1);
                        //Krylov evolution of two site tensor and effective Hamiltonian
                        if (b + 1 == N)
                        {
                            tstepq = 2 * tstep;
                        };
                        if (b + 1 < N)
                        {
                            tstepq = tstep;
                        };
                        applyExp(Heff, phi, -Cplx_i * tstepq / 2.);  

                        //moving gauge to the left
                        auto u = commonInds(psi(b), psi(b + 1));
                        //Print(1);
                        auto [A, D, B] = svd(phi, u, {"Truncate",false});
                        //Print(2);
                        psi.ref(b + 1) = B;
                        auto BondTensor = D * A;

                        //Effective Hamiltonian for site b
                        Heff.numCenter(0);
                        Heff.position(b+1, psi);

                        //Backward evolution of site b+1
                        //The ifs guarantee that we dont evolve backwards at the boundary
                        tstepq = tstep;
                        applyExp(Heff, BondTensor, +Cplx_i * tstep / 2.);
                        psi.ref(b) *= BondTensor;
                    }
                }
            } // for loop over b
        }
    }
return 0;
}
} // namespace itensor

#endif
