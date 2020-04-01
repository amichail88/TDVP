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
                    Heff.position(b, psi);

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
                    Heff.position(b, psi);

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
                    Heff.position(b, psi);

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
                    Heff.position(b, psi);

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
} // namespace itensor

#endif
