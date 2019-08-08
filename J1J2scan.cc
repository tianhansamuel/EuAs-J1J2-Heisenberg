#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <math.h>       /* exp */
using namespace std;
using namespace itensor;
int main()
    {

int Nx;

double J2;

double hh;

int LL=400;
cout << "Number of X unit cells Nx";
cin >> Nx;
cout << "The magnetic coupling J2";
cin >> J2;

ofstream myfile;
myfile.open ("Nx_"+std::to_string(Nx)+"_J2_"+std::to_string(J2)+"f.dat");

int N = 2*(Nx);  
auto sites = SpinHalf(N,{"ConserveQNs=",false});


for (int ip=1; ip<=LL;++ip)
{
hh=0.6*ip/LL;
//Create 1d Heisenberg Hamiltonian

auto ampo = AutoMPO(sites);
for(int jj = 1; jj <=Nx-1; ++jj)
{
    int jj1=2*jj-1;
    int jj2=2*jj;
    int jj3=2*jj+1;
    int jj4=2*jj+2;

    println(jj1,jj2,jj3,jj4);
    ampo += 0.5*hh,"S-",jj1;
    ampo += 0.5*hh,"S+",jj1;
    ampo += 0.5*hh,"S-",jj2;
    ampo += 0.5*hh,"S+",jj2;

    ampo += 0.5,"S+",jj1,"S-",jj2;
    ampo += 0.5,"S-",jj1,"S+",jj2;
    ampo += "Sz",jj1,"Sz",jj2;

    ampo += 0.5,"S+",jj2,"S-",jj3;
    ampo += 0.5,"S-",jj2,"S+",jj3;
    ampo += "Sz",jj2,"Sz",jj3;

    ampo += 0.5*J2,"S+",jj1,"S-",jj3;
    ampo += 0.5*J2,"S-",jj1,"S+",jj3;
    ampo += J2,"Sz",jj1,"Sz",jj3;

    ampo += 0.5*J2,"S+",jj2,"S-",jj4;
    ampo += 0.5*J2,"S-",jj2,"S+",jj4;
    ampo += J2,"Sz",jj2,"Sz",jj4;

}

ampo += 0.5*hh,"S-",2*Nx-1;
ampo += 0.5*hh,"S+",2*Nx-1;
ampo += 0.5*hh,"S-",2*Nx;
ampo += 0.5*hh,"S+",2*Nx;

ampo += 0.5,"S+",2*Nx-1,"S-",2*Nx;
ampo += 0.5,"S-",2*Nx-1,"S+",2*Nx;
ampo += "Sz",2*Nx-1,"Sz",2*Nx;

ampo += 0.5,"S+",2*Nx,"S-",1;
ampo += 0.5,"S-",2*Nx,"S+",1;
ampo += "Sz",2*Nx,"Sz",1;

ampo += 0.5*J2,"S+",2*Nx-1,"S-",1;
ampo += 0.5*J2,"S-",2*Nx-1,"S+",1;
ampo += J2,"Sz",2*Nx-1,"Sz",1;

ampo += 0.5*J2,"S+",2*Nx,"S-",2;
ampo += 0.5*J2,"S-",2*Nx,"S+",2;
ampo += J2,"Sz",2*Nx,"Sz",2;

auto H = toMPO(ampo);

//Set up random initial wavefunction
auto state = InitState(sites);
for(int i = 1; i <= N; ++i) 
        {
    if (i%2==0) {state.set(i,"Up");}
    else {state.set(i,"Dn");}
        }
auto psi0 = MPS(state);

//Perform 5 sweeps of DMRG
auto sweeps = Sweeps(20);

sweeps.cutoff() = 1E-5,1E-6,1E-7,1E-8,1E-9,1E-10,1E-15;
//Specify max number of states kept each sweep
sweeps.maxdim() = 50, 50, 100, 100, 200,400,600,800,1000;

//Run the DMRG algorithm

auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
//Continue to analyze wavefunction afterward 


for(int j = 1; j <= N; ++j)
    {
    //Make site j the MPS "orthogonality center"
    psi.position(j);
    //Measure magnetization
    auto Szj =elt(dag(prime(psi(j),"Site"))*sites.op("Sz",j)*psi(j));
    println("Sz_",j," = ",Szj);

    auto Sxj =elt(dag(prime(psi(j),"Site"))*sites.op("Sx",j)*psi(j));
    println("Sx_",j," = ",Sxj);
myfile << hh << "\t";
myfile << j << "\t";
myfile << Sxj << "\t";
myfile << Szj << "\t";
myfile << energy << "\n";
    }
println(hh,"+",energy);
}
   return 0;
    }


