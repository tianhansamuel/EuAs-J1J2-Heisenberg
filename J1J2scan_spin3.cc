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


cout << "Number of X unit cells Nx";
cin >> Nx;
cout << "The magnetic coupling J2";
cin >> J2;

ofstream myfile;
myfile.open ("Nx_"+std::to_string(Nx)+"_J2_"+std::to_string(J2)+"spin3.dat");

int N = 2*(Nx);  

int LL=50;


//Create 1d Heisenberg Hamiltonian

for (int ii = 1; ii <=LL-1; ++ii)
{

auto sites = SpinThreeHalf(N,{"ConserveQNs=",false});

double hh=3.0*ii/LL;

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

auto sweeps = Sweeps(30);
sweeps.maxdim() = 50,100,200,300,400,500;
sweeps.cutoff() = 9E-6;

sweeps.niter() = 2;
sweeps.noise() = 1E-5,1E-6,1E-7,1E-8,0.0;


auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet",true});


myfile << hh << "\t";
    //
    // Measuring S.S
    //

    //Sum total S.S to check that it's 
    //equal to ground state energy
    Real totalSdS = 0.;

    println("\nj S.S = ");
    for( auto b : range1(N-1) ) 
        { 
        psi.position(b);

        auto bondket = psi(b)*psi(b+1);
        auto bondbra = dag(prime(bondket,"Site")); 

        auto zzop = op(sites,"Sz",b)*op(sites,"Sz",b+1); 
        auto pmop = 0.5*op(sites,"S+",b)*op(sites,"S-",b+1); 
        auto mpop = 0.5*op(sites,"S-",b)*op(sites,"S+",b+1); 

        auto zz = elt(bondbra*zzop*bondket);
        auto pm = elt(bondbra*pmop*bondket);
        auto mp = elt(bondbra*mpop*bondket);

        printfln("%d %.12f",b,zz+pm+mp);

        myfile << b << "\t";
        myfile << zz << "\t";
        myfile << pm << "\t";
        myfile << mp << "\t";
        myfile << zz+pm+mp << "\t";
        }

    for( auto b : range1(N-2) ) 
        { 
        psi.position(b);

auto op_i = op(sites,"Sz",b);
auto op_j = op(sites,"Sz",b+2);

auto op_ip = op(sites,"Sp",b);
auto op_jm = op(sites,"Sm",b+2);


auto op_im = op(sites,"Sm",b);
auto op_jp = op(sites,"Sp",b+2);

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(b); 

//Create the bra/dual version of the MPS psi
auto psidag = dag(psi);

//Prime the link indices to make them distinct from
//the original ket links
psidag.prime("Link");

//index linking i-1 to i:
auto li_1 = leftLinkIndex(psi,b);

auto C = prime(psi(b),li_1)*op_i;
C *= prime(psidag(b),"Site");

auto C1 = prime(psi(b),li_1)*op_ip;
C1 *= prime(psidag(b),"Site");

auto C2 = prime(psi(b),li_1)*op_im;
C2 *= prime(psidag(b),"Site");


for(int k = b+1; k < b+2; ++k)
    {
    C *= psi(k);
    C *= psidag(k);

    C1 *= psi(k);
    C1 *= psidag(k);

    C2 *= psi(k);
    C2 *= psidag(k);

    }
//index linking j to j+1:
auto lj = rightLinkIndex(psi,b+2);

C *= prime(psi(b+2),lj)*op_j;
C *= prime(psidag(b+2),"Site");

C1 *= prime(psi(b+2),lj)*op_jm;
C1 *= prime(psidag(b+2),"Site");

C2 *= prime(psi(b+2),lj)*op_jp;
C2 *= prime(psidag(b+2),"Site");

auto result = elt(C); //or eltC(C) if expecting complex

auto result1 = 0.5*elt(C1); //or eltC(C) if expecting complex

auto result2 = 0.5*elt(C2); //or eltC(C) if expecting complex


myfile << b << "\t";
myfile << result << "\t";
myfile << result1 << "\t";
myfile << result2 << "\t";
myfile << result+result1+result2 << "\t";


printfln("%d %.12f",b,result+result1+result2);

}

/////////////////////////////////////////////////////////


    for( auto b : range1(N-5) ) 
        { 
        psi.position(b);

auto op_i = op(sites,"Sz",b);
auto op_j = op(sites,"Sz",b+5);

auto op_ip = op(sites,"Sp",b);
auto op_jm = op(sites,"Sm",b+5);


auto op_im = op(sites,"Sm",b);
auto op_jp = op(sites,"Sp",b+5);

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(b); 

//Create the bra/dual version of the MPS psi
auto psidag = dag(psi);

//Prime the link indices to make them distinct from
//the original ket links
psidag.prime("Link");

//index linking i-1 to i:
auto li_1 = leftLinkIndex(psi,b);

auto C = prime(psi(b),li_1)*op_i;
C *= prime(psidag(b),"Site");

auto C1 = prime(psi(b),li_1)*op_ip;
C1 *= prime(psidag(b),"Site");

auto C2 = prime(psi(b),li_1)*op_im;
C2 *= prime(psidag(b),"Site");


for(int k = b+1; k < b+5; ++k)
    {
    C *= psi(k);
    C *= psidag(k);

    C1 *= psi(k);
    C1 *= psidag(k);

    C2 *= psi(k);
    C2 *= psidag(k);

    }
//index linking j to j+1:
auto lj = rightLinkIndex(psi,b+5);

C *= prime(psi(b+5),lj)*op_j;
C *= prime(psidag(b+5),"Site");

C1 *= prime(psi(b+5),lj)*op_jm;
C1 *= prime(psidag(b+5),"Site");

C2 *= prime(psi(b+5),lj)*op_jp;
C2 *= prime(psidag(b+5),"Site");

auto result = elt(C); //or eltC(C) if expecting complex

auto result1 = 0.5*elt(C1); //or eltC(C) if expecting complex

auto result2 = 0.5*elt(C2); //or eltC(C) if expecting complex


myfile << b << "\t";
myfile << result << "\t";
myfile << result1 << "\t";
myfile << result2 << "\t";
myfile << result+result1+result2 << "\t";


printfln("%d %.12f",b,result+result1+result2);

}
////////////////////////////////////////////////////////////////////////////////////////////


    for( auto b : range1(N-6) ) 
        { 
        psi.position(b);

auto op_i = op(sites,"Sz",b);
auto op_j = op(sites,"Sz",b+6);

auto op_ip = op(sites,"Sp",b);
auto op_jm = op(sites,"Sm",b+6);


auto op_im = op(sites,"Sm",b);
auto op_jp = op(sites,"Sp",b+6);

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(b); 

//Create the bra/dual version of the MPS psi
auto psidag = dag(psi);

//Prime the link indices to make them distinct from
//the original ket links
psidag.prime("Link");

//index linking i-1 to i:
auto li_1 = leftLinkIndex(psi,b);

auto C = prime(psi(b),li_1)*op_i;
C *= prime(psidag(b),"Site");

auto C1 = prime(psi(b),li_1)*op_ip;
C1 *= prime(psidag(b),"Site");

auto C2 = prime(psi(b),li_1)*op_im;
C2 *= prime(psidag(b),"Site");


for(int k = b+1; k < b+6; ++k)
    {
    C *= psi(k);
    C *= psidag(k);

    C1 *= psi(k);
    C1 *= psidag(k);

    C2 *= psi(k);
    C2 *= psidag(k);

    }
//index linking j to j+1:
auto lj = rightLinkIndex(psi,b+6);

C *= prime(psi(b+6),lj)*op_j;
C *= prime(psidag(b+6),"Site");

C1 *= prime(psi(b+6),lj)*op_jm;
C1 *= prime(psidag(b+6),"Site");

C2 *= prime(psi(b+6),lj)*op_jp;
C2 *= prime(psidag(b+6),"Site");

auto result = elt(C); //or eltC(C) if expecting complex

auto result1 = 0.5*elt(C1); //or eltC(C) if expecting complex

auto result2 = 0.5*elt(C2); //or eltC(C) if expecting complex


myfile << b << "\t";
myfile << result << "\t";
myfile << result1 << "\t";
myfile << result2 << "\t";
myfile << result+result1+result2 << "\t";


printfln("%d %.12f",b,result+result1+result2);

}


myfile << energy << "\n";
}
   return 0;
    }


