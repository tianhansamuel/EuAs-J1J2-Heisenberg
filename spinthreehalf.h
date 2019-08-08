//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_SPINTHREEHALF_H
#define __ITENSOR_SPINTHREEHALF_H
#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class SpinThreeHalfSite;

using SpinThreeHalf = BasicSiteSet<SpinThreeHalfSite>;

class SpinThreeHalfSite
    {

  Index s;
	public:

    SpinThreeHalfSite(Index I) : s(I) { }

    SpinThreeHalfSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,S=7/2");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);
        if(conserveSz)
            {
            s = Index{QN({"Sz",+7}),1,
                      QN({"Sz",+5}),1,
                      QN({"Sz",+3}),1,
                      QN({"Sz",1}),1,
                      QN({"Sz",-1}),1,
                      QN({"Sz",-3}),1,
                      QN({"Sz",-5}),1,
                      QN({"Sz",-7}),1,Out,ts};
            }
        else
            {
            s = Index{8,ts};
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "Up") 
            {
            return s(1);
            }
        else 
        if(state == "Dn") 
            {
            return s(2);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args = Args::global()) const
        {
        const Real val0 = std::sqrt(7.0);
        const Real val1 = std::sqrt(12.0);
        const Real val2 = std::sqrt(15.0);
        const Real val3 = std::sqrt(16.0);

        auto sP = prime(s);
        auto Up  = s(1);
        auto UpP = sP(1);
        auto Up1 = s(2);
        auto Up1P = sP(2);
        auto Up2 = s(3);
        auto Up2P = sP(3);
        auto Up3 = s(4);
        auto Up3P = sP(4);
        auto Dn3 = s(5);
        auto Dn3P = sP(5);
        auto Dn2 = s(6);
        auto Dn2P = sP(6);
        auto Dn1 = s(7);
        auto Dn1P = sP(7);
        auto Dn  = s(8);
        auto DnP = sP(8);

        auto Op = ITensor(dag(s),sP);

        if(opname == "Sz")
            {
            Op.set(Up,UpP,+3.5);
            Op.set(Up1,Up1P,+2.5);
            Op.set(Up2,Up2P,+1.5);
            Op.set(Up3,Up3P,+0.5);
            Op.set(Dn3,Dn3P,-0.5);
            Op.set(Dn2,Dn2P,-1.5);
            Op.set(Dn1,Dn1P,-2.5);
            Op.set(Dn,DnP,-3.5);

            }
        else

        if(opname == "Sp" || opname == "S+")
            {

            Op.set(Up1,UpP,val0);
            Op.set(Up2,Up1P,val1); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Up3,Up2P,val2); // val1 = sqrt(6)/2 = = 1.2247...

            Op.set(Dn3,Up3P,val3); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Dn2,Dn3P,val2); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Dn1,Dn2P,val1); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Dn,Dn1P,val0);

            }
        else
        if(opname == "Sm" || opname == "S-")
            {

            Op.set(Up,Up1P,val0);
            Op.set(Up1,Up2P,val1); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Up2,Up3P,val2); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Up3,Dn3P,val3); // val1 = sqrt(6)/2 = = 1.2247...

            Op.set(Dn2,Dn3P,val2); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Dn1,Dn2P,val1); // val1 = sqrt(6)/2 = = 1.2247...
            Op.set(Dn,Dn1P,val0);
            }
        else
        if(opname == "projUp")
            {
            Op.set(Up,UpP,1);
            }
        else
        if(opname == "projUp1")
            {
            Op.set(Up1,Up1P,1);
            }
        else
        if(opname == "projU2")
            {
            Op.set(Up2,Up2P,1);
            }
        else
        if(opname == "projUp3")
            {
            Op.set(Up3,Up3P,1);
            }
        else


        if(opname == "projDn3")
            {
            Op.set(Dn3,Dn3P,1);
            }
        else

        if(opname == "projDn2")
            {
            Op.set(Dn2,Dn2P,1);
            }
        else
        if(opname == "projDn1")
            {
            Op.set(Dn1,Dn1P,1);
            }
        else
        if(opname == "projDn")
            {
            Op.set(Dn,DnP,1);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    SpinThreeHalfSite(int n, Args const& args = Args::global())
        {
        *this = SpinThreeHalfSite({args,"SiteNumber=",n});
        }

    };

} //namespace itensor
#endif
