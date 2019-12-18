/*
 * @BEGIN LICENSE
 *
 * marte by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

//#include "aux.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/mintshelper.h"
#include "ambit/tensor.h"
#include "ambit/blocked_tensor.h"

namespace psi{ namespace marte{
 
extern "C" PSI_API
int read_options(std::string name, Options &options)
{
    if (name == "MARTE" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}


extern "C" PSI_API
SharedWavefunction marte(SharedWavefunction ref_wfn, Options& options)
{

// Read Psi4 Options

    int print = options.get_int("PRINT");

    // Have the reference (SCF) wavefunction, ref_wfn
    if(!ref_wfn) throw PSIEXCEPTION("SCF has not been run yet!");

    psi::outfile->Printf("Running g-CCSD\n");

// Get SCF data
    int nelec = ref_wfn->nalpha() + ref_wfn->nbeta();
    if(nelec%2 != 0) throw PSIEXCEPTION("Only closed-shells systems are allowed for MARTE now");
    int nmo = ref_wfn->nmo();
    int ndocc = nelec/2;
    int nvir = nmo - ndocc;
    SharedMatrix C = ref_wfn->Ca();

// Get integrals

    // Not pointers for now ??
    MintsHelper mints = MintsHelper(ref_wfn->basisset());
    SharedMatrix psi4TEI = mints.mo_eri(C, C, C, C);

// Build ambit tensors

    std::vector<int> occ_space(ndocc), vir_space(nvir);

    for(int i = 0; i < ndocc; i++) {
        occ_space[i] = i;
    }

//    gprint(occ_space);

//    std::cout << "[ ";
//    for(int i = 0; i < occ_space.size(); i++) {
//        std::cout << occ_space[i] << " ";
//    }
//    std::cout << "]" << std::endl;

//    gprint(occ_space);

    ambit::SpinType spin_free = ambit::SpinType::NoSpin;
    ambit::initialize();
    ambit::BlockedTensor::add_mo_space("o","i,j,k,l",{0,1,2,3,4},spin_free);
    ambit::BlockedTensor::add_mo_space("v","a,b,c,d",{5,6,7,8,9},spin_free);
    ambit::BlockedTensor::add_composite_mo_space("g","p,q,r,s,t",{"o","v"});
//    ambit::Tensor TEI = ambit::Tensor::build(ambit::CoreTensor, "TEI", {100,100});
//    ambit::helpers::psi4::convert(psi4TEI, TEI);

    return ref_wfn;
}

}} // End Namespaces
