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

#include "aux.h"
#include "helper_ccsd.h"
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
#include "ambit/helpers/psi4/convert.h"
#include <vector>

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
    SharedVector epsilon = ref_wfn->epsilon_a();

// Get integrals

    MintsHelper mints = MintsHelper(ref_wfn->basisset());
    SharedMatrix psi4TEI = mints.mo_eri(C, C, C, C);

// Build the six MO space cases for the Two-electron repulsion integral

    size_t s_o = (size_t) ndocc;
    size_t s_v = (size_t) nvir;

    ambit::initialize();

    // => Case 1: V(o,o,o,o)
    ambit::Tensor Voooo = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Two-electron Repulsion integral. MO space (o,o,o,o)", {s_o, s_o, s_o, s_o});
    Voooo.iterate([psi4TEI, nmo](const std::vector<size_t>& indices,double& value){
        // General MO indices
        size_t p = indices[0];
        size_t q = indices[1];
        size_t r = indices[2];
        size_t s = indices[3];

        // Compound indices used in the Psi4 Matrix object
        size_t x = p*nmo + q;
        size_t y = r*nmo + s;
        value = psi4TEI->get(x,y);
    });

    // => Case 2: V(v,v,v,v)
    ambit::Tensor Vvvvv = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Two-electron Repulsion integral. MO space (v,v,v,v)", {s_v, s_v, s_v, s_v});
    Vvvvv.iterate([psi4TEI, nmo, ndocc](const std::vector<size_t>& indices,double& value){
        // General MO indices. Virtuals need to be shifted by ndocc
        size_t p = indices[0] + ndocc;
        size_t q = indices[1] + ndocc;
        size_t r = indices[2] + ndocc;
        size_t s = indices[3] + ndocc;

        // Compound indices used in the Psi4 Matrix object
        size_t x = p*nmo + q;
        size_t y = r*nmo + s;
        value = psi4TEI->get(x,y);
    });

    // => Case 3: V(o,o,o,v)
    ambit::Tensor Vooov = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Two-electron Repulsion integral. MO space (o,o,o,v)", {s_o, s_o, s_o, s_v});
    Vooov.iterate([psi4TEI, nmo, ndocc](const std::vector<size_t>& indices,double& value){
        // General MO indices. Virtuals need to be shifted by ndocc
        size_t p = indices[0];
        size_t q = indices[1];
        size_t r = indices[2];
        size_t s = indices[3] + ndocc;

        // Compound indices used in the Psi4 Matrix object
        size_t x = p*nmo + q;
        size_t y = r*nmo + s;
        value = psi4TEI->get(x,y);
    });

    // => Case 4: V(o,o,v,v)
    ambit::Tensor Voovv = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Two-electron Repulsion integral. MO space (o,o,v,v)", {s_o, s_o, s_v, s_v});
    Voovv.iterate([psi4TEI, nmo, ndocc](const std::vector<size_t>& indices,double& value){
        // General MO indices. Virtuals need to be shifted by ndocc
        size_t p = indices[0];
        size_t q = indices[1];
        size_t r = indices[2] + ndocc;
        size_t s = indices[3] + ndocc;

        // Compound indices used in the Psi4 Matrix object
        size_t x = p*nmo + q;
        size_t y = r*nmo + s;
        value = psi4TEI->get(x,y);
    });

    // => Case 5: V(o,v,o,v)
    ambit::Tensor Vovov = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Two-electron Repulsion integral. MO space (o,v,o,v)", {s_o, s_v, s_o, s_v});
    Vovov.iterate([psi4TEI, nmo, ndocc](const std::vector<size_t>& indices,double& value){
        // General MO indices. Virtuals need to be shifted by ndocc
        size_t p = indices[0];
        size_t q = indices[1] + ndocc;
        size_t r = indices[2];
        size_t s = indices[3] + ndocc;

        // Compound indices used in the Psi4 Matrix object
        size_t x = p*nmo + q;
        size_t y = r*nmo + s;
        value = psi4TEI->get(x,y);
    });

    // => Case 6: V(o,v,v,v)
    ambit::Tensor Vovvv = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Two-electron Repulsion integral. MO space (o,v,v,v)", {s_o, s_v, s_v, s_v});
    Vovvv.iterate([psi4TEI, nmo, ndocc](const std::vector<size_t>& indices,double& value){
        // General MO indices. Virtuals need to be shifted by ndocc
        size_t p = indices[0];
        size_t q = indices[1] + ndocc;
        size_t r = indices[2] + ndocc;
        size_t s = indices[3] + ndocc;

        // Compound indices used in the Psi4 Matrix object
        size_t x = p*nmo + q;
        size_t y = r*nmo + s;
        value = psi4TEI->get(x,y);
    });

    ambit::Tensor T1 = ambit::Tensor::build(ambit::TensorType::CoreTensor, "T1 amplitudes", {s_o, s_v});
    ambit::Tensor D1 = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Inverse D1 auxiliar", {s_o, s_v});

    D1.iterate([epsilon, ndocc, nvir](const std::vector<size_t>& indices,double& value){
        size_t i = indices[0];
        size_t a = indices[1] + ndocc;
        value = 1.0/((*epsilon)[a] - (*epsilon)[i]);
    });

    ambit::Tensor T2 = ambit::Tensor::build(ambit::TensorType::CoreTensor, "T2 amplitudes", {s_o, s_o, s_v, s_v});
    ambit::Tensor D2 = ambit::Tensor::build(ambit::TensorType::CoreTensor, "Inverse D2 auxiliar", {s_o, s_o, s_v, s_v});

    D2.iterate([epsilon, ndocc, nvir](const std::vector<size_t>& indices,double& value){
        size_t i = indices[0];
        size_t j = indices[1];
        size_t a = indices[2] + ndocc;
        size_t b = indices[3] + ndocc;
        value = 1.0/((*epsilon)[i] + (*epsilon)[j] - (*epsilon)[a] - (*epsilon)[b]);
    });

// MP2 guess for T2

    T2("ijab") = Vovov("iajb")*D2("ijab");
    double Ecc = cc_energy(T1, T2, Vovov);

    std::cout << "CC MP2 Energy: " << Ecc << std::endl;
    std::vector<size_t> v = T2.dims();
    std::cout << v[0] << std::endl;

    

    ambit::finalize();
    return ref_wfn;

}

}} // End Namespaces
