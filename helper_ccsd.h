#include "ambit/tensor.h"
#include <vector>

double cc_energy(ambit::Tensor& T1, ambit::Tensor& T2, ambit::Tensor& Vovov);

double T1_residue(ambit::Tensor& Told, ambit::Tensor& Tnew);

double T2_residue(ambit::Tensor& Told, ambit::Tensor& Tnew);

ambit::Tensor build_tau(ambit::Tensor& T1, ambit::Tensor& T2);

ambit::Tensor build_Te(ambit::Tensor& T1, ambit::Tensor& T2);

ambit::Tensor build_A2p(ambit::Tensor& tau, ambit::Tensor& Voooo);

ambit::Tensor build_B2p(ambit::Tensor& tau, ambit::Tensor& Vvvvv);

ambit::Tensor build_C1(ambit::Tensor& T1, ambit::Tensor& Voovv);

ambit::Tensor build_C2(ambit::Tensor& T2, ambit::Tensor& Voovv);

ambit::Tensor build_C2p(ambit::Tensor& tau, ambit::Tensor& Voovv);

ambit::Tensor build_D1(ambit::Tensor& T1, ambit::Tensor& Vovov);

ambit::Tensor build_D2p(ambit::Tensor& tau, ambit::Tensor& Vovov);

ambit::Tensor build_D2ps(ambit::Tensor& tau, ambit::Tensor& Vovov);

ambit::Tensor build_D2a(ambit::Tensor& T2, ambit::Tensor& Vovov);

ambit::Tensor build_D2b(ambit::Tensor& T2, ambit::Tensor& Vovov);

ambit::Tensor build_D2c(ambit::Tensor& T2, ambit::Tensor& Vovov);

ambit::Tensor build_E1s(ambit::Tensor& T1, ambit::Tensor& Vooov);

ambit::Tensor build_E1(ambit::Tensor& T1, ambit::Tensor& Vooov);

ambit::Tensor build_E2a(ambit::Tensor& T2, ambit::Tensor& Vooov);

ambit::Tensor build_E2b(ambit::Tensor& T2, ambit::Tensor& Vooov);

ambit::Tensor build_E2c(ambit::Tensor& T2, ambit::Tensor& Vooov);

ambit::Tensor build_F11(ambit::Tensor& T1, ambit::Tensor& Vovvv);

ambit::Tensor build_F12(ambit::Tensor& T1, ambit::Tensor& Vovvv);

ambit::Tensor build_F1s(ambit::Tensor& T1, ambit::Tensor& Vovvv);

ambit::Tensor build_F2a(ambit::Tensor& T2, ambit::Tensor& Vovvv);

ambit::Tensor build_F2p(ambit::Tensor& tau, ambit::Tensor& Vovvv);

ambit::Tensor build_giu(ambit::Tensor& E1, ambit::Tensor& D2p);

ambit::Tensor build_gap(ambit::Tensor& F1s, ambit::Tensor& D2ps);

std::vector<double> cc_iterate(ambit::Tensor& T1, ambit::Tensor& T2, ambit::Tensor& auxD1, ambit::Tensor& auxD2, double& Ecc,
    ambit::Tensor& Voooo, ambit::Tensor& Vvvvv, ambit::Tensor& Vooov, ambit::Tensor& Voovv, ambit::Tensor& Vovov, ambit::Tensor& Vovvv);
