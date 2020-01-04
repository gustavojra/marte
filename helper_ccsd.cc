#include "helper_ccsd.h"
#include "ambit/tensor.h"
#include <iostream>
#include <vector>

using ambit::Tensor;
using ambit::TensorType;

double cc_energy(Tensor& T1, Tensor& T2, Tensor& Vovov) {
    double energy;
    Tensor tau = build_tau(T1, T2);
    energy = Vovov("iajb")*(2.0 * tau("ijab") - tau("jiab"));
    return energy;
}

double T1_residue(Tensor& Told, Tensor& Tnew) {
    int N = Told.dims()[0] * Told.dims()[1];
    Tensor Dif = Tensor::build(TensorType::CoreTensor, "Differences", Told.dims());
    Dif("up") = Told("up") - Tnew("up");
    return Dif.norm(2)/N;
}

double T2_residue(Tensor& Told, Tensor& Tnew) {
    int N = Told.dims()[0] * Told.dims()[2];
    Tensor Dif = Tensor::build(TensorType::CoreTensor, "Differences", Told.dims());
    Dif("uvpg") = Told("uvpg") - Tnew("uvpg");
    return Dif.norm(2)/(N*N);
}

Tensor build_tau(Tensor& T1, Tensor& T2) {
    Tensor tau = Tensor::build(TensorType::CoreTensor, "Intermediate tau array", T2.dims());
    tau("ijab") = T1("ia")*T1("jb");
    tau("ijab") += T2("ijab");
    return tau;
}

Tensor build_Te(Tensor& T1, Tensor& T2) {
    Tensor Te = Tensor::build(TensorType::CoreTensor, "Intermediate Te array", T2.dims());
    Te("ijab") = T1("ia")*T1("jb");
    Te("ijab") += 0.5*T2("ijab");
    return Te;
}
    
Tensor build_A2p(Tensor& tau, Tensor& Voooo) {
    Tensor A2p = Tensor::build(TensorType::CoreTensor, "Intermediate A2' array", tau.dims());
    A2p("uvpg") = Voooo("uivj")*tau("ijpg");
    return A2p;
}

Tensor build_B2p(Tensor& tau, Tensor& Vvvvv) {
    Tensor B2p = Tensor::build(TensorType::CoreTensor, "Intermediate B2' array", tau.dims());
    B2p("uvpg") = Vvvvv("apbg")*tau("uvab");
    return B2p;
}

Tensor build_C1(Tensor& T1, Tensor& Voovv) {
    size_t o = T1.dims()[0];
    size_t v = T1.dims()[1];
    Tensor C1 = Tensor::build(TensorType::CoreTensor, "Intermediate C1 array", {o, o, v});
    C1("uip") = Voovv("uiap")*T1("ia");
    return C1;
}

Tensor build_C2(Tensor& T2, Tensor& Voovv) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor C2 = Tensor::build(TensorType::CoreTensor, "Intermediate C2 array", {v,o,o,v});
    C2("pvug") = Voovv("uiap")*T2("viga");
    return C2;
}

Tensor build_C2p(Tensor& tau, Tensor& Voovv) {
    size_t o = tau.dims()[0];
    size_t v = tau.dims()[2];
    Tensor C2p = Tensor::build(TensorType::CoreTensor, "Intermediate C2' array", {v,o,o,v});
    C2p("pvug") = Voovv("iuag")*tau("ivpa");
    return C2p;
}

Tensor build_D1(Tensor& T1, Tensor& Vovov) {
    size_t o = T1.dims()[0];
    size_t v = T1.dims()[1];
    Tensor D1 = Tensor::build(TensorType::CoreTensor, "Intermediate D1 array", {o,o,v,o});
    //maybe double check this one
    D1("uvpi") = Vovov("upia")*T1("va");
    return D1;
}

Tensor build_D2p(Tensor& tau, Tensor& Vovov) {
    size_t o = tau.dims()[0];
    Tensor D2p = Tensor::build(TensorType::CoreTensor, "Intermediate D2' array", {o,o,o,o});
    D2p("uvij") = Vovov("iajb")*tau("uvab");
    return D2p;
}

Tensor build_D2ps(Tensor& tau, Tensor& Vovov) {
    size_t v = tau.dims()[2];
    Tensor D2ps = Tensor::build(TensorType::CoreTensor, "Intermediate D2*' array", {v,v,v,v});
    D2ps("acpb") = Vovov("iajc")*tau("ijpb");
    return D2ps;
}

Tensor build_D2a(Tensor& T2, Tensor& Vovov) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor D2a = Tensor::build(TensorType::CoreTensor, "Intermediate D2a array", {v,o,o,v});
    D2a("avig") = Vovov("jbia") * (2*T2("vjgb") - T2("vjbg"));
    return D2a;
}

Tensor build_D2b(Tensor& T2, Tensor& Vovov) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor D2b = Tensor::build(TensorType::CoreTensor, "Intermediate D2b array", {v,o,o,v});
    D2b("avig") = Vovov("ibja") * T2("vjgb");
    return D2b;
}

Tensor build_D2c(Tensor& T2, Tensor& Vovov) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor D2c = Tensor::build(TensorType::CoreTensor, "Intermediate D2c array", {v,o,o,v});
    D2c("avig") = Vovov("ibja") * T2("vjbg");
    return D2c;
}

Tensor build_E1s(Tensor& T1, Tensor& Vooov) {
    size_t o = T1.dims()[0];
    size_t v = T1.dims()[1];
    Tensor E1s = Tensor::build(TensorType::CoreTensor, "Intermediate E1* array", {o,o,v,v});
    E1s("uvpg") = Vooov("viup")*T1("ig");
    return E1s;
}

Tensor build_E1(Tensor& T1, Tensor& Vooov) {
    size_t o = T1.dims()[0];
    Tensor E1 = Tensor::build(TensorType::CoreTensor, "Intermediate E1 array", {o,o,o,o});
    E1("uvij") = Vooov("iuja")*T1("va");
    return E1;
}

Tensor build_E2a(Tensor& T2, Tensor& Vooov) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor E2a = Tensor::build(TensorType::CoreTensor, "Intermediate E2a array", {o,o,o,v});
    E2a("uvig") = Vooov("iujb")*(2*T2("vjgb")-T2("vjbg"));
    return E2a;
}

Tensor build_E2b(Tensor& T2, Tensor& Vooov) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor E2b = Tensor::build(TensorType::CoreTensor, "Intermediate E2c array", {o,o,o,v});
    E2b("uvig") = Vooov("juib")*T2("vjgb");
    return E2b;
}

Tensor build_E2c(Tensor& T2, Tensor& Vooov) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor E2c = Tensor::build(TensorType::CoreTensor, "Intermediate E2c array", {o,o,o,v});
    E2c("uvig") = Vooov("juib")*T2("vjbg");
    return E2c;
}

Tensor build_F11(Tensor& T1, Tensor& Vovvv) {
    size_t o = T1.dims()[0];
    size_t v = T1.dims()[1];
    Tensor F11 = Tensor::build(TensorType::CoreTensor, "Intermediate F11 array", {v,o,v,o});
    F11("bvpi") = Vovvv("iapb")*T1("va");
    return F11;
}

Tensor build_F12(Tensor& T1, Tensor& Vovvv) {
    size_t o = T1.dims()[0];
    size_t v = T1.dims()[1];
    Tensor F12 = Tensor::build(TensorType::CoreTensor, "Intermediate F12 array", {v,o,o,v});
    F12("bvip") = Vovvv("ibpa")*T1("va");
    return F12;
}

Tensor build_F1s(Tensor& T1, Tensor& Vovvv) {
    size_t v = T1.dims()[1];
    Tensor F1s = Tensor::build(TensorType::CoreTensor, "Intermediate F1* array", {v,v,v,v});
    F1s("acpb") = Vovvv("icpa")*T1("ib");
    return F1s;
}

Tensor build_F2a(Tensor& T2, Tensor& Vovvv) {
    size_t o = T2.dims()[0];
    size_t v = T2.dims()[2];
    Tensor F2a = Tensor::build(TensorType::CoreTensor, "Intermediate F2a array", {v,o,v});
    F2a("aup") = Vovvv("ibpa")*(2*T2("uiab")-T2("uiba"));
    return F2a;
}

Tensor build_F2p(Tensor& tau, Tensor& Vovvv) {
    size_t o = tau.dims()[0];
    size_t v = tau.dims()[2];
    Tensor F2p = Tensor::build(TensorType::CoreTensor, "Intermediate F2' array", {o,o,v,o});
    F2p("uvpi") = Vovvv("ibpa")*tau("uvab");
    return F2p;
}

Tensor build_giu(Tensor& E1, Tensor& D2p) {
    size_t o = E1.dims()[0];
    Tensor giu = Tensor::build(TensorType::CoreTensor, "Intermediate g(iu) array", {o,o});
    Tensor X = Tensor::build(TensorType::CoreTensor, "X Holding values", {o,o,o,o});
    Tensor Y = Tensor::build(TensorType::CoreTensor, "Y Holding values", {o,o,o,o});
    Y("pqrs") = E1("pqrs") + D2p("pqrs");
    X("pqrs") = 2*Y("pqrs") - Y("pqsr");
    Y.set(0);
    // Can this step be replaced? Contraction of type X('ujij') -> X('ui')
    Y.iterate([](const std::vector<size_t>& indices, double& value) {
        if(indices[1] == indices[3]){value = 1.0;};
    });
    giu("ui") = X("udij")*Y("udij");
    return giu;
}

Tensor build_gap(Tensor& F1s, Tensor& D2ps) {
    size_t v = F1s.dims()[0];
    Tensor gap = Tensor::build(TensorType::CoreTensor, "Intermediate g(ap) array", {v,v});
    Tensor X = Tensor::build(TensorType::CoreTensor, "X Holding values", {v,v,v,v});
    Tensor Y = Tensor::build(TensorType::CoreTensor, "Y Holding values", {v,v,v,v});
    Y("pqrs") = F1s("pqrs") - D2ps("pqrs");
    X("pqrs") = 2*Y("pqrs") - Y("qprs");
    Y.set(0);
    // Can this step be replaced? Contraction of type X('ujij') -> X('ui')
    Y.iterate([](const std::vector<size_t>& indices, double& value) {
        if(indices[1] == indices[3]){value = 1.0;};
    });
    gap("ap") = X("abpc")*Y("abpc");
    return gap;
}

std::vector<double> cc_iterate(Tensor& T1, Tensor& T2, Tensor& auxD1, Tensor& auxD2, double& Ecc,
    Tensor& Voooo, Tensor& Vvvvv, Tensor& Vooov, Tensor& Voovv, Tensor& Vovov, Tensor& Vovvv) {

    std::vector<double> out;
    Tensor T1new = Tensor::build(TensorType::CoreTensor, "Updated T1 amplitudes", T1.dims());
    Tensor T2new = Tensor::build(TensorType::CoreTensor, "Updated T2 amplitudes", T2.dims());
    Tensor tau  = build_tau(T1, T2);
    Tensor Te   = build_Te(T1, T2);
    Tensor A2p  = build_A2p(tau, Voooo);
    Tensor B2p  = build_B2p(tau, Vvvvv);
    Tensor C1   = build_C1(T1, Voovv);
    Tensor C2   = build_C2(T2, Voovv);
    Tensor C2p  = build_C2p(tau, Voovv);
    Tensor D1   = build_D1(T1, Vovov);
    Tensor D2p  = build_D2p(tau, Vovov);
    Tensor D2ps = build_D2ps(tau, Vovov);
    Tensor D2a  = build_D2a(T2, Vovov);
    Tensor D2b  = build_D2b(T2, Vovov);
    Tensor D2c  = build_D2c(T2, Vovov);
    Tensor E1s  = build_E1s(T1, Vooov);
    Tensor E1   = build_E1(T1, Vooov);
    Tensor E2a  = build_E2a(T2, Vooov);
    Tensor E2b  = build_E2b(T2, Vooov);
    Tensor E2c  = build_E2c(T2, Vooov);
    Tensor F11  = build_F11(T1, Vovvv);
    Tensor F12  = build_F12(T1, Vovvv);
    Tensor F1s  = build_F1s(T1, Vovvv);
    Tensor F2a  = build_F2a(T2, Vovvv);
    Tensor F2p  = build_F2p(tau, Vovvv);
    Tensor giu  = build_giu(E1, D2p);
    Tensor gap  = build_gap(F1s, D2ps);

    // Update T2

    Tensor J = Tensor::build(TensorType::CoreTensor, "J intermediate", T2.dims());
    Tensor S = Tensor::build(TensorType::CoreTensor, "S intermediate", T2.dims());
    
    J("uvpg") += gap("ag")*T2("uvpa");
    J("uvpg") -= giu("vi")*T2("uipg");

    // check this
    S("uvpg") = 0.5*(A2p("uvpg") + B2p("uvpg")) - E1s("uvpg") - C2("pvug") - C2p("pvug") + D2a("pvug") + F12("pvug");

    Tensor X = Tensor::build(TensorType::CoreTensor, "Holding", D2a.dims());
    Tensor Y = Tensor::build(TensorType::CoreTensor, "Holding", T2.dims());
    X("avig") = D2a("avig") - D2b("avig");
    Y("uipa") = T2("uipa") - Te("uiap");
    
    S("uvpg") += X("avig")*Y("uipa");
    S("uvpg") += 0.5*D2c("avig")*T2("uipa");
    S("uvpg") += D2c("auig")*Te("viap");

    S("uvpg") += (0.5*D2p("uvij") + E1("uvij"))*tau("ijpg");

    S("uvpg") -= (D1("uvpi") + F2p("uvpi"))*T1("ig");

    S("uvpg") -= (E2a("uvig") - E2b("uvig") - E2c("vuig"))*T1("ip");

    S("uvpg") -= F11("avgi")*T2("uipa"); 
    S("uvpg") -= F11("avpi")*T2("uiag");
    S("uvpg") += (2*T2("uipa") - T2("uiap"))*F12("avig");
    
    T2new("uvpg") += auxD2("uvpg")*(Vovov("upvg") + J("uvpg") + J("vugp") + S("uvpg") + S("vugp"));

    // Update T1

    Tensor K = Tensor::build(TensorType::CoreTensor, "Holding", T1.dims());

    K("up") += giu("ui")*T1("ip");
    K("up") -= gap("ap")*T1("ua");
    K("up") -= 2*D1("juai")*T1("ja")*T1("ip");
    K("up") += D1("iuaj")*T1("ja")*T1("ip");
    K("up") -= (2*(D2a("auip") - D2b("auip")) + D2c("auip"))*T1("ia");
    // Resolve this
    Tensor Onevov = Tensor::build(TensorType::CoreTensor, "Holding", F2a.dims());
    Onevov.set(1.0);
    K("up") -= F2a("aup")*Onevov("aup");
    
    Tensor Oneoov = Tensor::build(TensorType::CoreTensor, "Holding", C1.dims());
    Oneoov.set(1.0);

    // Can this step be replaced? Contraction of type X('ujij') -> X('ui')
    Tensor Oneooov = Tensor::build(TensorType::CoreTensor, "Holding", E2a.dims());
    Oneooov.iterate([](const std::vector<size_t>& indices, double& value) {
        if(indices[1] == indices[2]){value = 1.0;};
    });

    //careful with this
    Tensor Hooov = Tensor::build(TensorType::CoreTensor, "Holding", E2a.dims());
    Hooov("uisp") = 0.5*(E2a("uisp") - E2b("uisp")) + E2c("uisp") - 2*D1("uips");

    K("up") += Hooov("uisp")*Oneooov("uisp");
    K("up") += C1("uip")*Oneoov("uip");

    T1new("up") = auxD1("up")*K("up");

    // New energy

    double Ecc_new = cc_energy(T1new, T2new, Vovov);

    out.push_back(Ecc_new - Ecc);

    Ecc = Ecc_new;

    // T1 residue

    out.push_back(T1_residue(T1new, T1));

    // T2 residue

    out.push_back(T2_residue(T2new, T2));

    // Update amplitudes

    T1 = T1new;
    T2 = T2new;

    return out;
}
