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

Tensor build_tau(Tensor& T1, Tensor& T2) {
    Tensor tau = Tensor::build(TensorType::CoreTensor, "Intermediate tau array", T2.dims());
    tau("ijab") = T1("ia")*T1("jb");
    tau("ijab") += T2("ijab");
    return tau;
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
    C2p("pvug") = Voovv("iuag")*tau("ivga");
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
    size_t v = tau.dims()[2];
    Tensor D2p = Tensor::build(TensorType::CoreTensor, "Intermediate D2' array", {v,v,v,v});
    D2p("uvij") = Vovov("iajb")*tau("uvab");
    return D2p;
}

Tensor build_D2ps(Tensor& tau, Tensor& Vovov) {
    size_t o = tau.dims()[0];
    Tensor D2ps = Tensor::build(TensorType::CoreTensor, "Intermediate D2*' array", {o,o,o,o});
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

std::vector<double> cc_iterate(Tensor& T1, Tensor& T2, Tensor& auxD1, Tensor& auxD2, double& Ecc) {
    std::vector<double> = out;
    Tensor T1new = Tensor::build(TensorType::CoreTensor, "Updated T1 amplitudes", T1.dims());
    Tensor T2new = Tensor::build(TensorType::CoreTensor, "Updated T2 amplitudes", T2.dims());
    
    return out;
}
    
