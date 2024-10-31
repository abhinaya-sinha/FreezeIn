//This C++ program exposes the functions in the FreezeIn.h library to python
//using a light-weight header-only pybind11 library (included in this
//repository)

/********************/
/* FreezeIn Library */
/********************/

#include "FreezeIn.h"

/******************/
/* Pybind Library */
/******************/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

/************************************/
/* Exposing C++ functions to python */
/************************************/

PYBIND11_MODULE(FreezeIn, mod)
{
    //Read_gstar(choice, gstarpath)
    mod.def("Read_gstar", &Read_gstar, R"pbdoc(
    Read tabulated data for effective number of degrees of freedom from various
    .tab files in the gstar folder with three columns:
    {Temperature in GeV, g*S, g*}.

    Inputs
    ------

    choices: "standard": Gondolo-Gelmini (LambdaQCD = 150 MeV) (default)
             "HP_A": Hindmarsh-Philipsen equation of state A
             "HP_B": Hindmarsh-Philipsen equation of state B
             "HP_B2": Hindmarsh-Philipsen equation of state B2
             "HP_B3": Hindmarsh-Philipsen equation of state B3
             "HP_C": Hindmarsh-Philipsen equation of state C
             (These are taken directly from MicrOMEGAs package)
    
    gstarpath: By default, set to the path to the gstar folder provided with
               this package
    )pbdoc", py::arg("choice")="standard", py::arg("gstarpath")=GSTARPATH);
    
    //gstar(T)
    mod.def("gstar", &gstar, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV
    anom_mass: mass of additional anomalons 10TeV by default

    Returns
    -------

    The effective number of degrees of freedom for energy density g*

    (By default uses the standard Gondolo-Gelmini g*(T). To use other choices
    for g*: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"),py::arg("anom_mass")=10000.0);
    
    //gstarS(T)
    mod.def("gstarS", &gstarS, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV
    anom_mass: mass of additional anomalons 10TeV by default

    Returns
    -------

    The effective number of degrees of freedom for entropy density g*S

    (By default uses the standard Gondolo-Gelmini g*S(T). To use other choices
    for g*S: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"),py::arg("anom_mass")=10000.0);
    
    /*
    //RhoVisible(T)
    mod.def("RhoVisible", &RhoVisible, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Energy density in the visible sector
    )pbdoc", py::arg("T"));
    
    //EntropyVisible(T)
    mod.def("EntropyVisible", &EntropyVisible, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Comoving entropy in the visible sector
    )pbdoc", py::arg("T"));

    //Hubble(T)
    mod.def("Hubble", &Hubble, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Hubble rate
    )pbdoc", py::arg("T")); */

    //Yield_FreezeIn(mchi, Ve, Ae, Vu, Au, Vd, Ad, Vc, Ac, ma, anom_mass, LambdaQCD, Trh)
    mod.def("Yield_FreezeIn", &Yield_FreezeIn, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    Ve: Vector coupling of dark photon to leptons
    Ae: Axial coupling of dark photon to leptons
    Vu: Vector coupling of dark photon to up-type quarks
    Au: Axial coupling of dark photon to up-type quarks
    Vd: Vector coupling of dark photon to down-type quarks
    Ad: Axial coupling of dark photon to down-type quarks
    Vc: Vector coupling of dark photon to dark matter fermion
    Ac: Axial coupling of dark photon to dark matter fermion
    ma: Mass of dark photon
    anom_mass: Mass of anomalons in GeV. Set to 0 = no anomalons by default
    LambdaQCD: QCD confinement scale in GeV. Set to 0.15 GeV by default
    Trh: "instantaneous reheating temperature". Set to Infinity by default
    Vebewsb: Vector coupling of dark photon to leptons before EWSB
    Aebewsb: Axial coupling of dark photon to leptons before EWSB
    Vubewsb: Vector coupling of dark photon to up-type quarks before EWSB
    Aubewsb: Axial coupling of dark photon to up-type quarks before EWSB
    Vdbewsb: Vector coupling of dark photon to down-type quarks before EWSB
    Adbewsb: Axial coupling of dark photon to down-type quarks before EWSB
    Vcbewsb: Vector coupling of dark photon to dark matter fermion before EWSB
    Acbewsb: Axial coupling of dark photon to dark matter fermion before EWSB
    )pbdoc",py::arg("mchi"), py::arg("Ve"), py::arg("Ae"), py::arg("Vu"), py::arg("Au"), py::arg("Vd"), py::arg("Ad"), py::arg("Vc"), py::arg("Ac"), py::arg("ma"), py::arg("anom_mass")=0.0, py::arg("LambdaQCD")=0.15, py::arg("Trh")=0.0, py::arg("Vebewsb")=0.0, py::arg("Aebewsb")=0.0,py::arg("Vubewsb")=0.0,py::arg("Aubewsb")=0.0,py::arg("Vdbewsb")=0.0,py::arg("Adbewsb")=0.0,py::arg("Vcbewsb")=0.0,py::arg("Acbewsb")=0.0);

    //SigmaDDe(mchi, kappa)
    mod.def("SigmaDDe", &SigmaDDe, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    Ve: Vector coupling of dark photon to electrons
    Ae: Axial coupling of dark photon to electrions
    Vc: Vector coupling of dark photon to dark matter fermion
    Ac: Axial coupling of dark photon to dark matter fermion
    ma: mass of dark photon

    Returns
    -------

    Direct detection cross section (in cm^2) through the light dark photon
    mediator
    )pbdoc", py::arg("mchi"), py::arg("Ve"), py::arg("Ae"), py::arg("Vc"), py::arg("Ac"), py::arg("ma"));

};
