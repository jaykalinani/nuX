<!-- <img align="top" src="Docs/figures/nux.png" width="140"> -->

**nuX** is a GPU-accelerated neutrino-transport module for dynamical spacetimes, written in C++ and designed to run within the [Einstein Toolkit](https://einsteintoolkit.org/) using the [CarpetX](https://github.com/EinsteinToolkit/CarpetX) driver. **CarpetX** itself is built atop [AMReX](https://amrex-codes.github.io), a framework for block-structured adaptive mesh refinement (AMR).  
**nuX** provides a two-moment (M1) evolution of radiation variables with analytic closures and multi-species support, interfacing with our GRMHD code [AsterX](https://github.com/EinsteinToolkit/AsterX).

<!-- 
[![GitHub CI](https://github.com/jaykalinani/nuX/workflows/CI/badge.svg)](https://github.com/jaykalinani/nuX/actions)  
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](./LICENSE.md)
-->
---

## Overview

- Two-moment (**M1**) neutrino transport with analytic closures (e.g., Minerbo/Levermore)  
- Multi-species support (typically $ν_e$, $\barν_e$, $ν_x$)  
- Source terms for neutrino–matter coupling and stress–energy feedback ($T^{μν}_{\rm rad}$) coupled with **AsterX** 
- Robust floors, masking, and diagnostics for production-quality BNS/CCSN simulations

---

## Available Modules

- `nuX_M1` — Core M1 evolution: fluxes, closures, source terms, analysis/diagnostics, stress–energy output  
- `nuX_Rates` — Physical interaction rates & opacities via `bns_nurates`
- `nuX_FakeRates` — Lightweight fake/constant rates for testing 
- `nuX_RatesToy` — Toy module for testing 
- `nuX_Seeds` — Initial data for radiation fields as well as MHD variables
- `nuX_Utils` — Tensor utilities, metric helpers, math wrappers 

---

## Getting Started

* Instructions for downloading and building nuX with the Einstein Toolkit are available [here](https://github.com/EinsteinToolkit/CarpetX/wiki/Getting-Started).  
* Simfactory files for various clusters and setup instructions can be found [here](https://github.com/lwJi/ETK-Compile-Guides).    

## Useful Repositories

* [AsterX](https://github.com/EinsteinToolkit/AsterX) - GRMHD code
* [CarpetX](https://github.com/EinsteinToolkit/CarpetX) – Next-generation driver for the Einstein Toolkit  
* [SpacetimeX](https://github.com/EinsteinToolkit/SpacetimeX) – Modules for spacetime evolution 
* [BNSTools](https://github.com/jaykalinani/BNSTools) – Utilities supporting BNS merger simulations  
