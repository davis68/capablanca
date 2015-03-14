# Introduction #

capablanca is an MPI-based implementation of a first-nearest-neighbor dissolution and deposition model written in C++. An electrochemically-based model of crystal dissolution has been developed in which (electro-)chemical reactions are formulated as Monte Carlo rules. The electrochemical model used assumes bonding in the solid to be a function of first nearest neighbors only, although more general reactions are also supported.

capablanca is an implementation of a first-nearest-neighbor dissolution model used in my master's thesis (available at http://www.mendeley.com/profiles/neal-davis1/) to simulate crystal dissolution. Ultimately, that work was towards developing a first-principles-based model of used nuclear fuel dissolution, which would lead to improved understanding of chemical reaction mechanisms relevant to separations processes. This code represents an implementation of a piece of such a model, although currently experimental data to provide validation are lacking. The method, however, is general for metallically- or ionically-bonded systems.

The model developed in my thesis and implemented in capablanca successfully reproduces the mechanics of crystal dissolution in a variety of crystallographic orientations and kinetic scenarios.

# Publications #
### Conference Proceedings (1) ###

> Neal E Davis, Rizwan Uddin (2011) Capablanca: A parallel code for simulating electrochemical surface dissolution. In _International Conference on Mathematics and Computational Methods Applied to Nuclear Science and Engineering (M&C 2011)_.

### Thesis (1) ###

> Neal E Davis (2011) Modelling and simulation of the dissolution of a physical system with application to head-end operations in aqueous reprocessing. University of Illinois at Urbana--Champaign.