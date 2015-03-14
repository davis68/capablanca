capablanca is an MPI-based implementation of a first-nearest-neighbor dissolution and deposition model written in C++. An electrochemically-based model of crystal dissolution has been developed in which (electro-)chemical reactions are formulated as Monte Carlo rules. The electrochemical model used assumes bonding in the solid to be a function of first nearest neighbors only, although more general reactions are also supported.

capablanca is an implementation of a first-nearest-neighbor dissolution model used in my master's thesis (available at http://www.mendeley.com/profiles/neal-davis1/ ) to simulate crystal dissolution. Ultimately, that work was towards developing a first-principles-based model of used nuclear fuel dissolution, which would lead to improved understanding of chemical reaction mechanisms relevant to separations processes. This code represents an implementation of a piece of such a model, although currently experimental data to provide validation are lacking. The method, however, is general for metallically- or ionically-bonded systems.

The model developed in my thesis and implemented in capablanca successfully reproduces the mechanics of crystal dissolution in a variety of crystallographic orientations and kinetic scenarios.

Details for usage are found at http://code.google.com/p/capablanca/wiki/CommandLineUsageOfCapablanca and http://code.google.com/p/capablanca/wiki/CommandLineUsageOfTools.

For questions, suggestions, or support, email me at davis 68 (no space) at illinois dot edu.  If you find `capablanca` to be useful, please cite
```
Neal E Davis, Rizwan-uddin (2011) Capablanca: A parallel code for simulating electrochemical surface dissolution. In International Conference on Mathematics and Computational Methods Applied to Nuclear Science and Engineering (M&C 2011). 
```


The following pages reflect earlier versions of this code: https://launchpad.net/capablanca and http://www.xp-dev.com/sc/browse/70712/