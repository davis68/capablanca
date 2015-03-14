# Introduction #
`capablanca` is an MPI-based implementation of the dissolution and deposition model written in C++.  It models a system of particles which can evolve according to an arbitrary set of rules which fit the parameters of the dissolution and deposition equations according to the Monte Carlo acceptance criterion.

The program uses MPI and was tested in the MPICH-2 runtime environment with the Fedora 11 (GNU/Linux) operating system.  Scalability of the problem size and number of processors has not been studied at this time.

`capablanca` uses the Boost C++ uniform random number generator `uniform_01` with the Mersenne twister `mt19937` as the random number generator parameter (Maurer (2006)).  This Mersenne twister function provides random distribution in 623 dimensions with a cycle length of 2<sup>19937</sup>-1.  The seed is given by the system clock.  It is written in the output file with the simulation rate information to facilitate reproducibility of simulations.

# Running `capablanca` #
Usage is as `capablanca [OPTION] ...` where the options are shown in the following table.  If no option is specified for the rule set file or configuration file, the program loads the default files `rules.rs` and `standard.conf`.  Any options set in those files are superseded by the command-line options specified.  Several of these options are intended for debugging usage only.
| Option | Description |
|:-------|:------------|
| `-c CONFIG.CFG` | specify configuration file |
| `-C` | make problem cyclical in _x_- and _y_-directions |
| `-d` | suppress deposition |
| `-h` | display a help message |
| `-E `_P_ | specify potential _P_ |
| `-i `_T_ | specify _T_ steps between output files |
| `-N` | output XYZ file with neighbor count, not state |
| `-n` | force recalculation of nearest neighbors |
| `-p DATAFILE.XYZ` | specify particle position file |
| `-R `_s_ | specify random seed _s_ |
| `-r RULE_SET.RS` | specify rule set file |
| `-S` | force full surface calculation |
| `-s` | output surface at each time step as well |
| `-T `_t_ | specify temperature _t_ |
| `-t `_a_ | specify _a_ time steps to run |
| `-V` | use verbose output |
| `-v` | display version information |
| `-z `_Z_ | specify coordination number _Z_ |

The input file is in the `xyz` file format, which specifies that the number of particles in the file be on the first line with each tab-delimited succeeding line consisting of an identifier and the _x_-, _y_-, and _z_-coordinates of that particle.

The configuration file consists of four lines, respectively specifying the reduced temperature, the inner potential difference, the output interval, and the number of desired time steps.  The default configuration file name is `standard.conf`.

The rule set file defines the transitions possible to each particle in the system.  The first three white-space-delimited numbers specify the number of rules contained in the file, the total number of states, and the state dividing dissolved states from deposited states.  These are followed by each rule consisting of an explanatory comment and eight numbers:  the starting state for a reaction towards dissolution, the dissolution prefactor, the deposition prefactor, the transfer coefficient, the energy required to break each bond, the activation energy for a reaction, the number of electrons transferred in the reaction, and the state if the reaction is accepted.  (Refer to the master's thesis on the documentation page for the use of these variables.)  The default rule set file name is `rules.rs`.  An example file is given in below.  A limitation of `capablanca` is that it currently only permits one reaction per state.

A rule file is tied to a specific input file by the states, which are arbitrarily assigned numbers that describe the particle states in the rule-set file and internally to the program.  The identification of any oxidation state or species with any number is done by the user.  However, the dissolved states must all be assigned numbers greater than the dissolution number, which separates the dissolved states from the deposited states.
```
        2	3	2
        #0-->1
        0	1.0	1.0	0.4	1.0	0.0	0	1
        #1-->2
        1	1.0	1.0	0.5	0.0	1.0	0	2
```

There are three exponential terms in the dissolution equation used (see master's thesis for more details; sorry, Google Code wiki language doesn't support equations very well).  The first exponential term applies to nearest-neighbor-dependent dissolution.  The second describes an electrochemical reaction in which the activation energy is known and is not a function of the number of nearest neighbors.  The third pertains to a potential difference applied across the cell.  Phenomena represented by the first and second exponential terms are mutually exclusive, and so the energy factor in one of the exponential terms should always be set to zero.  Similarly, if a reaction cannot proceed in the reverse direction, the prefactor for that direction should be set to zero.

To follow surface evolution, the exterior surface of the crystal, consisting of all contiguous non-substrate atoms with less than a full complement of nearest neighbors (and thus a non-zero probability of dissolution), is detected.  An example of dissolution in a two-dimensional crystal (_Z_ = 4) is given.  This contiguity permits vacancy defects to be treated properly.  Only atoms on the surface are capable of reaction or dissolution, and so the main program iterates over the surface, evolving it as additional atoms dissolve, exposing the substrate to the dissolving medium.  Sites with a dissolution probability of less than 1 in 10<sup>15</sup> parts are omitted from the calculation in order to reduce the number of invocations of the random number generator function.

To balance the calculation load between processors, the internal representation of the crystal is sliced, or divided according to the _y_-coordinate between the separate processors.  In the case of dissolution from a single crystal face (bounded by inert material on five sides of the crystal and thus permitting the dissolution characteristics of a single face to be isolated), this keeps the processor loads approximately equal.  (Load imbalance may develop if a large crystal, being dissolved from all directions simultaneously, dissolves completely in the range of _y_-coordinates assigned to a specific processor.)  In this case, `capablanca` scales in a straightforward manner, since calculations continue with a similar approximate distribution between processors.

Some further remarks on the operation of `capablanca`:  The nearest-neighbor loading routines occasionally fail for an unknown reason; in this case, execution should force recalculation with the `-n` flag on every execution.  When loading an odd number of particles, a handful are occasionally dropped incorrectly when the particles are read.  If the cluster crashes, the calculated data (written in `xyz_temp` files) can be stitched together by using the associated program `collate`.