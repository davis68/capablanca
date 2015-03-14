# Introduction #
`ruylopez`, `petroff`, and `polgar` are sequential programs which provide crystal structures for use in `capablanca`.  They are written in C++ using OpenMP.

They have been tested with the Ubuntu (Lucid Lynx and later) and Fedora 11 GNU/Linux operating systems.

# Running `ruylopez` #
`ruylopez` generates a grain boundary between two crystals.  Usage is as `ruylopez T N [t p] i` where the options are shown in the following table.  The extension in each dimension is specificed in the configuration file `grain.conf`.
| Option | Description |
|:-------|:------------|
| `-c CONFIG.CFG` | specify configuration file |
| `T`  | orientation type (0 - fcc; 1 - bcc; 2 - sc) |
| `N`  | proportional to number of unit cells per grain |
| `-h` | display a help message |
| `t`  | angle from _yz_ plane for grain 1 in negative units of pi |
| `p`  | angle from _yz_ plane for grain 2 in positive units of pi |
| `i`  | relative twist orientation of grains 1 and 2 in units of pi |

The output file is in the `xyz` file format, which specifies that the number of particles in the file be on the first line with each tab-delimited succeeding line consisting of an identifier and the _x_-, _y_-, and _z_-coordinates of that particle.

The configuration file consists of three lines indicating the number of units extending in the _x_, _y_, and _z_ directions.

# Running `petroff` #
`petroff` generates a screw dislocation in a crystal (with as-yet-unsatisfactory results).  Usage is as `petroff T N` where the options are shown in the following table.  (The spacing in the screw portion is still undesirable.)
| Option | Description |
|:-------|:------------|
| `T`  | orientation type (0 - fcc; 1 - bcc; 2 - sc) |
| `N`  | number of layers to generate |
| `V`  | vacancy defect concentration in units of 1e6 |

# Running `polgar` #
`polgar` generates a crystal face surrounded by inert material on the edges, thereby exposing only a single face to dissolution and allowing face-specific rates to be studied.  Usage is as `polgar T x y z p t [V]` where the options are shown in the following table.
| Option | Description |
|:-------|:------------|
| `T`  | orientation type (0 - fcc; 1 - bcc; 2 - sc) |
| `x`  | extension in _x_ direction |
| `y`  | extension in _y_ direction |
| `z`  | extension in _z_ direction |
| `t`  | angle from _yz_ plane for grain in units of pi |
| `p`  | angle from _xy_ plane for grain in units of pi |
| `V`  | vacancy defect concentration in units of 1e6 |

# Alternative methods of generating data #
Aten is an atomic layout editor introduced in 2010 ("Aten - An application for the creation, editing, and visualization of coordinates for glasses, liquids, crystals, and molecules", T. G. A. Youngs, _J. Comp. Chem. 31_, 639--648, (2009)).  Using this editor's unit cell replication features, raw data for `capablanca` can be easily generated.  Aten can be downloaded at [Project Aten](http://www.projectaten.net/).  I'll post more information as I continue to work with it.