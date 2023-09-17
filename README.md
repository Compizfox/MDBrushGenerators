# MDBrushGenerators
![GitHub license](https://img.shields.io/github/license/Compizfox/MDBrushGenerators)
![Python version](https://img.shields.io/badge/Python-%3E3.6-orange)
![GitHub last commit](https://img.shields.io/github/last-commit/Compizfox/MDBrushGenerators)
[![DOI](https://zenodo.org/badge/279256604.svg)](https://zenodo.org/badge/latestdoi/279256604)

MDBrushGenerators contains scripts for generating LAMMPS initial data files comprising polymer brush systems (polymers
end-grafted to a surface).

![Snapshot](https://user-images.githubusercontent.com/7603719/87414213-b8f02400-c5cb-11ea-87a8-f4b63c076801.png)

The core of the repository is `BrushGenerator`: an abstract, extensible, model-agnostic class for building a polymer
brush system and writing the initial data file in the format LAMMPS accepts. It makes use of a Poisson-disk point set
generator (implemented in a separate class `PoissonDiskGenerator`) to generate coordinates of the grafting layer.

In the current version, the only finalised polymer model implementation is the Kremer-Grest model, found in
`KremerGrestBrushGenerator.py`. `CrosslinkGenerator` defines a variant of a Kremer-Grest model with single-particle side
chains and/or backbone heteroparticles.

## Installation
Simply clone the repository using Git:

```console
foo@bar:~$ git clone https://github.com/Compizfox/MDBrushGenerators.git
```

or download a release and extract it. The Git approach has the advantage that updating is as easy as `git pull`.

### Dependencies
MDBrushGenerators requires at least Python 3.6, Poetry, and a Rust compiler to build the `PoissonDiskGenerator`.

Install the dependencies using Poetry:

```console
$ poetry install
```

Then, cd into the `PoissonDiskGenerator` directory and build the Rust code using Maturin:

```console
$ cd PoissonDiskGenerator
$ maturin develop -r
```

The `PoissonDiskGenerator` module should then be installed in the Poetry virtualenv.

## Usage
### Kremer-Grest
A ready-to-run example script that uses `KremerGrestBrushGenerator` can be found in `kremergrest_brush_generator.py`. It
expects you to set the desired chain length, grafting density, and box size in the script and generates a polymer brush
system using those parameters. It checks for over-density (the condition when the requested grafting density exceeds the
maximal Poisson-disk point set density) and plots the coordinates of the grafting points using Matplotlib. Finally, it
writes the gzip-compressed initial data file.

## Developer documentation
### `BrushGenerator`
Model-specific classes subclass (inherit from) `BrushGenerator`, with the most important method being overridden being
`_build_bead()`. This method is called iteratively in `build()` and should be implemented according to the polymer
model. Its function is to 'build a bead' of the polymer chain and append the coordinates, bonds, angles, and dihedrals
(if applicable) to the object's non-final lists.

When `build()` is run, it calls `_build_bead()` for every bead for every chain and converts the non-final lists to
Pandas DataFrames afterwards. Then, in `write()`, the DataFrames are read and converted to a LAMMPS data file, together
with force field coefficients.

### `PoissonDiskGenerator`
The grafting coordinates are sampled from a Poisson-disk distribution. A Poisson-disk point set is a set of `n` points
randomly sampled from a uniform distribution, with the constraint that no pair of points is closer than a given distance
`r`.

The class `PoissonDiskGenerator` class generates Poisson-disk point sets using cell list accelerated dart throwing.

This algorithm functions in roughly the following way:

- Divide the domain into square cells of size `r/sqrt(2)` and put them in a (flat) list of active cells.
- While the number of points is less than the desired number of points:
	- Choose an active cell and 'throw a dart' in it: sample a point from an uniform distribution.
	- Check for overlap by checking the distance to all points in the 20 neighbouring cells (a 5x5 grid centred around
the current cell, excluding the current cell itself and the 4 corner cells) using the cell lookup list. If there is no
overlap, add the point to the total list and to the cell lookup list, and remove the current cell from the active cells
list.

Two data structures are maintained: the active cells list, which is a list containing tuples of cell indices of the
cells that do not contain a point. It starts out containing all the cells. When a point is added to a cell, the cell is
removed form the active cells list. The cell lookup list is a 2D matrix that starts out empty. When a point is added,
its coordinates get added to the respective cell in the list, in order to accelerate searching for overlaps.

This algorithm can end in one of three ways: either the desired number of points is reached, no active cells are left,
or there are active cells left but the point set is maximal and no extra points can be added to the set. In the first
case, the result (consisting of a set of point coordinates) is satisfactory. In the second case, the result is a valid
point set, but the density is lower than requested. In the third case, the algorithm doesn't halt intrinsically and an
iteration limit (of the outer loop) is required to halt the program after some time.

This algorithm appears (empirically) to be `O(n)`, i.e. linear in time with respect to `n` at constant `r` and point
density.

## License
This project is free software licensed under the GPL. See [LICENSE](LICENSE) for details.

