
# Brief Description

`BsaLib`, a Modern Fortran Library for the Bispectral Stochastic Analysis 
of structures under non-Gaussian stationary random actions.

> NOTE: currently, only wind action is included in the library, but other phenomena (waves for instance) can be easily integrated. 
> See [further developments](#furhter_developments).


# License

`BsaLib` is release under the **GNU Lesser General Public License v3.0**.
Visit the [GPL official website](https://www.gnu.org/licenses/gpl-3.0.html) for more information.



# Build system

`BsaLib` uses [CMake](https://cmake.org/) as a build system generator.


# Documentation

`BsaLib` support documentation generation from sources. It requires the following tools:

- [Python](https://www.python.org/)
- [FORD](https://forddocs.readthedocs.io/en/latest/user_guide/getting_started.html)

To build it:

```
cmake --build <build_dir> --target docs
```


# Code structure

There are two parts in this repository:

1. [`BsaLib`](#bsalib_core), the core library;
2. [`bsa`](#bsa_executable), the built-in executable program.



## <a id="bsalib_core">BsaLib: core library</a>

`BsaLib` is the main core of this repository, found under `./BsaLib/`
It consists of the main library and its API to which anyone could link to and interact with.
To use `BsaLib`, simply import the main `BsaLib` API module:

```Fortran
program test
    use, non_intrinsic :: BsaLib
    implicit none (type, external)
    ! your declarations here

    ! your logic here

    ! initialise BsaLib
    call bsa_Init()

    ! set BsaLib internal state through its API procedures

    ! once done, run BsaLib
    call bsa_Run( ... args ... )

    ! finally, release BsaLib memory
    call bsa_Finalise()
end program
```

Being designed as a *plug-in* library, it needs the hosting program/library 
to provide some data needed by `BsaLib` in order to function properly.
For more details on the public API, visit the [main documentation page](https://html-preview.github.io/?url=https://github.com/miEsMar/BsaLib/blob/main/doc/index.html).



## <a id="bsa_executable">bsa: executable program</a>

As a side project package, a single-source executable file is provided under `./bsa/`. 
It emulates what one would normally do when using `BsaLib` as a *plug-in* for its own 
library/program. 
On the other hand, this program is thought and provided for all those interested in using 
`BsaLib` but not having any hosting library/program. 
Nonetheless, even if this program is provided, the user would yet need to provide the data 
for `BsaLib` to function properly. If any of this data is not provided, the `BsaLib` runtime check 
would detect it and abort correct logic flow.

For that, the provided executable program relies on reading two input files:

1. `BsaLib` related settings (formatted file, named `bsa.bsadata`). 
For details, read the [dedicated section](#bsadata_input_file).
2. External data file (named `bsa.extdata`), in binary format, containing 8-byte floating point records 
(`real64` of the `iso_fortran_env` compiler intrinsic module). 



# Cross-platform support

`BsaLib` strives to be as cross-platform as possible, so that any user can use it regardless of 
the tooling availability.
Currently, the code has been compiled and tested under three different OS-compiler configurations:

1. `Windows OS Build 10.0.19045` - Intel Fortran Compilers (`ifort 2021.7.1-20221019_000000`, `ifx 2022.2.1-20221101`)
2. `MacOS` - GFortran 13.2
3. `Linux Centos Fedora 8.7` - Intel Fortran Compilers (`ifort 2021.10.0-20230609`, `ifx 2023.2.0-20230721`)


For the Proper Orthogonal Decomposition (POD) problem, two approaches have been tested, for configurations
1 and 3:

- Linkage to proprietary `Intel MKL` (Math Kernel Libraries)
- Linkage to [LAPACK](https://www.netlib.org/lapack/) native implementation 
(**NOTE**: from direct source build in configuration 1)


# <a id="furhter_developments">What's missing? Further developments</a>

Mathematical:

- [ ] Integrate models for other non-Gaussian actions (waves, for instance)
- [ ] Extend to non-diagonal modal Frequency Response Function (FRF) models
- [ ] Extend to frequency-dependent modal matrices (e.g. aeroelastic phenomena in Wind Engineering)
- [ ] Add library core to generate spectra (PSDs and Bispectra) from time series

Numerical:

- [ ] Adapt Classic approach to dump BFM info as in Mesher (easy)
- [ ] Improve internal policy mechanism and integration (WIP)
- [ ] Provide a general API for user defined models integration
- [ ] Complete full support for spatial (in-plane) symmetries of a real-valued bispectrum
- [ ] Compute nodal correlation internally (don't require it as user data)
- [ ] Add support for Mesher zones' interest modes
- [ ] Add a local caching system
- [ ] Add `MPI` support (for running `BsaLib` on multi-node clusters)
- [ ] Improve and extend GPU offloading capabilities
- [ ] Provide automation service to convert user-level model function into GPU-kernel code
- [ ] Enhance the mechanism that lets the user provide its own desired exporting function, 
so that `BsaLib` is not tight to any specific exporting format.
- [ ] Integrate a built-in bispectrum post-processing Visualiser (using [Vulkan](https://www.vulkan.org/), for optimal performances)


# Known Issues

There is one main known issue in the current version. 

> Using `OpenMP` parallelisation in the second Post-meshing phase, 
> execution time is higher compared to serialised version. 
> This is due to the necessary synchronisation between threads when accessing 
> shared file I/O when reading each zone's dumped data, causing the `critical` 
> section to be the main bottleneck in this part.
> For a proper use of `OpenMP` parallelisation in Post-meshing phase, a fundamental
> rethinking of the algorithmic structure needs to be done.
> A first, temporary, possible solution, could be the usage of thread-level private
> I/O with dedicated units, avoiding the need to synchronise when reading back 
> information in post-meshing phase. However, this might soon become inelegant solution 
> when the number of threads would start increasing considerably.
> For this, as a temporary solution, a conditional compilation flag 
> (`BSA_USE_POST_MESH_OMP`) can be used to control effective use of 
> `OpenMP` parallelisation in Post-meshing phase, **disabled by default**.




# Related Scientific Publications
-------------------------------

1. [Non-Gaussian buffeting analysis of large structures by means of a Proper Orthogonal Decomposition](https://www.sciencedirect.com/science/article/abs/pii/S0167610523002799?via%3Dihub)
2. [A multiple timescale approach of bispectral correlation](https://www.sciencedirect.com/science/article/abs/pii/S0167610522003786?via%3Dihub)
3. [On the background and biresonant components of the random response of single degree-of-freedom systems under non-Gaussian random loading](https://www.sciencedirect.com/science/article/abs/pii/S0141029611001507)





# Appendix


## <a id="bsadata_input_file">`bsa.bsadata` input file structure</a>

This file is required by the `bsa` executable to gather information 
related to `BsaLib` specific settings. 

There are 5 sections (**WARNING**: to be provided in the given order!):

1. [GENERAL](#general)
2. [CLASSIC](#classic)
3. [MESHER](#mesher)
4. [DIRECTIONS](#directions)
5. [TURBULENCE](#turbulence)


### General

Currently supported general settings (in order):

- **Analysis type ID** (`int`):

    Valid values:
    - `1` (classic type);
    - `2` (mesher type);
    - `3` (both).

    Default: `1`.


- **Version ID** (`int`) [UNUSED]

- **PSD (Power Spectral Density) scaling convention ID** (`int`)

    Valid values:
    - `1` (circular frequencies convention, $\sigma^2 = \int_{-\infty}^{\infty} S(\omega) d\omega$);
    - `2` (frequency convention, $\sigma^2 = \int_{0}^{\infty} S(f) df$).

    Default: `1`.


- **PSD computation switch** (`int`)

    Valid values:
    - `0` (OFF, no PSDs computed);
    - `1` (ON).

    Default: `1`.


- **Bispectra computation switch** (`int`)

    Valid values:
    - `0` (OFF, no bispectra computed);
    - `1` (ON).

    Default: `1`.


- **Only diagonal elements computation switch** (`int`)

    Controls whether all tensors elements are computed, 
    or computation limited to tensors elements whose indices are equal in all dimensions (i.e. $a_{iii}$).\
    Valid values:
    - `0` (OFF, all elements computed);
    - `1` (ON, only diagonal elements computed).

    Default: `0`.

    ---
    **NOTE**:\
    Indeed, `1` results in a much faster computation, specially for Bispectra computation. 
    It is however less precise due to the neglecting of outer-diagonal elements,
    which might play a crucial role in the final statistical moments estimates.\
    This is explained in detail in the [companion paper](https://doi.org/10.1016/j.jweia.2023.105576).

    ---


- **Bispectra spatial symmetry value** (`int`) [EXPERIMENTAL]

    Controls how much information is implicitly assumed 
    leveraging spatial (in-plane) symmetries of real part of bispectra.\
    Valid values:
    - `0` (FULL, all information computed);
    - `2` (HALF);
    - `4` (FOURTH).

    Default: `0`.

    ---
    **NOTE**:\
    at one point, default should become `2`, since this would be the best solution.

    ---


- **Tensor symmetry switch** (`int`) [EXPERIMENTAL]

    Allows to control whether symmetrical elements of a spectra tensor 
    (PSD or Bispectra) are computed (i.e. $a_{ijk}$ computed only once for all possible permutations of 
    indices $(i,j,k)$.) or not.\
    Valid values:
    - `0` (OFF, all tensor elements are computed);
    - `1` (ON, only unique symmetrical tensor elements are computed).

    Default: `0`.

    ---
    **NOTE**:\
    at one points, default should become `1`.

    ---

    ---
    **NOTE2**:\
    this flag has only sense if diagonal elements only flag is `OFF`.

    ---


- **Test flag** (`int`)

    If ON, disables some internal checks, specially at the level of discretisation refinement 
    compared to suggested values.\
    Valid values:
    - `0` (OFF, check enabled);
    - `1` (ON, checks disabled).

    Default: `0`.



### Classic

The classic group controls settings regarding the "classic" approach, which includes 
"conventional" spectral and bispectral analysis implementations.

There are only two main values to be set:

- **Number of discretisation frequencies** $\tt N_{freqs}$ (`int`)

    To be considered from `0` to $f_{max} \text{ [Hz]}$ (or $\omega_{max}\ [\frac{rad}{s}]$ if 
    circular frequency convetion used).\
    Valid values: `> 0`.

    Default: NONE, a value must be provided!


- **Delta frequency** (refinement) $\Delta f \text{ [Hz]}$ (`double`)

    After n. of frequencies is give, $\Delta f$ specifies the (regular) spacing between each frequency.\
    Valid values: `> 0.`.

    Default: NONE, a value must be provided!



### Mesher

- **POD (Proper Orthogonal Decomposition) Switch** (`int`)

    Controls whether POD is used for decomposing base wind flow field.
    See [reference](https://doi.org/10.1016/j.jweia.2023.105576) for more details.\
    Valid values:
    - `0` (OFF, no POD);
    - `1` (ON).

    Default: `1`.

- **Background zone base refinement** (`int`)

    Specifies the base number of meshing points (per side) to disretise the rectangular zone 
    covering the background peak at the origin $(0, 0)$.\
    Valid values: `> 0`.

    Default: NONE, must be provided!

- **Background zone area extension** (`real`)

    Specifies the integer multiplier of the base background zone extension.\
    Valid values: `> 0.`.

    Default: `1.` (NO extension).


- **Peak zone area extension** (`real`)

    Specifies the integer multiplier of a peak zone extension, where a peak zone is a meshing zone 
    covering a part of the 2D space where a (bi-)resonant peak is located.\
    Valid values: `> 0.`.

    Default: `1.` (NO extension).


- **Max area extension** (`real`)

    Specifies the integer multiplier of the max area extension.\
    Internaly, a maximum area (of interest) limit is determined.\
    This factor allows to extend this limit up to the desired value.\
    Valid values: `> 0.`.

    Default: `1,` (NO extension).


- **Full coverage switch** (`int`) [DEPRECATED]


- **Dump modal info** switch (`int`)

    Controls whether to include modal info in dump file (used in post-processing tasks), or not.\
    Valid values:
    - `0` (NO, do not include);
    - `1` include.

    Default: `1`.


- **Wind directions** (`int - int[]`)

    See [example](#bsadata_input_file_example).


- Turbulence (`int - int[]`)

    Specifies the (3) spatial turbulent component to be accounted in the determination of the 
    stochastic wind loading.\
    First entry is an `int` specifying the effective number of turbulent components.\
    Following, a list (one entry per line) of component indices,
    where each index should be in the set ${1, 2, 3}$ ($(x, y, z)$ components in an Euclidean 3D space).

    ---
    **NOTE**:\
    in the [example below](#bsadata_input_file_example),
    we declare to consider 2 turbulent components, component `1` (x) and `3` (z).

    ---




## <a id="bsadata_input_file_example">Example of working `bsa.bsadata` file</a>

```text
-GENERAL:
   2
   0
   0
   1
   1
   0
   0
   0
   1
-CLASSIC:
   500
   0.010
-MESHER:
   1
   50
   2
   2
   2
   1
   1
-DIRECTIONS:
   1
   1
-TURBULENCE:
   2
   1
   3
```


