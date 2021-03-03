# Angora: A free finite-difference time-domain (FDTD) electromagnetic simulation package

*Angora* is a free, open-source software package that computes numerical solutions to electromagnetic radiation and scattering problems. It is based on the [finite-difference time-domain (FDTD) method](http://en.wikipedia.org/wiki/Finite-difference_time-domain), which is one of the most popular approaches for solving Maxwell's equations of electrodynamics. It has been developed by [Ilker R. Capoglu](mailto:capoglu@angorafdtd.org) and [Di Zhang](mailto:dizhang2016@u.northwestern.edu), and is currently maintained by Aya Eid and Patrick (Yuanzhe) Su. <dfn>Angora</dfn> has been written mostly in C++, with a few external C functions.

## Main features of Angora

* Flexible and user-friendly configuration via text-based configuration files.

*   An automatic build and install mechanism for the GNU/Linux operating system.
*   **Full parallelizability** in three dimensions, based on the Message Passing Interface ([MPI](http://www.mcs.anl.gov/research/projects/mpi/)) library.
*   Support for **planar multilayered spaces.**
    *   Total-field/scattered-field (TF/SF) plane-wave source.
    *   Phasor-domain near-field-to-far-field transformation (NFFFT).
*   Complex optical beams (focused beam, Hermite-Gaussian beam)
*   Convolution perfectly-matched layer (**CPML**) absorbing boundaries.
*   Support for **[HDF5](http://www.hdfgroup.org/HDF5/)**, a portable file storage format.
*   Support for **dispersive media** with Drude/Lorentz pole pairs.
*   Support for external inputs of materials and shapes.
*   Field-value recording (2D, 1D, single point) and movie generation.
*   Generation of **random** permittivity, permeability, conductivity profiles.
*   Numerical **optical imaging** simulation.

## Documentation

Documentation is available in both html and pdf formats:

|                                                 |                                                    |
| ----------------------------------------------- | -------------------------------------------------- |
| [PDF](http://www.angorafdtd.org/doc/angora.pdf) | [ HTML](http://www.angorafdtd.org/doc/angora.html) |

If you would like to receive updates on the development and documentation of <dfn>Angora</dfn>, please send an email to [help@angorafdtd.org](mailto:help@angorafdtd.org). You will only be contacted for <dfn>Angora</dfn>-related communication, and your email address will **absolutely never** be shared with anyone.

## Downloads

* * *

<dfn>Angora</dfn> is currently only available for the **GNU/Linux** operating system. If you would like to port <dfn>Angora</dfn> to another operating system, please [contact us](mailto:capoglu@angorafdtd.org). Contributions are always welcome.

The latest version of <dfn>Angora</dfn> can be downloaded here:

|                                            |                                                              | Requirements                                                 |
| ------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| GNU/Linux x86 64-bit binary                | [angora-0.22.5-bin.tar.gz](http://www.angorafdtd.org/angora/angora-0.22.5-bin.tar.gz) | GNU C shared libraries (>2.4)<br/>GNU standard C++ shared libraries (>3.4) |
| GNU/Linux x86 64-bit binary (with OpenMPI) | **[angora-0.22.5-bin-openmpi.tar.gz](http://www.angorafdtd.org/angora/angora-0.22.5-bin-openmpi.tar.gz)** | GNU C shared libraries (>2.8)<br/>GNU standard C++ shared libraries (>3.4)<br/>OpenMPI shared libraries (>1.4) |
| Source code                                | **[angora-0.22.5.tar.gz](http://www.angorafdtd.org/angora/angora-0.22.5.tar.gz)** | Following libraries: [blitz++](http://sourceforge.net/projects/blitz/), [libconfig](http://www.hyperrealm.com/libconfig/), [hdf5](http://www.hdfgroup.org/HDF5/), [boost](http://www.boost.org/) |

A list of latest important updates can be found in the [ChangeLog](http://www.angorafdtd.org/angora/ChangeLog).

## External Scripts

Some useful scripts for reading and displaying <dfn>Angora</dfn> output can be found [here](http://www.angorafdtd.org/scripts.html).

## License

<dfn>Angora</dfn> is licensed under the [GNU Public License (GPL)](http://www.gnu.org/licenses/gpl.html). For inquiries on other licensing options, please contact the [authors](mailto:capoglu@angorafdtd.org).

## Acknowledgments

* * *

The <dfn>Angora</dfn> project owes a great deal to the invaluable support and guidance provided by (in chronological order) Glenn S. Smith (Emeritus Professor, Georgia Institute of Technology), Vadim Backman (Professor, Northwestern University), and Allen Taflove (Professor, Northwestern University). The <dfn>Angora</dfn> project has been funded by the following sources (in chronological order): John Pippin Chair in Electromagnetics within the School of Electrical and Computer Engineering in Georgia Institute of Technology; National Institutes of Health (NIH) grant R01EB003682.

## Support

For questions, suggestions, and comments, contact us at [help@angorafdtd.org](mailto:help@angorafdtd.org). Please send bug reports to [bugs@angorafdtd.org](mailto:bugs@angorafdtd.org)

* * *



##### Copyright (C) 2006-2017 Ilker R. Capoglu, Di Zhang