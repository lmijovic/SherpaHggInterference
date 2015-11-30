/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Architecture identified as Darwin MacOS */
/* #undef ARCH_DARWIN */

/* Architecture identified as Linux */
#define ARCH_LINUX "1"

/* Architecture identified as Unix */
/* #undef ARCH_UNIX */

/* BlackHat directory */
/* #undef BLACKHAT_PATH */

/* binreloc activation */
/* #undef ENABLE_BINRELOC */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `mkdir' function. */
#define HAVE_MKDIR 1

/* If available, contains the Python version number currently in use. */
/* #undef HAVE_PYTHON */

/* Have the SQLITE3 library */
#define HAVE_SQLITE3 /**/

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* HEPEVT common block size */
#define HEPEVT_CB_SIZE 10000

/* ld path name set to LD_LIBRARY_PATH */
#define LD_PATH_NAME "LD_LIBRARY_PATH"

/* LHAPDF directory */
/* #undef LHAPDF_PATH */

/* library suffix set to .so */
#define LIB_SUFFIX ".so"

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* Openloops installation prefix */
/* #undef OPENLOOPS_PREFIX */

/* Name of package */
#define PACKAGE "SHERPA-MC"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "sherpa@projects.hepforge.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "SHERPA MC"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "SHERPA MC 2.2.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "SHERPA-MC"

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.2.0"

/* Sherpa include directory */
#define SHERPA_INCLUDE_PATH "/afs/cern.ch/work/s/soffi/CMSSW720-HggInt/src/hintmcnlo/include/SHERPA-MC"

/* Sherpa library directory */
#define SHERPA_LIBRARY_PATH "/afs/cern.ch/work/s/soffi/CMSSW720-HggInt/src/hintmcnlo/lib/SHERPA-MC"

/* Sherpa version name */
#define SHERPA_NAME "Cho Oyu"

/* Sherpa installation prefix */
#define SHERPA_PREFIX "/afs/cern.ch/work/s/soffi/CMSSW720-HggInt/src/hintmcnlo"

/* Sherpa data directory */
#define SHERPA_SHARE_PATH "/afs/cern.ch/work/s/soffi/CMSSW720-HggInt/src/hintmcnlo/share/SHERPA-MC"

/* Sherpa subversion */
#define SHERPA_SUBVERSION "2.0"

/* Sherpa version */
#define SHERPA_VERSION "2"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Using BLACKHAT */
/* #undef USING__BLACKHAT */

/* cernlib found */
/* #undef USING__CERNLIB */

/* Using colour */
#define USING__COLOUR "1"

/* using delphes */
/* #undef USING__DELPHES */

/* Using FASTJET */
/* #undef USING__FASTJET */

/* Using FASTJET 3 */
/* #undef USING__FASTJET__3 */

/* using gzip */
/* #undef USING__GZIP */

/* Using HEPMC2 */
#define USING__HEPMC2 "1"

/* HepMCDefs.h available */
#define USING__HEPMC2__DEFS "1"

/* HepMC::IO_GenEvent available */
#define USING__HEPMC2__IOGENEVENT "1"

/* HepMC::Units available */
#define USING__HEPMC2__UNITS "1"

/* hztool found */
/* #undef USING__HZTOOL */

/* using LHAPDF */
/* #undef USING__LHAPDF */

/* using LHAPDF6 */
/* #undef USING__LHAPDF6 */

/* Using MCFM */
/* #undef USING__MCFM */

/* using MPI */
/* #undef USING__MPI */

/* Pythia interface enabled */
/* #undef USING__PYTHIA */

/* using Rivet */
/* #undef USING__RIVET */

/* setSumOfWeights function available in Rivet */
/* #undef USING__RIVET__SETSOW */

/* Rivet uses YODA as its histogramming backend */
/* #undef USING__RIVET__YODA */

/* using ROOT */
/* #undef USING__ROOT */

/* Version number of package */
#define VERSION "2.2.0"

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
