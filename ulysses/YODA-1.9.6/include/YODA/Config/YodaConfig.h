/* include/YODA/Config/YodaConfig.h.  Generated from YodaConfig.h.in by configure.  */
#ifndef YODA_YODACONFIG_H
#define YODA_YODACONFIG_H


/* Define to the address where bug reports for this package should be sent. */
#define YODA_BUGREPORT "yoda@projects.hepforge.org"

/* Define to the full name of this package. */
#define YODA_NAME "YODA"

/* Define to the full name and version of this package. */
#define YODA_STRING "YODA 1.9.6"

/* Define to the one symbol short name of this package. */
#define YODA_TARNAME "YODA"

/* Define to the version of this package. */
#define YODA_VERSION "1.9.6"


#include <string>
namespace YODA {
  /// Namespaced version string function
  inline std::string version() {
    return YODA_VERSION;
  }
}


#endif
