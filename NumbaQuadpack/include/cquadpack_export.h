#ifndef CQUADPACK_EXPORT_H
#define CQUADPACK_EXPORT_H

/* Static replacement for the cmake-generated export header.
 * All C sources are compiled into a single Python extension,
 * so no import/export decoration is required beyond basic
 * symbol visibility. */

#if defined(_WIN32) || defined(__CYGWIN__)
#  define CQUADPACK_EXPORT __declspec(dllexport)
#  define CQUADPACK_NO_EXPORT
#  define CQUADPACK_DEPRECATED __declspec(deprecated)
#elif defined(__GNUC__) || defined(__clang__)
#  define CQUADPACK_EXPORT __attribute__((visibility("default")))
#  define CQUADPACK_NO_EXPORT __attribute__((visibility("hidden")))
#  define CQUADPACK_DEPRECATED __attribute__((__deprecated__))
#else
#  define CQUADPACK_EXPORT
#  define CQUADPACK_NO_EXPORT
#  define CQUADPACK_DEPRECATED
#endif

#define CQUADPACK_DEPRECATED_EXPORT    CQUADPACK_EXPORT    CQUADPACK_DEPRECATED
#define CQUADPACK_DEPRECATED_NO_EXPORT CQUADPACK_NO_EXPORT CQUADPACK_DEPRECATED

#endif /* CQUADPACK_EXPORT_H */
