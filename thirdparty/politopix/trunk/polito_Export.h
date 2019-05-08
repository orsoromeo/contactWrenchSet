
#ifndef polito_EXPORT_H
#define polito_EXPORT_H

#ifdef polito_BUILT_AS_STATIC
#  define polito_EXPORT
#  define POLITO_NO_EXPORT
#else
#  ifndef polito_EXPORT
#    ifdef polito_EXPORTS
        /* We are building this library */
#      define polito_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define polito_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef POLITO_NO_EXPORT
#    define POLITO_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef POLITO_DEPRECATED
#  define POLITO_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef POLITO_DEPRECATED_EXPORT
#  define POLITO_DEPRECATED_EXPORT polito_EXPORT POLITO_DEPRECATED
#endif

#ifndef POLITO_DEPRECATED_NO_EXPORT
#  define POLITO_DEPRECATED_NO_EXPORT POLITO_NO_EXPORT POLITO_DEPRECATED
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define POLITO_NO_DEPRECATED
#endif

#endif
