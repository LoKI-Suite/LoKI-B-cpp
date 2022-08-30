#ifndef H_LOKI_B_EXPORTS_H
#define H_LOKI_B_EXPORTS_H

/** \file Define macro lokib_export for dll/so imports, exports
 *
 *  Macro lokib_export must be placed before the return type of a
 *  free function declaration and before the name of a struct/class
 *  declaration when that function/class is not defined inline.
 *
 *  When building the loki-b library, it tells the compiler that the
 *  function or class must be exported from the dll/so that is being
 *  created. This is achieved by defining LOKIB_BUILDING_LIBLOKIB
 *  when the files that constitute the dll/so are being compiled.
 *
 *  When the library interfaces (header files) are used by other code,
 *  so LOKIB_BUILDING_LIBLOKIB is not defined, the lokib_export macro
 *  is defined to an instruction to the compiler that said function,
 *  class must be imported from the dll/so.
 *
 *  The precise values of lokib_export in these cases also depends on
 *  the compiler and host system.
 */

#if defined(LOKIB_BUILDING_LIBLOKIB)

// when building a shared object:

#if defined(_MSC_VER)
    //  Microsoft
    #define lokib_export __declspec(dllexport)
#elif defined(__GNUC__)
    //  GCC & CLANG
    #define lokib_export __attribute__((visibility("default")))
#else
    //  do nothing and hope for the best?
    #define lokib_export
#endif

#else

// when using code from a shared object:

#if defined(_MSC_VER)
    //  Microsoft
    #define lokib_export __declspec(dllimport)
#elif defined(__GNUC__)
    //  GCC & CLANG
    #define lokib_export
#else
    #define lokib_export
#endif

#endif // defined(LOKIB_BUILDING_LIBLOKIB)


#endif // H_LOKI_B_EXPORTS_H
