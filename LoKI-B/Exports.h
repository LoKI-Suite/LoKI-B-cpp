#ifndef H_LOKI_B_EXPORTS_H
#define H_LOKI_B_EXPORTS_H

// For now we want to enable this only for MSVC builds
// We could also add something like 'if doing a shared build',
// at present we do that unconditionally.

#ifdef _MSC_VER

#define LOKIB_USE_DLL_IMPORTS

#if defined(LOKIB_BUILDING_LIBLOKIB)
#define lokib_export __declspec( dllexport )
#elif defined(LOKIB_USE_DLL_IMPORTS)
#define lokib_export __declspec( dllimport )
#else
#define lokib_export
#endif

#else

#define lokib_export

#endif // _MSC_VER

#endif // H_LOKI_B_EXPORTS_H
