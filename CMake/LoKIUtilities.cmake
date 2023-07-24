# This part contains content from PLASIMO's CMake/PlasimoUtilities.cmake,
# 'plasimo' has been renamed into 'loki' in function and macro names.

#####################################################
#Supported C++ compilers: GCC, MSVC, Clang, Intel
# Applies, at a directory scope, CXX flag(s) for a specific compiler
# E.g. to apply C++ flags for GCC only, a dedicated function is available:
# loki_gcc_cxxflags(-O2 -pedantic)
# Compiler agnostic version: loki_cxxflags
#####################################################

#GCC
function(loki_gcc_cxxflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:${option}>)
  endforeach()
endfunction(loki_gcc_cxxflags)

#MSVC
function(loki_msvc_cxxflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:${option}>)
  endforeach()
endfunction(loki_msvc_cxxflags)

#Intel
function(loki_intel_cxxflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Intel>>:${option}>)
  endforeach()
endfunction(loki_intel_cxxflags)

#Clang
function(loki_clang_cxxflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>>:${option}>)
  endforeach()
endfunction(loki_clang_cxxflags)

#Compiler agnostic
function(loki_cxxflags)
  foreach(option ${ARGV})
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${option}>)
  endforeach()
endfunction(loki_cxxflags)

#####################################################
#Supported C compilers: GCC, MSVC, Clang, Intel
# Applies, at a directory scope, C flag(s) for a specific compiler
# E.g. to apply C flags for GCC only, a dedicated function is available:
# loki_gcc_cflags(-O2 -pedantic)
#####################################################

#GCC
function(loki_gcc_cflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:GNU>>:${option}>)
  endforeach()
endfunction(loki_gcc_cflags)

#MSVC
function(loki_msvc_cflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:MSVC>>:${option}>)
  endforeach()
endfunction(loki_msvc_cflags)

#Intel
function(loki_intel_cflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:Intel>>:${option}>)
  endforeach()
endfunction(loki_intel_cflags)

#Clang
function(loki_clang_cflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:Clang>>:${option}>)
  endforeach()
endfunction(loki_clang_cflags)

#Compiler agnostic
function(loki_cflags)
  foreach(option ${ARGV})
    add_compile_options($<$<COMPILE_LANGUAGE:C>:${option}>)
  endforeach()
endfunction(loki_cflags)

#####################################################
#Supported C/C++ compilers: GCC, MSVC, Clang, Intel
# Applies, at a directory scope, C AND C++ flag(s) for a specific compiler
# E.g. to apply C flags for GCC only, a dedicated function is available:
# loki_gcc_ccxxflags(-O2 -pedantic)
#####################################################

#GCC
function(loki_gcc_ccxxflags)
  loki_gcc_cflags(${ARGV})
  loki_gcc_cxxflags(${ARGV})
endfunction(loki_gcc_ccxxflags)

#MSVC
function(loki_msvc_ccxxflags)
  loki_msvc_cflags(${ARGV})
  loki_msvc_cxxflags(${ARGV})
endfunction(loki_msvc_ccxxflags)

#Intel
function(loki_intel_ccxxflags)
  loki_intel_cflags(${ARGV})
  loki_intel_cxxflags(${ARGV})
endfunction(loki_intel_ccxxflags)

#Clang
function(loki_clang_ccxxflags)
  loki_clang_cflags(${ARGV})
  loki_clang_cxxflags(${ARGV})
endfunction(loki_clang_ccxxflags)

#Compiler agnostic
function(loki_ccxxflags)
  loki_cflags(${ARGV})
  loki_cxxflags(${ARGV})
endfunction(loki_ccxxflags)


#####################################################
#Supported Fortran compilers: GCC, Intel
# Applies, at a directory scope, CXX flag(s) for a specific compiler
# E.g. to apply C++ flags for GCC only, a dedicated function is available:
# loki_gcc_cxxflags(-O2 -pedantic)
#####################################################

#GCC
function(loki_gcc_fflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:${option}>)
  endforeach()
endfunction(loki_gcc_fflags)

#Intel
function(loki_intel_fflags)
  foreach(option ${ARGV})
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:Intel>>:${option}>)
  endforeach()
endfunction(loki_intel_fflags)

#Compiler agnostic
function(loki_fflags)
  foreach(option ${ARGV})
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:${option}>)
  endforeach()
endfunction(loki_fflags)


