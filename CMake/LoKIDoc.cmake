find_program(LOKI_BUILD_HTML_DOCS tex4ht)

# Initially, these targets do nothing. doc-projects can add a dependency of
# these targets on the targets that are defined for those projects.
add_custom_target(doc_clean)
add_custom_target(doc_pdf)

set(DOC_ALL_DEPS doc_pdf)

if(LOKI_BUILD_HTML_DOCS)
  message(STATUS "Building HTML documentation")

  add_custom_target(doc_html)
  list(APPEND DOC_ALL_DEPS doc_html)
endif()

add_custom_target(doc_all COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
                                  --target ${DOC_ALL_DEPS})

include(CMake/UseLATEX.cmake)

# CMAKE commands to create a target that produces output files (pictures,
# typically) by running a Gnuplot script.
#
# Function lokib_add_gnuplot_pictures accepts the target name as argument.
# Options are: GNUFILE: the name of the script to run. OUTFILES: the names of
# the files that are created by the script DEPFILES: the names of the files on
# which the output targets depend. Normally this is the set of data files (if
# any) that are plotted. The GNUFILE is added automatically as a dependency.
#
# As an example, assume that a GNUPLOT script md2d_pictures.gnu in some
# directory creates two files efield.eps and potential.pdf from two data files
# efield_cv.dat and potential_np.dat. This can be achieved by the following text
# in the CMakeLists.txt in that directory:
#
# lokib_add_gnuplot_pictures(doc_md2d_figures_em GNUFILE md2d_pictures.gnu
# OUTFILES efield.eps potential.pdf DEPFILES efield_cv.dat potential_np.dat )
#
# The target name doc_md2d_figures_em can be chosen arbitrarily, but must be
# unique throughout the project. We suggest the following naming scheme: if the
# pictures reside in a subdirectory figures/em of a directory where UseLATEX is
# used to create a document with target name doc_md2d, use that target name and
# append the subdirectory name, with slashes changed into underscores. The latex
# document should be made to depend on the Gnuplot target, so the pictures will
# be created/updated before latex is run. This means, for this example, that
# doc_md2d_figures_em should be added to the DEPENDS option of the latex target,
# as in:
#
# add_latex_document( doc_md2d ... DEPENDS doc_md2d_figures_em ... )
#
# TODO: in principle DEPFILES and OUTFILES can be extracted from the gnuplot
# file by scanning for content like plot 'file' and set output 'file'. But we
# have to be careful; such scan must be repeated when the gnuplot file changes.
# That means that cmake has to be run again.
#
function(lokib_add_gnuplot_pictures tgt_name)

  set(options)
  set(oneValueArgs GNUFILE)
  set(multiValueArgs OUTFILES DEPFILES)
  cmake_parse_arguments(PIC "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})
  add_custom_command(
    OUTPUT ${tgt_name}.tag
    # PRODUCTS ${PIC_OUTFILES}
    COMMAND ${CMAKE_COMMAND} -E echo
            "Running ${GNUPLOT_EXECUTABLE} ${PIC_GNUFILE}"
    COMMAND ${GNUPLOT_EXECUTABLE} ${PIC_GNUFILE}
    COMMAND ${CMAKE_COMMAND} -E touch ${tgt_name}.tag
    DEPENDS ${PIC_DEPFILES})
  add_custom_target(${tgt_name} DEPENDS ${tgt_name}.tag)
endfunction()

macro(lokib_prepare_latex_doc_dir tgt_name)

  set(LATEX_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
  set(DOC_TGT_ALL_DEPS doc_${tgt_name}_pdf)

  if(LOKI_BUILD_HTML_DOCS)
    list(APPEND DOC_TGT_ALL_DEPS doc_${tgt_name}_html)
  endif()

  add_custom_target(
    doc_${tgt_name}_all
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target
            ${DOC_TGT_ALL_DEPS}
    COMMENT "Building ${tgt_name} documentation files.")

  # -E rm is introduced in version 3.17, in favor of -E remove_directory
  add_custom_target(
    doc_${tgt_name}_clean
    COMMAND ${CMAKE_COMMAND} -E rm -f ${tgt_name}*.*
    COMMAND ${CMAKE_COMMAND} -E rm -f *.png
    COMMENT "Cleaning ${tgt_name} files."
    DEPENDS doc_${tgt_name}_auxclean)

  add_dependencies(doc_clean doc_${tgt_name}_clean)
  add_dependencies(doc_pdf doc_${tgt_name}_pdf)

  if(LOKI_BUILD_HTML_DOCS)
    add_dependencies(doc_html doc_${tgt_name}_html)
  endif()

endmacro(lokib_prepare_latex_doc_dir)

# Install the files that are produced by LoKI-B latex target 'doc_${tgt_name}'
# to directory dst_dir relative to the installation prefix. Here dst_dir is an
# optional second argument of this function; If dst_dir is not provided, it is
# recalculated as the directory of the CMakeLists.txt file from which this
# function is called, relative to the toplevel source directory. The html files
# (and related files, such as pictures, and the css) are copied to subdirectory
# html/, the PDF file to the destination directory itself.
#
macro(lokib_install_latex_doc_dir tgt_name)

  # The directory we wish to install
  set(tgt_dir ${CMAKE_CURRENT_BINARY_DIR})
  # Calculate the (relative) destination directory. If dst_dir is provided by
  # the caller (as second argument), we use that directory unmodified.
  if(${ARGC} EQUAL 2)
    set(dst_dir ${ARGV1})
  else()
    # Absolute source directory of the cmake target currently being built
    set(abs_dir ${CMAKE_CURRENT_SOURCE_DIR})
    # Calculate the path of the parent dir, relative to CMAKE_SOURCE_DIR
    file(RELATIVE_PATH dst_dir ${CMAKE_SOURCE_DIR} ${abs_dir})
  endif()
  # message("tgt: ${tgt_dir}") message("dst: ${dst_dir}") trailing slash:
  # install the *contents* of ${tgt_dir}
  install(
    DIRECTORY ${tgt_dir}/
    DESTINATION ${dst_dir}/html
    FILES_MATCHING
    PATTERN "${tgt_name}*.html"
    PATTERN "${tgt_name}*.xml"
    PATTERN "${tgt_name}*.png"
    PATTERN "${tgt_name}*.css")
  install(
    FILES ${tgt_dir}/${tgt_name}.pdf
    DESTINATION ${dst_dir}
    OPTIONAL)

endmacro(lokib_install_latex_doc_dir tgt_name)

macro(lokib_finalize_latex_doc_dir tgt_name)

  # nothing, at the moment

endmacro(lokib_finalize_latex_doc_dir)
