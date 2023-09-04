lokib_prepare_latex_doc_dir(cpp_notes)
add_latex_document(cpp_notes.tex
    TARGET_NAME
        doc_cpp_notes
    BIBFILES
        ${CMAKE_SOURCE_DIR}/doc/reflist.bib
    IMAGE_DIRS
#        pictures
    INCLUDE_DIRECTORIES
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}
    INPUTS
        cpp_notes.tex
    EXCLUDE_FROM_ALL
)
lokib_finalize_latex_doc_dir(cpp_notes)

#lokib_install_latex_doc_dir(cpp_notes)
