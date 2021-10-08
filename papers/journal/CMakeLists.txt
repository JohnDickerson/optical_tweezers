CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT( GPU_TWEEZERS NONE)
SET( BIBTEX_COMPILER bibtex )

# The UseLatex.cmake file is up one level
INCLUDE(UseLATEX.cmake)

#SET(LATEX_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
SET(INPUT_FILES)# introduction.tex related.tex)

ADD_LATEX_DOCUMENT( gpu_tweezers.tex
  INPUTS ${INPUT_FILES}
  BIBFILES gpu_tweezers.bib
  IMAGE_DIRS images figures
  DEFAULT_PDF
  )
