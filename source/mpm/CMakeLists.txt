## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set(_src
  material_point.cc
  mpm_handler.cc
  constitutive_model.cc
  )

set(_inst
  material_point.inst.in
  mpm_handler.inst.in
  constitutive_model.inst.in
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/mpm/*.h
  )

define_object_library(object_mpm OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_mpm "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
