---
project: BsaLib
summary: BsaLib, a Modern Fortran Library for the Bispectral Stochastic Analysis.
project_github: https://github.com/miEsMar/BsaLib
author: Michele Esposito Marzino
email: michele.espositomarzino@gmail.com
source: false
src_dir: ./BsaLib/src
output_dir: ./docs
exclude_dir: BsaLib/src/*
preprocess: true
preprocessor: ifort /E /I:.\BsaLib\src
predocmark: >
docmark_alt: #
predocmark_alt: <
display: public
         protected
graph: false
coloured_edges: true
search: true
sort: alpha
license: LGPL v3
print_creation_date: false
creation_date: %Y-%m-%d %H:%M %z
dbg: true
extra_mods: iso_fortran_env:https://fortranwiki.org/fortran/show/iso_fortran_env
md_extensions: markdown.extensions.toc
---

[TOC]

{!README.md!}

<!-- # `BsaLib` documentation generation file -->

<!-- This file is read by [FORD](https://forddocs.readthedocs.io/en/latest/user_guide/getting_started.html) -->
<!-- to generate HTML documentation for `BsaLib`. -->
