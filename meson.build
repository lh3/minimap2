project(
  'minimap2',
  'c',
  version : '2.9.0',
  default_options : [
    'buildtype=release',
    'default_library=shared',
    'warning_level=3',
    'c_std=gnu11'],
  license : 'MIT',
  meson_version : '>= 0.41.1')

##########
# CFLAGS #
##########

mm2_flags = ['-DHAVE_KALLOC']
c = meson.get_compiler('c')
foreach cflag: [
  '-Wc++-compat',
  '-Wno-unused-parameter',
  '-Wno-unused-variable',
  '-Wno-sign-compare']
    if c.has_argument(cflag)
      mm2_flags += cflag
    endif
endforeach

if get_option('enable-build-sse4')
  mm2_flags += ['-msse4']
endif

################
# dependencies #
################

# threads
mm2_thread_dep = dependency('threads', required : true)

# zlib
mm2_zlib_dep = dependency('zlib', required : true)

# math
mm2_math_dep = c.find_library('m')

mm2_deps = [mm2_thread_dep, mm2_zlib_dep, mm2_math_dep]

###########
# headers #
###########

if not meson.is_subproject()
  install_headers(
    files([
      'bseq.h',
      'getopt.h',
      'kalloc.h',
      'kdq.h',
      'khash.h',
      'kseq.h',
      'ksort.h',
      'ksw2.h',
      'kthread.h',
      'kvec.h',
      'minimap.h',
      'mmpriv.h',
      'sdust.h']),
    subdir : 'minimap2')
endif

mm2_include_directories = include_directories('.')

#############
# libraries #
#############

# install library if
# - either running as a proper project
# - or using shared libraries
mm2_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

mm2_cpp_sources = files([
  'align.c',
  'bseq.c',
  'chain.c',
  'esterr.c',
  'format.c',
  'getopt.c',
  'hit.c',
  'index.c',
  'kalloc.c',
  'ksw2_extd2_sse.c',
  'ksw2_exts2_sse.c',
  'ksw2_extz2_sse.c',
  'ksw2_ll_sse.c',
  'kthread.c',
  'map.c',
  'misc.c',
  'options.c',
  'pe.c',
  'sdust.c',
  'sketch.c'])

mm2_lib = library(
  'minimap2',
  mm2_cpp_sources,
  install : mm2_lib_install,
  # use boost SONAME practice:
  #   cause ld startup issues before
  #   you even have the chance of running
  #   into ABI issues.
  soversion : meson.project_version(),
  version : meson.project_version(),
  dependencies : mm2_deps,
  include_directories : mm2_include_directories,
  c_args : mm2_flags)

###############
# executables #
###############

# minimap2
if not meson.is_subproject()
  mm2_main = executable(
    'minimap2', files([
      'main.c']),
    install : true,
    dependencies : mm2_deps,
    include_directories : mm2_include_directories,
    link_with : mm2_lib,
    c_args : mm2_flags)
endif

###################
# dependency info #
###################

if not meson.is_subproject()
  import('pkgconfig').generate(
    libraries : mm2_lib,
    version : meson.project_version(),
    name : 'minimap2',
    filebase : 'minimap2',
    subdirs : 'minimap2',
    description : 'A versatile pairwise aligner for genomic and spliced nucleotide sequences')
endif

mm2_dep = declare_dependency(
  include_directories : mm2_include_directories,
  link_with : mm2_lib,
  dependencies : mm2_deps,
  version : meson.project_version())
