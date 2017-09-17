try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension

cmdclass = {}

try:
	from Cython.Build import build_ext
except ImportError: # without Cython
	module_src = 'python/minimap2.c'
else: # with Cython
	module_src = 'python/minimap2.pyx'
	cmdclass['build_ext'] = build_ext

import sys
sys.path.append('python')

setup(
	name = 'minimap2',
	version = '2.2rc',
	url = 'https://github.com/lh3/minimap2',
	description = 'Minimap2 python binding',
	author = 'Heng Li',
	author_email = 'lh3@me.com',
	license = 'MIT',
	keywords = ['bioinformatics', 'sequence-alignment'],
    ext_modules = [Extension('minimap2',
		sources = [module_src, 'align.c', 'bseq.c', 'chain.c', 'format.c', 'hit.c', 'index.c',
				   'ksw2_extd2_sse.c', 'ksw2_exts2_sse.c', 'ksw2_extz2_sse.c', 'ksw2_ll_sse.c',
				   'kalloc.c', 'kthread.c', 'map.c', 'misc.c', 'sdust.c', 'sketch.c'],
		depends = ['minimap.h', 'bseq.h', 'kalloc.h', 'kdq.h', 'khash.h', 'kseq.h', 'ksort.h',
				   'ksw2.h', 'kthread.h', 'kvec.h', 'mmpriv.h', 'sdust.h',
				   'python/cminimap2.h', 'python/cminimap2.pxd'],
		extra_compile_args = ['-msse4'], # WARNING: ancient x86_64 CPUs don't have SSE4
		include_dirs = ['.'],
		libraries = ['z', 'm', 'pthread'])],
	cmdclass = cmdclass)
