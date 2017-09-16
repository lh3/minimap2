try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension

from Cython.Build import cythonize
import sys

sys.path.append('python')

setup(
	name = 'minimap2',
	version = '2.1.1',
	url = 'https://github.com/lh3/minimap2',
	description = 'Minimap2 python binding',
	author = 'Heng Li',
	author_email = 'lh3@me.com',
	keywords = ['bioinformatics', 'sequence-alignment'],
    ext_modules = cythonize([Extension('minimap2',
		['python/minimap2.pyx', 'align.c', 'bseq.c', 'chain.c', 'format.c', 'hit.c', 'index.c', 'kalloc.c',
		 'ksw2_extd2_sse.c', 'ksw2_exts2_sse.c', 'ksw2_extz2_sse.c', 'ksw2_ll_sse.c', 'kthread.c', 'map.c',
		 'misc.c', 'sdust.c', 'sketch.c'],
		extra_compile_args = ['-msse4'], # WARNING: ancient x86_64 CPUs don't have SSE4
		include_dirs = ['.'])])
)
