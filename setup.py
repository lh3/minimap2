try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension

import sys, platform, os, re, sysconfig

sys.path.append('python')

extra_compile_args = ['-DHAVE_KALLOC']
include_dirs = ["."]

def simd_args(machine):
	if machine.startswith(('aarch64', 'arm64')):
		return ['-ftree-vectorize', '-DKSW_SSE2_ONLY', '-D__SSE2__', '-Isse2neon/']
	return ['-msse4.1'] # WARNING: ancient x86_64 CPUs don't have SSE4

# on macOS ARCHFLAGS/CFLAGS say what we are compiling for, which is not always the host
archs = []
if sys.platform == 'darwin':
	archs = sorted(set(re.findall(r'-arch\s+(\S+)',
		os.environ.get('ARCHFLAGS', sysconfig.get_config_var('CFLAGS') or ''))))
if len(archs) > 1: # universal build: give each slice its own SIMD flags
	for arch in archs:
		for flag in simd_args(arch):
			extra_compile_args.extend(['-Xarch_' + arch, flag])
else:
	extra_compile_args.extend(simd_args(archs[0] if archs else platform.machine()))

def readme():
	with open('python/README.rst') as f:
		return f.read()

setup(
	name = 'mappy',
	version = '2.31',
	url = 'https://github.com/lh3/minimap2',
	description = 'Minimap2 python binding',
	long_description = readme(),
	author = 'Heng Li',
	author_email = 'lh3@me.com',
	license = 'MIT',
	keywords = 'sequence-alignment',
	scripts = ['python/minimap2.py'],
	ext_modules = [Extension('mappy',
		sources = ['python/mappy.pyx', 'align.c', 'bseq.c', 'lchain.c', 'seed.c', 'format.c', 'hit.c', 'index.c', 'pe.c', 'jump.c', 'options.c',
				   'ksw2_extd2_sse.c', 'ksw2_exts2_sse.c', 'ksw2_extz2_sse.c', 'ksw2_ll_sse.c',
				   'kalloc.c', 'kthread.c', 'map.c', 'misc.c', 'sdust.c', 'sketch.c', 'esterr.c', 'splitidx.c'],
		depends = ['minimap.h', 'bseq.h', 'kalloc.h', 'kdq.h', 'khash.h', 'kseq.h', 'ksort.h',
				   'ksw2.h', 'kthread.h', 'kvec.h', 'mmpriv.h', 'sdust.h',
				   'python/cmappy.h', 'python/cmappy.pxd'],
		extra_compile_args = extra_compile_args,
		include_dirs = include_dirs,
		libraries = ['z', 'm', 'pthread'])],
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'License :: OSI Approved :: MIT License',
		'Operating System :: POSIX',
		'Programming Language :: C',
		'Programming Language :: Cython',
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics'],
	setup_requires=["cython"])
