try:
	from setuptools import setup, Extension
	from setuptools.command.build_ext import build_ext
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension
	from distutils.command.build_ext import build_ext

import sys, platform, subprocess


def readme():
	with open('python/README.rst') as f:
		return f.read()


class LibMM2Build(build_ext):
	# Uses Makefile to build library, avoids duplicating logic
	# determining which objects to compile but does require
	# end users to have Make (since precompiled wheels are not
	# distributed on PyPI).
	def run(self):
		def compile_libminimap2(*args, **kwargs):
			cmd = ['make', 'libminimap2.a'] + list(args)
			subprocess.check_call(cmd)
		options = []
		if platform.machine() in ["aarch64", "arm64"]:
			options = ["arm_neon=1", "aarch64=1"]
		self.execute(
			compile_libminimap2, options,
			'Compiling libminimap2 using Makefile')
		build_ext.run(self)


setup(
	name = 'mappy',
	version = '2.24',
	url = 'https://github.com/lh3/minimap2',
	description = 'Minimap2 python binding',
	long_description = readme(),
	author = 'Heng Li',
	author_email = 'lh3@me.com',
	license = 'MIT',
	keywords = 'sequence-alignment',
	scripts = ['python/minimap2.py'],
	cmdclass = {'build_ext': LibMM2Build},
	ext_modules = [
		Extension(
			'mappy',
			sources = ['python/mappy.pyx'],
			depends = ['python/cmappy.h', 'python/cmappy.pxd'],
			include_dirs = ['.'],
			extra_objects = ['libminimap2.a'],
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
