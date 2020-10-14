from setuptools import setup, find_packages, Extension


libans = Extension('libans',
					language='c++',
					libraries=['stdc++'],
					sources = ['c/binaryANS.cpp'],
					depends=['c/ANStoolkit.cpp'])

setup(name = 'anspy',
		version = '1.0',
		description = 'ANS toolkit in Python',
		author = 'Luca Colagrande',
		author_email = 'luca.colagrande3@gmail.com',
		packages=find_packages(),
		install_requires=[
							'bitstring',
							'numpy',
							'cffi',
							'pathlib'
						],
		ext_modules = [libans],
		python_requires='>=3.6'
		)
