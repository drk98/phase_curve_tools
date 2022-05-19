#!/usr/bin/env python

"""The setup script."""


from setuptools import setup, find_packages, dist
dist.Distribution().fetch_build_eggs(['Cython>=0.29.0', 'numpy>=1.22'])

from Cython.Build import cythonize
import numpy as np

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy>=1.22']

test_requirements = ['pytest>=3', 'numpy>=1.22']

setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        'setuptools>=18.0',
        'cython',
        'numpy'
    ],
    author="Daniel Kramer",
    author_email='drk98@nau.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A Python package for asteroid phase curves",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='phase_curve_tools',
    name='phase_curve_tools',
    packages=find_packages(include=['phase_curve_tools', 'phase_curve_tools.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/drk98/phase_curve_tools',
    version='0.1.0',
    zip_safe=False,
    ext_modules = cythonize("phase_curve_tools/*.pyx"),
    package_data={"phase_curve_tools": ["phase_curve_tools/*.h", "phase_curve_tools/*.pyx"]},
    include_dirs = [np.get_include()],
)