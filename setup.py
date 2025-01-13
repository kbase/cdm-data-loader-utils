# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.md") as f:
    _readme = f.read()

with open("LICENSE") as f:
    _license = f.read()

setup(
    name="CDM-Data-Loader-Utils",
    version="0.0.1",
    description="CDM-Data-Loader-Utils",
    long_description_content_type="text/x-rst",
    long_description=_readme,
    author="KBase developers",
    author_email="KBase developers",
    url="https://github.com/kbase/cdm-data-loader-utils",
    license=_license,
    packages=find_packages(exclude=("docs")),
    package_data={

    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Natural Language :: English",
    ],
    install_requires=[
        "pandas >= 2.2.2",
    ],
    tests_require=[
        "pytest",
    ],
    project_urls={
        "Documentation": "",
        "Issues": "",
    },
)
