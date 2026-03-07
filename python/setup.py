#!/usr/bin/env python3
"""FlexAID∆S: Entropy-driven molecular docking with flexible side-chains.

Setup script for building Python bindings via pybind11.
"""

import sys
import os
from pathlib import Path
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

# Read version from package
version_file = Path(__file__).parent / "flexaidds" / "__version__.py"
exec(version_file.read_text())

# Read README for long description
readme_file = Path(__file__).parent.parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""


class CMakeExtension(Extension):
    """Extension that uses CMake to build pybind11 module."""
    
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """Custom build_ext that invokes CMake for compilation."""
    
    def build_extension(self, ext):
        if not isinstance(ext, CMakeExtension):
            return super().build_extension(ext)
        
        import subprocess
        from setuptools import distutils
        
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            "-DBUILD_PYTHON_BINDINGS=ON",
            f"-DCMAKE_BUILD_TYPE={'Debug' if self.debug else 'Release'}",
        ]
        
        build_args = ["--config", "Debug" if self.debug else "Release"]
        
        # Platform-specific build configuration
        if sys.platform.startswith("darwin"):
            cmake_args += ["-DCMAKE_OSX_DEPLOYMENT_TARGET=10.14"]
        
        # Parallel build
        build_args += ["--", "-j4"]
        
        # Create build directory
        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        
        # Run CMake configure
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=build_temp
        )
        
        # Run CMake build
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=build_temp
        )


setup(
    name="flexaidds",
    version=__version__,  # noqa: F821 (defined in __version__.py exec)
    author="Louis-Philippe Morency",
    author_email="lp.morency@umontreal.ca",
    description="Entropy-driven molecular docking with flexible side-chains",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lmorency/FlexAIDdS",
    project_urls={
        "Documentation": "https://flexaidds.readthedocs.io",
        "Source": "https://github.com/lmorency/FlexAIDdS",
        "Tracker": "https://github.com/lmorency/FlexAIDdS/issues",
    },
    packages=find_packages(),
    ext_modules=[CMakeExtension("flexaidds._core", sourcedir="..")],
    cmdclass={"build_ext": CMakeBuild},
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
    ],
    extras_require={
        "visualization": ["pymol>=2.5.0"],
        "analysis": ["pandas>=1.3.0", "matplotlib>=3.4.0"],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "mypy>=0.950",
            "sphinx>=4.5.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: C++",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ],
    keywords="molecular-docking drug-discovery thermodynamics entropy cheminformatics",
    zip_safe=False,
)
