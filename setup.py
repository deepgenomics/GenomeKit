# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.

from setuptools import setup, find_packages
from setuptools.command.egg_info import egg_info
from setup import c_ext

COPYRIGHT_FILE = "COPYRIGHT.txt"
LICENSE_FILE = "LICENSE"

tests_require = [
    "twobitreader>=3.1",
]

version = "4.3.0"

# See https://stackoverflow.com/questions/9977889/how-to-include-license-file-in-setup-py-script/66443941#66443941
class egg_info_ex(egg_info):
    """Includes license file into `.egg-info` folder."""

    def run(self):
        # don't duplicate license into `.egg-info` when building a distribution
        if not self.distribution.have_run.get('install', True):
            # `install` command is in progress, copy license
            self.mkpath(self.egg_info)
            self.copy_file(COPYRIGHT_FILE, self.egg_info)
            self.copy_file(LICENSE_FILE, self.egg_info)

        egg_info.run(self)

if __name__ == "__main__":
    setup(
        author="Deep Genomics",
        author_email="info@deepgenomics.com",
        python_requires=">=3.8, <4",
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: Apache Software License",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
        ],
        description="GenomeKit is a Python library for fast and easy access to genomic resources such as sequence, data tracks, and annotations.",
        license="Apache License 2.0",
        license_files=(COPYRIGHT_FILE, LICENSE_FILE,),
        name="genomekit",
        packages=find_packages(include=["genome_kit"]),
        project_urls={
            "Documentation": "https://deepgenomics.github.io/GenomeKit"
        },
        cmdclass={
            'build_ext': c_ext.NoCWarningsBuildExt,
            'egg_info': egg_info_ex
        },
        ext_modules=[c_ext.extension],
        test_suite="tests",
        tests_require=tests_require,
        url="https://github.com/deepgenomics/GenomeKit",
        version=version,
        zip_safe=False,
    )
