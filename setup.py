# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
import os
from pathlib import Path

from setuptools import setup, find_packages
from setuptools.command.egg_info import egg_info
from setup import c_ext

COPYRIGHT_FILE = "COPYRIGHT.txt"
LICENSE_FILE = "LICENSE"

tests_require = [
    "twobitreader>=3.1",
]

version = "7.1.0"

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
    install_requires = []
    if os.environ.get("GK_BUILD_WHEELS", None) is not None:
        install_requires = [
            "appdirs>=1.4.0",
            "numpy<2.0dev0",
            "google-cloud-storage>=2.10.0",
            "boto3",
            "tqdm",
            "setuptools",
            "importlib-metadata",
            "typing-extensions",
        ]
    setup(
        author="Deep Genomics",
        author_email="info@deepgenomics.com",
        python_requires=">=3.9, <4",
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: Apache Software License",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
        ],
        description="GenomeKit is a Python library for fast and easy access to genomic resources such as sequence, data tracks, and annotations.",
        long_description=(Path(__file__).parent / "README.md").read_text(),
        long_description_content_type='text/markdown',
        install_requires=install_requires,
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
        url=f"https://github.com/deepgenomics/GenomeKit",
        version=version,
        zip_safe=False,
    )
