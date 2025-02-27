# Changelog

## [6.4.0](https://github.com/deepgenomics/GenomeKit/compare/v6.3.0...v6.4.0) (2025-02-27)


### Features

* add appris data ([#132](https://github.com/deepgenomics/GenomeKit/issues/132)) ([6914b05](https://github.com/deepgenomics/GenomeKit/commit/6914b055d348892070be40c153f9d3e77adbc530))
* support gencode v47 ([#131](https://github.com/deepgenomics/GenomeKit/issues/131)) ([210ec32](https://github.com/deepgenomics/GenomeKit/commit/210ec329e487949b431f0262bdad9de1ecbcbf22))

## [6.3.0](https://github.com/deepgenomics/GenomeKit/compare/v6.2.0...v6.3.0) (2025-01-23)


### Features

* support MANE Select and Clinical Plus transcripts ([#127](https://github.com/deepgenomics/GenomeKit/issues/127)) ([98d0757](https://github.com/deepgenomics/GenomeKit/commit/98d075746d3c8a8937149e498d8fef85e1579db1))
* support refseq anno GCF_000001405.40-RS_2024_08 ([#125](https://github.com/deepgenomics/GenomeKit/issues/125)) ([b8e1ca3](https://github.com/deepgenomics/GenomeKit/commit/b8e1ca35116b4d6b4ac04f5d3526f12786a02b05)), closes [#124](https://github.com/deepgenomics/GenomeKit/issues/124)

## [6.2.0](https://github.com/deepgenomics/GenomeKit/compare/v6.1.1...v6.2.0) (2024-12-20)


### Features

* S3 data manager ([#121](https://github.com/deepgenomics/GenomeKit/issues/121)) ([85f1cbd](https://github.com/deepgenomics/GenomeKit/commit/85f1cbd93f17235e195a026e19db6509d43c44a2))

## [6.1.1](https://github.com/deepgenomics/GenomeKit/compare/v6.1.0...v6.1.1) (2024-12-16)


### Documentation

* simplify strandedness instructions for bedgraph/wig ([#118](https://github.com/deepgenomics/GenomeKit/issues/118)) ([60c44dd](https://github.com/deepgenomics/GenomeKit/commit/60c44ddedd09fdd713ace12e18d5fa34ce844e15))


### Miscellaneous Chores

* release 6.1.1 ([#119](https://github.com/deepgenomics/GenomeKit/issues/119)) ([e081b55](https://github.com/deepgenomics/GenomeKit/commit/e081b55cf3cd29b306510d69dc5db7a2fc0cbb67))

## [6.1.0](https://github.com/deepgenomics/GenomeKit/compare/v6.0.3...v6.1.0) (2024-12-03)


### Features

* deprecate py38, support py11 py12 ([#101](https://github.com/deepgenomics/GenomeKit/issues/101)) ([cd94123](https://github.com/deepgenomics/GenomeKit/commit/cd941236dd57abbb9f477e8fec89513c88775617))


### Documentation

* prepare required documentation for AWS Open Data ([#116](https://github.com/deepgenomics/GenomeKit/issues/116)) ([364d723](https://github.com/deepgenomics/GenomeKit/commit/364d72323eaee22147064f7ad6d7027dbafea384))
* starter_build script download ([#96](https://github.com/deepgenomics/GenomeKit/issues/96)) ([1245558](https://github.com/deepgenomics/GenomeKit/commit/12455581374025313e767a6fdfef96b23a77987e))

## [6.0.3](https://github.com/deepgenomics/GenomeKit/compare/v6.0.2...v6.0.3) (2024-11-06)


### Miscellaneous Chores

* fix windows build ([#106](https://github.com/deepgenomics/GenomeKit/issues/108)) ([5fe3814](https://github.com/deepgenomics/GenomeKit/commit/5fe3814f7c7ec0654aabbe2d49064b98b834e564))

## [6.0.2](https://github.com/deepgenomics/GenomeKit/compare/v6.0.1...v6.0.2) (2024-11-05)


### Miscellaneous Chores

* fix windows build ([#106](https://github.com/deepgenomics/GenomeKit/issues/106)) ([e1754c7](https://github.com/deepgenomics/GenomeKit/commit/e1754c71a7b1787e12a3c88ef98a65e1106b0094))

## [6.0.1](https://github.com/deepgenomics/GenomeKit/compare/v6.0.0...v6.0.1) (2024-10-18)


### Bug Fixes

* **data:** exception raising was dependent on an envvar before ([#99](https://github.com/deepgenomics/GenomeKit/issues/99)) ([3d461c1](https://github.com/deepgenomics/GenomeKit/commit/3d461c184d71eaca48bebc60f1f40dcba95fcafe))
* **dna:** consistent allow_outside_chromosome defaults ([#98](https://github.com/deepgenomics/GenomeKit/issues/98)) ([611970d](https://github.com/deepgenomics/GenomeKit/commit/611970d9b1c6cfaca2e2dd7938ec9d3d46cef697))
* **track:** OOB issue with negative strands on u1 tracks ([#104](https://github.com/deepgenomics/GenomeKit/issues/104)) ([53ae2f5](https://github.com/deepgenomics/GenomeKit/commit/53ae2f5881c2c97b572cb787098544f59e2e2d14))

## [6.0.0](https://github.com/deepgenomics/GenomeKit/compare/v5.2.5...v6.0.0) (2024-09-17)


### ⚠ BREAKING CHANGES

* **vcfbin:** VCFTable.GT_HETEROZYGOUS is replaced with VCFTable.GT_HETEROZYGOUS_UNPHASED

### Features

* **vcfbin:** support phased VCF ([#92](https://github.com/deepgenomics/GenomeKit/issues/92)) ([a1f0b2a](https://github.com/deepgenomics/GenomeKit/commit/a1f0b2ad8a86b21dcf68134e9e8c87ce9f79d665))

## [5.2.5](https://github.com/deepgenomics/GenomeKit/compare/v5.2.4...v5.2.5) (2024-08-27)


### Bug Fixes

* **track:** sparse encoding broken ([#88](https://github.com/deepgenomics/GenomeKit/issues/88)) ([c8f6732](https://github.com/deepgenomics/GenomeKit/commit/c8f6732e1be4325e56675c958c852a24a652056b))

## [5.2.4](https://github.com/deepgenomics/GenomeKit/compare/v5.2.3...v5.2.4) (2024-08-22)


### Bug Fixes

* **track:** handle res&gt;1 for track interval interrogation ([#84](https://github.com/deepgenomics/GenomeKit/issues/84)) ([246a6dc](https://github.com/deepgenomics/GenomeKit/commit/246a6dc26f524aa70cb9c21c0c5ab9c2d0caeff5))

## [5.2.3](https://github.com/deepgenomics/GenomeKit/compare/v5.2.2...v5.2.3) (2024-08-16)


### Bug Fixes

* assert that appris transcripts match the requested gene ([#70](https://github.com/deepgenomics/GenomeKit/issues/70)) ([fd93609](https://github.com/deepgenomics/GenomeKit/commit/fd93609f13d67d3eac84aae3591d30a833108c4f))

## [5.2.2](https://github.com/deepgenomics/GenomeKit/compare/v5.2.1...v5.2.2) (2024-08-13)


### Bug Fixes

* **track:** check inplace dtype ([#79](https://github.com/deepgenomics/GenomeKit/issues/79)) ([5a2cb38](https://github.com/deepgenomics/GenomeKit/commit/5a2cb3895565846ae0440d59da120ee7ca8bc6e0))

## [5.2.1](https://github.com/deepgenomics/GenomeKit/compare/v5.2.0...v5.2.1) (2024-08-12)


### Bug Fixes

* **track:** segfault when decoding large dims ([#75](https://github.com/deepgenomics/GenomeKit/issues/75)) ([6f5c29e](https://github.com/deepgenomics/GenomeKit/commit/6f5c29ec610018994f0d451624edad0053013c15))
* **track:** strided decoding ([#77](https://github.com/deepgenomics/GenomeKit/issues/77)) ([1fc1b37](https://github.com/deepgenomics/GenomeKit/commit/1fc1b378875a0376c7f7d41a1335a8cccc1b8325))

## [5.2.0](https://github.com/deepgenomics/GenomeKit/compare/v5.1.1...v5.2.0) (2024-08-01)


### Features

* **gtrack:** inplace output ([#73](https://github.com/deepgenomics/GenomeKit/issues/73)) ([17838b7](https://github.com/deepgenomics/GenomeKit/commit/17838b7a7b78aa067eecf817193bee1ec07f9245))
* **gtrack:** interval iteration ([#71](https://github.com/deepgenomics/GenomeKit/issues/71)) ([1e2b672](https://github.com/deepgenomics/GenomeKit/commit/1e2b67222fd3504dfd112925893f633bf665f11c))

## [5.1.1](https://github.com/deepgenomics/GenomeKit/compare/v5.1.0...v5.1.1) (2024-06-22)


### Bug Fixes

* broken installation on osx-arm64 ([#65](https://github.com/deepgenomics/GenomeKit/issues/65)) ([6bd68c2](https://github.com/deepgenomics/GenomeKit/commit/6bd68c2ae79ac37cca1b01a9dc3cde582c92bfa2))
* check genome.dna(Interval) arg type ([#68](https://github.com/deepgenomics/GenomeKit/issues/68)) ([f7b6e5c](https://github.com/deepgenomics/GenomeKit/commit/f7b6e5c998b8492cb5142b4fb83a5c3b86c4c49f)), closes [#67](https://github.com/deepgenomics/GenomeKit/issues/67)

## [5.1.0](https://github.com/deepgenomics/GenomeKit/compare/v5.0.2...v5.1.0) (2024-06-13)


### Features

* read/write files for refg name hash lookups ([#62](https://github.com/deepgenomics/GenomeKit/issues/62)) ([8074701](https://github.com/deepgenomics/GenomeKit/commit/80747011a83d731b7e1f594e98f1996ac69085fd))


### Bug Fixes

* avoid accessing released memory ([#63](https://github.com/deepgenomics/GenomeKit/issues/63)) ([1ab7456](https://github.com/deepgenomics/GenomeKit/commit/1ab74561c5a80ed5a4c47c02c8cbdeb55040260c))
* ensure correctness of signed/unsigned comparisons ([#55](https://github.com/deepgenomics/GenomeKit/issues/55)) ([23c1d03](https://github.com/deepgenomics/GenomeKit/commit/23c1d0356a582933c8971f3eef408562727c18a6))

## [5.0.2](https://github.com/deepgenomics/GenomeKit/compare/v5.0.1...v5.0.2) (2024-05-09)


### Bug Fixes

* avoid errors when blob metadata is None ([#54](https://github.com/deepgenomics/GenomeKit/issues/54)) ([e850c54](https://github.com/deepgenomics/GenomeKit/commit/e850c5417f285e2f4f82d30cff2555c12e590be9))
* avoid user confusion when loading annotation ([#53](https://github.com/deepgenomics/GenomeKit/issues/53)) ([e556030](https://github.com/deepgenomics/GenomeKit/commit/e556030bebd2c214ada791909925fa151a016017))
* provide a better hint in case of permissions error ([#42](https://github.com/deepgenomics/GenomeKit/issues/42)) ([cb21280](https://github.com/deepgenomics/GenomeKit/commit/cb212800dd9017b3627ce1d6093da9e0fa3c4d3f))


### Documentation

* acknowledge past contributors ([#43](https://github.com/deepgenomics/GenomeKit/issues/43)) ([d0dcfd5](https://github.com/deepgenomics/GenomeKit/commit/d0dcfd5df09b1f36e0627031bad6461766b24b12))
* restore and document link to versionized docs ([#51](https://github.com/deepgenomics/GenomeKit/issues/51)) ([cc22955](https://github.com/deepgenomics/GenomeKit/commit/cc22955ec125d7f802ff1a2d42eec20dc16f8db8)), closes [#50](https://github.com/deepgenomics/GenomeKit/issues/50)

## [5.0.1](https://github.com/deepgenomics/GenomeKit/compare/v5.0.0...v5.0.1) (2024-04-08)


### Bug Fixes

* allow local-only files to be loaded ([#38](https://github.com/deepgenomics/GenomeKit/issues/38)) ([a56a5dc](https://github.com/deepgenomics/GenomeKit/commit/a56a5dc690d507a04b429956b5b3fb96e64839d8)), closes [#37](https://github.com/deepgenomics/GenomeKit/issues/37)


### Documentation

* provide a data starter script ([#39](https://github.com/deepgenomics/GenomeKit/issues/39)) ([a257da1](https://github.com/deepgenomics/GenomeKit/commit/a257da1954e750209a0df9848b5d998901dc3d7f))

## [5.0.0](https://github.com/deepgenomics/GenomeKit/compare/v4.3.0...v5.0.0) (2024-04-04)


### Features

* use new dedicated requester-pays public bucket ([#33](https://github.com/deepgenomics/GenomeKit/issues/33)) ([e4fb7a1](https://github.com/deepgenomics/GenomeKit/commit/e4fb7a196e0201d653e7a563ca2150c393a0e749))


### Miscellaneous Chores

* release 5.0.0 ([#36](https://github.com/deepgenomics/GenomeKit/issues/36)) ([cc939ef](https://github.com/deepgenomics/GenomeKit/commit/cc939ef03a039703a0138d60d8ce1ed5709c8a26))

## [4.3.0](https://github.com/deepgenomics/GenomeKit/compare/v4.2.1...v4.3.0) (2024-03-27)


### Features

* compatibility with python 3.10 ([#28](https://github.com/deepgenomics/GenomeKit/issues/28)) ([1d51397](https://github.com/deepgenomics/GenomeKit/commit/1d51397db0f84cb5d543f6dc9de331871ff1c6ab))


### Bug Fixes

* missing mock for dna boundary ([#25](https://github.com/deepgenomics/GenomeKit/issues/25)) ([b058b53](https://github.com/deepgenomics/GenomeKit/commit/b058b5380a98f0b7766f4fa9de880085d0935a1d))

## [4.2.1](https://github.com/deepgenomics/GenomeKit/compare/v4.2.0...v4.2.1) (2024-03-14)


### Miscellaneous Chores

* skip gcs-dependent tests on conda-forge ([#23](https://github.com/deepgenomics/GenomeKit/issues/23)) ([1595cba](https://github.com/deepgenomics/GenomeKit/commit/1595cba8c0071b6ae127aff482cbc9ca5c96ff92))

## [4.2.0](https://github.com/deepgenomics/GenomeKit/compare/v4.1.8...v4.2.0) (2024-03-14)


### Features

* list available genome configs ([#549](https://github.com/deepgenomics/GenomeKit/issues/549)) ([#19](https://github.com/deepgenomics/GenomeKit/issues/19)) ([fb9776c](https://github.com/deepgenomics/GenomeKit/commit/fb9776cd78b296cdb7cdd903d052a0ed6c68178a))


### Documentation

* more standard cmake instructions ([#21](https://github.com/deepgenomics/GenomeKit/issues/21)) ([55a73d8](https://github.com/deepgenomics/GenomeKit/commit/55a73d86301b9514f1e5df5c98428e4f1eb3d464))

## [4.1.8](https://github.com/deepgenomics/GenomeKit/compare/v4.1.7...v4.1.8) (2024-03-14)


### Miscellaneous Chores

* skip bioconda-dependant tests on conda-forge ([#17](https://github.com/deepgenomics/GenomeKit/issues/17)) ([cad0584](https://github.com/deepgenomics/GenomeKit/commit/cad0584f84cf494a1bf6a17efdbabda1c26ec8f0))

## [4.1.7](https://github.com/deepgenomics/GenomeKit/compare/v4.1.6...v4.1.7) (2024-03-13)


### Miscellaneous Chores

* avoid import on Windows ([#15](https://github.com/deepgenomics/GenomeKit/issues/15)) ([0aef946](https://github.com/deepgenomics/GenomeKit/commit/0aef9466fbcea725f16f731e8ce7e499b78d0e32))

## [4.1.6](https://github.com/deepgenomics/GenomeKit/compare/v4.1.5...v4.1.6) (2024-03-13)


### Miscellaneous Chores

* reenable tests on windows ([#13](https://github.com/deepgenomics/GenomeKit/issues/13)) ([fb90358](https://github.com/deepgenomics/GenomeKit/commit/fb9035829806da8ebc2fbfbd627e75c1803f4472))

## [4.1.5](https://github.com/deepgenomics/GenomeKit/compare/v4.1.4...v4.1.5) (2024-03-12)


### Miscellaneous Chores

* revert the windows test change ([#11](https://github.com/deepgenomics/GenomeKit/issues/11)) ([879adb5](https://github.com/deepgenomics/GenomeKit/commit/879adb54963ba7788d681108ca8f494f27637790))

## [4.1.4](https://github.com/deepgenomics/GenomeKit/compare/v4.1.3...v4.1.4) (2024-03-11)


### Miscellaneous Chores

* skip windows tests correctly ([#9](https://github.com/deepgenomics/GenomeKit/issues/9)) ([ce7f1b3](https://github.com/deepgenomics/GenomeKit/commit/ce7f1b3d728e71acf910930c3e900b6899966bf0))

## [4.1.3](https://github.com/deepgenomics/GenomeKit/compare/v4.1.2...v4.1.3) (2024-03-11)


### Miscellaneous Chores

* avoid some tests on windows ([#7](https://github.com/deepgenomics/GenomeKit/issues/7)) ([7359a45](https://github.com/deepgenomics/GenomeKit/commit/7359a450de76cefba51a14dc74c42f659e230ce5))

## [4.1.2](https://github.com/deepgenomics/GenomeKit/compare/v4.1.1...v4.1.2) (2024-03-11)


### Miscellaneous Chores

* avoid compile error on mac ([#5](https://github.com/deepgenomics/GenomeKit/issues/5)) ([1f9d59e](https://github.com/deepgenomics/GenomeKit/commit/1f9d59e35bc66e6c6ae87f7cc12aa6db06dc72a6))

## [4.1.1](https://github.com/deepgenomics/GenomeKit/compare/v4.1.0...v4.1.1) (2024-03-11)


### Miscellaneous Chores

* avoid CI failure on conda-forge ([#3](https://github.com/deepgenomics/GenomeKit/issues/3)) ([468272e](https://github.com/deepgenomics/GenomeKit/commit/468272eb8a0debaaa245771ba64b8017b8dfd247))

## 4.1.0 (2024-03-06)


### Miscellaneous Chores

* move releases to conda-forge ([#1](https://github.com/deepgenomics/GenomeKit/issues/1)) ([c2ddf93](https://github.com/deepgenomics/GenomeKit/commit/c2ddf932ec83acd3483d472dda6eddbdf3b4d4ea))
