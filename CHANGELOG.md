# Changelog

## [3.1.1](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.1.0...v3.1.1) (2025-10-29)


### Bug Fixes

* when all samples in an analysis had no 2nd skani species hit, it would fail, now handles this gracefully with empty cols. ([6e79375](https://github.com/RIVM-bioinformatics/juno-assembly/commit/6e793759617a837131d2dcd2c0fb8ae55e5e0bcc))

## [3.1.0](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.12...v3.1.0) (2025-09-24)


### Features

* added skani species labelling and new species warnings for ANI &lt;95% ([#50](https://github.com/RIVM-bioinformatics/juno-assembly/issues/50)) ([5b66f70](https://github.com/RIVM-bioinformatics/juno-assembly/commit/5b66f70cd20269061f3beaa931e80745b4991f4d))


### Bug Fixes

* Ensured log files are generated for all processes (see [#48](https://github.com/RIVM-bioinformatics/juno-assembly/issues/48)) ([434727e](https://github.com/RIVM-bioinformatics/juno-assembly/commit/434727ef64ddbca0b77f4353eb3799cf34eecc57))
* Made `subsample.py` more robust to errors (see [#49](https://github.com/RIVM-bioinformatics/juno-assembly/issues/49)) ([434727e](https://github.com/RIVM-bioinformatics/juno-assembly/commit/434727ef64ddbca0b77f4353eb3799cf34eecc57))

## [3.0.12](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.11...v3.0.12) (2025-09-11)


### Bug Fixes

* Increase multiqc memory to handle large datasets ([280ff94](https://github.com/RIVM-bioinformatics/juno-assembly/commit/280ff94a05ec14a8cd105b9e74ebf616d0dc905b))

## [3.0.11](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.10...v3.0.11) (2025-07-01)


### Bug Fixes

* Add -Xmx100g flag to sort_paired_fastq rule ([ff7e42b](https://github.com/RIVM-bioinformatics/juno-assembly/commit/ff7e42b3924b23fdb12e38a7c4116447e09baaa9))

## [3.0.10](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.9...v3.0.10) (2025-05-12)


### Bug Fixes

* Add flag to pileup rule to increase memory usage ([0ea1c91](https://github.com/RIVM-bioinformatics/juno-assembly/commit/0ea1c9133c74c5df7c63ee37615691254ffcb953))

## [3.0.9](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.8...v3.0.9) (2025-04-15)


### Bug Fixes

* updated juno library version to latest ([8307927](https://github.com/RIVM-bioinformatics/juno-assembly/commit/8307927f1d44b999b930b354a4ce2ad1d140fd93))

## [3.0.8](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.7...v3.0.8) (2025-04-08)


### Bug Fixes

* increase memory limit for kraken in order to run properly on cluster ([5c34d5c](https://github.com/RIVM-bioinformatics/juno-assembly/commit/5c34d5ca6a210f9551c9be4964e698dfa2357988))

## [3.0.7](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.6...v3.0.7) (2024-07-12)


### Bug Fixes

* added fix for bug that caused some samples with results missing in qc report ([e0d5685](https://github.com/RIVM-bioinformatics/juno-assembly/commit/e0d568550757ab859611ed1927ae03dca3c817a3))

## [3.0.6](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.5...v3.0.6) (2024-06-24)


### Bug Fixes

* increase sort mem ([f359b22](https://github.com/RIVM-bioinformatics/juno-assembly/commit/f359b228ca7394f41011cb85d21f64bcda5d797d))

## [3.0.5](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.4...v3.0.5) (2024-05-21)


### Dependencies

* fix pulp version ([0669908](https://github.com/RIVM-bioinformatics/juno-assembly/commit/066990861d48e0ec511361aa1924eb1201747c0a))

## [3.0.4](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.3...v3.0.4) (2024-05-15)


### Bug Fixes

* dont save k2 classification per read ([e7f3733](https://github.com/RIVM-bioinformatics/juno-assembly/commit/e7f373388b55c9e0026c35e80afcac9abcdee751))
* increase max mem for bbmap ([2137cd0](https://github.com/RIVM-bioinformatics/juno-assembly/commit/2137cd0eccb2fc7a09a1d90f594b5004c3c98d67))


### Dependencies

* prevent installation of python 1.6 in env ([023ebaf](https://github.com/RIVM-bioinformatics/juno-assembly/commit/023ebaf9988624ac10dc47e540b8bb44aa134455))

## [3.0.3](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.2...v3.0.3) (2023-11-13)


### Bug Fixes

* update conda envs ([094d0f1](https://github.com/RIVM-bioinformatics/juno-assembly/commit/094d0f17329d20ee8991d6d235bdc77508961c10))

## [3.0.2](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.1...v3.0.2) (2023-08-10)


### Bug Fixes

* sort reads before subsample for repeatability ([c14a50f](https://github.com/RIVM-bioinformatics/juno-assembly/commit/c14a50f3d286bc3f2e4410ba3208b96337cc8211))

## [3.0.1](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v3.0.0...v3.0.1) (2023-07-19)


### Bug Fixes

* conda env file path ([6b03b3c](https://github.com/RIVM-bioinformatics/juno-assembly/commit/6b03b3c39c1a3a6c0283a9cc5dd11133a98db1d1))


### Dependencies

* remove anaconda ([c9440f8](https://github.com/RIVM-bioinformatics/juno-assembly/commit/c9440f8de666bd88ac9e89b41c002b20a2462561))

## [2.2.0](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v2.1.3...v2.2.0) (2023-06-30)


### Features

* select checkm genus ([53d7f2b](https://github.com/RIVM-bioinformatics/juno-assembly/commit/53d7f2b257e9da013b83567ecfff747b5805ccd2))


### Bug Fixes

* capitalise input genus ([6a3f78b](https://github.com/RIVM-bioinformatics/juno-assembly/commit/6a3f78b651e950a2a7bd9af7c6ffcca6670e0a96))

## [2.1.3](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v2.1.2...v2.1.3) (2023-06-09)


### Bug Fixes

* update run pipeline to circumvent a bug in mamba installation procedure within irods ([d97d625](https://github.com/RIVM-bioinformatics/juno-assembly/commit/d97d625eca859ebe734df9e1c85a5e244b3b2804))

## [2.1.2](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v2.1.1...v2.1.2) (2023-04-20)


### Miscellaneous Chores

* release 2.1.2 ([e8e2a73](https://github.com/RIVM-bioinformatics/juno-assembly/commit/e8e2a73b86f8dbd2b7736e68f71c820e0bf9f2e6))

## [2.1.1](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v2.1.0...v2.1.1) (2023-03-06)


### Bug Fixes

* Fix singularity prefix to be str in stead of Path ([3949e4d](https://github.com/RIVM-bioinformatics/juno-assembly/commit/3949e4d28a9a8d5a8f6c1130c0b750f24a4cae52))
* Unset prefix if not using containers ([9907dd4](https://github.com/RIVM-bioinformatics/juno-assembly/commit/9907dd460e377b39911ed4b89061ba1eea67634e))


### Dependencies

* Remove defaults channel ([b54fc55](https://github.com/RIVM-bioinformatics/juno-assembly/commit/b54fc555730825adac2040f2f83c3eeba110c5e7))

## [2.1.0](https://github.com/RIVM-bioinformatics/juno-assembly/compare/v2.0.8...v2.1.0) (2023-03-02)


### Features

* Add qc report ([6f56ccc](https://github.com/RIVM-bioinformatics/juno-assembly/commit/6f56cccefcd2f287bc501bbd99208edb00783359))


### Bug Fixes

* Added units to columns ([20bfc70](https://github.com/RIVM-bioinformatics/juno-assembly/commit/20bfc70df191f773804c2bfcbe681d133580b740))
* Fix github action using singularity ([9580933](https://github.com/RIVM-bioinformatics/juno-assembly/commit/9580933cda27216ede7b62e41236d4931d003fee))
* Set sample column datatype to string ([d794089](https://github.com/RIVM-bioinformatics/juno-assembly/commit/d7940898d27753fdb33d23bc305a91ba7efd1824))


### Dependencies

* Add openpyxl as depedency ([249e5aa](https://github.com/RIVM-bioinformatics/juno-assembly/commit/249e5aa63bc4a6c9f4c30f885a89512af5749a2e))
