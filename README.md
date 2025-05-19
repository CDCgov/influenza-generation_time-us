# Estimating the generation time for influenza transmission using household data in the United States

[![License](https://img.shields.io/badge/license-Apache_2.0_license-brightgreen)](LICENSE)
[![GitHub stars](https://img.shields.io/github/stars/CDCgov/influenza-generation_time-us)](https://github.com/CDCgov/influenza-generation_time-us/stargazers)

## Overview

This repository contains the code required to reproduce the results from the study [Chan et al. 2024](https://doi.org/10.1016/j.epidem.2025.100815).
The original [code](https://github.com/will-s-hart/UK-generation-times), developed by [Hart et al. 2022](https://doi.org/10.7554/elife.70767), was translated from MATLAB to R.

For this demonstration, we are not sharing the US household data. Instead, we used the UK household data from the original study. This demonstration was conducted using 1,000 iterations for simplification. 
In the folder `R/Results`, we have uploaded the posterior distributions generated from 1,000,000 iterations, as presented using the US household data.

### Running the code

`R/main.R` contains the functions to estimate generation time across all data stratifications in parallel, while `R/main_scenario.R` demonstrates estimation for each data stratification. The estimation process consists of *parameter assumption*, *Bayesian data augmentation Markov Chain Monte Carlo (MCMC)*, and *plotting posterior distribution*.

### Data sharing statement

The household data are available upon reasonable request and upon completion of required approvals. The R code for estimating the generation time is available at <https://github.com/CDCgov/influenza-generation_time-us>.

### Disclaimer

The conclusions, findings, and opinions expressed by authors contributing to this article do not necessarily reflect the official position of the U.S. Department of Health and Human Services, the Public Health Service, the Centers for Disease Control and Prevention, or the authors' affiliated institutions.

### Contact information

**Louis Yat Hin Chan, PhD, MSc**  
CDC Steven M. Teutsch Prevention Effectiveness (PE) Fellow – Analytics and Modeling Track, Class of 2023  
Applied Research and Modeling (ARM) Team  
Epidemiology and Prevention Branch (EPB)  
Influenza division (ID)  
National Center for Immunization and Respiratory Diseases (NCIRD)  
Centers for Disease Control and Prevention (CDC)  

Work address: 1600 Clifton Road, Atlanta, GA 30329  
Work email: <LouisChan@cdc.gov>  

## Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
