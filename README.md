# Prognostics Algorithm Library

The Prognostics Algorithm Library is a suite of algorithms implemented in the MATLAB programming language for model-based prognostics (remaining life computation). It includes algorithms for state estimation and prediction, including uncertainty propagation. The algorithms take as inputs component models developed in Matlab, and perform estimation and prediction functions. The library allows the rapid development of prognostics solutions for given models of components and systems. Different algorithms can be easily swapped to do comparative studies and evaluations of different algorithms to select the best for the application at hand.

## Citation

Publications making use of software products obtained from this repository are requested to acknowledge the assistance received by using this repository. Please cite: "M. Daigle; Prognostics Algorithm Library [Computer software]. (2016). Retrieved from https://github.com/nasa/PrognosticsAlgorithmLibrary".

## Installation

Installation can be done in one of two ways. Either (1) use the MATLAB toolbox installer provided in the install folder, which will install the toolbox in your local toolboxes folder and add the folder to your MATLAB path, or (2) copy the source from the MATLAB folder to any desired directory, and add that directory to your MATLAB path. Do not add the subdirectories (the package directories) to your MATLAB path. If the first option is used, then the MATLAB add-on manager can be used to uninstall the package; otherwise, the installation can be removed manually by removing the directory from your MATLAB path and deleting the source.

## User Manual

Instructions for using the software can be found in the following sources:
- [Wiki](https://github.com/nasa/PrognosticsAlgorithmLibrary/wiki)
- [User Manual](https://github.com/nasa/PrognosticsAlgorithmLibrary/blob/master/docs/PrognosticsAlgorithmLibrary-UserManual.pdf)

## Dependencies

Some modules included in this library are dependent on the [Prognostics Model Library](https://github.com/nasa/PrognosticsModelLibrary).

## Compatibility

The PrognosticsModelLibrary has been tested with Matlab R2016a, but should work with older versions, down to at least R2012a.

## Contributions

All contributions are welcome. Issues may be opened using GitHub. To contribute directly, open a pull request against the "develop" branch. Pull requests will be evaluated and integrated into the next official release.

## License

This software is released under the [NASA Open Source Agreement Version 1.3](https://github.com/nasa/PrognosticsAlgorithmLibrary/blob/master/LICENSE.pdf).

## Notices

Copyright © 2016 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.

### Disclaimers
```
No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
```
