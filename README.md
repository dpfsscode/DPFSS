# DPFSS

### Introduction for the scheme

Our paper poposes a DPFSS scheme in the extended model of DPSS to pass functions of the shared secrets to the new committee.
And in our DPFSS construction, we apply Shamir's secret sharing scheme to share and reconstruct a secret.
We provide the code of our scheme (dpfss.h, dpfss.cpp), where the Shamir's secret sharing is additionally prepared in util.h, util.cpp.

### Introduction for the code
This repo consists of several files for running the dpfss scheme, which is detailed as follow:
- `dpfss.h`: header file of DPFSS, which clarifies the class definition and interface DPFSS.
- `util.h`: header file for useful utils such as Shamir's secret sharing, which clarifies the main process to share  and reconstruct a secret.
- `dpfss.cpp`: implementation DPFSS protocol, which consist of the share generation, function handoff and reconstrction.
- `util.cpp`: implementation for utils decleared in `util.h`.
- `test1.cpp`: a sample test file to test to run the procedure of DPFSS scheme.

### Dependencies:

Our code was built on the library of flint 2.8.0, which support large integer operation based on GMP 6.2.1, and use relic 0.5.0 to implement the bilinear pairing. Our code successfully excecuted on a computer with Ubuntu 20.04.5 LTS and the compiler is g++ with version 9.4.0.

### How to run the code

Open the termial and type the following command:

    g++ -o test test1.cpp dpfss.cpp util.cpp -lgmp -lflint -lrelic

This will produce a executable file called `test`, then run `./test` to execute the file.

