# Copyright (c) 2017 Yusuke Nomura, Terumasa Tadano 
# Usage of qe2respack (an interface from Quantum Espresso to RESPACK)

1. Edit the Makefile (edit "FC", "FFLAGS", "LINKER", "LDFLAGS", "LAPACK", and "QEROOT" according to your environment). 
   Please edit "QEROOT" to make it possible to use Quantum Espresso iotk toolkit. 
   One also needs to give a link to a LAPACK library (edit "LAPACK" and "LDFLAGS").
   If one compiles Quantum Espresso in MPI environment, one needs a MPI compiler as "LINKER" to build qe2respack.

2. Compile qe2respack.
   Type "make". 

3. Use qe2respack.sh to create input files of RESPACK.
   Edit qe2respack_root in qe2respack.sh. 
   Usage example: "bash qe2respack.sh ./outdir/prefix.save/"
   ./outdir/prefix.save/ is the Quantum Espresso output directory. 
   qe2respack.sh will create dir-wfn directory under a directory where qe2respack.sh is called. 
