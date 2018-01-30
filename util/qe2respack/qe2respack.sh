#!/bin/sh
# Copyright (c) 2017 Yusuke Nomura, Terumasa Tadano 

if [ $# -ne 1 ]; then
  echo "Please specify the output directory where data-file.xml exists" 1>&2 
  echo "Example: bash qe2respack.sh ./outdir/prefix.save/ " 1>&2 
  exit 1
fi

 qe2respack_root="/home/xxx/RESPACK-20170820/util/qe2respack"               #your machine 
#qe2respack_root="/home/r0088/r008806/src/RESPACK-20170820/util/qe2respack" #sysmteB@issp 

if [ ! -e dir-wfn ]; then
  mkdir dir-wfn
  $qe2respack_root/qe2respack $1
else
  mv dir-wfn dir-wfn_original
  mkdir dir-wfn
  $qe2respack_root/qe2respack $1
  echo "dir-wfn exists. Original directory is saved as dir-wfn_original." 
  echo "If it is not necessary, please remove it (rm -r dir-wfn_original)." 
fi

