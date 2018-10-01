#!/bin/sh
# Copyright (c) 2017,2018 Yoshihide Yoshimoto

set -- `getopt b:s: "$@"`
if [ $? != 0 ]; then
    echo "Usage: $0 [-b path_to_wfn2respack ] [-s path_to_strconv ] wfn_base_name" 1>&2
    exit 1
fi

WFN2RESPACK=wfn2respack
STRCONV=strconv

while [ $# -gt 0 ]; do
    case $1 in
	--) shift; break ;;
	-b) WFN2RESPACK=`readlink -f $2`; shift;;
	-s) STRCONV=`readlink -f $2`; shift;;
    esac
    shift
done

which $WFN2RESPACK > /dev/null
if  [ $? = 1 ]; then
    echo "Error: $WFN2RESPACK was not found in your PATH"
    echo "specify the path to $WFN2RESPACK by -b option"
    exit 1
fi
which $STRCONV > /dev/null
if  [ $? = 1 ]; then
    echo "Error: $STRCONV was not found in your PATH"
    echo "specify the path to $STRCONV by -s option"
    exit 1
fi

target=$1

error=0
if [ ! -f $target.wfn ]; then
    echo "$target.wfn not found"
    error=1
fi
if [ ! -f $target.str ]; then
   echo "$target.str not found"
   error=1
fi
if [ $error = 1 ]; then
    exit 1
fi

if [ ! -e dir-wfn ]; then
    mkdir dir-wfn
fi
( cd dir-wfn; $WFN2RESPACK ../$target.wfn 2> /dev/null ; exit $? )
status=$?
if [ $status -ne 0 ]; then
    echo "wfn2respack failed"
    rmdir dir-wfn > /dev/null 2>&1
    exit 1
fi

( cd dir-wfn; $STRCONV respack ../$target.str 2> /dev/null ; exit $? )
if [ $status -ne 0 ]; then
    echo "strconv failed"
    rmdir dir-wfn > /dev/null 2>&1
    exit 1
fi

