#!/bin/csh

##########################################################################
#
#  Peter MU Ung @ MSSM
#
#  Run fp_comparison for a list of libraries against a target template
#  basically a wrapper to run multiple fp_comparison with one script
#
##########################################################################

if ($#argv != 4) then
  echo ''
  echo '    Usage: x.csh [Template Mol]'
  echo '                 [Library List]'
  echo '                 [Fingerprint for Ranking] (ecfp4, dl, maccs, total)'
  echo '                 [no. of CPU]'
  echo '    e.g.   x.csh templ.smi lib.list ecfp4 4'
  echo ''
  exit
endif

set mol  = $argv[1]
set list = $argv[2]
set fp   = $argv[3]
set cpu  = $argv[4]

##########################################################################

set moln = `basename $mol .smi`

foreach dbase (`cat $list`)

  set name = `basename $dbase .sdf.gz`
  
  if (! -e $moln.$name.fp.txt) then
    time ./1_fp_2_compare.mpi.py \
      $mol $dbase $fp $cpu $moln.$name
  else
    echo $moln.$name.fp.txt exists
  endif

end
