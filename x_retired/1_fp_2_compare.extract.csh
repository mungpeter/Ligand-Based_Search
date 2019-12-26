#!/bin/csh

if ($argv != 4) then
  echo ''
  echo "    Usage: x.csh [list name] [column to sort]"
  echo "                 [print top hits] [output prefix]"
  echo ''
  exit
endif

set lst = $argv[1]
set col = $argv[2]
set top = $argv[3]
set out = $argv[4]

foreach s (`cat $lst`)
#  grep '#' $s > $out.txt
  grep -v '#' $s | sort -gk$col >> $out.txt
  sed -n "1,$top p" $out.txt > temp

  ~/Dropbox/9_scripts/3_program/scripts/extract_column.pl \
    temp 6 1 $out.smi
end
