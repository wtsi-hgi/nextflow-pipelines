#!/usr/bin/env bash
echo running del_work_dirs_failed.sh script...
workdir=$1

if [ -z "$workdir" ]
then
      echo script first argument workdir is empty.
else
      echo script first argument workdir is "$workdir"
      echo will now find and remove work dirs with .exitcode files with non-0 exitcode...
      find "$workdir" -maxdepth 3 -mindepth 3 \
	   -name .exitcode \
	   -exec cat {} \; -exec echo ,{} \; | \
	  sed s'/\/.exitcode /\n/'g | sed s'/\/.exitcode//'g | \
	  grep -v '^0,' | sed s'/^[0-9]*,//'g | \
	  xargs --no-run-if-empty rm -r
      echo find and remove done.
fi
