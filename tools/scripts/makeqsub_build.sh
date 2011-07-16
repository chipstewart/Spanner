#!/bin/bash

program="Spanner"
args="-d 30 "
inbase="/d1/data/1000G/20101123/phase1.mosaik2/25.samples/pipeline/mosaik"
outbase="/d1/data/1000G/20101123/phase1.mosaik2/25.samples/pipeline/Spanner"
here=$PWD


subscript="$here/build.sh"

echo "#!/bin/bash" >  $subscript

i=1;


for S in $inbase/*
do
	  SAMPLE=$(basename $S)
	  echo $SAMPLE

	  A1=$inbase/$SAMPLE/bam
	  
	  AWC="$A1/*.PAIRED.20101123.bam"
	  #echo "$A1/*.PAIRED.20101123.bam"
	  #echo $AWC
	  
	  
	  nfiles=$(ls $AWC 2> /dev/null | wc -l)
	  echo $nfiles
	  if [ "$nfiles" != "0" ]
	  then
			echo "find $A1"
			echo $nfiles
	  else
			echo "no $A1"
			A1=$inbase/$SAMPLE/not_recal
			AWC="$A1/*.PAIRED.20101123.bam"
			nfiles=$(ls $AWC  2> /dev/null | wc -l)
			if [ "$nfiles" != "0" ]
			then
				echo "find $A1"
				echo $nfiles
			else
				echo "no $A1"
			fi
	  fi

	  #echo $A1

	  builddir=$outbase/$SAMPLE

	  for file in $A1/*.PAIRED.20101123.bam
	  do

		if [[ "$file" == *SINGLE* ]] 
		then
			echo "skip " $file
			continue
		fi

		i=$(($i+1))

		if [ -f $file ]
		then
			f1=$(basename $file)
			echo "find " $f1
		else
			echo "no " $file	
			continue
		fi

		# mkdir -p $builddir

		f2=${f1%.PAIRED.20101123.bam*}

		bamfile=$A1/${f2}.PAIRED.20101123.bam
		fragfile=$A1/${f2}.PAIRED.20101123.stat
		specialbamfile=$A1/${f2}.PAIRED.20101123.special.bam

		inarg="-i $bamfile "
		outarg="-o $builddir"
	 
		build="$program -b $args -i $bamfile -o $builddir -f $fragfile"
		specialbuild="$program -b $args -i $specialbamfile -o $builddir -f $fragfile"

		script="$here/jobs/Spanner_Build_${SAMPLE}_$f2.sh"

		echo "#PBS -V" > $script
		echo "#PBS -l walltime=23:00:00" >> $script
		echo "#PBS -j oe " >> $script
		echo "#PBS -o log/SpannerBuild_${SAMPLE}_$f2.log" >> $script   
		echo "#PBS -l nodes=1:ppn=2" >> $script
		echo "echo \$(date)" >> $script 
		echo "mkdir -p $builddir" >> $script
		echo "cd  $here" >> $script
		echo "pwd" >> $script
		echo "	" >> $script   
		echo "ls -latr $file" >> $script
		echo "echo \$(date)" >> $script
		echo $build >> $script
		echo "	" >> $script
		echo "ls -latr $specialbamfile" >> $script		   
		echo "echo \$(date)" >> $script
		echo $specialbuild >> $script
		echo "	" >> $script
		echo "echo \$(date)" >> $script
		echo "qsub -q low -N ${SAMPLE}_${f2}_build $script" >>	$subscript
		echo "sleep 0.01" >> $subscript
		#break
	  done
	  #break
 done

