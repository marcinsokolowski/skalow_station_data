#!/bin/bash

datadir=/data_archive/
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir=$1
fi

template="2021*pulsar*"
if [[ -n "$2" && "$2" != "-" ]]; then
   template=$2
fi

freq_start=410
if [[ -n "$3" && "$3" != "-" ]]; then
   freq_start=$3
fi

objects="J* B*"
if [[ -n "$4" && "$4" != "-" ]]; then
   objects="$4"
fi

n_cc_channels=1
if [[ -n "$5" && "$5" != "-" ]]; then
   n_cc_channels=$5
fi

force=1
if [[ -n "$6" && "$6" != "-" ]]; then
   force=$6
fi

conversion_options=""
if [[ -n "$7" && "$7" != "-" ]]; then
  conversion_options="$7"   
fi

n_channels=128
if [[ -n "$8" && "$8" != "-" ]]; then
   n_channels=$8
fi

period=0.0894189988
if [[ -n "$9" && "$9" != "-" ]]; then
   period=$9
fi

skip_bytes=4096
ext="dada"
if [[ -n "${10}" && "${10}" != "-" ]]; then
   ext=${10}
fi
if [[ ${ext} != "dada" ]]; then
   skip_bytes=0
fi

file_template="*"
if [[ -n "${11}" && "${11}" != "-" ]]; then
   file_template="${11}"
fi

options=""
if [[ -n "${12}" && "${12}" != "-" ]]; then
   options=${12}
fi

channel_to_process=0
if [[ -n "${13}" && "${13}" != "-" ]]; then
   channel_to_process=${13}
fi
freq_channel_processed=$(($freq_start+$channel_to_process))

outdir="ch${freq_channel_processed}"
if [[ -n "${14}" && "${14}" != "-" ]]; then
   outdir="${14}"
fi
mkdir -p ${outdir}

transposed=0
if [[ -n "${15}" && "${15}" != "-" ]]; then
   transposed=${15}
fi


export PATH=$HOME/github/hdf5_correlator/scripts/:$PATH

done_file=processed.txt


echo "############################################"
echo "Script version 2022-11-28"
echo "PARAMETERS:"
echo "############################################"
echo "datadir       = $datadir"
echo "template      = $template"
echo "file template = $file_template"
echo "freq_start    = $freq_start -> freq_channel_processed = $freq_channel_processed"
echo "objects       = $objects"
echo "n_cc_channels = $n_cc_channels"
echo "channel_to_process = $channel_to_process"
echo "n_channels    = $n_channels"
echo "force         = $force"
echo "period        = $period"
echo "file format   = $ext (skip $skip_bytes [bytes])"
echo "options       = $options"
echo "outdir        = $outdir"
echo "transposed    = $transposed"
echo "############################################"
date

cd $datadir
pwd

for dir in `ls -d ${template}`
do
   echo
   echo "Processing $dir"
   pwd
   cd ${dir}
   pwd
   
#   if [[ -s ${done_file} ]]; then
#      echo "\t$dir already processed"
#   else
   if [[ 1 -gt 0 ]]; then
      pwd
      echo "ls -d ${objects}"
      ls -d ${objects}
      
      for subdir in `ls -d ${objects} 2>/dev/null`
      do
         object=`echo ${subdir} | cut -b 1-10`
         echo
         echo "INFO : processing subdirectory $subdir -> object = $object"
         
         if [[ -d ${subdir} ]]; then
            echo "cd ${subdir}/"
            cd ${subdir}/ # was also ${freq_start}/
            pwd
            ls -al
         
            if [[ -s ${done_file} && $force -le 0 ]]; then
               echo "${subdir} / object ${object} already processed (remove file ${done_file} to repeat processing)"
            else         
               echo "INFO : processing ${subdir} / object ${object} ..."
               for channel in `ls -d ??? ?? 2>/dev/null`
               do
                  cd $channel
                  pwd
                  echo "Check input files:"
                  ls ${file_template}.${ext}
                  for dada_file in `ls ${file_template}.${ext} 2>/dev/null`
                  do
                     # channel_1_1_1659060532.876270.dada
                     ch=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ch;}'`
                     if [[ $ch -gt 0 ]]; then
                        # for the case of processing single BIG file use channel from file name
                        # otherwise use parameter channel_to_process to calculate the index of the channel to 
                        echo "DEBUG : ch = $ch > 0 from the file name -> processing this channel in a big file"
                        channel_to_process=$ch
                     else
                        # channel = 0 obtained from the .dat/.dada file name should not overwrite the parameter because it most likely means that either it is the first processed channel 
                        # or we are processing single files per frequency channel -> the parameter passed from the script above will have the correct value of the processed channel (see the loop in skalow_station_fold_allcc.sh)
                        # $ch <= 0 -> probably =0 which means we need to keep channel_to_process as obtained from the parameters                         
                        echo "DEBUG : ch = $ch -> keeping channel_to_process = $channel_to_process"
                     fi
                     channel_total=`echo "$channel $channel_to_process" | awk '{printf("%d\n",($1+$2));}'`
                     freq_mhz=`echo "$channel $channel_to_process" | awk '{printf("%.6f\n",($1+$2)*(400.00/512.00));}'`
                     ux=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ux;}'`
                     utc=`date -u -d "1970-01-01 UTC $ux seconds" +"%Y%m%dT%H%M%S"`
                     outfile=${utc}_ch${channel_total}
                     outdir=ch${channel_total}
                     
                     echo ".dada file = $dada_file"
                     echo "ch = $channel + $channel_to_process = $channel_total -> freq = $freq_mhz [MHz]"
                     echo "ux = $ux -> utc  = $utc"
                     echo "outfile = $outfile"
                     echo "outdir  = $outdir"
                  
                     processed_file=${dada_file%%dada}processed
            
                     if [[ -s $processed_file && $force -le 0 ]]; then
                        echo "File $dada_file already processed, in order to re-process remove file $processed_file"
                     else
                         # skalow_spectrometer test.dada -f test -p 0 -C 1 -c 0 -s 4096 -Z  -m -1 -F 410 -N 128 -O dynspec -a 7 -P 0.0894189988 -D 2
                         mkdir -p ${outdir}
                         echo "skalow_spectrometer $dada_file -f test -p 0 -C ${n_cc_channels} -c ${channel_to_process} -s ${skip_bytes} -Z  -m -1 -F ${channel} -N ${n_channels} -O dynspec -a 7 -P ${period} -D 2 -A ${outdir} ${options}"
                         skalow_spectrometer $dada_file -f test -p 0 -C ${n_cc_channels} -c ${channel_to_process} -s ${skip_bytes} -Z  -m -1 -F ${channel} -N ${n_channels} -O dynspec -a 7 -P ${period} -D 2 -A ${outdir} ${options}
                         
                         # calculate mean vs. freq. :
                         cd ${outdir}
                         # calculate stokes I and Q :
                         file_x=`ls dynspec_*_x.fits | tail -1`
                         file_y=${file_x%%_x.fits}_y.fits
                         xy_phase=${file_x%%_x.fits}_xyphase.fits
                         stokes_i=${file_x%%_x.fits}_i.fits
                         stokes_q=${file_x%%_x.fits}_q.fits
                         stokes_u=${file_x%%_x.fits}_u.fits
                         stokes_v=${file_x%%_x.fits}_v.fits
                         echo "calcfits_bg ${file_x} + ${file_y} ${stokes_i}"
                         calcfits_bg ${file_x} + ${file_y} ${stokes_i}

                         echo "calcfits_bg ${file_x} - ${file_y} ${stokes_q}"
                         calcfits_bg ${file_x} - ${file_y} ${stokes_q}
                         
                         if [[ $transposed -le 0 ]]; then
                            echo "calcfits_bg ${xy_phase} s"
                            calcfits_bg ${xy_phase} s
                         
                            echo "calcfits_bg ${stokes_i} s"
                            calcfits_bg ${stokes_i} s
                         
                            echo "calcfits_bg ${stokes_q} s"
                            calcfits_bg ${stokes_q} s
                         
                            echo "calcfits_bg ${stokes_u} s"
                            calcfits_bg ${stokes_u} s
                         
                            echo "calcfits_bg ${stokes_v} s"
                            calcfits_bg ${stokes_v} s
                         else
                            echo "WARNING : calculating mean spectrum for transposed dynamic spectrum not implemented yet"
                         fi   
                         cd ..
                     fi
                  done
                  cd ..
               done
            fi
         
            echo "date > ${done_file}"
            date > ${done_file}         
            cd ../ # was also ..
         else
            echo "WARNING : subdirectory $subdir does not exist -> skipped"
         fi
      done
   else
      # echo "process_all_objects.sh \"${objects}\" ${conjugate}"
      # process_all_objects.sh "${objects}" ${conjugate}   
      echo "ERROR : not expected option -> needs to be fixed"
   fi
      
   echo "date > ${done_file}"
   date > ${done_file}
#   fi
   cd ..
done

