#!/bin/bash

datadir=/data_archive/
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir=$1
fi

template="2021*pulsar*"
if [[ -n "$2" && "$2" != "-" ]]; then
   template=$2
fi

objects="J* B*"
if [[ -n "$3" && "$3" != "-" ]]; then
   objects="$3"
fi

options=""
if [[ -n "$4" && "$4" != "-" ]]; then
   options=$4
fi

fits2fil=0
if [[ -n "$5" && "$5" != "-" ]]; then
   fits2fil=$5
fi

transposed=0
if [[ -n "$6" && "$6" != "-" ]]; then
   transposed=$6
fi

n_avg=7
if [[ -n "$7" && "$7" != "-" ]]; then
   n_avg=$7
fi


export PATH=$HOME/github/hdf5_correlator/scripts/:$PATH

done_file=processed.txt

echo "############################################"
echo "PARAMETERS:"
echo "############################################"
echo "datadir       = $datadir"
echo "template      = $template"
echo "objects       = $objects"
echo "options       = $options"
echo "fits2fil      = $fits2fil"
echo "transposed    = $transposed"
echo "n_avg         = $n_avg"
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
            cd ${subdir}/ # was also ${freq_start}/
            pwd
            ls -al
         
            echo "INFO : processing ${subdir} / object ${object} ..."
            for channel in `ls -d ??? ?? 2>/dev/null`
            do
               cd $channel
               pwd
               echo "Check input files:"

               # added now stitching :
               ls ch???/dynspec_avg${n_avg}_i.fits > fits_list
               ls ch???/dynspec_avg${n_avg}_q.fits > fits_list_q
               ls ch???/dynspec_avg${n_avg}_u.fits > fits_list_u
               ls ch???/dynspec_avg${n_avg}_v.fits > fits_list_v
               ls ch???/dynspec_avg${n_avg}_x.fits > fits_list_x
               ls ch???/dynspec_avg${n_avg}_y.fits > fits_list_y
               
               progname=skalow_stitch_data
               if [[ $transposed -gt 0 ]]; then
                  progname=skalow_stitch_data_transposed
               fi

               echo "${progname} fits_list all_channels_ch${channel}_i.fits"
               ${progname} fits_list all_channels_ch${channel}_i.fits

               echo "${progname} fits_list_q all_channels_ch${channel}_q.fits"
               ${progname} fits_list_q all_channels_ch${channel}_q.fits

               echo "${progname} fits_list_u all_channels_ch${channel}_u.fits"
               ${progname} fits_list_u all_channels_ch${channel}_u.fits

               echo "${progname} fits_list_v all_channels_ch${channel}_v.fits"
               ${progname} fits_list_v all_channels_ch${channel}_v.fits

               echo "${progname} fits_list_x all_channels_ch${channel}_x.fits"
               ${progname} fits_list_x all_channels_ch${channel}_x.fits

               echo "${progname} fits_list_y all_channels_ch${channel}_y.fits"
               ${progname} fits_list_y all_channels_ch${channel}_y.fits

               if [[ $transposed -gt 0 ]]; then
                  echo "WARNING : total power not yet implemented for transposed version"
               else
                  echo "total_power all_channels_ch${channel}_i.fits all_channels_ch${channel}_i_total_power.txt -c 12"
                  total_power all_channels_ch${channel}_i.fits all_channels_ch${channel}_i_total_power.txt -c 12

                  echo "total_power all_channels_ch${channel}_x.fits all_channels_ch${channel}_x_total_power.txt -c 12"
                  total_power all_channels_ch${channel}_x.fits all_channels_ch${channel}_x_total_power.txt -c 12

                  echo "total_power all_channels_ch${channel}_y.fits all_channels_ch${channel}_y_total_power.txt -c 12"
                  total_power all_channels_ch${channel}_y.fits all_channels_ch${channel}_y_total_power.txt -c 12
               fi

               if [[ $fits2fil -gt 0 ]]; then
                  echo "fits2fil2 all_channels_ch${channel}_i.fits  all_channels_ch${channel}_i.fil -u"
                  fits2fil2 all_channels_ch${channel}_i.fits  all_channels_ch${channel}_i.fil -u
               else
                  echo "WARNING : fits2fil conversion is not required"
               fi

               cd ..
            done
         
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

