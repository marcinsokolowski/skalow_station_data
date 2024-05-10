 /*
 * SigprocFile.cpp
 *
 *  Created on: 27 Oct 2016
 *      Author: ban115
 */

#include "SigprocFile.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <sys/types.h>
#include <stdint.h>
#include <math.h>
#include "InvalidSourceFormat.h"
#include <fcntl.h> // posix_fadvise

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <libnova_interface.h>


int gDebugLevel = 0;

/* Same as strstr but goes through *all* the string - even if it contains nulls
 *
 */
char* mystrnstr(const char* s1, const char* s2, size_t n)
{
	size_t s2_len = strlen(s2);
	for (size_t i = 0; i < n - s2_len; i++) {
		if(strncmp(s1+i, s2, s2_len) == 0) {
			return (char*) ( s1 + i);
		}
	}

	return NULL;
}

SigprocFile::SigprocFile()
: m_file(NULL), m_filename(NULL), m_samples_read(0), m_fd(0), m_nifs(0), m_nchans(0), m_nbits(0), m_src_raj(0.00), m_src_dej(0.00) 
{
}

static char* find_string_in_header( char* header, int header_len, const char* template_str )
{
   int start = 0;
   int template_len = strlen(template_str);
   
   while( start < (header_len-template_len) ){
      int match=1;
      for(int i=0;i<template_len;i++){
         if( header[start+i] != template_str[i] ){            
            if( gDebugLevel > 0 ) {printf("start = %d stopped at %c != %c\n",start,header[start+i],template_str[i]);}
            match=0;
            break;
         }
      }
      
      if( match ){
         if( gDebugLevel > 0 ) { printf("FOUND !!!\n"); }
         return (header+start);
      }
      
      start++;
   }
   
   if( gDebugLevel > 0 ) { printf("NOT FOUND !!!\n"); }

   return NULL;
}

char SigprocFile::get_keyword_type( const char* keyword )
{
// f - flag
// i - integer
// s - string
// d - double
// u - unknown 

    if( strcmp(keyword,"HEADER_START") == 0 ){ return  'f'; }
    if( strcmp(keyword,"telescope_id") == 0 ){ return  'i'; }
    if( strcmp(keyword,"machine_id") == 0 ){ return  'i'; }
    if( strcmp(keyword,"data_type") == 0 ){ return  'i'; }
    if( strcmp(keyword,"rawdatafile") == 0 ){ return  's'; }
    if( strcmp(keyword,"source_name") == 0 ){ return  's'; }
    if( strcmp(keyword,"barycentric") == 0 ){ return  'i'; }
    if( strcmp(keyword,"pulsarcentric") == 0 ){ return  'i'; }
    if( strcmp(keyword,"az_start") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"za_start") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"src_raj") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"src_dej") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"tstart") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"tstart_utc") == 0 ){ return  's'; } 
    if( strcmp(keyword,"tsamp") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"nbits") == 0 ){ return  'i'; }
    if( strcmp(keyword,"signed") == 0 ){ return  'b'; }
    if( strcmp(keyword,"nsamples") == 0 ){ return  'i'; }
    if( strcmp(keyword,"nbeams") == 0 ){ return  'i'; }
    if( strcmp(keyword,"ibeam") == 0 ){ return  'i'; }
    if( strcmp(keyword,"fch1") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"foff") == 0 ){ return  'd'; }
    if( strcmp(keyword,"FREQUENCY_START") == 0 ){ return  'f'; }
    if( strcmp(keyword,"fchannel") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"FREQUENCY_END") == 0 ){ return  'f'; }
    if( strcmp(keyword,"nchans") == 0 ){ return  'i'; }
    if( strcmp(keyword,"nifs") == 0 ){ return  'i'; }
    if( strcmp(keyword,"refdm") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"period") == 0 ){ return  'd'; } 
    if( strcmp(keyword,"npuls") == 0 ){ return  'q'; }
    if( strcmp(keyword,"nbins") == 0 ){ return  'i'; }
    if( strcmp(keyword,"HEADER_END") == 0 ){ return  'f'; }

    printf("ERROR : unknown type of keyword = %s\n",keyword);

    return 'u';
}


int SigprocFile::SetHeaderValue( const char* keyword, double value )
{
    return SetHeaderValue( m_hdr, m_hdr_nbytes, keyword, value );
}

int SigprocFile::SetHeaderValue( const char* keyword, int value )
{
    return SetHeaderValue( m_hdr, m_hdr_nbytes, keyword, value, value );
}

int SigprocFile::SetHeaderValue( char* pHeader, int header_len, const char* keyword, double value, int value_int, const char* value_str, int start_index )
{
   char* p_foff_ptr = find_string_in_header( pHeader, header_len, keyword  );
   int len = strlen(keyword); // was just hardcoded 4 
   char keyword_type = get_keyword_type(  keyword );
   char szValue[16];
   memset(szValue,'\0',16);

   
   if( p_foff_ptr ){
      printf("keyword %s found in header OK\n",keyword);           
      char* p_value = (p_foff_ptr + len);
      
      double foff_value, minus_foff;
      int foff_value_int, minus_foff_int;
      
      switch( keyword_type ){
         case 'd' :
           memcpy(szValue,p_value,sizeof(double));
         
           foff_value = *((double*)(szValue));
           printf("CHECK : current header value of %s  = %.8f, changing to %.8f\n",keyword,foff_value,value);                   
            
           minus_foff = value;
           memcpy(p_value, &minus_foff, sizeof(double) );
           break;
           
         case 'i' :
           memcpy(szValue,p_value,sizeof(int));
         
           foff_value_int = *((int*)(szValue));
           printf("CHECK : current header value of %s  = %d, changing to %d\n",keyword,foff_value_int,value_int);                   
            
           minus_foff_int = value_int;
           memcpy(p_value, &minus_foff_int, sizeof(int) );
           
           break;
           
         default :
            printf("ERROR : unsupported type of keyword %s\n",keyword);
            break;
           
      }
      
      return 1;
   }else{
      if( start_index >= 0 ){
         // write length of keyword and keyword :
         memcpy( pHeader + start_index, &len, sizeof(len) );
         start_index += sizeof(len);
         
         memcpy( pHeader + start_index, keyword, len );
         start_index += len;
         int len_string = 0;

         switch( keyword_type ){
            case 'd' :
               memcpy(pHeader + start_index, &value, sizeof(double));         
               start_index += sizeof(double);
               break;
           
             case 'i' :
               memcpy( pHeader + start_index, &value_int, sizeof(int) );
               start_index += sizeof(int);         
               break;
               
             case 's' :
               len_string = strlen(value_str);
               memcpy( pHeader + start_index , &len_string, sizeof(int) );
               start_index += sizeof(int);

             
	       memcpy(pHeader + start_index, value_str, len_string);
               start_index += len_string;
               break;

             case 'f' :
//	       memcpy(pHeader + start_index, keyword, len);
//               start_index += len;
               // already written above
               break;

           
             default :
                printf("ERROR : unsupported type of keyword %s\n",keyword);
                break;
           
          }         
      }else{
         printf("ERROR : could not find keyword %s in header\n",keyword);
      }
   }
   
   return start_index;
}

int SigprocFile::FillHeader( bool recalc_tstart, bool fill_radec )
{
   bool b_tstart_utc = false;
   
   m_hdr_nbytes = MAX_HDR_SIZE;
   memset( m_hdr, '\0', m_hdr_nbytes );
   

   double mjd = m_tstart;
   if( recalc_tstart ){
      mjd = ux2mjd( int(m_tstart), (m_tstart-int(m_tstart))*1000000.00 ); 
//   void get_ymd_hms_ut( time_t ut_time, int& year, int& month, int& day,
//                  int& hour, int& minute, double& sec )
   }

   // Only for test/presentation reasons :
   int year, month, day, hour, minute;
   double sec;
   get_ymd_hms_ut( m_tstart, year, month, day, hour, minute, sec );
   sec += (m_tstart-int(m_tstart)); // add fractional part 
   char szUTC[128]; // 2019-07-18T14:53:13.920
   sprintf(szUTC,"%04d-%02d-%02dT%02d:%02d:%.3f",year, month, day, hour, minute, sec );   
   printf("DEBUG : CFilFile::WriteHeader ux = %.6f -> mjd = %.8f and %s\n",m_tstart,mjd,szUTC);
   
   int idx=0;
   // int SigprocFile::SetHeaderValue( char* pHeader, int header_len, const char* keyword, double value, int value_int, const char* value_str, int start_index )
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "HEADER_START", 0.00, 0 , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "nchans", 0.00      , m_nchans, NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "nifs"  , 0.00      , m_nifs  , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "fch1"  , m_fch1    , 0       , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "foff"  , m_foff    , 0       , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "tstart", mjd       , 0       , NULL, idx );
   if( b_tstart_utc ){
      idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "tstart_utc", 0.00  , 0       , szUTC, idx );
   }
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "tsamp" , m_tsamp   , 0       , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "nbits" , 0.00      , m_nbits , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "src_raj" , m_src_raj   , 0     , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "src_dej" , m_src_dej  , 0     , NULL, idx );
   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "telescope_id" , 0.00, m_telescope_id, NULL, idx ); 
   printf("DEBUG : added src_raj, src_decj and telescope_id to the header [%.8f,%.8f,%d]\n",m_src_raj,m_src_dej,m_telescope_id);

   idx = SetHeaderValue( m_hdr, m_hdr_nbytes, "HEADER_END" , 0.00 , 0       , NULL, idx);
   m_hdr_nbytes = idx;
   
   printf("DEBUG : set nbits := %d\n",m_nbits);
      
   return idx;         
}

int SigprocFile::Write( const char* filename, SigprocFile& right, int bFlipFreq, bool bSetBits )
{
   if(!m_file){
      m_file = fopen( filename , "wb" );
      m_fd = fileno(m_file);

      int written = 0;
      char* pHeader = new char[right.m_hdr_nbytes];
      memcpy(pHeader,right.m_hdr,right.m_hdr_nbytes);
      
//      if( bSetBits ){
//         SetHeaderValue( pHeader, right.m_hdr_nbytes, "nbits" , m_nbits );
//         printf("nbits set to %d\n",m_nbits);
//      }

      
      if( bFlipFreq ){
         double foff = right.foff();
         double fch1 = right.fch1();
         
         printf("ORIGNAL foff = %.4f , fch1 = %.4f\n",right.foff(),right.fch1());
         if( foff > 0 ){
            // means it was >0 and now we are changing it to <0 :
            // WAS foff > 0 -> fch1 lower end ->
            // NOW fch1 is upper end :
            fch1 =  fch1 + right.nchans()*fabs(foff);
            foff = -foff;
         }else{
            // foff < 0 
            // means it is now going to be >0  and was <0 which means fch1 - was the upper end of frequency 
            // WAS foff < 0 -> fch1 upper end ->
            // NOW fch1 is lower end :
            fch1 =  fch1 - right.nchans()*fabs(foff);
            foff = -foff;
         }         
         printf("CHANGED TO : fch1 = %.4f , foff = %.4f\n",fch1,foff);
         
         SetHeaderValue( pHeader, right.m_hdr_nbytes, "foff", fabs(foff) );
         SetHeaderValue( pHeader, right.m_hdr_nbytes, "fch1", fch1 ); // was 140.1560
         printf("foff set to %.4f\n",foff);         
         
         written = fwrite( pHeader, sizeof(char), right.m_hdr_nbytes, m_file );
         printf("Written %d bytes to flipped header\n",written);         
      }else{    
         written = fwrite( pHeader, sizeof(char), right.m_hdr_nbytes, m_file );
      }
      
      delete [] pHeader;
      
      if( written != right.m_hdr_nbytes ){
         printf("ERROR : could not write header written %d bytes, should be %d bytes !!!\n",written,(int)m_hdr_nbytes);
         fclose(m_file);
         m_file = NULL;
         m_fd   = 0;
         return -m_hdr_nbytes;
      }
   }else{
      printf("File already opened -> nothing written !\n");
   }
   
   return 0;
}

int SigprocFile::CopyFilFile( const char* filename, SigprocFile& infile )
{
   float* data = new float[10240472];
   int n_read = fread(data, sizeof(uint8_t), 10240472, infile.m_file);
   printf("Read %d bytes\n",n_read);

   // close if previously opened :
   name( filename );
   if( m_file ){ 
      fclose(m_file); 
   }
   
   // create a new file :
   m_file = fopen( m_filename, "wb" );
   m_fd = fileno(m_file);
   int written = fwrite( infile.m_hdr, sizeof(char), infile.m_hdr_nbytes, m_file );
   
   int total_written = 0;
   total_written += WriteData( data, n_read/sizeof(float) );

   // read and write the remaining data - usually the files can be longer than just 10240472 bytes ( = 10.24 MB )   
   while( (n_read = fread(data, sizeof(uint8_t), 10240472, infile.m_file)) >0 ){
       total_written += WriteData( data, n_read/sizeof(float) );
   }
   
   delete [] data;
   
   printf("SigprocFile::CopyFilFile( SigprocFile& infile ) : written %d bytes\n",total_written);
   
   return 0;

}

int SigprocFile::Write( const char* filename )
{
   printf("ERROR : this function should not be used as this is wrong !!! it uses m_file for both read and write tile\n");
   printf("ERROR : this means that the entire file cannot be read/written !!!\n");
   printf("ERROR : instead use Write( SigprocFile& infile )\n");
   exit(-1);

   float* data = new float[10240472];
   int n_read = fread(data, sizeof(uint8_t), 10240472, m_file);
   printf("Read %d bytes\n",n_read);

   m_file = fopen( filename , "wb" );
   m_fd = fileno(m_file);
   int written = fwrite( m_hdr, sizeof(char), m_hdr_nbytes, m_file );
   WriteData( data, n_read/sizeof(float) );

   // read and write the remaining data - usually the files can be longer than just 10240472 bytes ( = 10.24 MB )   
   while( (n_read = fread(data, sizeof(uint8_t), 10240472, m_file)) >0 ){
       WriteData( data, n_read/sizeof(float) );
   }
   
   delete [] data;
   
   return 0;
}

int SigprocFile::WriteData( float* buffer, int n_channels )
{
    if( m_file ){
        int written = fwrite( buffer, sizeof(float), n_channels, m_file );
        if( written != n_channels ){
           printf("ERROR : could not write header written %d bytes, should be %d bytes !!!\n",written,n_channels);
           fclose(m_file);
           m_file = NULL;
           m_fd   = 0;
           return -n_channels;
        }
        
        return written*sizeof(float);
    }
    
    return 0;
}

int SigprocFile::WriteData( unsigned char* buffer, int n_channels )
{
    if( m_file ){
        int written = fwrite( buffer, sizeof(unsigned char), n_channels, m_file );
        if( written != n_channels ){
           printf("ERROR : could not write header written %d bytes, should be %d bytes !!!\n",written,n_channels);
           fclose(m_file);
           m_file = NULL;
           m_fd   = 0;
           return -n_channels;
        }
        
        return written;
    }
    
    return 0;
}

bool SigprocFile::DoesFileExist(const char* fname)
{
        bool bRet=false;
        if(fname && fname[0] ){
                // mystring szFile=fname;
                // szFile.env2str();

                if( access( fname, F_OK ) == 0 ) {
                        bRet = true;
                }
        }
        return bRet;
}


int SigprocFile::GetFileSize( const char* filename )
{
        if( !DoesFileExist( filename ) ){
                return -1;
        }

        struct stat buf;
        stat( filename, &buf );
        return (int)buf.st_size;
}


int  SigprocFile::MergeCoarseChannels( std::vector<string>& fil_file_list, const char* out_file , double*& avg_spectrum, int foff_sign )
{
   int max_filfiles = fil_file_list.size();
   SigprocFile* infiles[fil_file_list.size()]; // maximum 24 
   
   if( fil_file_list.size() > max_filfiles || fil_file_list.size() <= 0 ){
      printf("ERROR : maximum number of fil files for SigprocFile::MergeCoarseChannels exceeded (%d given , %d is limit , no files is not allowed either)\n",int(fil_file_list.size()),max_filfiles);
      return -1;
   }

   // read zero file into output buffer :   
   SigprocFile outfil_file( fil_file_list[0].c_str() );
   // SigprocFile( int nbits, int nifs, int nchans, double fch1, double foff, double tstart, double tsamp );
   // SigprocFile filfile_norm( filfile.nbits(), filfile.nifs(), filfile.nchans(), filfile.fch1(), filfile.foff(), filfile.tstart(), filfile.tsamp() );
   outfil_file.nbits( 32 ); // always output file is 32 bits float (no matter input)
   outfil_file.Close();
   
   int n_channels_file0 = outfil_file.nchans();
   int n_out_channels = 0;
   float foff = 0.00;
   for(int i=0;i<fil_file_list.size();i++){
      infiles[i] = new SigprocFile( fil_file_list[i].c_str() );
      foff = (infiles[i])->foff();
 
      n_out_channels += (infiles[i])->nchans();     
      if( (infiles[i])->nchans() != n_channels_file0 ){
         printf("WARNING : fil file %s has %d (first one has %d)\n",fil_file_list[i].c_str(),(infiles[i])->nchans() ,n_channels_file0);
      }
      printf("DEBUG : added %d channels -> total %d channels\n",(infiles[i])->nchans(),n_out_channels);
   }
   printf("DEBUG : number of output channels = %d\n",n_out_channels);

   if( !avg_spectrum ){
      avg_spectrum = new double[n_out_channels];
      printf("INFO : allocated buffer for average spectrum of %d channels\n",n_out_channels);
   }   
   if( avg_spectrum ){
      for(int j=0;j<n_out_channels;j++){
         avg_spectrum[j] = 0;
      }
   }

   
   outfil_file.SetHeaderValue( "nchans" , n_out_channels );
   // outfil_file.SetHeaderValue( "foff", foff*foff_sign );
   outfil_file.foff( foff*foff_sign );
   outfil_file.FillHeader( false, true );
   outfil_file.SetHeaderValue( "nchans" , n_out_channels );
   outfil_file.WriteHeader( out_file , false /* do not close */ , true /* new file -> set m_file := out_f */ );
   float* out_spectrum = new float[n_out_channels];
   unsigned char*  char_buffer = new unsigned char[n_out_channels];
   bool all_ok = true;
   int n_out_spectra = 0;
   int sizeof_sample = sizeof(float);
   int out_bytes_count = n_out_channels*sizeof(float);
   
   while( all_ok ){
       float* p_out_ptr = out_spectrum;
       memset( out_spectrum, '\0', out_bytes_count );
       int n_read_total = 0;
       
       // read 1 spectrum from each file and put into out_spectrum :
       for(int i=0;i<fil_file_list.size();i++){
           // int read( void* buffer, int size );
           int in_channels = (infiles[i])->nchans();
           int bytes_to_read = in_channels*sizeof(float);
           int n_read = 0;
           
           if( (infiles[i])->nbits() >= 32 ){ 
              n_read = (infiles[i])->read( p_out_ptr , bytes_to_read );
           }else{
              // the case of merging 1-byte .fil files:
              if( (infiles[i])->nbits() == 8 ){ 
                 sizeof_sample = sizeof(unsigned char);
                 bytes_to_read = in_channels*sizeof_sample;
                 n_read = (infiles[i])->read( char_buffer, bytes_to_read );
                 for(int ch=0;ch<in_channels;ch++){
                    p_out_ptr[ch] = char_buffer[ch];
                 }
              }
           }
           if( n_read != bytes_to_read ){
              printf("END OF FILE : read %d bytes expected %d bytes -> cannot continue\n",n_read,bytes_to_read);
              all_ok = false;
              break;
           }
           if( n_out_spectra <= 0 ){
              printf("INFO : auto-detected sizeof sample %d bits\n",(infiles[i])->nbits()); 
           }
           p_out_ptr += in_channels;
           n_read_total += n_read;
       }           
       
       if( avg_spectrum ){
          for(int j=0;j<n_out_channels;j++){
             avg_spectrum[j] += out_spectrum[j];                          
          }
          
          if( n_out_spectra == 0 ){
             printf("DEBUG spectra :\n");
             for(int j=0;j<n_out_channels;j++){
                printf("\t%d %.8f\n",j,avg_spectrum[j]);
             }
          }
       }
       
       if( n_read_total == n_out_channels*sizeof_sample ){
       // int SigprocFile::WriteData( float* buffer, int n_channels )
          outfil_file.WriteData( out_spectrum , n_out_channels );
          n_out_spectra++;
       }else{
          printf("END OF file : number bytes to be written = %d vs. %d actually read\n",out_bytes_count,n_read_total);
          all_ok = false;
          printf("Written %d spectra in total\n",n_out_spectra);
       }
   }
   delete [] out_spectrum;
   delete [] char_buffer;
    
   outfil_file.Close();
   
   if( avg_spectrum ){
      for(int j=0;j<n_out_channels;j++){
         avg_spectrum[j] = avg_spectrum[j] / n_out_spectra;
      }
   }
   
   printf("------------------------------------------------\n");
   printf("STATISTICS :\n");
   printf("------------------------------------------------\n");
   printf("Spectra count = %d\n",n_out_spectra);
   printf("------------------------------------------------\n");

   return n_out_channels;
}

int SigprocFile::WriteAveragedChannels( const char* out_file, int n_avg_factor )
{
//   SigprocFile infile( in_file );
   int old_channels = nchans();
   int new_channels = old_channels / n_avg_factor;
   
   SetHeaderValue( "nchans" , new_channels );
   SetHeaderValue( "foff", foff()*n_avg_factor );
//   WriteHeader( out_file );   
   FILE* out_f = fopen(out_file,"w");
   int ret = fwrite( m_hdr, sizeof(char), m_hdr_nbytes, out_f );
   
   char* data_buffer = new char[m_data_nbytes];
   int read_bytes = fread(data_buffer, 1, m_data_nbytes, m_file);
   
   float* data_float = (float*)data_buffer;
   int n_floats = m_data_nbytes / int(sizeof(float));
   printf("Data size = %d bytes -> %d floats\n",int(m_data_nbytes),n_floats);
   int n_spectra = n_floats / nchans();
   
   float* out_spectrum = new float[new_channels]; 
   float* in_spectrum  = new float[old_channels];
   
   for(int i=0;i<n_spectra;i++){
      memcpy( in_spectrum , data_float + i*nchans(), sizeof(float)*nchans() );
      
      int out_idx = 0;
      int in_idx  = 0;
      int cnt=0;
      double sum = 0.00;
      while( in_idx < nchans() ){
         sum += in_spectrum[in_idx];
         cnt++;
         
         if( cnt == n_avg_factor ){
            out_spectrum[out_idx] = (sum/cnt);
            cnt = 0;
            sum = 0.00;
            out_idx++;            
         }
         in_idx++;         
      }
      
      if( out_idx == new_channels ){      
          fwrite( out_spectrum, sizeof(float), out_idx, out_f );      
      }else{
          printf("ERROR : number of channels after averaging = %d != expected %d\n",out_idx,new_channels);
      }
   }
 
   delete [] out_spectrum;   
   delete [] data_buffer;
   fclose( out_f );
   
   return 0;
}

int SigprocFile::WriteHeader( const char* filename , bool bClose, bool bNewFile )
{
   FILE* out_f = fopen(filename,"w");
   int ret = fwrite( m_hdr, sizeof(char), m_hdr_nbytes, out_f );
   
   if( bClose ){
      fclose(out_f);
   }
   
   if( bNewFile ){
      m_file = out_f;      
      m_fd = fileno(m_file);      
   }
   
   return ret;
}

/*void SigprocFile::rewind()
{
   if( m_file ){
      fclose( m_file );
      m_file = fopen( m_filename , "r" );
      printf("DEBUG : re-opened file %s\n",m_filename);

      if (! m_file) {
         printf("SigprocFile: could not open file: %s - %s\n",m_filename, strerror(errno));
	 throw InvalidSourceFormat();
      }
   }
}*/

SigprocFile::SigprocFile( int nbits, int nifs, int nchans, double fch1, double foff, double tstart, double tsamp )
: m_nbits(nbits), m_nifs(nifs), m_nchans(nchans), m_fch1(fch1), m_foff(foff), m_tstart(tstart), m_tsamp(tsamp), m_file(NULL), m_filename(NULL),m_hdr_nbytes(MAX_HDR_SIZE), m_src_raj(0.00), m_src_dej(0.00) 
{}

void SigprocFile::name( const char* filename )
{
   if( m_filename ){
      delete [] m_filename;
   }
   m_filename = new char[strlen(filename) + 1];
   m_filename = strcpy(m_filename, filename);
}

SigprocFile::SigprocFile(const char* filename, int header_size /* = -1 */) 
: m_file(NULL), m_filename(NULL), m_samples_read(0), m_fd(0), m_nifs(0), m_nchans(0), m_nbits(0), m_src_raj(0.00), m_src_dej(0.00) 
{
        m_hdr_nbytes = MAX_HDR_SIZE;
        m_file_nbytes = GetFileSize( filename );

	m_file = fopen(filename, "r");
	m_filename = new char[strlen(filename) + 1];
	m_filename = strcpy(m_filename, filename);
	m_samples_read = 0;

	if (! m_file) {
		printf("SigprocFile: could not open file: %s - %s\n",filename, strerror(errno));
		throw InvalidSourceFormat();
	}

        if( header_size < 0 ){
            // 472
            header_size = MAX_HDR_SIZE;
        }

	// find end of header
	size_t size = fread(m_hdr, sizeof(char), MAX_HDR_SIZE, m_file);
	char* hdr_end = mystrnstr(m_hdr, "HEADER_END", MAX_HDR_SIZE);
	if (hdr_end == NULL) {
		printf("SigprocFile: File %s does not contain HEADER_END\n", filename);
		throw InvalidSourceFormat();
	}
	// TODO: Check it starts with HEADER_START
	m_hdr_nbytes = (size_t)((hdr_end - m_hdr) + strlen("HEADER_END"));
	m_data_nbytes = m_file_nbytes - m_hdr_nbytes;
	assert(m_hdr_nbytes < MAX_HDR_SIZE);
	m_hdr[m_hdr_nbytes] = 0;
	seek_sample(0); // ,m_hdr_nbytes);
	printf("File %s has header %d bytes long. i.e. 0x%x\n", filename, int(m_hdr_nbytes), int(m_hdr_nbytes));
	m_nbits = header_int("nbits");
	m_nifs = header_int("nifs");
	m_nchans = header_int("nchans");
	m_fch1 = header_double("fch1");
	m_foff = header_double("foff");
	m_tstart = header_double("tstart");
	m_tsamp = header_double("tsamp");
	m_src_raj = header_double("src_raj");
	m_src_dej = header_double("src_dej");
	
	printf("DEBUG : m_tstart = %.8f (SigprocFile::SigprocFile)\n",m_tstart);

	// tell_linux we'll be reading sequentailly
	m_fd = fileno(m_file);

	if (posix_fadvise(m_fd, 0, 0, POSIX_FADV_SEQUENTIAL) != 0) {
	  perror("Could not set advice\n");
	  exit(EXIT_FAILURE);
	}

	// block size must be a multiple of 8 bits
	if (nbits()*nifs()*nchans() % 8  != 0) {
		printf("SigprocFile: File %s nbits=%d not supported\n", filename, nbits());
		throw InvalidSourceFormat();
	}
}

void SigprocFile::Close()
{
   if( m_file ){
      fclose(m_file);
      m_file = NULL;
   }
   
}

SigprocFile::~SigprocFile() {
    Close();
    if( m_filename ){
       delete[] m_filename;
    }        
}

const char* SigprocFile::header_find(const char* hname) const
{
	char* hstart = mystrnstr(m_hdr, hname, MAX_HDR_SIZE);
	if (hstart == NULL) {
		printf("Could not find header %s in file %s\n", hname, m_filename);
		exit(EXIT_FAILURE);
	}

	int hlen = *((int*)(hstart - sizeof(int)));
	if (hlen != strlen(hname)) {
		printf("Could not find header %s n file %s (but found a substring)\n", hname, m_filename);
		exit(EXIT_FAILURE);
	}

	char* dstart = (hstart + strlen(hname));
	return dstart;
}
double SigprocFile::header_double(const char* hname) const {

	double* vstart = (double*) header_find(hname);
	double value = *vstart;
	return value;
}

int SigprocFile::header_int(const char* hname) const {
	int* vstart = (int*) header_find(hname);
	int value = *vstart;
	return value;
}

size_t SigprocFile::seek_sample(size_t t ) // , size_t hdr_size )
{
	size_t boff = t*nifs()*nchans()*nbits()/8 + m_hdr_nbytes;
	printf("DEBUG : t = %ld , m_hdr_nbytes = %ld -> boff = %ld, file = %ld\n",t,m_hdr_nbytes,boff,(long int)m_file);
	if(fseek(m_file, boff, SEEK_SET) < 0) {
		printf("SigprocFile: Could not seek to offset of file %s\n. Error: %s", m_filename, strerror(errno));
		assert(0);
	}

	m_current_sample = t;
	return boff;
}

void SigprocFile::advise_block(off_t bytes_per_block)
{
  int nblocks = 16;
  off_t offset = ftell(m_file);
  if (posix_fadvise(m_fd, offset, bytes_per_block*nblocks, POSIX_FADV_WILLNEED) != 0) {
    perror("Couln't set advice for next block\n");
    exit(EXIT_FAILURE);
  }

  // tell linux we don't need the stuff we've read
  if (posix_fadvise(m_fd, 0, offset, POSIX_FADV_DONTNEED) != 0) {
    perror("Coulnt set advise for previous data\n");
    exit(EXIT_FAILURE);
  }
}

 // int bytes_read = fread( buffer, 1, file_size,  filfile.
int SigprocFile::read( void* buffer, int size )
{
   int bytes_read = fread( buffer, 1, size, m_file );
   return bytes_read;
}


size_t SigprocFile::read_samples_uint8(size_t nt, uint8_t* output, bool bDoAssert /*=true*/ )
{
	// RETURNS TBF ordering. WARNING: This will deeply confuse the output if nbeams != 1
	if ( bDoAssert ){
   	   assert(nifs() == 1); // Otherwise users will be confused. TODO: get sources to tell user what data order is
   	}
	size_t nelements = nt*m_nifs*m_nchans;
	size_t nreq_bytes = nelements*m_nbits/8;
	size_t nbytes_read = fread(output, sizeof(uint8_t), nreq_bytes, m_file);
	size_t output_nsamp = nbytes_read * 8 / m_nbits / m_nifs / m_nchans;
	m_samples_read += output_nsamp;
	m_current_sample += nt;
	advise_block(sizeof(uint8_t)*nreq_bytes);
	return output_nsamp;
}


size_t SigprocFile::read_samples(void** output)
{
	// Usually we use the buffer in the FileSet, so we'll just ignore this for now.
	assert(0);
	return 0;
}

double SigprocFile::last_sample_elapsed_seconds()
{
	double toff_sec = ((double) m_samples_read) * m_tsamp;
	return toff_sec;
}

double SigprocFile::last_sample_mjd()
{
	double mjd = m_tstart + last_sample_elapsed_seconds()/86400.0;
	return mjd;
}
