/*
 * SigprocFile.h
 *
 *  Created on: 27 Oct 2016
 *      Author: ban115
 */

#ifndef SIGPROCFILE_H_
#define SIGPROCFILE_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <string>
#include "DataSource.h"

// const size_t MAX_HDR_SIZE = 4096;
#define MAX_HDR_SIZE 4096

/**
 * Reads a sigproc file.
 * Sigproc order is TIF order -
 * T = time, I = IF, F= channel. F fastest
 *
 * @see http://sigproc.sourceforge.net/sigproc.pdf
 */
class SigprocFile : public DataSource {
public:
	SigprocFile(const char* filename, int header_size=-1);
	SigprocFile( int nbits, int nifs, int nchans, double fch1, double foff, double tstart, double tsamp );
	SigprocFile();
	virtual ~SigprocFile();
	void Close();
	const char* header_find(const char* hname) const;
	int header_int(const char* hname) const;
	double header_double(const char* hname) const;
	size_t read_samples_uint8(size_t nt, uint8_t* output, bool bDoAssert=true);
	size_t read_samples(void** output);
	double last_sample_elapsed_seconds();
	double last_sample_mjd();
	char get_keyword_type( const char* keyword );
	
	// int bytes_read = fread( buffer, 1, file_size,  filfile.
	int read( void* buffer, int size );
	
	// write :
	int Write( const char* filename, SigprocFile& right, int bFlipFreq=0, bool bSetBits=false );
	int Write( const char* filename );
	int CopyFilFile( const char* filename, SigprocFile& infile );
	int WriteHeader( const char* filename , bool bClose=true, bool bNewFile=false );
	int WriteData( float* buffer, int n_channels );
	int WriteData( unsigned char* buffer, int n_channels );
	int FillHeader( bool recalc_tstart=true, bool fill_radec=false );
	int SetHeaderValue( char* pHeader, int header_len, const char* keyword, double value, int value_int=0, const char* value_str=NULL, int start_index=-1 );
	int SetHeaderValue( const char* keyword, double value );
	int SetHeaderValue( const char* keyword, int value );
	
	int WriteAveragedChannels( const char* out_file, int n_avg_factor );
	
	// output is always nbits=32 (float) .fil file (no matter input - 1 byte char or 4 bytes float) 
	static int MergeCoarseChannels( std::vector<string>& fil_file_list, const char* out_file, double*& avg_spectrum, int foff_sign=1 );

	int nifs() {
		return m_nifs;
	}
	int nbits() {
		return m_nbits;
	}
	int nbits( int nbits ) {
		m_nbits = nbits;
		return m_nbits;
	}
	int nbeams() {
		return nifs();
	}
	int nchans() {
		return m_nchans;
	}
	int npols() {
		return 1;
	}
	int nants() {
		return 1;
	}
	size_t samples_read() {
		return m_samples_read;
	}
	size_t current_sample() {
		return m_current_sample;
	}
	double fch1() {
		return m_fch1;
	}
	double foff() {
		return m_foff;
	}
	double foff( double _foff) {
		m_foff = _foff;
		return m_foff;
	}
	double tstart() {
		return m_tstart;
	}
	double tsamp() {
		return m_tsamp;
	}
	DataOrder data_order() {
		return DataOrder_TFBP;
	}

	char* name() {
		return m_filename;
	}
	
	void name( const char* filename );
	
	char* gethdr() {
	        return m_hdr;
	}
	
	int gethdrbytes(){
	    return m_hdr_nbytes;
	}
	
	int getdatabytes(){
            return m_data_nbytes;   	   
	}	

	int getfilebytes(){
            return m_file_nbytes;   	   
	}	
	
	const char* sourcename()
	{
	   return m_sourcename.c_str();
	}

	void sourcename( const char* src )
	{
	   m_sourcename = src;
	}
	
	void src_raj( double _src_raj )
	{
	   m_src_raj = _src_raj;
	}
	
	double src_raj()
	{
	   return m_src_raj;
	}
	
	void src_dej( double _src_dej )
	{
	   m_src_dej = _src_dej;
	}
	
	double src_dej()
	{
	   return m_src_dej;
	}
	
	void telescope_id( int tele_id ){
	   m_telescope_id = tele_id;	   	   
	}
	
	int telescope_id()
	{
	   return m_telescope_id;
	}
	
	
	void rewind()
	{
	   if( m_file ){
	      printf("DEBUG : rewind\n");
	      ::rewind( m_file );
	      m_samples_read = 0;
	      seek_sample(0);
              // ::fseek( m_file, 0L, SEEK_SET);
	   }
	}

	void advise_block(off_t nt);

	size_t seek_sample( size_t t );

// static functions :
        static bool DoesFileExist(const char* fname);
        static int GetFileSize( const char* filename );        

private:
	double m_fch1;
	double m_foff;
	double m_tstart; // unixtime
	double m_tsamp;
	int m_nifs;
	int m_nchans;
	int m_nbits;
	std::string m_sourcename;
	int m_telescope_id;
	double m_src_raj;
	double m_src_dej;

	size_t m_samples_read;
	size_t m_current_sample;

	FILE* m_file;
	int m_fd;
	char* m_filename;
	char m_hdr[MAX_HDR_SIZE];
	size_t m_hdr_nbytes;
	size_t m_file_nbytes;
	size_t m_data_nbytes;

};

#endif /* SIGPROCFILE_H_ */
