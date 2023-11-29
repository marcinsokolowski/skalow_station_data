#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <complex>

// 8192 -> 8200 (for timestamp)
#define N_CHANNELS (24*32*4) // 10kHz channels 24*32*4 = 3072
#define MAX_N_SAMPLES 2*N_CHANNELS // Re/Im of 24*32 fine channels 
#define DATA_TYPE     short

int gVerb=0;
int gSnapshotIdx=-1;
int gDumpChannel=0;
char filename[1024];
char gOutFile[1024],gOutFileAvg[1024],gOutFileDat[1024];
int gAddIndex=0;
int gAverageNSpectra=0;
int gPol=1; // 0 - X , 1 - Y 
int gNChannels=1;
int gSaveNChannels=1;

#define BYTES_PER_SAMPLE 4 // this is re/im for X and re/im for Y -> 2 x 2 bytes = 4 bytes 

int gNBytesPerSample = BYTES_PER_SAMPLE;
int n_samples = 1024 * BYTES_PER_SAMPLE; // number of samples to read, used to be N_CHANNELS; 
int gStartIndex=-1;
int gSkipNBytes=-1;
int gProcessNSamples = -1;
int gMultiplier=1;
int gCheckZeros=0;

double gCorrWithInvJones=-1; // normalisation factor
bool gXYPhase = false;

// # /home/msok/Desktop/Curtin/students/2022/Chris_Lee_Honours_Project/logbook/20220412_beam-correcting_VELA_data.odt
// std::complex<double> gInvJones00( -0.19391154 , -0.00342333 ); // -0.133340316282896 , 0.039036116382880 );
// std::complex<double> gInvJones01( 1.1539054   , +0.37828104 ); // -0.860154030652396 , 0.187848858551529 );
// std::complex<double> gInvJones10( -1.07909468 , -0.24393281 ); // 0.764011070552511 , -0.244298897390238 );
// std::complex<double> gInvJones11( -0.18581259 , -0.04685302 ); // -0.122434290299445 , 0.069153015019698 );

// # page 5 of /home/msok/Desktop/Curtin/students/2022/Chris_Lee_Honours_Project/logbook/20220509_VELA_obs_amps_and_InvAmps_BeamCorr.odt
std::complex<double> gInvJones00( 0.01824534 , 0.00622549 );
std::complex<double> gInvJones01( 0.07101837 , 0.02065822 );
std::complex<double> gInvJones10( -0.06648148 , -0.01180349 );
std::complex<double> gInvJones11(  0.01690944 , 0.0054325 );


void usage()
{
   printf("read_binary_station_beam_test_order1_2pol voltages_r.dat -d -u -c -f outfile -i integration_index -v -a AVERAGE_N_POINTS -p POL -s START_INDEX -n PROCESS_TIMESTAMPS -M multiplier -Z -S SAVE_N_CHANNELS\n");
   printf("\n");
   printf("Options :\n");
//   printf("\t-u - assumes that UNIXTIME timestamp is added to line as extra record of DATA_TYPE structure\n");
   printf("\t-d - increases debug level\n");
   printf("\t-v - increases debug level\n");   
   printf("\t-c CHANNEL : specifies channel to be dumped [default 0]\n");      
   printf("\t-C N_CHANNELS : number of channels [default %d]\n",gNChannels);
   printf("\t-f OUTFILENAME : name of file to which data is dumped\n");
   printf("\t-i INTEGRATION_INDEX : index of integration (= adcv dump) to be dumped, by default all are dumped\n");   
   printf("\t-a AVERAGE_N_POINTS : average N data points [default %d], <=1 -> do not average\n",gAverageNSpectra);
   printf("\t-p POL : POL = 0 = X , POL = 1 = Y , POL = -1 X+Y (power)\n");
   printf("\t-s START_INDEX : where to set start index [default %d] , <0 -> starting from the start of the file (and save these many bytes to file header_saved.txt\n",gStartIndex);   
   printf("\t-S SAVE_N_CHANNELS : how many channels to save into the output file [default %d]\n",gSaveNChannels);
   printf("\t-n PROCESS_TIMESTAMPS : number of timestamps to process [default %d and <0 is all]\n",gProcessNSamples);
   printf("\t-M multiplier : multiply timestep index (for example 128 fine channels * 7 averages to get ~1ms)\n");
   printf("\t-Z : check statistics of ZEROS\n");
   printf("\t-X N_BYTES : skip these many bytes [default %d]\n",gSkipNBytes);
   printf("\t-I : if >0 correct with inverse Jones and multiply by this factor [default %.4f]\n",gCorrWithInvJones);
   printf("\t-Y : enable XY phase correction\n");
   
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "vudc:C:f:i:r:m:p:s:n:a:M:ZS:X:I:Y";
   int opt,opt_param,i;
   
//   strcpy(filename,"");
   strcpy(gOutFile,"out_file");
      
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'c':
            if( optarg ){   
               if( atol(optarg)>=0 && atol(optarg)<=512 ){
                 gDumpChannel = atol(optarg);
               }else{
                 printf("ERROR : wrong number of channels %d , not in range [0-512]!\n",(int)atol(optarg));
               }
            }
            break;

         case 'f':
            if( optarg ){   
               strcpy(gOutFile,optarg);
            }
            break;

         case 'i':
            if( optarg ){   
               gSnapshotIdx = atol(optarg);
            }
            break;

         case 'I':
            gCorrWithInvJones = atof(optarg);
            break;

         case 'a':
            if( optarg ){   
               gAverageNSpectra = atol(optarg);
            }
            break;

         case 'M':
            if( optarg ){   
               gMultiplier = atol(optarg);
            }
            break;

         case 'C':
            if( optarg ){   
               gNChannels = atol(optarg);
            }
            break;

         case 'p':
            if( optarg ){   
               gPol = atol(optarg);
               if( gPol >= 2 ){
                  gPol = 1;
               }
               if( gPol <= -1 ){
                  gPol = 0;
               }
            }
            break;
            
         case 'd':
            gVerb++;
            break;

         case 's':
            if( optarg ){
               gStartIndex = atol( optarg );
            }
            break;
            
         case 'S':
            if( optarg ){
               gSaveNChannels = atol( optarg );
            }
            break;
            
         case 'n':
            if( optarg ){
               gProcessNSamples = atol( optarg );
            }
            break;
            
         case 'v':
            gVerb++;
            break;
            
         case 'X':
            if( optarg ){
               gSkipNBytes = atol( optarg );
            }
            break;
            
         case 'Z':
            gCheckZeros = 1;
            break;

         case 'Y':
            gXYPhase = true;
            break;
            
         default:
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();                                          
      }
   }
   
/*   if( strlen(filename) == 0 ){
      if( gStartIndex > 0 && gProcessNSamples > 0 ){
         sprintf(filename,"%ld_%ld.dat",gStartIndex,(gStartIndex+gProcessNSamples));
      }else{
         strcpy(filename,"test.dat");
      }
   }*/

   sprintf(gOutFileAvg,"%s_avg.txt",gOutFile);
   sprintf(gOutFileDat,"%s.dat",gOutFile);

   if( strlen(gOutFile) == 0 ){
      if( gStartIndex > 0 && gProcessNSamples > 0 ){
         sprintf(gOutFile,"%d_%d.txt",gStartIndex,(gStartIndex+gProcessNSamples));
      }else{
         strcpy(gOutFile,"test.txt");
      }
      
      printf("DEBUG : output file name was empty -> set to %s\n",gOutFile);
   }   
   
   gStartIndex = gStartIndex * gMultiplier;
   gProcessNSamples = gProcessNSamples * gMultiplier;
}                                                                           

void printf_parameters()
{
  printf("#####################################\n");
  printf("PARAMETERS :\n");
  printf("#####################################\n");
  printf("Input file      = %s\n",filename);
  if( strlen(gOutFile) ){
     printf("Output file        = %s\n",gOutFile);
  }else{
     printf("Output file        = -\n");
  }
  printf("Dump channeld idx  = %d\n",gDumpChannel);
  printf("Start Index        = %d (save to header file -s option)\n",gStartIndex);
  printf("Skip N bytes       = %d (no save to header file -X option)\n",gSkipNBytes);
  printf("# of timesteps to process = %d\n",gProcessNSamples);
  printf("# averaged points  = %d\n",gAverageNSpectra);
  printf("Multiplier         = %d\n",gMultiplier);
  printf("Verb level         = %d\n",gVerb);
  printf("Check zeros        = %d\n",gCheckZeros);
  printf("Number of channels = %d\n",gNChannels);
  printf("Save number of channels = %d\n",gSaveNChannels);
  printf("Correct with inverse Jones = %.4f\n",gCorrWithInvJones);
  printf("XY phase correction = %d\n",gXYPhase);
  printf("#####################################\n");       
}

int main(int argc,char* argv[])
{
  if( argc<2 || (argc==2 && strncmp(argv[1],"-h",2)==0) ){
    usage();
  }

  if( argc>=2 ){
    strcpy(filename,argv[1]);
  }
  
  // parsing and printing paramters :
  parse_cmdline(argc-1,argv+1);
  printf_parameters();

  n_samples=gNBytesPerSample*gNChannels*1024;  
  FILE *f=NULL;  
  FILE* outf=NULL;
  char* buffer = new char[n_samples];
  int line_size = n_samples;
  int avg_count=0;
  int n = n_samples;

  
  if( gVerb ){
    printf("Reading file %s\n",filename);fflush(stdout);
  }
  f = fopen(filename, "rb");  
  if( strlen(gOutFileAvg) ){
     outf = fopen(gOutFileAvg,"w");
  }
  
  FILE* out_dat_f = fopen( gOutFileDat , "w" );
  
  if (f)
  {
      int i=0,line=0;
      int line_raw=0;
      long int total_processed_samples = 0;
      long int zeros_count = 0;
      
      double accum = 0;
      int accum_count = 0;
      long int avg_index = 0;
      
      if( gSkipNBytes > 0 ){
         fseek( f , gSkipNBytes*1, SEEK_SET ); // was *BYTES_PER_SAMPLE
         printf("DEBUG : set file position to index = %d (option -X)\n",gSkipNBytes*1); // was *BYTES_PER_SAMPLE
      }else{      
         printf("INFO : skiping of bytes with fseek (-X option) is not required, checking if -s to skip and save (.dada header is required)\n");
         // add fseek to find some specific pulses at a given index :                  
         if ( gStartIndex > 0 ){
              printf("INFO skip and save (option -c) is required (gStartIndex = %d bytes)\n",gStartIndex);
              char* buffer_tmp = new char[gStartIndex+10];
              int n_tmp = fread(buffer_tmp, 1, gStartIndex, f);
              if ( n_tmp != gStartIndex ){
                 printf("ERROR : could not read %d bytes of header only read %d\n",gStartIndex,n_tmp);
              }else{
                 printf("INFO : read %d bytes, same as expected %d bytes\n",n_tmp,gStartIndex);
              }
              FILE* out_header_f = fopen("header_saved.txt","wb");
//              fprintf(out_header_f,"%s",buffer_tmp);
              int fd = fileno(out_header_f);
              size_t written = write(fd, buffer_tmp, n_tmp);
              fclose(out_header_f);
              delete [] buffer_tmp;
         }else{
            printf("INFO skip and save (option -c) is not required (gStartIndex = %d bytes)\n",gStartIndex);
         }
      }
      printf("PROGRESS : starting to read %d bytes in each iteration\n",gNBytesPerSample*line_size);fflush(stdout);
      
      char* out_buffer = new char[line_size/gNChannels];

      // gNBytesPerSample = BYTES_PER_SAMPLE = 4        
      bool bContinue = true;
      while( (n = fread(buffer, 1, n_samples, f)) > 0 && bContinue ){ // reads a line of all fine channels 
      
        int timesample=0;
        int bytes_written = 0;
        for(int i=(gDumpChannel*gNBytesPerSample);i<line_size;i+=(gNBytesPerSample*gNChannels)){
           // ordering 1 : X Y | X Y | X Y | X Y ...
           //            CH1    CH2  CH1   CH2 ...
           char x_re = buffer[i];
           char x_im = buffer[i+1];
           char y_re = buffer[i+2];        
           char y_im = buffer[i+3];
           
           if( gCorrWithInvJones > 0 ){
              std::complex<double> vx( x_re, x_im );
              std::complex<double> vy( y_re, y_im );
              
              std::complex<double> e_theta = (vx*gInvJones00 + vy*gInvJones01)*gCorrWithInvJones;
              std::complex<double> e_phi   = (vx*gInvJones10 + vy*gInvJones11)*gCorrWithInvJones;
              
              x_re = e_theta.real();
              x_im = e_theta.imag();
              y_re = e_phi.real();
              y_im = e_phi.imag();
              
              buffer[i] = x_re;
              buffer[i+1] = x_im;
              buffer[i+2] = y_re;
              buffer[i+3] = y_im;
           }
           
/*           if( gXYPhase ){
              double diff = (freq_mhz-323.75950000);
              double xy_phase_rad = 2.00*M_PI*( diff*0.10675035 + diff*diff*(-0.00649167) + diff*diff*diff*0.00001584442515310631 -0.00010545467580092109 ) + M_PI/2.00;
              std::complex<double> xy_phase = tan( xy_phase_rad );
           }*/
           
           
           if( total_processed_samples <= 100 ){
              printf("Line %d : X = %d / %d , Y =  %d / %d\n",line,x_re,x_im,y_re,y_im);
           }

           // just 1 pol. output :
           int write_n_bytes = gSaveNChannels*4; // 2 polarisations * 2 for real/imag = 4 bytes per channel
           memcpy( &(out_buffer[bytes_written]), &(buffer[i]), write_n_bytes ); // 2 for 1 polarisation instead of gNBytesPerSample
           bytes_written += write_n_bytes;
           
           // 2 pols output (not looking good):
           // memcpy( &(out_buffer[bytes_written]), &(buffer[i]), 4 );
           // bytes_written += 4;
        
           double x_power = x_re*x_re + x_im*x_im;
           double y_power = y_re*y_re + y_im*y_im;
           double i_power = ( x_power + y_power )*0.5; // same as option B in calcfits/main.cpp : left.AddImages( right , 0.5 );
           
           if( gPol == 0 ){
              i_power = x_power;
           }
           if( gPol == 1 ){
              i_power = y_power;
           }
           
           if( gVerb >= 3 ){
              printf("%ld : %.4f\n",total_processed_samples,i_power);
           }

           total_processed_samples++;           
           if( gCheckZeros ){
              if( i_power <= 0 ){
                 zeros_count++;
              }
           }

//           if( gAverageNSpectra > 1 ){
/*              if( accum_count <= 0 ){
                 printf("DEBUG : first sample = %.8f\n",i_power);
              }*/
/*              accum += i_power;
              accum_count++;
              
              if( accum_count == gAverageNSpectra ){
                 fprintf(outf,"%ld %.8f\n",avg_index,accum/accum_count);
                 accum_count = 0;
                 accum = 0;
                 avg_index++;
              }*/
//           }else{           
//              fprintf(outf,"%ld %.8f\n",total_processed_samples,i_power);
//           }
           
//           total_processed_samples++;
           
/*           if( gProcessNSamples > 0 ){
              if( total_processed_samples >= gProcessNSamples ){
                 printf("PROGRESS : processed %ld samples (required %ld) -> exiting the data processing loop\n",total_processed_samples,gProcessNSamples);                 
                 bContinue = false;
                 break;
              }
           }*/
           
           timesample++;
        }
        
        if( gVerb>2){ printf("timesample = %d vs. expected %d\n",timesample,(line_size/(gNBytesPerSample*gNChannels))); }
        if( timesample == (line_size/(gNBytesPerSample*gNChannels)) ){        
//           fwrite( out_buffer, timesample, gNBytesPerSample, out_dat_f );
           fwrite( out_buffer, bytes_written, 1, out_dat_f );
        }else{
           printf("ERROR : expected %d time samples, got %d = %d/(%d*%d)\n",(line_size/(gNBytesPerSample*gNChannels)),timesample,line_size,gNBytesPerSample,gNChannels);
           exit(-1);
        }

        if( gVerb > 0 ){
           printf("read line = %d\n",line);
        }
        line++;
      }      
      
      if( gCheckZeros ){
         printf("################################## STATISTICS ##################################\n");
         printf("Number of zeros          = %ld ( %.3f %% )\n",zeros_count,(double(zeros_count)/double(total_processed_samples))*100.00);
         printf("Total number of samples  = %ld\n",total_processed_samples);
         printf("################################################################################\n");
      }
  }else{
     printf("ERROR : could not open file %s\n",filename);
  }

  fclose(f);
  
  if( out_dat_f ){
     fclose( out_dat_f );
  }
  
  if( outf ){
     fclose(outf);
  }
}

