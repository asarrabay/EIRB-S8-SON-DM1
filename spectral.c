#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <complex.h>
#include <sndfile.h>
#include <fftw3.h>

#include <math.h>

#include "gnuplot_i.h"



#define EPSILON 30
#define NUM_SIZE 20

/* taille de la fenetre */
// #define	FRAME_SIZE 4096
#define	FRAME_SIZE 1760 // 40ms à 44100 Hz
/* avancement */
// #define HOP_SIZE 1024
#define HOP_SIZE 1760 // 30ms à 44100 Hz

#define CONVERT(freqA, freqB, num1, num2, num3)\
    if (abs(freqB - 1209) < EPSILON ) {\
        return num1;\
    } else if(abs(freqB - 1336) < EPSILON){\
        return num2;\
    } else if(abs(freqB - 1477) < EPSILON){\
        return num3;\
    } else{\
        fprintf(stderr, "%d/%d Hz ne correspond à aucun numéro\n", freqA, freqB);\
        return 'X';\
    }\

int PLOT = 0;

static gnuplot_ctrl *h;


void dft (double s[FRAME_SIZE], complex S[FRAME_SIZE]){
    for (size_t m = 0; m < FRAME_SIZE; m++) {
        S[m] = 0;
        for (size_t n = 0; n < FRAME_SIZE; n++) {
            S[m] += s[n] * cexp(-I * 2 * M_PI / FRAME_SIZE * n * m);
        }
    }
}



static void
print_usage (char *progname)
{
    fprintf (stderr, "\nUsage : %s <input file> [plot]\n", progname) ;
    puts ("\n");

}
static void
fill_buffer(double *buffer, double *new_buffer)
{
    int i;
    double tmp[FRAME_SIZE-HOP_SIZE];

    /* save */
    for (i=0;i<FRAME_SIZE-HOP_SIZE;i++)
    tmp[i] = buffer[i+HOP_SIZE];

    /* save offset */
    for (i=0;i<(FRAME_SIZE-HOP_SIZE);i++)
    {
        buffer[i] = tmp[i];
    }

    for (i=0;i<HOP_SIZE;i++)
    {
        buffer[FRAME_SIZE-HOP_SIZE+i] = new_buffer[i];
    }
}

static int
read_n_samples (SNDFILE * infile, double * buffer, int channels, int n)
{

    if (channels == 1)
    {
        /* MONO */
        int readcount ;

        readcount = sf_readf_double (infile, buffer, n);

        return readcount==n;
    }
    else if (channels == 2)
    {
        /* STEREO */
        double buf [2 * n] ;
        int readcount, k ;
        readcount = sf_readf_double (infile, buf, n);
        for (k = 0 ; k < readcount ; k++)
        buffer[k] = (buf [k * 2]+buf [k * 2+1])/2.0 ;

        return readcount==n;
    }
    else
    {
        /* FORMAT ERROR */
        fprintf (stderr, "Channel format error.\n");
    }

    return 0;
}

static int
read_samples (SNDFILE * infile, double * buffer, int channels)
{
    return read_n_samples (infile, buffer, channels, HOP_SIZE);
}

char getNumber(int freqA, int freqB){
    if (freqA > freqB) {
        int tmp = freqA;
        freqA = freqB;
        freqB = tmp;
    }
    if (abs(freqA - 710) < EPSILON) {
        CONVERT(freqA, freqB, '1', '2', '3');
    } else if(abs(freqA - 780) < EPSILON) {
        CONVERT(freqA, freqB, '4', '5', '6');
    } else if(abs(freqA - 860) < EPSILON){
        CONVERT(freqA, freqB, '7', '8', '9');
    } else if(abs(freqA - 950) < EPSILON){
        CONVERT(freqA, freqB, '*', '0', '#');
    } else{
        fprintf(stderr, "%d/%d Hz ne correspond à aucun numéro\n", freqA, freqB);
        return 'X';
    }
}

//Exercice 2
float rms(double *buffer, int size)
{
  float s = 0;
  for (int k =0; k<size; k++)
  {
    s += buffer[k] * buffer[k];
  }
  return s/size;
}

int
main (int argc, char * argv [])
{	char 		*progname, *infilename;
    SNDFILE	 	*infile = NULL ;
    SF_INFO	 	sfinfo ;

    progname = strrchr (argv [0], '/') ;
    progname = progname ? progname + 1 : argv [0] ;

    if (argc > 2) {
        PLOT = 1;
    }
    if (argc > 3)
    {	print_usage (progname) ;
        return 1 ;
    } ;

    infilename = argv [1] ;

    if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
    {	fprintf (stderr, "Not able to open input file %s.\n", infilename) ;
    puts (sf_strerror (NULL)) ;
    return 1 ;
} ;


/* Read WAV */
int nb_frames = 0;
double new_buffer[HOP_SIZE];
double buffer[FRAME_SIZE];

/* Plot Init */
h = gnuplot_init();
gnuplot_setstyle(h, "lines");

int i;
for (i=0;i<(FRAME_SIZE/HOP_SIZE-1);i++)
{
    if (read_samples (infile, new_buffer, sfinfo.channels)==1)
    fill_buffer(buffer, new_buffer);
    else
    {
        fprintf(stderr, "not enough samples !!\n");
        return 1;
    }
}

/* Info file */
fprintf(stderr, "sample rate %d\n", sfinfo.samplerate);
fprintf(stderr, "channels %d\n", sfinfo.channels);
fprintf(stderr, "size %d\n", (int)sfinfo.frames);

static fftw_plan plan;
fftw_complex data_in[sfinfo.samplerate];
fftw_complex data_out[sfinfo.samplerate];
for (size_t i = 0; i < sfinfo.samplerate; i++) {
    data_in[i] = 0;
    data_out[i] = 0;
}
plan = fftw_plan_dft_1d(sfinfo.samplerate, data_in, data_out, FFTW_FORWARD, FFTW_ESTIMATE);



// fftw_real s[FRAME_SIZE]; /* domaine temporel */



for (int i = FRAME_SIZE; i < sfinfo.samplerate; i++) {
    data_in[i] = 0;
}
int afterSilence = 1;
int cpt = 0;
char num[NUM_SIZE];
while (read_samples (infile, new_buffer, sfinfo.channels)==1)
{
    /* Process Samples */
    fprintf(stderr, "Processing frame %d\n",nb_frames);

    /* hop size */
    fill_buffer(buffer, new_buffer);
    fprintf(stderr,"RMS %f\n",rms(buffer,FRAME_SIZE));

    for (i=0; i<FRAME_SIZE; i++){
        // data_in[i] = buffer[i];
        // Fenetre de Hann :
        data_in[i] = buffer[i] * (0.5 - 0.5 * cos(2.0* M_PI * i/FRAME_SIZE));
    }
    fftw_execute (plan);


    double amp[FRAME_SIZE];
    int ind_max = -1;
    int ind_max2 = -1;
    for (size_t i = 0; i < FRAME_SIZE; i++) {
        amp[i] = cabs(data_out[i]);
    }




    for (int i = 0; i < FRAME_SIZE; i++) {
        if (amp[i] > 30.0 && amp[i-1] <= amp[i] && amp[i] > amp[i+1]) {
            fprintf(stderr, "FREQUENCE ============================= %d : %lf\n", i, amp[i]);
            if (ind_max == -1) {
                ind_max = i;
            } else{
                ind_max2 = i;
            }
        }
    }

    fprintf(stderr, "frequence = %d / %d\n", ind_max, ind_max2);

    if (ind_max == -1 || ind_max2 == -1) {
        afterSilence = 1;
        fprintf(stderr, "SILENCE\n");
    } else {
        char c_num = getNumber(ind_max, ind_max2);
        fprintf(stderr, "NUM : %c\n", c_num);
        if (afterSilence) {
            afterSilence = 0;
            num[cpt] = c_num;
            cpt++;
            num[cpt] = '\0';
            fprintf(stderr, "\t%s\n", num);
        }
    }

    if (PLOT) {
        gnuplot_resetplot(h);
        gnuplot_plot_x(h,amp,FRAME_SIZE,"temporal frame");
        sleep(1);
    }

    nb_frames++;
}

num[cpt] = '\0';
printf("numero : %s\n", num);

fftw_destroy_plan (plan);

sf_close (infile) ;

return 0 ;
} /* main */
