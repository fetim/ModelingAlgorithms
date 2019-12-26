#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int larger(int a ,int b)
{
    if (a>b)
        return a;
    return b;
}

float* ricker(int n, float fcorte, float tlag, float dt)
{
	float pi = 4 * atan(1);
	float fc = fcorte / (3 * sqrt(pi));
    float aux1 = 0.0;
    float aux2 = 0.0;
	float* output = (float*)malloc(n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		float td = i * dt - tlag;
        output[i] = (1 - 2 * pi*(pi*fc*td)*(pi*fc*td))*exp(-pi * (pi*fc*td)*(pi*fc*td));
    }
	return(output);
}

float* ricker_short(float fcut, float dt)
{
	float pi      = 4 * atan(1);
	float fc      = fcut / (3 * sqrt(pi));
	float t_src   = 0.0;
	float t0_src  = 4*sqrt(pi)/fcut;	
    float aux     = 0.0;
	int Nt_src    = 2*t0_src/dt + 1;

	float* output = (float*)malloc(Nt_src * sizeof(float));
	for (int i = 0; i < Nt_src; i++)
	{
		t_src = i*dt-t0_src;
		aux = pi*(pi*fc*t_src)*(pi*fc*t_src);
		output[i] = (2*aux-1)*exp(-aux);
		printf("output[%i] = %f \n",i, output[i]);
    }
	printf("Source samples = %d \n", Nt_src);
	return(output);
}

void export_float32(char* name, int N_POINTS, float* vector)
{
	FILE* fp;
	fp = fopen(name, "wb");
	if (fp != NULL)
	{
		fwrite((char*)vector, N_POINTS * sizeof(float), 1, fp);
	}
	fclose(fp);
	
	printf("%s successfully written. \n", name);

}

void export_ascii(char* name, int N_lines, float* vector)
{
	FILE * fp;
	fp = fopen (name,"w");
	for(int i = 0; i<N_lines;i++)
	{		
		fprintf(fp, "%f \n", vector[i]);
	}
	printf("%s successfully written. \n", name);
	fclose(fp);
}


float* import_float32(char* name, int N_POINTS)
{
	FILE * fp;
	size_t result;

	fp = fopen(name, "rb");
	if (fp == NULL) { fputs("File error", stderr); exit(1); }

	float* buffer = (float*)malloc(N_POINTS * sizeof(float));

	if (buffer == NULL) { fputs("Memory error", stderr); exit(2); }

	result = fread(buffer, sizeof(float), N_POINTS, fp);
	if (result != N_POINTS) { fputs("Reading error", stderr); exit(3); }

	fclose(fp);
	return (buffer);
}


float* import_ascii(char* name,int N_lines)
{
	FILE * fp;
	char str[100];
	int count = 0;
	float* output = (float*)malloc(N_lines * sizeof(float));
	fp = fopen (name, "r");  

	if(fp == NULL) {
      perror("Error opening file");
    }
    else{
		for(int i = 0; i<N_lines;i++)
		{
			fscanf(fp, "%s", str);
			output[i] = atoi(str);
			count++;
		}
		printf("\n\n");
		printf("Reading %s.  Number of lines = %i \n",name, count);
		fclose (fp);
		return(output);
	}
}

int main()
{
    /* Model parameters*/
    int SIZE_P        = (int)1000;
	int SIZE_L        = (int)1000;
	float dP          = (float)10;
	float dL          = (float)10;

	/* Time parameters*/
	int Nt            = (int)2001;
	float dt          = (float)1.0e-3;

	 /* Source parameters*/
	 int sx           = 200;
	 int sz           = 20;
	 float fcut       = 30;

	 /* Receiver position */
	 int rz           = 20;

    /* Velocity model */
    float* CP = (float*)malloc(SIZE_P*SIZE_L * sizeof(float));
    for (int i = 0; i <= SIZE_P*SIZE_L; i++)
    {
        CP[i] = 1500;
    }
	
    /* Source wavelet */
    // float* wavelet = ricker(Nt, 30, Nt*dt/5, dt);	

	int src_samples;

    float* wavelet = ricker_short(30, dt);
		
	export_float32("waveletricker.bin", src_samples, wavelet);

}