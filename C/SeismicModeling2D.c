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
	float* output = (float*)malloc(n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		float td = i * dt - tlag;
        output[i] = (1 - 2 * pi*(pi*fc*td)*(pi*fc*td))*exp(-pi * (pi*fc*td)*(pi*fc*td));
    }
	return(output);
}

float* ricker_short(int n,float fcut, float dt)
{
	float pi      = 4 * atan(1);
	float fc      = fcut / (3 * sqrt(pi));
	float t_src   = 0.0;
	float t0_src  = 4*sqrt(pi)/fcut;	
    float aux     = 0.0;
	int Nt_src    = 2*t0_src/dt + 1;

	float* output = (float*)malloc(n * sizeof(float));for (int i=0; i < n;i++) output[i]=0;
	for (int i = 0; i < Nt_src; i++)
	{
		t_src = i*dt-t0_src;
		aux = pi*(pi*fc*t_src)*(pi*fc*t_src);
		output[i] = -(2*aux-1)*exp(-aux);
    }
	// printf("Source samples = %d \n", Nt_src);
	return(output);
}

float ricker_short_Nsample(float fcut, float dt)
{
	float pi      = 4 * atan(1);
	float fc      = fcut / (3 * sqrt(pi));
	float t0_src  = 4*sqrt(pi)/fcut;	
	float Nt_src    = 2*t0_src/dt + 1;

	// printf("Source samples = %f \n", Nt_src);
	return(Nt_src);
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
	
	printf("%s written successfully. Number of lines = %d \n", name,N_POINTS);

}

void export_ascii(char* name, int N_lines, float* vector)
{
	FILE * fp;
	fp = fopen (name,"w");
	for(int i = 0; i<N_lines;i++)
	{		
		fprintf(fp, "%f \n", vector[i]);
	}
	printf("%s written successfully. \n", name);
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
    int Nx        = (int)1000;
	int Nz        = (int)1000;
	float dx      = (float)10;
	float dz      = (float)10;

	int ini_x   = (int)4;
	int end_x   = (int)Nx-4;
	int ini_z   = (int)4;
	int end_z   = (int)Nz-4;

	/* Time parameters*/
	int Nt            = (int)2001;
	float dt          = (float)1.0e-3;
	int Nsnap         = (int)10;
	int snaptime      = Nt/Nsnap;

	/* Source parameters*/
	int sx           = Nx/2;
	int sz           = 20;
	float fcut       = 30;

	/* Receiver position */
	int rz           = 20;
	int Nchannel     = Nx;
	int Nsamples     = Nt;

	int count = 0;

    /* Velocity model */
    float* VP = (float*)malloc(Nx*Nz*sizeof(float));
    for (int index = 0; index <= Nx*Nz; index++)
    {
		int i = index / Nz;
		int j = index - i * Nz;
		if (j <= 100){
        	VP[j + Nz*i] = 1500;
		}
		else{
			VP[j + Nz*i] = 1800;
		}
    }
	
    /* Source wavelet */
    //float* wavelet = ricker(Nt, 30, Nt*dt/5, dt);	

	/* Source wavelet */
    float* wavelet = ricker_short(Nt,30, dt);	
	float src_samples = ricker_short_Nsample(30, dt);
	
	/* Allocate arrays */
	float* P1 = (float*)malloc(Nx*Nz*sizeof(float)); for (int i=0; i < Nx*Nz;i++) P1[i]=0;
	float* P2 = (float*)malloc(Nx*Nz*sizeof(float)); for (int i=0; i < Nx*Nz;i++) P2[i]=0;
	float* P3 = (float*)malloc(Nx*Nz*sizeof(float)); for (int i=0; i < Nx*Nz;i++) P3[i]=0;
	
	float* Seismogram = (float*)malloc(Nx*Nt*sizeof(float)); for (int i=0; i < Nx*Nt;i++) Seismogram[i]=0;

	float* snapshot = (float*)malloc(Nsnap*Nx*Nz*sizeof(float));
	
	// export_float32("waveletricker.bin", src_samples, wavelet);
	for (int n=0; n<Nt;n++)
	{		
		/* Injecting the source*/
		P2[sz + sx*Nz] = P2[sz + sx*Nz] - wavelet[n];
	
		/* Solve wave equation */
		for (int index = 4; index < Nx*Nz-4;index++)
		{
			int i = index / Nz;
			int j = index - i * Nz;
			if ( (i > ini_x) && (i < end_x) && (j > ini_z) && (j < end_z) )
			{
				float p_zz = (-1*P2[(j-2) + Nz*i]+16*P2[(j-1) + Nz*i]-30*P2[j + Nz*i]+16*P2[(j+1) + Nz*i]-1*P2[(j+2) + Nz*i])/(12*dz*dz);
				float p_xx = (-1*P2[j + Nz*(i-2)]+16*P2[j + Nz*(i-1)]-30*P2[j + Nz*i]+16*P2[j + Nz*(i+1)]-1*P2[j + Nz*(i+2)])/(12*dx*dx);

				P3[index] = 2*P2[j + Nz*i] - P1[j + Nz*i] + (dt*dt)*(VP[j + Nz*i]*VP[j + Nz*i])*(p_xx + p_zz);	
			}
		}

		/* Registering seismogram*/
		for (int rx = 0; rx < Nchannel;rx++)
		{
			Seismogram[n + rx*Nt] = P3[rz + rx*Nz] ;
		}
		/*Update fields*/
		for (int index = 4; index < Nx*Nz-4;index++)
		{
			P1[index] = P2[index];
			P2[index] = P3[index];			
		}
		if (n % snaptime == 0){			
			
			// /*Registering Snap shot*/
			for (int i=0; i < Nx*Nz;i++) snapshot[i + count*(Nx*Nz)]=P3[i]+1.0e-3*VP[i];
			count = count + 1;
			printf("Propagation time = %f. Registering snapshot %d. \n", dt*n,count);
		}
	}
	/*Save seismogram in disk*/
	export_float32("seismogram.bin", Nchannel*Nt, Seismogram);	

	/*Writting Snapshot in disk */
	export_float32("snapshot.bin", Nx*Nz*Nsnap, snapshot);

}