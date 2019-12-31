#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#ifdef _OPENMP
	#include<omp.h>
#endif

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
	if (n<Nt_src)
	{
		printf("Error! Nt (%d) < source samples (%d)\n",n,Nt_src);
		exit(0);
	}

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
      perror("Error opening file\n");
    }
    else{
		for(int i = 0; i<N_lines;i++)
		{
			fscanf(fp, "%s", str);
			output[i] = atof(str);
			count++;
		}
		printf("\n\n");
		printf("Reading %s.  Number of lines = %i \n",name, count);
		fclose (fp);
		return(output);
	}
}

float max(int N_lines, float* vector) {
    float value;    
    value = vector[0];
    for(int i=0; i<N_lines; i++) {
            if(value < vector[i]) {
			value = vector[i];            
    		}
	}
    // printf("Maximum speed in the model = %.0f m/s\n",value);
    return value;
}

float min(int N_lines, float* vector) {
    float value;    
    value = vector[0];
    for(int i=0; i<N_lines; i++) {
            if(value > vector[i]) {
			value = vector[i];            
    		}
	}
    // printf("Mininum speed in the model = %.0f m/s\n",value);
    return value;
}

void check_acoustic_stability(float dt,float dh,int N_lines,float vpmax, float vpmin, float fcut, int T_order, int S_order)
{		
	if ((T_order == 2) && (S_order==4))
	{
		float alpha  = 5;		
		float beta   = 4;
		float dh_max = vpmin/(alpha*fcut);
		float dt_max = dh_max/(beta*vpmax);		
			
		if ( (dt <= dt_max) && (dh <= dh_max ))
		{
			printf("Dispersion and stability conditions ... OK!\n");			
		}
		else
		{
			printf("\nRead time and space parameters\n");
			printf("Time                        = %f\n", dt);
			printf("Space                       = %f\n", dh);					
			
			printf("\nRecommended time spacing and time interval\n");
			printf("Time  - Stability limit  - dt < %f\n", dh/(beta*vpmax));
			printf("Time  - Stability limit  - dt < %f\n", dt_max);
			printf("Space - Dispersion limit - dh < %f\n", dh_max);		
			printf("The program will be stopped! \n \n");	
			exit(0);
		}
		
	}
	else
	{
		printf("Select the right order \n");
		printf(" Time derivative order =%d Space derivative order = %d, \n", T_order, S_order);
		
	}
}
int main()
{
    /* Model parameters*/
    float* parameters = import_ascii("../../parameters/2D_acoustic_modeling.dat",7);

  	int Nx        = (int)parameters[0];
	int Nz        = (int)parameters[1];
	float dx      = (float)parameters[2];
	float dz      = (float)parameters[3];

	printf("Nx     = %d \n", Nx);
	printf("Nz     = %d \n", Nz);
	printf("dx     = %f \n", dx);
	printf("dz     = %f \n", dz);

	int ini_x   = (int)4;
	int end_x   = (int)Nx-4;
	int ini_z   = (int)4;
	int end_z   = (int)Nz-4;

	/* Time parameters*/
	int Nt            = (int)parameters[4];
	float dt          = (float)parameters[5];
	
	printf("Nt     = %d \n", Nt);
	printf("dt     = %f \n", dt);	

	/* Snapshots */
	int reg_snapshot  = 1;
	int Nsnap         = (int)10;
	int snaptime      = Nt/Nsnap;

	/* Source parameters*/
int Nshot        = 3;

	int* sx          = (int*)malloc(Nshot*sizeof(int));for (int i=0; i < Nshot;i++) sx[i]=Nx/2 + 5*i  ;
	int* sz          = (int*)malloc(Nshot*sizeof(int));for (int i=0; i < Nshot;i++) sz[i]=20;

	float fcut       = (float)parameters[6];
	printf("fcut   = %f \n",fcut);

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
        	VP[j + Nz*i] = 1600;
		}
		else{
			VP[j + Nz*i] = 1800;
		}
    }
	
	float vpmax = max(Nx*Nz,VP);
	float vpmin = min(Nx*Nz,VP);
	check_acoustic_stability(dt,dz,Nx*Nz,vpmax,vpmin,fcut,2,4);

    /* Source wavelet */
    //float* wavelet = ricker(Nt, 30, Nt*dt/5, dt);	

	/* Source wavelet */
    float* wavelet = ricker_short(Nt,30, dt);	
	int src_samples = ricker_short_Nsample(30, dt);
	
	/* Allocate arrays */
	float* P1 = (float*)malloc(Nx*Nz*sizeof(float)); for (int i=0; i < Nx*Nz;i++) P1[i]=0;
	float* P2 = (float*)malloc(Nx*Nz*sizeof(float)); for (int i=0; i < Nx*Nz;i++) P2[i]=0;
	float* P3 = (float*)malloc(Nx*Nz*sizeof(float)); for (int i=0; i < Nx*Nz;i++) P3[i]=0;
	
	float* Seismogram = (float*)malloc(Nchannel*Nsamples*Nshot*sizeof(float)); for (int i=0; i < Nchannel*Nsamples;i++) Seismogram[i]=0;

	float* snapshot = (float*)malloc((Nsnap+1)*Nx*Nz*sizeof(float));
	
	// export_float32("waveletricker.bin", src_samples, wavelet);
	for (int shot=0; shot<Nshot;shot++)
	{
		printf("Running shot %d \n", shot+1);
		for (int n=0; n<Nt;n++)
		{		
			/* Injecting the source*/
			P2[sz[shot] + sx[shot]*Nz] = P2[sz[shot] + sx[shot]*Nz] - wavelet[n];
		
			/* Solve wave equation */
			#pragma omp parallel shared(P1,P2,P3,Seismogram,dt,dx,dz,Nx,Nz,ini_x,ini_z,end_x,end_z,Nt,rz)
			{
				#pragma omp for
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
				#pragma omp for nowait
				for (int rx = 0; rx < Nchannel;rx++)
				{
						Seismogram[n + (rx*Nsamples) + shot*(Nchannel*Nsamples)] = P3[rz + rx*Nz] ;
				}
				/* Update fields */			
				#pragma omp for nowait
				for (int index = 4; index < Nx*Nz-4;index++)
				{
					P1[index] = P2[index];
					P2[index] = P3[index];			
				}
			}

			/* Registering Snap shot */
			if ((n % snaptime == 0) && reg_snapshot==shot+1)
			{			
				for (int i=0; i < Nx*Nz;i++) snapshot[i + count*(Nx*Nz)]=P3[i]+1.0e-3*VP[i];
				count = count + 1;
				printf("Propagation time = %f. Registering snapshot %d. \n", dt*n,count);
			}
		}
		// Restarting fields // 
		for (int i=0; i < Nx*Nz;i++) P1[i]=0;
		for (int i=0; i < Nx*Nz;i++) P2[i]=0;
		for (int i=0; i < Nx*Nz;i++) P3[i]=0;
	}
	/*Save seismogram in disk*/
	export_float32("seismogram.bin", Nchannel*Nsamples*Nshot, Seismogram);	

	/*Writting Snapshot in disk */
	if (reg_snapshot){export_float32("snapshot.bin", Nx*Nz*Nsnap, snapshot);}

}