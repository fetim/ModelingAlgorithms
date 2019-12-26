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

void export_float32(char* name, int N_POINTS, float* vector)
{
	FILE* fp;
	fp = fopen(name, "wb");
	if (fp != NULL)
	{
		fwrite((char*)vector, N_POINTS * sizeof(float), 1, fp);
	}
	fclose(fp);
	
	printf("%s successfully complete a written \n", name);

}

void export_ascii(char* name, int N_lines, float* vector)
{
	FILE * fp;
	fp = fopen (name,"w");
	for(int i = 0; i<N_lines;i++)
	{		
		fprintf(fp, "%f \n", vector[i]);
	}
	printf("%s successfully complete a written, \n", name);
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
		printf("\n The content of the file %s is  :\n",name);
		for(int i = 0; i<N_lines;i++)
		{
			fscanf(fp, "%s", str);
			output[i] = atoi(str);
			printf("%f input \n", output[i]);
			count++;
		}
		printf("\n\n");
		printf("Number of lines = %i \n", count);
		fclose (fp);
		return(output);
	}
}

int main()
{
    /* Main parameters*/
    int TAM_P         = (int)1000;
	int TAM_L         = (int)1000;
	float dP          = (float)10;
	float dL          = (float)10;
	int N_SOURCES     = (int)1;
	int N_ITERACAO    = (int)8001;
	float dt          = (float)4.0e-4;
    float ratio       = (float)1;
	int BORDA         = (int)100;
	int N_REC_VP      = (int)1000;
	
    int N_POINTS = TAM_P*TAM_L;
    
    /* Velocity model */
    float* CP = (float*)malloc(N_POINTS * sizeof(float));
    for (int i = 0; i <= N_POINTS; i++)
    {
        CP[i] = 1500;
    }

	
	// char fname[20]="666";
	// char nametest[20];
	// // float recebe;

	// printf(" %s \n",fname);
	
	// sprintf(nametest, "%f", dt);

	// printf(" %s \n",nametest);
		
	// recebe = atoi(fname);
	
	// printf(" %f \n",recebe/4);
    /* Source wavelet */
    float* wavelet = ricker(N_ITERACAO, 30, N_ITERACAO*dt/2, dt);

	float* arraytest = (float*)malloc(N_ITERACAO * sizeof(float));

	//export_ascii("asciifile.dat",N_ITERACAO, wavelet);


	arraytest = import_ascii("asciifile.dat",N_ITERACAO);

	/* Write a binary in disk */
    // export_float32("waveletricker2.bin", N_ITERACAO, wavelet);

    
    for (int i=0;i < 4; i++)
    {
        printf("%i sample %f \n",i, arraytest[i]);
    }
    // return 0;


}