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
	float* saida = (float*)malloc(n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		float td = i * dt - tlag;
        saida[i] = (1 - 2 * pi*(pi*fc*td)*(pi*fc*td))*exp(-pi * (pi*fc*td)*(pi*fc*td));
    }
	return(saida);
}

void export_float32(char* nome, int N_PONTOS, float* vetor)
{
	FILE * fp;
	fp = fopen(nome, "wb");
	if (fp != NULL)
	{
		fwrite((char*)vetor, N_PONTOS * sizeof(float), 1, fp);
	}
	fclose(fp);
	//	printf("%s escrito com sucesso!\n", nome);

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
	
    int N_PONTOS = TAM_P*TAM_L;
    
    /* Velocity model */
    float* CP = (float*)malloc(N_PONTOS * sizeof(float));
    for (int i = 0; i <= N_PONTOS; i++)
    {
        CP[i] = 1500;
    }

    /* Source wavelet */
    float* wavelet = ricker(N_ITERACAO, 30, N_ITERACAO*dt/2, dt);
    
    export_float32("waveletricker.bin", N_ITERACAO, wavelet);

    
    // for (int i=0;i<= N_ITERACAO; i++)
    // {
    //     printf("%i sample %f \n",i, wavelet[i]);
    // }
    // return 0;

}