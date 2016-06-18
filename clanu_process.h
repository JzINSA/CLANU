// Clanu 2015 - 2016
// T. Grenier : thomas.grenier@insa-lyon.fr

#ifndef CLANU_PROCESS_H
#define CLANU_PROCESS_H

// THESE 3 DECLARATIONS MUST NOT BE MODIFIED
void Question1(float **Rout, float **Gout, float **Bout, float **Rin, float **Gin, float **Bin, float **Mask, int width, int height, double param);
void InpaintingBW(float **Iout, float **Iin, float **Mask, int width, int height, double param);
void InpaintingColor(float **Rout, float **Gout, float **Bout, float **Rin, float **Gin, float **Bin, float **Mask, int width, int height, double param);

//Add your own functions' declaration below

float **matriceA(float **I, float **Mask, float **temp, int w , int h);
float **subtractMatrix(float **A, float **B, float **temp, int width, int height);
float **addMatrix(float **A, float **B, float **temp, int width, int height);
float **multiplyScalaire(float **A, float alpha, float **temp, int width, int height);
float psMatrix(float **A,float **B, int width, int height);
void equalMatrix(float **A,float **B, int width, int height);
void divide255(float **I,int width, int height);
void multiply255(float **I,int width, int height);
void complementMask(float **Mask,int width, int height);

// nothing after this line
#endif 

