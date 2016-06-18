// Clanu 2015 - 2016
// T. Grenier : thomas.grenier@insa-lyon.fr

#include "clanu_process.h"
#include "imageconvert.h"
#include <algorithm>
#include <stdio.h>
using namespace std;
#include <iostream>
#include <cmath>
// to complete for Q1
void Question1(float **Rout, float **Gout, float **Bout, float **Rin, float **Gin, float **Bin, float **Mask, int width, int height, double param)
{
    for(int x=0; x<width; x++)
        for(int y=0;y<height; y++)
            {

           if(Mask[x][y]>0){
                Rout[x][y] = 255 -  Rin[x][y];
                if(Gin != 0) Gout[x][y] = 255 -  Gin[x][y];
                if(Bin != 0) Bout[x][y] = 255 - Bin[x][y];
            }else{
                Rout[x][y] =Rin[x][y];
                if(Gin != 0) Gout[x][y] = Gin[x][y];
                if(Bin != 0) Bout[x][y] = Bin[x][y];
            }
            }
}

float **matriceA(float **I, float **Mask,float **temp, int w ,int h)
{

   for(int l=0;l<w;l++)
   {
       for(int m=0;m<h;m++)
       {
           temp[l][m]=0;
       }
   }

    int iplusb= 0;
    int  imoinsb=0;
    int  iplus=0;
    int  imoins=0;

    int  jplusb=0;
    int  jmoinsb=0;    
    int  jplus=0;
    int  jmoins=0;


  for(int i=0;i<h;i++)
  {
      if( i<h-1)  { iplusb=1;  }        else{ iplusb=0; }

      if( i+1<h-1)  { iplus=i+1; }    else{ iplus=h-1; }

      if(i>0){  imoinsb=1; }            else { imoinsb=0;}

      if( i>1){ imoins=i-1; }            else{ imoins=0;}


      for(int j=0;j<w;j++)
      {
          if( j<w-1) {jplusb=1; }         else {jplusb=0;}

          if( j+1<w-1){jplus=j+1; }     else {jplus=w-1; }

          if(j>0){ jmoinsb=1;}            else{ jmoinsb=0;}

          if( j>1) { jmoins=j-1; }         else  {  jmoins=0; }

        if(Mask[j][i]>0)
        {
            temp[j][i] =(1/36)*(16*I[j][i] +4*(iplusb*I[j][iplus] + imoinsb*I[j][imoins] + jplusb*I[jplus][i] +jmoinsb*I[jmoins][i])
                                + ((iplusb)*(jplusb)*I[jplus][iplus] + (imoinsb)*(jplusb)*I[jplus][imoins] + (imoinsb)*(jmoinsb)*I[jmoins][imoins])
                                +(iplusb)*(jmoinsb)*I[jmoins][iplus])  ;

        }
        else
        {
             temp[j][i]=-w*h*(-8*I[j][i] + 1*(iplusb*I[j][iplus] + imoinsb*I[j][imoins] + jplusb*I[jplus][i] + jmoinsb*I[jmoins][i])
                            + (jplusb*iplusb*I[jplus][iplus] + jplusb*imoinsb*I[jplus][imoins]
                               + jmoinsb*imoinsb*I[jmoins][imoins] + jmoinsb*iplusb*I[jmoins][iplus]));
        }
      }
  }
      return temp;
}

float **subtractMatrix(float **A,float **B,float **temp, int width, int height)
{
    for(int l=0;l<width;l++)
    {
        for(int m=0;m<height;m++)
        {
            temp[l][m]=0;
        }
    }
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            temp[i][j]=A[i][j]-B[i][j];
        }
    }
    return temp;
}

float **addMatrix(float **A,float **B,float **temp, int width, int height)
{
    for(int l=0;l<width;l++)
    {
        for(int m=0;m<height;m++)
        {
            temp[l][m]=0;
        }
    }
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            temp[i][j]=A[i][j]+B[i][j];
        }
    }
    return temp;
}

float **multiplyScalaire(float **A, float alpha,float **temp,int width, int height)
{
    for(int l=0;l<width;l++)
    {
        for(int m=0;m<height;m++)
        {
            temp[l][m]=0;
        }
    }

    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            temp[i][j]=alpha*A[i][j];

        }
    }
   return temp;

}

float psMatrix(float **A,float **B, int width, int height)
{
    float ps=0;
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            ps=ps+A[i][j]*B[i][j];
        }
    }
    return ps;
}

void equalMatrix(float **A,float **B, int width, int height)
{
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
           A[i][j]=B[i][j];
        }
    }

}

void divide255(float **I,int width, int height)
{
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
           I[i][j]=I[i][j]/255;
        }
    }
}

void multiply255(float **I,int width, int height)
{
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
           I[i][j]=I[i][j]*255;
        }
    }
}

void complementMask(float **Mask,int width, int height)
{
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            if(Mask[i][j]==0){
           Mask[i][j]=1;
            }else
            {
                Mask[i][j]=0;
            }
        }
    }
}


// to do ... IF Q2
void InpaintingBW(float **Iout, float **Iin, float **Mask, int width, int height, double param)
{

   Question1(Iout,0,0,Iin,0,0,Mask,width,height,param);
   divide255(Iout,width,height);
   complementMask(Mask,width,height);

   float res=100;
   float alphak=0;
   float Bk=0;
    float **xk = 0;
    float **rk =0;
    float **rki = 0;
    float **dk = 0;
    float **temp=0;
    float **temp2=0;
    float **temp3=0;
    float **ones=0;
     xk  = AllocateFloatArray( width,height);
     rk = AllocateFloatArray( width,height);
     rki  = AllocateFloatArray( width,height);
     dk = AllocateFloatArray( width,height);
     temp = AllocateFloatArray( width,height);
     temp2 = AllocateFloatArray( width,height);
     temp3 = AllocateFloatArray( width,height);

    equalMatrix(xk,Iout,width,height);

   temp2=matriceA(xk,Mask,temp2,width,height);
   temp=subtractMatrix(Iout,temp2,temp,width,height);
   equalMatrix(rk,temp,width,height);


    equalMatrix(dk,rk,width,height);


   while(res>14)
    {

         temp2=matriceA(dk,Mask,temp2,width,height);
         alphak=psMatrix(rk,rk,width,height) / psMatrix(temp2,dk,width,height);

         temp2=multiplyScalaire(dk,alphak,temp2,width,height);
         temp=addMatrix(xk,temp2,temp,width,height);
         equalMatrix(xk,temp,width,height);

         equalMatrix(rki,rk,width,height);

         temp2=matriceA(dk,Mask,temp2,width,height);
         temp3=multiplyScalaire(temp2,alphak,temp3,width,height);
         temp=subtractMatrix(rk,temp3,temp,width,height);
         equalMatrix(rk,temp,width,height);

         Bk=psMatrix(rk,rk,width,height) / psMatrix(rki,rki,width,height);
         temp2=multiplyScalaire(dk,Bk,temp2,width,height);
         temp=addMatrix(temp2,rk,temp,width,height);
         equalMatrix(dk,temp,width,height);

     res=psMatrix(rk,rk,width,height) / abs(width*height);

    }



equalMatrix(Iout,xk,width,height);
multiply255(Iout,width,height);

    DestroyFloatArray(&xk, width);
    DestroyFloatArray(&rk, width);
    DestroyFloatArray( &rki, width);
    DestroyFloatArray( &dk, width);
    DestroyFloatArray( &temp, width);
    DestroyFloatArray( &temp2, width);

}


// to do ... IF Q3
void InpaintingColor(float **Rout, float **Gout, float **Bout, float **Rin, float **Gin, float **Bin, float **Mask, int width, int height, double param)
{

   Question1(Rout,Gout,Bout,Rin,Gin,Bin,Mask,width,height,param);
    float ***Iout=new float**[3];
   Iout[0]=Rout;
   Iout[1]=Gout;
   Iout[2]=Bout;

   for(int i=0;i<3;i++){
   divide255(Iout[i],width,height);
   }
   complementMask(Mask,width,height);

    float **xk = 0;
    float **rk =0;
    float **rki = 0;
    float **dk = 0;
    float **temp=0;
    float **temp2=0;
    float **temp3=0;

     xk  = AllocateFloatArray( width,height);
     rk = AllocateFloatArray( width,height);
     rki  = AllocateFloatArray( width,height);
     dk = AllocateFloatArray( width,height);
     temp = AllocateFloatArray( width,height);
     temp2 = AllocateFloatArray( width,height);
     temp3 = AllocateFloatArray( width,height);

for(int k=0;k<3;k++){
    float res=100;
    float alphak=0;
    float Bk=0;
     equalMatrix(xk,Iout[k],width,height);

    temp2=matriceA(xk,Mask,temp2,width,height);
    temp=subtractMatrix(Iout[k],temp2,temp,width,height);
    equalMatrix(rk,temp,width,height);
    equalMatrix(dk,rk,width,height);

    while(res>14)
     {

          temp2=matriceA(dk,Mask,temp2,width,height);
          alphak=psMatrix(rk,rk,width,height) / psMatrix(temp2,dk,width,height);

          temp2=multiplyScalaire(dk,alphak,temp2,width,height);
          temp=addMatrix(xk,temp2,temp,width,height);
          equalMatrix(xk,temp,width,height);

          equalMatrix(rki,rk,width,height);

          temp2=matriceA(dk,Mask,temp2,width,height);
          temp3=multiplyScalaire(temp2,alphak,temp3,width,height);
          temp=subtractMatrix(rk,temp3,temp,width,height);
          equalMatrix(rk,temp,width,height);

          Bk=psMatrix(rk,rk,width,height) / psMatrix(rki,rki,width,height);
          temp2=multiplyScalaire(dk,Bk,temp2,width,height);
          temp=addMatrix(temp2,rk,temp,width,height);
          equalMatrix(dk,temp,width,height);


      res=psMatrix(rk,rk,width,height) / abs(width*height);



     }

    equalMatrix(Iout[k],xk,width,height);
}



    for(int i=0;i<3;i++){
    multiply255(Iout[i],width,height);
    }
    Rout=Iout[0];
   Gout= Iout[1];
    Bout=Iout[2];


     DestroyFloatArray(&xk, width);
     DestroyFloatArray(&rk, width);
     DestroyFloatArray( &rki, width);
     DestroyFloatArray( &dk, width);
     DestroyFloatArray( &temp, width);
     DestroyFloatArray( &temp2, width);

}

