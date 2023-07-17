#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
 {
    if(x<0)
    x=0;
    else if(x>=im.w)
    x=im.w-1;
    if(y<0)
    y=0;
    else if(y>=im.h)
    y=im.h-1;
   int p= x + y*im.w +c*im.w*im.h;
   return im.data[p];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if(x<0||x>im.w||y<0||y>im.h||c<0||c>im.c)
    return;
    int k= x + y*im.w + c*im.w*im.h;
      im.data[k]=v;
    
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    int k=im.w*im.h*im.c;
    for(int i=0;i<k;i++)
    copy.data[i]=im.data[i];
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    for(int i=0;i<im.w*im.h;i++){
        float r=im.data[i];
        float g=im.data[i+im.h*im.w];
        float b=im.data[i+2*im.h*im.w];
        float y= (0.299*r + 0.587*g + 0.114*b);
         gray.data[i]=y;
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    for(int i=0;i<im.w;i++){
        for(int j=0;j<im.h;j++){
            int k= i + j*im.w + c*im.h*im.w;
            im.data[k]+=v;
        }
    }
}

void clamp_image(image im)
{
    for(int i=0;i<im.h*im.w*im.c;i++){
        if(im.data[i]<0)
        im.data[i]=0;
        else if(im.data[i]>1)
        im.data[i]=1;
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for(int i=0;i<im.h*im.w;i++)
    {
        float r=im.data[i];
        float g=im.data[i+im.h*im.w];
        float b=im.data[i+2*im.h*im.w];
        float V= three_way_max(r,g,b);
        float m= three_way_min(r,g,b);
        float C= V-m;
        float S,H,H1=0;
        if(r==0&&g==0&&b==0)
          S=0;
          else
          S=C/V;

          if(C==0)
           H1=0;
           else if(V==r)
           H1=(g-b)/C;
           else if(V==g)
           H1=(b-r)/C +2;
           else if(V==b)
           H1=(r-g)/C +4;

           if(H1<0)
           H=(H1/6)+1;
           else
           H=(H1/6);

           im.data[i]=H;
           im.data[i+im.h*im.w]=S;
           im.data[i+2*im.h*im.w]=V; 

    }
}

void hsv_to_rgb(image im)
{
    for(int i=0;i<im.w*im.h;i++)
    {
        float H=im.data[i];
        float S=im.data[i+im.h*im.w];
        float V=im.data[i+2*im.h*im.w];
        float C= V*S;
        float H1=H*6;
        float m=V-C;
        float j=floor(H1);
        float f= H1-j;
        float n=V-C*f;
        float n1=V-C*(1-f);
        float r,g,b;
         if(H1>=0&&H1<1){
            r=V;
            g=n1;
            b=m;
         }
        else if(H1>=1&&H1<2){
            r=n;
            g=V;
            b=m;
         }
        else if(H1>=2&&H1<3){
            r=m;
            g=V;
            b=n1;
         }
        else if(H1>=3&&H1<4){
            r=m;
            g=n;
            b=V;
         }
        else if(H1>=4&&H1<5){
            r=n1;
            g=m;
            b=V;
         }
        else {
            r=V;
            g=m;
            b=n;
         }
        im.data[i]=r;
        im.data[i+im.h*im.w]=g;
        im.data[i+2*im.h*im.w]=b;
}
}


void scale_image(image im, int c, float v)
{
    for(int i=0;i<im.w;i++){
        for(int j=0;j<im.h;j++){
            int k= i + j*im.w + c*im.h*im.w;
            im.data[k]*=v;
        }
    }
}

void rgb_to_hcl(image im)
{
     for(int i=0;i<im.h*im.w;i++)
    {
        float r=im.data[i];
        float g=im.data[i+im.h*im.w];
        float b=im.data[i+2*im.h*im.w];
           //Gamma Decompression
           if (r <= 0.04045)
           r = r/ 12.92;
           else
           r = pow(((r + 0.055) / 1.055) , 2.4);
           if (g <= 0.04045)
           g = g / 12.92;
           else
           g = pow(((g + 0.055) / 1.055), 2.4); 
           if (b <= 0.04045)
           b = b / 12.92;
           else
           b = pow(((b + 0.055) / 1.055),2.4) ; 
             //RGB TO CIEXYZ
            float X = r * 0.4124564f + g * 0.3575761f + b * 0.1804375f;
            float Y = r * 0.2126729f + g * 0.7151522f + b * 0.0721750f;
            float Z = r * 0.0193339f + g * 0.1191920f + b * 0.9503041f;
              //CIEXYZ TO CIELUV
        float Yn = 1.00000;
        float u1 =(4*X)/(X + 15*Y + 3*Z);
        float v1 =(9*Y)/(X + 15*Y + 3*Z);
        float L;
        if(Y/Yn<=pow(6.0/29,3))
        L=pow(29.0/3,3)*Y/Yn;
        else
        L=116*pow(Y/Yn,1.0/3)-16;
        float u=13*L*(u1-0.2009);
        float v=13*L*(v1-0.4610);
        
        //CIELUV TO HCL
        float H = atan2(v,u)*(180.0/ 3.1416);
        if(H<0)
        H+=360.0;
        float C = sqrt(u * u + v * v);
        L=L;
         
        im.data[i]=H;
        im.data[i+im.h*im.w]=C;
        im.data[i+2*im.h*im.w]=L;
        
        
    }
}

void hcl_to_rgb(image im)
{
        for(int i=0;i<im.h*im.w;i++)
    {
        float H=im.data[i];
        float C=im.data[i+im.h*im.w];
        float L=im.data[i+2*im.h*im.w];
        
        //HCL TO CIELUV
        L=L;
        float p=H*(3.14159265358/180.0);
        float U=cos(p)*C;
        float V=sin(p)*C;
         
        //CIELUV TO CIEXYZ
         float u1=U/(13*L) + 0.2009;
         float v1=V/(13*L) + 0.4610;
         float X,Y,Z;
         float Yn=1.0000;
         if(L<=8)
         Y=Yn*L*pow(3.0/29,3);
         else
         Y=Yn*powf((L+16)/116.0,3);
         X=Y*((9*u1)/(4*v1));
         Z=Y*((12-3*u1-20*v1)/(4*v1));
        
         //CIEXYZ TO RGB
           float r1 = 3.2404542f * X - 1.5371385f * Y - 0.4985314f * Z;
           float g1 = -0.9692660f * X + 1.8760108f * Y + 0.0415560f * Z;
           float b1 = 0.0556434f * X - 0.2040259f * Y + 1.0572252f * Z;
         //Gamma Compression
         float r,g,b;
         if (r1<= 0.0031308)
         r= r1 * 12.92;
         else
         r = pow(r1, 1.0 / 2.4) * 1.055 - 0.055;
        if (g1<= 0.0031308)
         g= g1 * 12.92;
         else
         g = pow(g1, 1.0 / 2.4) * 1.055 - 0.055;
         if (b1<= 0.0031308)
         b= b1 * 12.92;
         else
         b = pow(b1, 1.0 / 2.4) * 1.055 - 0.055;

          im.data[i]=r;
          im.data[i+im.h*im.w]=g;
          im.data[i+2*im.h*im.w]=b;
    }
}