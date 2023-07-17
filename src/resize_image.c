#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    int x1 = round(x);
    int y1 = round(y);
    float pixel = get_pixel(im,x1,y1,c);
    return pixel;
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
   image new_im = make_image(w, h, im.c);
    float a = (float)im.w / w;
    float b = (float)im.h / h;

    for (int k = 0; k < new_im.c; k++)
    {
        for (int j = 0; j < new_im.h; j++)
        {
            for (int i = 0; i < new_im.w; i++)
            {
                float x = (i + 0.5) * a - 0.5;
                float y = (j + 0.5) * b - 0.5;
                new_im.data[i + j*w + k*w*h] = nn_interpolate(im, x, y, k);   
            }
        }
    }
    return new_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    float left,right,top,bottom;
    left = floor(x);
    right = ceil(x);
    top = floor(y);
    bottom = ceil(y);
    float d1,d2,d3,d4;
    d1 = get_pixel(im,left,top,c);
    d2 = get_pixel(im,right,top,c);
    d3 = get_pixel(im,left,bottom,c);
    d4 = get_pixel(im,right,bottom,c);
    float e1,e2,e3,e4;
    e1 = bottom - y;
    e2 = y - top;
    e3 = x - left;
    e4 = right - x;
    float f1,f2,f3;
    f1 = d1*e1 + d3*e2;
    f2 = d2*e1 + d4*e2;
    f3 = f2*e3 + f1*e4;
    return f3;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
     image new_im = make_image(w, h, im.c);
    float a = (float)im.w / w;
    float b = (float)im.h / h;
    for (int k = 0; k < new_im.c; k++)
    {
        for (int j = 0; j < new_im.h; j++)
        {
            for (int i = 0; i < new_im.w; i++)
            {
                float x = (i + 0.5) * a - 0.5;
                float y = (j + 0.5) * b - 0.5;
                new_im.data[i + j*w + k*w*h] = bilinear_interpolate(im, x, y, k);   
            }
        }
    }
    return new_im;
}

