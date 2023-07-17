#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    for(int i = 0; i < im.c; i++) 
    {
        float s = 0;
        for(int j = 0; j < im.h; j++) 
        {
            for(int k = 0; k < im.w; k++) 
            {
                s += get_pixel(im, k, j, i);
            }
        }
        for(int j = 0; j < im.h; j++) 
        {
            for(int k = 0; k < im.w; k++) 
            {
                float pixel = get_pixel(im, k, j, i) / s;
                set_pixel(im, k, j, i, pixel);
            }
        }
    }
}

image make_box_filter(int w)
{
    // TODO
    image new_im =make_image(w,w,1);
    for(int i=0;i<w;i++)
    {
        for(int j=0;j<w;j++)
        {
            set_pixel(new_im, j, i, 0, 1.0);
        }
    }
    l1_normalize(new_im);
    return new_im;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    assert(filter.c == im.c || filter.c == 1);
  image new_im = make_image(im.w, im.h, im.c);
  for (int i = 0; i < im.w; i++) 
  {
    for (int j = 0; j < im.h; j++) 
    {
      for (int k = 0; k < im.c; k++) 
      {
        float q = 0;
        int f;
        if (filter.c != 1) 
        {
          f = k;
        }
        else 
        {
          f = 0;
        }
        for (int fx = 0; fx < filter.w; fx++) 
        {
          for (int fy = 0; fy < filter.h; fy++) 
          {
            float value = filter.data[fx + fy*filter.w + filter.w*filter.h*f];
            int f1 = i - filter.w / 2 + fx;
            int f2 = j - filter.h / 2 + fy;
            q += get_pixel(im, f1, f2, k) * value;
          }
        }
        new_im.data[i + j*im.w + k*im.w*im.h] = q;
      }
    }
  }
  if (preserve == 0) 
  {
    image im2 = make_image(im.w, im.h, 1);
    for (int i = 0; i < im.w; i++) 
    {
      for (int j = 0; j < im.h; j++) 
      {
        float q = 0;
        for (int k = 0; k < im.c; k++) 
        {
          q += new_im.data[i + j*im.w + im.w*im.h*k];
        }
        im2.data[i + j*im.w] = q;
      }
    }
    return im2;
  }
  else {
    return new_im;
  }
}

image make_highpass_filter()
{
    // TODO
     image new_im = make_box_filter(3);
      int a1[9] = {0,-1,0,-1,4,-1,0,-1,0};
    for(int i=0;i<9;i++)
    {
        new_im.data[i] = a1[i];
    }
    return new_im;
}

image make_sharpen_filter()
{
    // TODO
   image new_im = make_box_filter(3);
    int a1[9] = {0,-1,0,-1,5,-1,0,-1,0};
    for(int i=0;i<9;i++)
    {
        new_im.data[i] = a1[i];
    }
    return new_im;
}

image make_emboss_filter()
{
    // TODO
    image new_im = make_box_filter(3);
     int a1[9] = {-2,-1,0,-1,1,1,0,1,2};
    for(int i=0;i<9;i++)
    {
        new_im.data[i] = a1[i];
    }
    return new_im;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We should use highpass and sharpen filters and avoid emboss filter to preserve when we run our convolution.
//         Highpass and sharpen filters enhances edges and fine details and preserve the original details.
//         Emboss filter alters the image and introduce the stylized effect which makes it less appropriate for preserving original image details. 

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    int a = 6*sigma;
    if(a%2 == 0)
    a += 1;
    image new_im = make_box_filter(a);
    int c = (a-1)/ 2;
    for(int i = 0; i<a ; i++)
    {
        for(int j = 0; j<a ; j++)
        {
            float x = j - c;
            float y = i - c;
            float g = (1.0/(TWOPI * pow(sigma,2))) * exp((-1) * (pow(x,2) + pow(y,2)) / (2 * pow(sigma,2)));
            new_im.data[j + i*a] = g;
        }
    }
    return new_im;
}

image add_image(image a, image b)
{
    // TODO
    image new_im = make_image(a.w, a.h, a.c);
    for(int i=0;i<new_im.c*new_im.h*new_im.w;i++)
    {
        new_im.data[i] = a.data[i] + b.data[i];
    }
    return new_im;
}

image sub_image(image a, image b)
{
    // TODO
    image new_im = make_image(a.w, a.h, a.c);
    for(int i=0;i<new_im.c*new_im.h*new_im.w;i++)
    {
        new_im.data[i] = a.data[i] - b.data[i];
    }
    return new_im;
}

image make_gx_filter()
{
    // TODO
    image new_im = make_box_filter(3);
     int a1[9] = {-1,0,1,-2,0,2,-1,0,1};
    for(int i=0;i<9;i++)
    {
        new_im.data[i] = a1[i];
    }
    return new_im;
}

image make_gy_filter()
{
    // TODO
    image new_im = make_box_filter(3);
    int a1[9] = {-1,-2,-1,0,0,0,1,2,1};
    for(int i=0;i<9;i++)
    {
        new_im.data[i] = a1[i];
    }
    return new_im;
}

void feature_normalize(image im)
{
    // TODO
    float max = -1.0;
  float min = INFINITY;
  for (int i = 0; i < im.w; i++) 
  {
    for (int j = 0; j < im.h; j++) 
    {
      for (int k = 0; k < im.c; k++) 
      {
        if (im.data[i + j*im.w + k*im.w*im.h] > max) {
          max = im.data[i + j*im.w + k*im.w*im.h];
        }
        if (im.data[i + j*im.w + k*im.w*im.h] < min) {
          min = im.data[i + j*im.w + k*im.w*im.h];
        }
      }
    }
  }
  if (max - min != 0) 
  {
    for (int i = 0; i < im.w; i++) 
    {
      for (int j = 0; j < im.h; j++) 
      {
        for (int k = 0; k < im.c; k++) 
        {
          im.data[i + j*im.w + k*im.w*im.h] = (im.data[i + j*im.w + k*im.w*im.h] - min) / (max - min);
        }
      }
    }
  }
  else {
    for (int i = 0; i < im.w; i++) 
    {
      for (int j = 0; j < im.h; j++) 
      {
        for (int k = 0; k < im.c; k++) 
        {
          im.data[i + j*im.w + k*im.w*im.h] = 0;
        }
      }
    }
  }
}

image *sobel_image(image im)
{
    // TODO
  image *new_im = calloc(2, sizeof(image));
  image f_gx = make_gx_filter();
  image f_gy = make_gy_filter();
  image gx = convolve_image(im, f_gx, 0);
  image gy = convolve_image(im, f_gy, 0);
  new_im[0] = make_image(im.w, im.h, 1);
  new_im[1] = make_image(im.w, im.h, 1);
  for (int i = 0; i < im.w; i++) 
  {
    for (int j = 0; j < im.h; j++) 
    {
      float mag = sqrtf(gx.data[i + j*im.w]*gx.data[i + j*im.w] + gy.data[i + j*im.w]*gy.data[i + j*im.w]);
      float grad = atan2(gy.data[i + j*im.w], gx.data[i + j*im.w]);
      new_im[0].data[i + j*im.w] = mag;
      new_im[1].data[i + j*im.w] = grad;
    }
  }
  return new_im;
}

image colorize_sobel(image im)
{
    // TODO
    image new_im = make_image(im.w,im.h,im.c);
    image im1;
    im1 = make_gx_filter();
    image gx = convolve_image(im,im1,1);
    im1 = make_gy_filter();
    image gy = convolve_image(im,im1,1);

    for(int i=0;i<im.c;i++)
    {
        for(int j=0;j<im.h;j++)
        {
            for(int k=0;k<im.w;k++)
            {
                float gx1 = get_pixel(gx, j, i, 0);
                float gy1 = get_pixel(gy, j, i, 0);
                float v = sqrt(pow(gx1,2) + pow(gy1,2) + atan(gy1/gx1) + 1);
                set_pixel(new_im, k, j, i, v);
            }
        }
    }
    image p = make_gaussian_filter(2);
    new_im = convolve_image(new_im, p, 1);
    return new_im;
}
