//author:  Ben Hixon
//date:    10-12-2010
//course:  493.69 (Computer Vision)
//description:  HW3, program 1: First we find the edge points of the image,
//		so that we can run the Hough transform on them and find lines.
//		Edges are the locations of high contrast, so we want points
//		where the derivative is high.
//
//		The purpose of this program is to create an 'edge image' of
//		a grayscale image.  The edge will be assigned a magnitude, and
//		the intensity of the points on this edge will be some integer
//		from 0 to 255 proportional to the magnitude.
//		We can use the squared gradient, which is
//		  S(x,y) = (I_x)^2 + (I_y)^2 
//			 = [(B2-B1)*d(x*sin(theta)-y*cos(theta)+rho)]^2
//		where I_x and I_y are the first derivatives, and d() is the dirac.
//		Then the magnitude is the square root of S, sqrt(S(x,y)).
//
//		To estimate the squared gradient, we use the mask by Sobel, which
//		combines via convolution the mean filter and first derivative filter:
//		  delta1 = gradx = {(-1,0,1),(-2,0,2),(-1,0,1)}
//		  delta2 = grady = {(1,2,1),(0,0,0),(-1,-2,-1)}.
//		so that if sqrt(gradx*gradx + grady*grady) > threshold, we have an edge.
//		These two 3x3 matrices are convolution filters representing the gradients
//		in the x and y directions, so for each point x,y in the image, we apply
//		the convolution filters, get the results, square them, add them together,
//		and ask if they're bigger than a threshold.  If so, then we save that point.
//
//		Since it's a 3x3 mask, we can't operate on the borders of the image, but they're
//		won't be an edge there anyway.  Convolution is flip, shift, multiply, but since
//		these filters are symmetric, flipping doesn't do anything, so it's just multiply
//		for every point (x,y) on the image.
//
//		There's a function for convolution that takes a 3x3 mask, an Image, and an (x,y),
//		and returns the convolution.


#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Image.h"
using namespace std;

int convolve(int mask[][3], Image* I, int x, int y);
//convolves at the indicated pixel.
//assumes a 3x3 mask.

void findEdges(Image* im, Image* edgesIm);
//takes two images and stores in the second the edge intensity of each pixel in the first.
//uses squared gradient to find edge intensity by convolving with the sobel mask.
//actual value stored is (magnitude/max_magnitude)*255, where magnitude is the sqrt
//	of the squared gradient.

int main(int argc, char ** argv){

	char* in_filename;
	char* out_filename;

	if(argc!=3){
		cout<<"USAGE: ./h1 input_image.pgm output_image.pgm\n";
		exit(1);
	}else{
				
		in_filename = *(argv+1);
		out_filename = *(argv+2);
	}

	Image *im = new Image;
	if ( readImage(im,in_filename)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	Image* edgesIm = new Image;	//we'll store our edge magnitudes here
	readImage(edgesIm,in_filename);

	findEdges(im, edgesIm);

	if ( writeImage(edgesIm,out_filename)!=0 ){
		cout<<"ERROR: Can't write to output\n";
		exit(1);
	}
}

void findEdges(Image* im, Image* edgesIm)
//takes two images and stores in the second the edge intensity of each pixel in the first.
//uses squared gradient to find edge intensity by convolving with the sobel mask.
//actual value stored is (magnitude/max_magnitude)*255, where magnitude is the sqrt
//	of the squared gradient.
//we do not modify pixel values for the borders of the image.
{
	int gradx=0,grady=0,squared_grad=0,magnitude=0,max_magnitude=0;


	int S_x[3][3] = {{-1, 0, 1 },	//sobel operator for gradx
			 {-2, 0, 2 },
			 {-1, 0, 1 }};

	int S_y[3][3] = {{ 1, 2, 1 },	//sobel operator for grady
			 { 0, 0, 0 },
			 {-1,-2,-1 }};


	int ncols=im->getNCols();
	int nrows=im->getNRows();

	//we'll call convolve(S_x,im,x,y) and convolve(S_y,im,x,y) on
	//every (x,y) pair, and sum their squares to get the squared grad.
	//then we'll set the new pixel values proportional to these magnitudes.

	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
		{
			gradx=convolve(S_x,im,i,j);
			grady=convolve(S_y,im,i,j);
			squared_grad = gradx*gradx+grady*grady;
			magnitude=int( sqrt(squared_grad) );
			if(magnitude>max_magnitude){ max_magnitude=magnitude; }	 
				//we'll set max mag to 255 and make all other proportional in next loop	
			
			edgesIm->setPixel( i, j, magnitude );	//set pixel i,j in edgesIm to sqrt of its squared gradient
		}

	//second loop through to set proportional with max_magnitude == 255
	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
		{
			magnitude = edgesIm->getPixel(i,j);		
			int propMag = int( float(magnitude)/float(max_magnitude) * 255);	
		
			//int threshold = 75;		//could have put program #2 in these five lines.
			//if(propMag > threshold) 
			//	edgesIm->setPixel( i,j,propMag );
			//else 
			//	edgesIm->setPixel( i,j,0 );	

			edgesIm->setPixel( i,j,propMag );
		}

}


int convolve(int mask[][3], Image* I, int x, int y)
//convolution is just multiplying the values of each element in the mask with the
//corresponding neighbor to x,y and then adding up those values.
//So the result is
//	m[0][0]*I[x-1][y-1] + m[0][1]*I[x-1][y] + m[0][2]*I[x-1][y+1] + ....
//where the x values are the row values and the y values are the column values.
{
	int c=0;

	for(int i=-1;i<2; i++)
		for(int j=-1; j<2; j++)
			c = c + mask[i+1][j+1] * I->getPixel(x+i,y+j);
	
	return c;
}
