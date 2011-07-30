//author:  Ben Hixon
//date:    9-28-2010
//course:  493.69 (Computer Vision)
//description:  HW2, program 1: convert a grayscale image to a binary image.
//		The program accepts three arguments: an input image, a gray-level threshold, 
//		and an output image.  It scans each pixel, and if less than the threshold 
//		sets its color to 0 (black), else sets it to 255 (white). 
//		The threshold that seems to work best for me is 125.
//use:  to compile: g++ -o p1 hw2p1.cpp Image.cpp Pgm.cpp
//	to run: ./p1 two_objects.pgm 125 binary.pgm

#include <iostream>
#include <cstdlib>
#include "Image.h"
using namespace std;

void makeBinary(Image* im, int t);

int main(int argc, char ** argv){

	char* in_filename;
	int threshold_val;
	char* out_filename;

	if(argc!=4){
		cout<<"USAGE: ./p1 grayscale.pgm threshold_val binary.pgm\n";
		exit(1);
	}else{
				
		in_filename = *(argv+1);
		threshold_val = atoi ( *(argv+2) );
		out_filename = *(argv+3);
	}

	Image *im = new Image;
	if ( readImage(im,in_filename)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	makeBinary(im, threshold_val);

	if ( writeImage(im,out_filename)!=0 ){
		cout<<"ERROR: Can't write to output\n";
		exit(1);
	}

}

void makeBinary(Image* im, int t)
{
	int ncols=im->getNCols();
	int nrows=im->getNRows();

	for(int i=0; i<nrows; i++)
	{
		for(int j=0; j<ncols; j++)
		{
			
			if( im->getPixel(i,j)<t )
				im->setPixel(i,j,0);
			else
				im->setPixel(i,j,1);
		}
	}

	im->setColors(1);	//set number of colors.  there's just one color here, white.

}




