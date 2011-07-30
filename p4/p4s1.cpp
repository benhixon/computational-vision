//author:  Ben Hixon
//course:  493.69
//program:  HW4, P1
//description:	Take an image of a sphere and find its center and radius
//		We find its center by thresholding and computing the centroid.
//		We find the radius by first letting 
//			diameter = (topmost-bottommost)+(rightmost-leftmost)/2
//		Then
//			radius = diameter/2
//		The output file contains a single line that looks like "x_coord y_coord radius"
//usage:
//	to compile:	g++ -o s1 hw4p1.cpp Image.cpp Pgm.cpp
//	to run:		./s1 input.pgm t parameters.txt

#include <iostream>
//#include <sstream>
#include <fstream>
#include <cstdlib>
#include "Image.h"
//#include <vector>
//#include <cmath>

using namespace std;

void get_parameters(Image* im, char* out_filename);
void make_binary(Image* im, int t);

int main(int argc, char ** argv){

	char* in_filename;
	int threshold;
	char* out_filename;

	if(argc!=4){
		cout<<"USAGE: ./s1 input.pgm t parameters.txt\n";
		exit(1);
	}else{
				
		in_filename = *(argv+1);
		threshold = atoi( *(argv+2) );
		out_filename = *(argv+3);
	}

	Image *im = new Image;
	if ( readImage(im,in_filename)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	make_binary(im, threshold);
	get_parameters(im, out_filename);
}

void make_binary(Image* im, int t)
{
	int ncols=im->getNCols();
	int nrows=im->getNRows();

	for(int i=0; i<nrows; i++)
		for(int j=0; j<ncols; j++)
		{
			
			if( im->getPixel(i,j)<t )
				im->setPixel(i,j,0);
			else
				im->setPixel(i,j,255);
		}
}

void get_parameters(Image* im, char* out_filename)
//(1) compute centroid
//(2) get radius
//(3) print centroid and radius to out_filename as "x0 y0 r"
//NOTE: to keep x and i together, we let the x-axis be the left-hand side of the image.
{

	int ncols=im->getNCols();

	int nrows=im->getNRows();
	double x_bar=0.0, y_bar=0.0, radius=0.0;
	int A=0;

	//(1) get area A by counting total number of non-zero pixels
	//(2) get x_coord by counting total x_coords and dividing by A
	//(3) get y_coord by counting total y_coords and dividing by A

	for(int i=1; i<nrows-1; i++)

		for(int j=1; j<ncols-1; j++)

			if( im->getPixel(i,j)>0 )
			{
				A++;
				x_bar = x_bar + i;

				y_bar = y_bar + j;
			}
	
	x_bar = x_bar / float(A);	//divide the integral by the area to get the coordinate

	y_bar = y_bar / float(A);

	int min_x=nrows;	//top = lowest row containing a nonzero pixel)
	int max_x=0;		//bottom = highest row containing a nonzero pixel
	int min_y=ncols;	//left = lowest column containing a nonzero pixel
	int max_y=0;		//right = highest column containing a nonzero pixel

	for(int i=1; i<nrows-1; i++)

		for(int j=1; j<ncols-1; j++)

			if( im->getPixel(i,j)>0 )
			{	
				if(j<min_y){ min_y=j; }
				if(j>max_y){ max_y=j; }
				if(i<min_x){ min_x=i; }
				if(i>max_x){ max_x=i; }
			}

	radius = 0.25 * (max_x-min_x + max_y-min_y);

	ofstream ofile(out_filename);

	if( ofile.is_open() )

	{

		ofile<<x_bar<<' '<<y_bar<<' '<<radius;
		ofile.close();
	
	}
	
	else cout<<"Error, can't open file"<<out_filename<<endl;
}
