//author:  Ben Hixon
//course:  CSCI 493.69
//program:  HW4, P2
//description:  Compute directions and intensities of light sources.
//		First we need to get the surface normal for each point (x,y) on the image.
//		We derive the following formula:
//			A sphere is given by the equation (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2.
//			We know that the surface normal at a point z(x,y) is 
//				N = (p,q,1), where p and q are the partial derivatives of x and y.
//			Therefore we have:
//
//			(1)	(x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2
//			(2)	z = sqrt ( r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) ) + z0
//			(3)	dzdx = -1 * (x-x0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) ); 
//			(4)     dzdy = -1 * (y-y0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) );
//			(5)	N = (p,q,1) = (dzdx,dzdy,1);
//
//			For a Lambertian surface, we can assume that the brightest point is the direction
//			of the light source, so we find the highest-magnitude pixel in the image and
//			compute the normal at that point to get the light-source direction.
//
//			We scale the direction vector so that its length is equal to the magnitude of that
//			pixel by normalizing and multiplying each component of the resulting 
//			unit vector by the pixel's magnitude.
//
//			The output file consists of three lines, each of the form "x1 y1 z1" for images 1, 2, and 3.//
//usage:
//			to compile: g++ -o s2 hw4p1.cpp Image.cpp Pgm.cpp
//			to run: ./s2 parameters.txt img1.pgm img2.pgm img3.pgm directions.txt

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Image.h"#include <vector>
#include <cmath>

using namespace std;

void get_directions(Image* im1, Image* im2, Image* im3, char* param_filename, char* out_filename);

int main(int argc, char ** argv){

	char* param_filename;
	char* im1_filename;
	char* im2_filename;
	char* im3_filename;
	char* directs_filename;

	if(argc!=6){
		cout<<"USAGE: ./s2 parameters.txt img1.pgm img2.pgm img3.pgm directions.txt\n";
		exit(1);
	}else{
				
		param_filename = *(argv+1);
		im1_filename = *(argv+2);
		im2_filename = *(argv+3);
		im3_filename = *(argv+4);
		directs_filename = *(argv+5);	
	}

	Image *im1 = new Image;
	if ( readImage(im1,im1_filename)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	Image *im2 = new Image;
	if ( readImage(im2,im2_filename)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	Image *im3 = new Image;
	if ( readImage(im3,im3_filename)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	get_directions(im1, im2, im3, param_filename, directs_filename);
}

void get_directions(Image* im1, Image* im2, Image* im3, char* param_filename, char* out_filename)
//(1) get x0,y0,r from param_filename
//(2) find x and y of highest pixel-value coordinate for each im
//(3) set x, y, and z for each im using 
//(1)	(x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2
//(2)	z = sqrt ( r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) ) + z0
//(3)	dzdx = -1 * (x-x0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) ); 
//(4)     dzdy = -1 * (y-y0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) );
//(5)	n = (p,q,1) = (dzdx,dzdy,1);
//(3) print to output file three lines each of form "nx ny nz"
{
	int x=1, y=1;
	double nx=0.0, ny=0.0, nz=0.0;
	double x0,y0,r;
	//double dzdx=0.0, dzdy=0.0;
	vector<double> coords;

	//read in x0,y0,r;
	ifstream infile(param_filename);
	infile>>x0>>y0>>r;

	int ncols=im1->getNCols();
	int nrows=im1->getNRows();

	//get max x,y for im1
	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
			if( im1->getPixel(i,j) > im1->getPixel(x,y) )
			{
				x=i;
				y=j;
			}

	//cout<<"im1's highest-mag pixel is at (x,y) = "<<x<<','<<y<<'\n';

	nx = -1 * (x-x0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) ); 
	ny = -1 * (y-y0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) );
	nz = 1;

	double mag = sqrt( nx*nx + ny*ny + nz*nz );
	double pVal = im1->getPixel(x,y);

	//normalize and scale to pixel value:
	nx = nx * pVal / mag;
	ny = ny * pVal / mag;
	nz = nz * pVal / mag;

	coords.push_back( nx );
	coords.push_back( ny );
	coords.push_back( nz );

	nx=ny=nz=0.0;
	x=y=1;

	//get max x,y for im2
	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
			if( im2->getPixel(i,j) > im2->getPixel(x,y) )
			{
				x=i;
				y=j;
			}

	//cout<<"im2's highest-mag pixel is at (x,y) = "<<x<<','<<y<<'\n';

	nx = -1 * (x-x0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) ); 
	ny = -1 * (y-y0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) );
	nz = 1;

	mag = sqrt( nx*nx + ny*ny + nz*nz );
	pVal = im2->getPixel(x,y);

	//normalize and scale to pixel value:
	nx = nx * pVal / mag;
	ny = ny * pVal / mag;
	nz = nz * pVal / mag;

	coords.push_back( nx );
	coords.push_back( ny );
	coords.push_back( nz );

	nx=ny=nz=0.0;
	x=y=1;

	//get max x,y for im3
	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
			if( im3->getPixel(i,j) > im3->getPixel(x,y) )
			{
				x=i;
				y=j;
			}

	//cout<<"im3's highest-mag pixel is at (x,y) = "<<x<<','<<y<<'\n';

	nx = -1 * (x-x0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) ); 
	ny = -1 * (y-y0) / sqrt(r*r - (x-x0)*(x-x0) - (y-y0)*(y-y0) );
	nz = 1;

	mag = sqrt( nx*nx + ny*ny + nz*nz );
	pVal = im3->getPixel(x,y);

	//normalize and scale to pixel value:
	nx = nx * pVal / mag;
	ny = ny * pVal / mag;
	nz = nz * pVal / mag;

	coords.push_back( nx );
	coords.push_back( ny );
	coords.push_back( nz );

	//print to file
	ofstream ofile(out_filename);
	if( ofile.is_open() )
	{		
		for(int i=0; i<9; i++)
		{
			if( (i+1)%3 == 0)
				ofile<<coords[i]<<'\n';
			else
				ofile<<coords[i]<<' ';
		}
		ofile.close();
	}
	
	else cout<<"Error, can't open file"<<out_filename<<endl;
}
