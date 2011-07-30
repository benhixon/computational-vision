//author:  Ben Hixon
//course:  CSCI 493.69
//program:  HW4, P4
//description:  Compute the albedo for every pixel visible from all three light sources.
//		We use the same technique used in P3 to find rho*n, where rho is the
//		albedo, except this time we do not have a step size so we compute it for
//		every pixel above the given threshold, and we scale the albedos to fit
//		within the range 0 to 255 and show that value in the output image.
//
//	Basically everything is the same as in s3 except for draw_albedo.
//
//usage:
//		to compile: g++ -o s4 hw4p1.cpp Image.cpp Pgm.cpp
//		to run: ./s4 directions.txt img1.pgm img2.pgm img3.pgm threshold albedomap.pgm


#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "Image.h"
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> vec;
typedef vector<vec> mat;

//to make these functions really simple I'll assume vectors are 3x1 and matrices are 3x3.
mat invert(mat A);
vec matTimesVec(mat A, vec v);

void get_albedos(Image* im1, Image* im2, Image* im3, char* directs_filename, int threshold);
void draw_albedos(Image* im, vector<vec> albedos);

int main(int argc, char ** argv){

	char* directs_filename;
	char* im1_filename;
	char* im2_filename;
	char* im3_filename;
	//int step;
	int threshold;
	char* out_filename;

	if(argc!=7){
		cout<<"USAGE: ./s3 directions.txt img1.pgm img2.pgm img3.pgm threshold albedo_output.pgm\n";
		exit(1);
	}else{
				
		directs_filename = *(argv+1);
		im1_filename = *(argv+2);
		im2_filename = *(argv+3);
		im3_filename = *(argv+4);
		//step = atoi( *(argv+5) );
		threshold = atoi( *(argv+5) );
		out_filename = *(argv+6);	
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
	get_albedos(im1, im2, im3, directs_filename, threshold);

	if ( writeImage(im1,out_filename)!=0 ){
		//im1 has needles on it, write out as output.pgm
		cout<<"ERROR: Can't write to output\n";
 
		exit(1);
	
	}
}

void get_albedos(Image* im1, Image* im2, Image* im3, char* directs_filename, int threshold)
//	(1) First get the vectors s1,s2,s3 from directs_filename.  
//	(2) Normalize them (set to unit vectors of themselves.)
//	(3) Then for every pixel above the threshold:
//
//	(A) set E1,E2,E3 to the pixel values in im1,im2,im3 respectively.
//	(B) find rho_times_n at the given pixel using the formula:
//		rho*n = S^-1 * E, where S is a matrix with rows s1, s2, s3, and E is the vector with components E1,E2,E3
//	(C) Having found rho*n, get length and set rho to it
//	(D) push (x_coord, y_coord), rho onto albedo vector
//
//	(4) After all pixels are processed, send the albedo vector to draw_albedo()
{

	vector<vec> albedos;	//every two entries are <x.y> followed by rho
	vec s1 (3,0);		//light source directions.  initialize to (0,0,0)
	vec s2 (3,0);
	vec s3 (3,0);
	mat S;		//matrix whose rows are s1,s2,s3
	int E1,E2,E3;	//pixel magnitudes
	vec E;		//vector whose components are E1,E2,E3
	vec rho_n;  	//albedo times surface normal: for each pixel, rho*n = S^-1 * E
	//vec unit_n;	//surface normal
	double rho;	//albedo = length of rho_n

	int ncols=im1->getNCols();
	int nrows=im1->getNRows();
	
	//(1) First get the vectors s1,s2,s3 from directs_filename

	ifstream infile(directs_filename);
	infile>>s1[0]>>s1[1]>>s1[2];
	infile>>s2[0]>>s2[1]>>s2[2];
	infile>>s3[0]>>s3[1]>>s3[2];

	
	//(2) normalize.  only three vectors so I just wrote it out without loops
	/*
	double mag = sqrt(s1[0]*s1[0] + s1[1]*s1[1] + s1[2]*s1[2]);
	s1[0]=s1[0]/mag;
	s1[1]=s1[1]/mag;
	s1[2]=s1[2]/mag;
	mag = sqrt(s2[0]*s2[0] + s2[1]*s2[1] + s2[2]*s2[2]);
	s2[0]=s2[0]/mag;
	s2[1]=s2[1]/mag;
	s2[2]=s2[2]/mag;
	mag = sqrt(s3[0]*s3[0] + s3[1]*s3[1] + s3[2]*s3[2]);
	s3[0]=s3[0]/mag;
	s3[1]=s3[1]/mag;
	s3[2]=s3[2]/mag;
	*/

	//add light source vectors to matrix S so we can invert it
	S.push_back(s1);
	S.push_back(s2);
	S.push_back(s3);

	mat S_inverse = invert(S);

	//(3) Then for every pixel above threshold:
	//(A) set E1,E2,E3 to the pixel values in im1,im2,im3 respectively.
	//(B) find rho_times_n at the given pixel using the formula:
	//	rho*n = S^-1 * E, where S is a matrix with rows s1, s2, s3, and E is the vector with components E1,E2,E3
	//(C) Having found rho*n, get length and set rho to it
	//(D) push (x_coord, y_coord), rho onto albedo vector
	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
		{
			E1 = im1->getPixel(i,j);
			E2 = im2->getPixel(i,j);
			E3 = im3->getPixel(i,j);

			if(E1>=threshold && E2>=threshold && E3>=threshold)
			{
				E.clear();
				E.push_back(E1);
				E.push_back(E2);
				E.push_back(E3);

				rho_n = matTimesVec( S_inverse,E );
				
				rho = sqrt( rho_n[0]*rho_n[0] + rho_n[1]*rho_n[1] + rho_n[2]*rho_n[2] );

				vec rhoVec;
				rhoVec.push_back(rho);				

				vector<double> coords;
				coords.push_back(i);
				coords.push_back(j);
				albedos.push_back(coords);
				albedos.push_back(rhoVec);
			}
		}
		
	draw_albedos(im1,albedos);
}

mat invert(mat A)
//| a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
//| a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
//| a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |
//	with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)
//The above equation is a common linear algebra formula, but I copied it directly from
//	http://www.dr-lex.be/random/matrix_inv.html
//I coded it myself though so any errors are mine.
{
	mat B;
	double a11=A[0][0],a12=A[0][1],a13=A[0][2];
	double a21=A[1][0],a22=A[1][1],a23=A[1][2];
	double a31=A[2][0],a32=A[2][1],a33=A[2][2];

	double det = a11*(a33*a22-a32*a23)-a21*(a33*a12-a32*a13)+a31*(a23*a12-a22*a13);

	vec row;
	row.push_back( (1/det)*(a33*a22-a32*a23) );
	row.push_back( (1/det)*-1*(a33*a12-a32*a13) );
	row.push_back( (1/det)*(a23*a12-a22*a13) );
	B.push_back(row);
	row.clear();
	row.push_back( (1/det)*-1*(a33*a21-a31*a23) );
	row.push_back( (1/det)*(a33*a11-a31*a13) );
	row.push_back( (1/det)*-1*(a23*a11-a21*a13) );
	B.push_back(row);
	row.clear();
	row.push_back( (1/det)*(a32*a21-a31*a22) );
	row.push_back( (1/det)*-1*(a32*a11-a31*a12) );
	row.push_back( (1/det)*(a22*a11-a21*a12) );
	B.push_back(row);

	return B;
}

vec matTimesVec(mat A, vec v)
//assume A=3x3, v=3x1, though this is obviously easily generalizable
{
	vec x;
	x.push_back( A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2] );
	x.push_back( A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2] );
	x.push_back( A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2] );

	return x;
}

void draw_albedos(Image* im1, vector<vec> albedos)//get max_rho.  for every rho, (255/max_rho)*rho = new pixel value
{
	int ncols=im1->getNCols();
	int nrows=im1->getNRows();
	int max_rho=0;
	double scale_factor = 0.0;

	for(int i=0; i<nrows; i++)
		for(int j=0; j<ncols; j++)
			im1->setPixel(i,j,0);	//make everything black

	int n = albedos.size() / 2;

	for(int i=0; i<n; i++)
		if( int(albedos[i*2+1][0]) > max_rho)
			max_rho = int(albedos[i*2+1][0]);

	scale_factor = 255/double(max_rho);


	for(int i=0; i<n; i++)
	{
		int x = int(albedos[i*2][0]);
		int y =	int(albedos[i*2][1]);

		int pixel_val = int( albedos[i*2+1][0]*scale_factor );
		im1->setPixel(x,y,pixel_val);


	}	
}
	


