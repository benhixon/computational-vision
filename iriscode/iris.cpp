//author: Ben Hixon  ( shixon@hunter.cuny.edu, ben.hixon@gmail.com )
//date 12/21/2010
//course: CSCI 493.69 (Computer Vision)
//instructor: Professor Stamos
//To compile:	g++ -o iris iris.cpp Image.cpp Pgm.cpp
//To run:	./iris iris1.pgm iris2.pgm


#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Image.h"
#include <cmath>
#include <complex>
#include <vector>
#include <string>
using namespace std;




void get_irisCode(Image* im, vector<int> &ic);
//Finds irisCode of im by sampling im to get r0,th0 and then calling Gabor() for this r0,th0 and varying values of
//omega. Adds a bit pair to the iriscode vector ic depending on which quadrant the phasor returned by Gabor() is in.

complex<double> Gabor(Image* im, double theta0, double r0, const vector<double> &center, double omega, double alpha, double beta);
//Given the image im, the patch sample (theta0,r0), the center, and wavelet parameters omega,alpha, and beta, find the Gabor transform,
//which is defined to be the integral of function given by the pixel value times the complex sinusoid 'carrier' function times the 
//gaussian 'envelope,' integrated all over the image.
//The Gaussian envelope depends on parameters alpha and beta, which therefore act like the window size of the Gaussian, while
//the wavelet depends on the value of omega, which varies from 0 to 2*3.14.
//I've hard-coded alpha and beta instead of varying them as Daugman does, since I don't really understand how alpha and beta inform
//the envelope and Daugman doesn't say what values he uses.  I do vary the wavelet's omega parameter though.

complex<double> gabor2(int ival,double om, double al, double be, double th0, double r0, double rho, double phi);
//Helper function for Gabor().  Returns the integrand of the Gabor function.

vector<double> get_center(Image* im);
//returns a vector with first value x0, second value y0, where (x0,y0) is the centroid of the pupil.

double get_HD(const vector<int> &ic1, const vector<int> &ic2);
//returns fractional Hamming Distance

int main(int argc, char ** argv){

	char* im1_fname;
	char* im2_fname;

	if(argc!=3){
		cout<<"USAGE: ./iris iris1.pgm iris2.pgm\n";
		exit(1);
	}else{
				
		im1_fname = *(argv+1);
		im2_fname = *(argv+2);
	}

	Image* im1 = new Image;
	if ( readImage(im1,im1_fname)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	Image* im2 = new Image;
	if ( readImage(im2,im2_fname)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	cout<<"Read "<<im1_fname<<" and "<<im2_fname<<endl;

	cout<<"Calculating IrisCode for "<<im1_fname<<endl;
	vector<int> ic1;
	get_irisCode(im1,ic1);

	//print ic1 to file
	string out_filename = string(im1_fname) + "iriscode.txt";
	ofstream ofile( out_filename.c_str() );
	if( ofile.is_open() )
	{	
		for(int i=0; i<ic1.size(); i++){ ofile<<ic1[i]; }
		ofile.close();
	}
	else cout<<"Error, can't open file"<<out_filename<<endl;


	cout<<"Calculating IrisCode for "<<im1_fname<<endl;
	vector<int> ic2;
	get_irisCode(im2,ic2);

	//print ic2 to file
	out_filename = string(im2_fname) + ".iriscode.txt";
	ofstream ofile2( out_filename.c_str() );
	if( ofile2.is_open() )
	{	
		for(int i=0; i<ic2.size(); i++){ ofile2<<ic2[i]; }
		ofile2.close();
	}
	else cout<<"Error, can't open file"<<out_filename<<endl;

	double HD = get_HD(ic1,ic2);
	cout<<"\nFractional Hamming Distance of "<<im1_fname<<" and "<<im2_fname<<" is: "<<HD<<'\n';	
}





void get_irisCode(Image* im, vector<int> &ic)
//
{
	vector<double> center = get_center(im);	
	cout<<"Center is "<<center[0]<<','<<center[1]<<endl;
	ic.clear();

	double theta0, r0,omega,alpha,beta;
	int i_xy;

	complex<double> G(0.0,0.0);

	int ncols=im->getNCols();
	int nrows=im->getNRows();

	int c=0;
	for(int i=0; i<nrows; i=i+15)
		for(int j=0; j<ncols; j=j+15)
		{
			i_xy = im->getPixel(i,j);
			if(i_xy > 50)
			{
				r0 = sqrt( (center[0]-i)*(center[0]-i) + (center[1]-j)*(center[1]-j) );
				theta0 = atan2(center[0]-i,center[1]-j);

				alpha = beta = 21;	//more gabor parameters (size of the wavelet)
							//i don't really know what happens if i vary these; they influence the size of the gaussian window
							//21 chosen randomly
				
				for(omega = 0.0; omega < 6.28; omega = omega + 0.75)
				//8 of these loops
				{
					c++;
					//cout<<"Getting bit pair "<<c<<"..."<<endl;
									
					G = Gabor(im, theta0, r0, center,omega,alpha,beta);
					if(c%500 == 0){ 
						cout<<"Gabor of bitpair "<<c<<" at patch (r0,th0,omega)=("<<r0<<','<<theta0<<','<<
								omega<<") is "<<G<<endl; 
					}
					//cout<<"Gabor of bitpair "<<c<<" at patch (r0,th0)=("<<r0<<','<<theta0<<") is "<<G<<endl;

					if( real(G) > 0 ){ ic.push_back(1); }
					else{ ic.push_back(0); }

					if( imag(G) > 0 ){ ic.push_back(1); }
					else{ ic.push_back(0); }
				}
			}
		}
}


complex<double> Gabor(Image* im, double theta0, double r0, const vector<double> &center, double omega, double alpha, double beta)		
//
{
	
	int ncols=im->getNCols();
	int nrows=im->getNRows();
	int i_xy = 0;
	complex<double> gabor (0.0,0.0);
	double rho, phi;
		
	for(int i=0; i<nrows; i++)
		for(int j=0; j<ncols; j++)
		{
			i_xy = im->getPixel(i,j);
			if(i_xy > 50)
			//threshold at 50 so as not to get pupil
			{
				//convert to polar coords
				rho = sqrt( (center[0]-i)*(center[0]-i) + (center[1]-j)*(center[1]-j) );
				phi = atan2(center[0]-i,center[1]-j);
				//cout<<"getting g with (r0,th0,omega,rho,phi)="<<r0<<','<<theta0<<','<<omega<<','<<rho<<','<<phi<<endl;
				gabor = gabor + gabor2(i_xy,omega,alpha,beta,theta0,r0,rho,phi);  //do i need rho drho dphi?
			}
		}


	return gabor;

}

complex<double> gabor2(int ival,double om, double al, double be, double th0, double r0, double rho, double phi)
//returns integrand of gabor transform 
//carrier is a complex exponential, for which e^(i*t) = cos(t) + i * sin(t)
{
	//for some reason these lines didn't work:	
	//complex<double> pwr_carr (0.0, -1 * om * (th0 - phi ) );
	//double pwr_env = -1 * (r0-rho)*(r0-rho)/(al*al) + (th0-phi)*(th0-phi)/(be*be);
	//return ival * exp(pwr_carr) * exp(pwr_env);

	//neither did:
	//double env = exp(-1 * (r0-rho)*(r0-rho)/(al*al) - (th0-phi)*(th0-phi)/(be*be) );
	//double t = -1 * om * (th0 - phi ) ;
	//complex<double> carr ( cos(t), sin(t) );
	//g_integrand = ival * carr * env;
	
	//complex is supposed to be overloaded for multiplication, but maybe not on our version of the standard library

	double pwr_env = -1 * (r0-rho)*(r0-rho)/(al*al) - (th0-phi)*(th0-phi)/(be*be);
	double env = exp(pwr_env);	//gaussian envelope

	double t = -1 * om * (th0 - phi ) ;
	complex<double> g ( ival * env * cos(t), ival * env * sin(t) );	//complex carrier: e^(i*t)
	
	return g;
}


double get_HD(const vector<int> &ic1, const vector<int> &ic2)
//returns fractional Hamming Distance, for the smaller vector size.
//just xor ic1 and ic2 together
//my iriscodes are different sized vectors so I just use the smaller length and stop, then normalize to it
//They would be the same size vector if I was more careful in my sampling and first normalized all images to be the 
//same sized pgm.
{
	int L1 = ic1.size();
	int L2 = ic2.size();
	int L = min(L1,L2);
	//cout<<"L1 is "<<L1<<", L2 is "<<L2<<", and L is "<<L<<endl;
	int HD=0;
	double HDf = 0.0;	//fractional HD

	for(int i=0; i<L; i++)
		if(ic1[i] != ic2[i])
			HD++;

	HDf = sqrt(HD)/sqrt(L);

	cout<<"int HD = "<<HD<<endl;
	cout<<"HDf = "<<HDf<<endl;

	return HDf;
}


vector<double> get_center(Image* im)
//returns a vector with first value x0, second value y0, where (x0,y0) is the center.
{
	int A = 0;
	double x=0;
	double y=0;
	vector<double> coords;
	int i_xy = 0;

	int ncols=im->getNCols();
	int nrows=im->getNRows();
	
	for(int i=0; i<nrows; i++)
		for(int j=0; j<ncols; j++)
		{
			i_xy = im->getPixel(i,j);
			if(i_xy < 50)
			{
				A++;
				x=x+i;
				y=y+j;
			}
		}

	x=x/double(A);
	y=y/double(A);	
		
	coords.push_back(x);
	coords.push_back(y);
	
	return coords;
}

