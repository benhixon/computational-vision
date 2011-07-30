//author:  Ben Hixon
//date:    9-28-2010
//course:  493.69 (Computer Vision)
//description:  HW2, program 3: takes a labeled image and for each object in
//		the image, finds its center of mass, minimum moment of inertia,
//		and orientation.  Store these values in a file, and also output
//		an image displaying for each object in the input image a dot for 
//		its center of mass and a line for its orientation.
//use:  to compile: g++ -o p3 hw2p3.cpp Image.cpp Pgm.cpp Line.cpp
//	to run: ./p3 labeled.pgm database.txt objects.pgm

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "Image.h"
#include <vector>
#include <cmath>
using namespace std;

void makeDatabase(Image* im, char* db_filename);  //makes vector of strings each containing the image info
void writeToFile(vector<string> strings, char* filename);  //writes the strings to a file

int main(int argc, char ** argv){

	char* in_filename;
	char* db_filename;
	char* out_filename;

	if(argc!=4){
		cout<<"USAGE: ./p3 labeled.pgm database.txt objects.pgm\n";
		exit(1);
	}else{
				
		in_filename = *(argv+1);
		db_filename = *(argv+2);
		out_filename = *(argv+3);
	}

	Image *im = new Image;
	if ( readImage(im,in_filename)!=0 ){
		cout<<"ERROR: Can't read from input\n";
		exit(1);
	}

	makeDatabase(im, db_filename);

	if ( writeImage(im,out_filename)!=0 ){
		cout<<"ERROR: Can't write to output\n";
		exit(1);
	}
}

void makeDatabase(Image* im, char* db_filename)
/* First we count the objects by counting the different pixel values in the image.
Then for each object:
	(1) Get area A by counting total number of pixels with that label (all pixels 
		with the same label are connected, and  the values for each pixel are the
		same so let b(x,y)=1.)
	(2) Get x value of centroid by adding up all x-vals of pixels in the image and dividing by A.
	(3) Get y value of centroid by addding up all y-vals of pixels in the image and dividing by A.
	(4) Get second moment of inertia by getting values a, b, and c.

Then draw centroids and orientation lines on each object.  We use atan2(a, b-c) to get theta, and 
then use x*sin(theta)-y*cos(theta)+rho=0 to get a second coordinate for the Line() function.
*/
{
	vector<string> database;  //each string has format "label x_bar y_bar a b c theta"
	vector<int> labels;	  //holds label values of the different objects.

	int ncols=im->getNCols();
	int nrows=im->getNRows();

	//get number and value of labels for each object (I don't assume they're labeled with the counting numbers)
	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
		{
			int p=im->getPixel(i,j);
			
			if( p!= 0)
			{			
				bool found=false;				
				//check to see if p's label is already in our labels vector.
				for(int k=0; k<labels.size(); k++)
					if(labels[k]==p){ found=true; }
				
				if(!found)	//it wasn't in our labels vector, so add it
					labels.push_back(p);
			}
		}
	//cout<<"There are "<<labels.size()<<" objects in the picture."<<endl;

	for(int k=0; k<labels.size(); k++)	
	/*for each label, raster scan the picture and accumulate totals.
		area += 1
		first moment x += j
		first moment y += i
		second moments: +=ii, +=ij, and +=jj 
	*/	
	{
		long int area = 0;
		double x_bar = 0.0, y_bar = 0.0;	//centroid
		double x2=0.0, y2=0.0;		//for a second point on the line of least inertia
		double a=0.0,b=0.0,c=0.0;	//second moments
		double theta=0.0;	//we know that tan(2*theta) == b/(a-c)
		double rho=0.0;		//we know that x*sin(theta)-y*cos(theta)+rho = 0	
		vector<int> xVals;	//to hold the x and y coordinate values of each point in each object.
		vector<int> yVals;
		
		//Raster scan: update area, xbar, ybar, and store values for getting a,b, and c later
		//There's probably a way to get a,b, and c in this loop as well but I couldn't quite see how to do it,
		//	since we need xbar and ybar in order to calculate a,b, and c.
		for(int i=1; i<nrows-1; i++)
			for(int j=1; j<ncols-1; j++)
			{	
				int p=im->getPixel(i,j);
				if( p==labels[k] )
				{
					area++;
					xVals.push_back(i);
					yVals.push_back(j);

					x_bar = x_bar + i;
					y_bar = y_bar + j;
				}
			}

		x_bar = x_bar / float(area);	//divide the integral by the area to get the coordinate
		y_bar = y_bar / float(area);

		//get a,b, and c using previously found x_bar, y_bar, and the stored points
		//We subtract x_bar from xVals[i] to correct for the change of variables as
		//described in lecture.
		for(int i=0; i<xVals.size(); i++)
		{
			a = a + ( xVals[i] - x_bar )*( xVals[i] - x_bar );
			b = b + 2*( xVals[i] - x_bar )*( yVals[i] - y_bar );
			c = c + ( yVals[i] - y_bar)*(yVals[i] - y_bar );
		}

		//from a, b, and c, get theta, using tan(2*theta) == b/(a-c)
		theta = (0.5) * atan2( b, a-c );	//atan2 will give us our theta taking the quadrant into consideration

		//our equation of the line is x*sin(theta)-y*cos(theta)+rho = 0.  So use x_bar, y_bar, and thta to get rho.
		//	Then to get another point near the centroid, just choose x=x_bar+dx and solve for y.		
		rho = y_bar*cos(theta) - x_bar*sin(theta);
		
		x2 = x_bar+50;
		y2 = (x2*sin(theta) + rho)/cos(theta);

		//draw dot as a mxm black pixel square around x_bar,y_bar
		for(int i=int(x_bar)-4; i<int(x_bar)+4; i++)
			for(int j=int(y_bar)-4; j<int(y_bar)+4; j++)
				im->setPixel(i,j,0);

		//draw line from centroid to x2,y2
		line(im, int(x_bar), int(y_bar), int(x2), int(y2), 255);
		
		//add info to string using the stringstream class
		string data;
		ostringstream str_stream;
		str_stream << labels[k]<<" " << x_bar<<" " << y_bar<< " " << a<< " " << b<< " " << c<< " " << theta <<'\n';
		data=str_stream.str();
		database.push_back(data);
	}

	//write the database vector to the file db_filename
	writeToFile(database, db_filename);
}

void writeToFile(vector<string> strings, char* filename)
{
	ofstream ofile(filename);
	if( ofile.is_open() )
	{
		for(int i=0; i<strings.size(); i++)
			ofile << strings[i];
		ofile.close();
	}
	else cout<<"Error, can't open file"<<filename<<endl;
}




