//author:  Ben Hixon
//date:    9-28-2010
//course:  493.69 (Computer Vision)
//description:  HW2, program 2: takes a binary image and labels its objects
//		The program accepts two arguments, an input image file and 
//		an output file.  It uses the sequential labeling algorithm
//		to label the regions. 
//use:  to compile: g++ -o p2 hw2p2.cpp Image.cpp Pgm.cpp
//	to run: ./p2 binary.pgm labeled.pgm

#include <iostream>
#include <cstdlib>
#include "Image.h"
#include <vector>
using namespace std;

void label(Image* im);	//runs the raster scan
void addToEquiv(int E[], int N, int W);	//used in raster scan to add label W to
					//label N's equivalence class.

int main(int argc, char ** argv){

	char* in_filename;
	char* out_filename;

	if(argc!=3){
		cout<<"USAGE: ./p2 binary.pgm labeled.pgm\n";
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

	label(im);

	if ( writeImage(im,out_filename)!=0 ){
		cout<<"ERROR: Can't write to output\n";
		exit(1);
	}
}

void label(Image* im)
/*Take a binary image and label its distinct objects different colors of gray.

We use the sequential labeling algorithm: scan the pixels from left to right 
starting from the top row.  For each pixel p, we know that its neighbors above, 
to the top-left, and to the left (represented by N, NW, and W respectively) are
already scanned and labeled.  

At the end, I relabel all the n labels with the first n natural numbers, and set 
the number of colors in im to the number of labels, to get a nice contrast of colors, 
as recommended in the assignment.

I have an array of equivalencies.  For each label n, if it's equivalent to something else
I'll assign that equivalency to E[n].  To set the representative equivalency values, every
time I encounter a pixel whose N and W neighbors are both labeled but have different 
labels, I know they're part of the same equivalence class so I look through E[W], whose 
value is W's label, and find every pixel with that same label and set it to E[N].*/
{
	int n=0;  //this is our label.  will be incremented every time a new region is encountered.

	int E[1024];	//this is the array of equivalence class labels.  for every label i, E[i] is the representative
			//	label of the equiv. class to which i belongs.
			//this could have been a vector instead of the ugly c-style array.

	for(int i=0; i<1024; i++){ E[i]=0; }	//initialize array of equivalence labels

	int ncols=im->getNCols();
	int nrows=im->getNRows();

	//To make things easier, I'll set every pixel on the border to black.  This won't change the 
	//labeling but will make the code easier, since I won't have to check if it's on the border.
	//For our default image, the border is black anyway so the following four lines do nothing.

	for(int i=0; i<ncols; i++){ im->setPixel(0,i,0); }  //set top row to black
	for(int i=0; i<ncols; i++){ im->setPixel(nrows-1,i,0); }  //set bottom row to black
	for(int i=1; i<nrows-1; i++){ im->setPixel(i,0,0); }  //set left column to black
	for(int i=1; i<nrows-1; i++){ im->setPixel(i,ncols-1,0); }  //set right column to black

	int p=0, N=0, NW=0, W=0;	//the pixel value and its north, northwest, and west neighbors

	//sequential label raster scan:
	for(int i=1; i<nrows-1; i++)
	{
		for(int j=1; j<ncols-1; j++)
		{
			p=im->getPixel(i,j);
			N=im->getPixel(i,j-1);
			NW=im->getPixel(i-1,j-1);
			W=im->getPixel(i-1,j);

			if(p!=0)  //p is part of an object, so label it
			{
				if(NW!=0) //if NW is labeled, set to NW's label
					im->setPixel(i,j,E[NW]);
				else if(N==0 && W!=0) //if only W is labeled, set to W's label
					im->setPixel(i,j,E[W]);
				else if(N!=0 && W==0) //if only N is labeled, set to N's label
					im->setPixel(i,j,E[N]);
				else if(N==0 && W==0)	//neither is labeled, make new label
				{
					n++;	//we've encountered a new label, so increment label numbers
					im->setPixel(i,j,n);
					E[n]=n;	//Since this is a new label, it doesn't yet belong to an equivalence class.
						//So make it its new equivalence class's representative label
				}
				else	
				//N and W are both labeled.
				//Either N==W or N!=W.  In either case, set p to N's label.
				//If N==W, no problem we're done.
				//If N!=W, they have different labels but are actually the same object, so
				//	add W to N's equivalence class.
				{
					im->setPixel(i,j,E[N]);

					if(E[N]!=E[W])	//different labels, so set every element in W's equivalence class to N's label
						addToEquiv(E,N,W);
				}
			}
		}
	}

	//second scan to set equivalence classes to the same label, and also to renumber the n labels to the first n natural numbers:

	vector<int> labels;	//to relabel the n labels to the first n natural numbers, we store the n labels in a vector so that the new
				//label is one plus the old label's position in the vector.  this makes the new labels the first n natural numbers.
				//labels.size() is the number of objects.

	for(int i=1; i<nrows-1; i++)
		for(int j=1; j<ncols-1; j++)
		{
			p=im->getPixel(i,j);
			
			if( p!= 0)
			{
				p=E[p];		//set p to its label			
				bool found=false;				
				//check to see if p's label is in our labels vector.
				for(int k=0; k<labels.size(); k++)
					if(labels[k]==p)
					{
						im->setPixel(i,j,k+1);	//set p to its label, which is its position offset by 1
						found=true;
					}
				
				if(!found)
				{
					//cout<<"New Label!  Value is: "<<p<<".  Object is #"<<labels.size()+1<<endl;
					labels.push_back(p);	//p is a new label, add to the labels vector.
					im->setPixel(i,j,labels.size());	//the new size of the vector is p's new label.
				
				}
			}
		}

	im->setColors(labels.size());	//sets number of colors to number of different objects. 			
}


void addToEquiv(int E[], int N, int W)
//for every n, E[n]=representative equivalence class label of n.
//precondition: E[N]!=0, E[W]!=0, and E[N]!=E[W]
//Action: get the representative label for N and for W, and find every n sharing
//	W's representative label and set it to N's label.
{
	int e=E[N];	//e is N's equiv. label
	int f=E[W];	//f is W's equiv. label
	
	for(int i=0; i<1024; i++)
		if(E[i]==f)
			E[i]=e; 	//if i is labeled with f, then label it with e instead, b/c they're in the same class.
}









