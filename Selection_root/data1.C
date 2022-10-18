//creating random cells distributed as a gauss
#include "Riostream.h"
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "TRandom3.h"


using namespace std;
void data1(){
    int i=0,n;
    
    n=20000;
    // creating random numbers for histogram
    TRandom3* k = new TRandom3(0);
    TRandom3* g = new TRandom3(0);
    TRandom3* z = new TRandom3(0);
    
    float x,y,ra,x2,y2;
    Float_t c[i];
    
    ra=(z->Rndm())*360;
    //cout<<ra<<endl;
    //inserting semi minor and major axis
    // loop for random nummbers and writing to file
    ofstream data;
    data.open("proef.txt");
    data<<n<<endl;
    for( i=0;i<n;i++)
    {
        y= k ->Gaus(0,10);
        x= g ->Gaus(0,30);
        x2=x*cos(ra)-y*sin(ra);
        y2=y*cos(ra)+x*sin(ra);
        // printing the data within 3 significant numbers
        
        data<<setprecision(3);
        data<<x2<<'\t'<<y2<<endl;
        //cout<<x2<<'\t'<<'\t'<<y2<<endl;
    }
    //cout<<ra<<endl;
    data.close();
}
