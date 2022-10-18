#include "TCanvas.h"
#include "TMath.h"
#include "Riostream.h"
#include <fstream>
#include <string>
#include "Tgaxis.h"
#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TEllipse.h"
#include <stdlib.h>
#include <stdio.h>

void data(){
  
    //variables
    Int_t n,i=0,bad=0,badc=0;
    fstream data;
   
    Float_t x[n],y[n],c1;
    Double_t mx,my,angle_deg,angle_rad,re[n],theta[n];
    Double_t p,pi,p_semi[n],p_major[n],p_axis,rd[n],c[n],px,py;
    
    // enter semi major and minor axis, where major>minor
    Double_t axis[3]= {120,90,-15};
    
    p_axis  = pow((axis[0]*axis[1]),2);
    c1=pow(axis[0],2)-pow(axis[1],2);
    px=1-pow((axis[1]/axis[0]),2);
    py=pow((axis[0]/axis[1]),2)-1;
    data.open("proef.txt");
    //getting data from file
    data>>n;
    i=0;
    while(!data.eof()){
            data>>x[i]>>y[i];
        cout<<x[i]<<'\t'<<y[i];
        i++;
    }
    TCanvas *g1 = new TCanvas("c1","histogram",600,500);
    g1->Divide(2,1);
    //properties histogram and making one
    TH2F* hist = new TH2F("title","original", 100,-200,200,100,-200,200);
    TH2F* hist2= new TH2F("title","rotated around O", 100,-200,200,100,-200,200);
    for( i=0;i<n;i++)
    {
        hist ->Fill(x[i],y[i]);
    }
    //getting mean values x and y
    mx      = hist->GetMean(1);
    my      = hist->GetMean(2);
    //finding angle of elliptical histogram
    
    
    //drawing ellips and histogram
    g1->cd(2);
    hist ->Draw("colz");
    TF1* f = new TF1("f","tan([1])*x",100);
    
    //setting canvas
    //fitting line
    f->SetLineWidth(1);
    f->SetLineStyle(9);
    f->SetLineColor(8);
    hist->Fit("f");
    p = f->GetParameter(1);
    pi= TMath::Pi();
    angle_deg =(p*180)/pi;
    angle_rad =(angle_deg)*(pi/180);
    cout<<"Here....D"<<endl;
    TEllipse* e = new TEllipse(mx,my,axis[0],axis[1],360,0,angle_deg);
    TEllipse* e2= new TEllipse(0,0,axis[0],axis[1],360,0,0);
    //setting properties ellips
    e->SetLineColor(1);
    e->SetLineStyle(9);
    e->SetLineWidth(2);
    e->SetFillStyle(0);
    //drawing ellips on same canvas as histogram
    e->Draw("same");
    
    // calculating distance center to data points
    for(i=0;i<n;i++){
        x[i]=x[i]-mx;
        y[i]=y[i]-my;
        x[i]=(x[i]*cos(-angle_rad))-(y[i]*sin(-angle_rad));
        y[i]=(x[i]*sin(-angle_rad))+(y[i]*cos(-angle_rad));
        hist2->Fill(x[i],y[i]);
        
        rd[i]=pow(x[i],2)+pow(y[i],2);
        c[i]=pow(x[i],2)*px+pow(y[i],2)*py;
        theta[i]=atan(y[i]/x[i]);
        p_semi[i]  = pow((axis[1]*cos(theta[i])),2);
        p_major[i] = pow((axis[0]*sin(theta[i])),2);
        
        
        
        re[i]=p_axis/(p_semi[i]+p_major[i]);
        //determining number of bad cells (outside ellips)
        if (rd[i]>=re[i] && c[i]>=c1 ) {
            bad++;
            
        }
        
        //cout<<re[i]<<'\t'<<rd[i]<<endl;
        
    }
    g1->cd(1);
    hist2->Draw("colz");
    e2->SetLineColor(1);
    e2->SetLineStyle(9);
    e2->SetLineWidth(2);
    e2->SetFillStyle(0);
    e2->Draw("same");
    g1->SaveAs("picture.eps");
   
    cout<<"number of entries outside ellips"<<endl;
    cout<<"number of cells:"<<'\t'<<bad<<endl;
    cout<<"done"<<endl;
}

