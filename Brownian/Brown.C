#include "TCanvas.h"
#include "TMath.h"
#include "TGaxis.h"
#include <stdlib.h>
#include "TLatex.h"
#include "TRandom3.h"
#include "TAxis.h"
#include "TColor.h"
#include "TGraph.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TPad.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;


void Brown(){
    
    
    Float_t k=1;
    Float_t T=300;
    Float_t tf,s,C1,g,y,tau;        y=1;        tau=1/y;             tf=22;      s=0.15625;
                                    g=2*y*k*T;  C1=(g*tau)/2;
    Float_t v_0,x_0,X;              v_0=0;      x_0=0;
    Int_t n,t0,r,x_max;             t0=0;       r=round((tf/s)+1);   n=2000;    x_max=1e3;
    
    Float_t v[r],v_mean[r],x[r],x_mean[r],v_a[r],X_i[r],t[r];
    Float_t x_low[r],x_high[r],v_comp[r],sum[r],sum1,x_damp[r],x_damp_mean[r],x_damp_ana[r];
    
    v[0]=v_0; x[0]=x_0; v_mean[0]=v_0; x_mean[0]=x_0;   v_comp[0]=1;    x_damp[0]=x_0;    x_damp_mean[0]=x_0;
    sum1=0;
    Float_t Sg,Sg_temp;          Sg=0; Sg_temp=0;
    Double_t p,tf2,a;            a=2;        tf2=tf;
    
    TRandom3 *c = new TRandom3();
    
    gROOT->SetBatch(1);
    TCanvas *c1;
    TCanvas *c2;
    
    TGraph  *gv;
    TGraph  *gx;
    TGraph  *gX;
    TGraph  *gx_damp;
    
    TH2F    *h2 = new TH2F("D(t)","p(x)",r,-x_max,x_max,r,-1,tf2+1);
    TH1F    *h1 = new TH1F("p(x)","p(x)",r,-x_max,x_max);
    h2->GetXaxis()->SetTitle("x");
    h2->GetYaxis()->SetTitle("t");

    TMultiGraph *mgv = new TMultiGraph();
    TMultiGraph *mgx = new TMultiGraph();
    TMultiGraph *mgX = new TMultiGraph();
    TMultiGraph *mgx_damp = new TMultiGraph();
    
    TLegend *l;
    
    
    for(Int_t j=0;j<n;j++){
        
        // creating the data
        
    for(Int_t i=0;i<r+1;i++){
        
        t[i]=s*i;
        
        X_i[i]=c->Gaus(0,1);
        
        v[i+1]=v[i]-(y*s*v[i])+(X_i[i]*sqrt(g*s));
        x[i+1]=x[i]+(s*v[i+1]);
        x_damp[i+1]=x_damp[i]+(X_i[i]/y)*sqrt(g*s);
        p=(1/y)*(X_i[i]*sqrt(g*s)-2*a*x_damp[i])+x_damp[i];

        cout<<v[i]<<'\t'<<x[i]<<'\t'<<X_i[i]<<endl;
        h2->Fill(p,t[i]);
        h1->Fill(p);
    }
        cout<<"---------------------------------"<<endl;
        
        gv = new TGraph(r,t,v);
        gx = new TGraph(r,t,x);
        gX = new TGraph(r,t,X_i);
        gx_damp = new TGraph(r,t,x_damp);
        mgv->Add(gv);
        mgx->Add(gx);
        mgX->Add(gX);
        mgx_damp->Add(gx_damp);
        // calculating the mean v_squared and mean x_squared
        if(j>0){
        for(Int_t i=0;i<r+1;i++){
            v_mean[i+1] = ((j-1)*v_mean[i+1]+pow(v[i+1],2))/j;
            x_mean[i+1] = ((j-1)*x_mean[i+1]+pow((x[i+1]-x[0]),2))/j;
            x_damp_mean[i+1]=((j-1)*x_damp_mean[i+1]+pow((x_damp[i+1]-x_damp[0]),2))/j;
            v_a[i] = (pow(v_0,2)-C1)*exp(-2*t[i]/tau)+C1;
            x_damp_ana[i]=(g/pow(y,2))*t[i];
            v_comp[i+1]=x_damp_mean[i+1]/x_damp_ana[i+1];
            cout<<v_mean[i]<<'\t'<<x_mean[i]<<'\t'<<x_damp_mean[i]<<endl;
            
            
        
        }
        cout<<"---------------------------------"<<endl;
        }
        
    }
    
    c1 = new TCanvas("c1","random_paths_and_velocities",300,300,400,300);
    c1->Divide(3,1);
    
    c1->cd(1);
    mgv->SetTitle("velocity"); mgv->GetXaxis()->SetTitle("t");
    mgv->Draw("AL");
    c1->cd(2);
    mgx->SetTitle("Position");  mgx->GetXaxis()->SetTitle("t");
    mgx->Draw("AL");
    c1->cd(3);
    mgX->SetTitle("X"); mgX->GetXaxis()->SetTitle("t");
    mgX->Draw("AL");
   
    
    c2 = new TCanvas("c2","random_damped_motion",300,300,400,300);
    c2->cd();
    mgx_damp->GetXaxis()->SetTitle("t");
    mgx_damp->Draw("AL");
    
    c1->SaveAs("v(t)_x(t)_X(t).pdf");
    c2->SaveAs("x(t)_damped.pdf");
    
    
    Int_t b; b=x_max;
    Double_t Vx[2*b],z[2*b];
    for (Int_t i=0;i<(2*b+1);i++){
            z[i]=-b+i;
        Vx[i]=a*pow((z[i]),2);
        
    }
   
    TGraph *gV = new TGraph(2*b,z,Vx);
    c1 = new TCanvas("c1","position Histogram",300,300,400,300);
    c1->Divide(2,1);
    c1->cd(1);
    h2->Draw("COLZ");
    
    c1->cd(2);
    h1->GetXaxis()->SetTitle("x");
    h1->GetYaxis()->SetTitle("N");
    h1->Draw("");
    c1->SaveAs("Potential.pdf");
    
    cout<<"---------------------------------"<<endl;
    cout<<"Plots saved as v(t)_x(t)_X(t).pdf"<<endl;
    cout<<"Damped motion saved as x_damp.pdf"<<endl;
    cout<<"Histograms saved as Potential.pdf"<<endl;
    cout<<"---------------------------------"<<endl;
    
    c1 = new TCanvas("c1","mean_squared_velocity_position",300,300,400,300);
    c1->Divide(2,1);
    c1->cd(1);
    gv = new TGraph(r,t,v_mean);
    gv->SetTitle("< v^{2}>"); //gv->SetMarkerStyle(kCircle); gv->SetMarkerSize(0.5);
    gv->SetLineColor(kRed);   gv->GetXaxis()->SetTitle("t");     gv->Draw("APL");
    gX = new TGraph(r,t,v_a);
    gX->SetTitle("<v_{analytic}^{2}>");    //gX->SetMarkerStyle(kStar); gX->SetMarkerSize(0.5);
    gX->SetLineColor(kGreen);    gX->GetXaxis()->SetTitle("t");   gX->Draw("SAME");
    
    l = new TLegend(0.1,0.1,0.3,0.2);
    l->AddEntry(gv,"Data");
    l->AddEntry(gX,"Analytic");
    l->Draw();
    
    c1->cd(2);
    gx = new TGraph(r,t,x_mean);
    gx->SetTitle("< x^{2}>");    //gx->SetMarkerStyle(kPlus); gx->SetMarkerSize(0.5);
    gx->SetLineColor(kBlue);    gx->GetXaxis()->SetTitle("t");  gx->Draw("APL");
    
    
    c1->SaveAs("v_mean_x_mean_v_analytic.pdf");
    
    c2 = new TCanvas("c2","mean_damped_motion",300,300,400,300);
    c2 -> Divide(2,1);
    c2->cd(1);
    gx_damp = new TGraph(r,t,x_damp_mean);
    gx_damp->SetTitle("<x_{damped}^{2}>"); gx_damp->SetMarkerStyle(kPlus); gx_damp->SetMarkerSize(0.5); gx_damp->SetLineColor(kRed);
    gx_damp->GetXaxis()->SetTitle("t"); gPad->SetLogy();
    gPad->SetLogx();    gx_damp->Draw("AP");
    
    gX = new TGraph(r,t,x_damp_ana);
    gX->SetTitle("<x_{damped}^{2}>"); gX->SetLineColor(kRed);
    gX->GetXaxis()->SetTitle("t"); gPad->SetLogy();
    gPad->SetLogx();    gX->Draw("SAME");
    
    l = new TLegend(0.1,0.75,0.5,0.9);
    l->AddEntry(gx_damp,"Data");
    l->AddEntry(gX,"Analytic fit");
    l->Draw();
    
    c2->cd(2);
    gv = new TGraph(r,t,v_comp);
    gv->SetTitle("<x_{data}^{2}>/<x_{an.}^{2}>"); gv->SetMarkerStyle(kStar); gv->SetMarkerSize(0.5);   gv->GetXaxis()->SetTitle("t");     gv->Draw("AP");
    
    
    c2->SaveAs("Mean_dampted_position.pdf");
    
    cout<<"---------------------------------"<<endl;
    cout<<"Mean Plots saved as v_mean_x_mean_v_analytic.pdf"<<endl;
    cout<<"Mean Damped motion saved as Mean_dampted_position.pdf"<<endl;
    cout<<"---------------------------------"<<endl;
    
    
    
    for(Int_t i=0;i<r+1;i++){
        
        v_comp[i+1] = v_mean[i+1]/v_a[i+1];
        sum1=sum1+v_comp[i];
        
        x_low[i] = k*T*pow(t[i],2);
        x_high[i] = 2*k*T*t[i];
      
        
        //cout<<v_comp[i]<<'\t'<<x_low[i]<<'\t'<<x_high[i]<<endl;
    }
    for(Int_t i=0;i<r;i++){
    sum[i]=(sum1/(r));
        Sg_temp=pow((v_comp[i]-sum[0]),2)+Sg_temp;
        Sg=sqrt(Sg_temp/(r*(r-1)));
     }
    
    
    c1 = new TCanvas("c1","comparison_plot",300,300,400,300);
    c1->Divide(2,1);
    c1->cd(1);
    gv = new TGraph(r,t,v_comp);
    gv->SetTitle("< v^{2}>/<v_{analytic}^{2}>"); gv->SetMarkerStyle(kStar);    gv->SetMarkerSize(0.5);    gv->GetXaxis()->SetTitle("t");   gv->Draw("AP");
    gX = new TGraph(r,t,sum);   gX->SetLineColor(kBlue);
    gX->GetXaxis()->SetTitle("t");  gX->Draw("SAME");
    
    l = new TLegend(0.1,0.1,0.3,0.2);
    l->AddEntry(gv,"v^{2}/v_{analytic}^{2}");
    l->AddEntry(gX,"Mean");
    l->Draw();
    
    c1->cd(2);
    
    //-------------------------------------------------------
    //high range t
    
    gv = new TGraph(r,t,x_mean);
    gPad->SetLogy();
    gPad->SetLogx();
    gv->SetLineColor(kBlue);    gv->SetTitle("<x^{2}>");    gv->GetXaxis()->SetTitle("t");
    gv->Draw("APL");
    
    gx = new TGraph(r,t,x_low);
    gPad->SetLogy();
    gPad->SetLogx();
    gx->SetLineColor(kRed); gx->GetXaxis()->SetTitle("t");
    gx->Draw("SAME");
    
    gX = new TGraph(r,t,x_high);
    gPad->SetLogy();
    gPad->SetLogx();
    gX->SetLineColor(kGreen);   gX->GetXaxis()->SetTitle("t");
    gX->Draw("SAME");
    
    l = new TLegend(0.1,0.75,0.5,0.9);
    l->AddEntry(gv,"<x^{2}>");
    l->AddEntry(gx,"t#rightarrow 0");
    l->AddEntry(gX,"t#rightarrow#infty");
    l->Draw();
    
    
    
    c1->SaveAs("Comparison_plot.pdf");
    
    cout<<"---------------------------------"<<endl;
    cout<<"mean_v_comp = "<<'\t'<<sum[0]<<" +/- "<<Sg<<endl;
    cout<<"Plots saved as comparison.pdf"<<endl;
    cout<<"---------------------------------"<<endl;

    
    
}



