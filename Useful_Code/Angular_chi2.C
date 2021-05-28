#include <stdio.h>
#include <map>
#include <fstream>

#include <TGraph.h>

//a2E is the experimentally fitted value of the data, error, A0 is the L=0 solution to eq3.9, Qk is a geometric attenuation factor, Bk is a Clsheb Gordon coeff., Rk1-3 are Racah Coeffients, dmin/max are the range of  
struct transition {

  double a2E;
  double a2Error;

  double a4E;
  double a4Error;

  double A0;

  double Bk;
  double Rk1;
  double Rk2;
  double Rk3;

  double Bk2; 
  double Rk12; 
  double Rk22; 
  double Rk32; 

  transition( double a2E_f, double a2Error_f, double a4E_f, double a4Error_f, double A0_f, double Bk_f, double Rk1_f,	double Rk2_f,	double Rk3_f,double Bk2_f, double Rk12_f, double Rk22_f, double Rk32_f ){

    a2E =     a2E_f;
    a2Error = a2Error_f;   
    a4E =     a4E_f; 
    a4Error = a4Error_f; 

    A0 =      A0_f;  

    Bk =      Bk_f; 
    Rk1 =     Rk1_f; 
    Rk2 =     Rk2_f; 
    Rk3 =     Rk3_f; 

    Bk2  =    Bk2_f; 
    Rk12 =    Rk12_f; 
    Rk22 =    Rk22_f; 
    Rk32 =    Rk32_f; 
  }

  /*
     transition( double a2E_f, double a2Error_f, double A0_f, double Bk_f, double Rk1_f, double Rk2_f, double Rk3_f ){

     a2E =     a2E_f;
     a2Error = a2Error_f;
     A0 =      A0_f;
     Bk =      Bk_f;
     Rk1 =     Rk1_f;
     Rk2 =     Rk2_f;
     Rk3 =     Rk3_f;
     }

*/



}; 

TGraph* ADBestFit(double step, transition gamma){

  double Qk = 0.03;


  double a2E = gamma.a2E;
  double a2Error = gamma.a2Error;
  double a4E = gamma.a4E;
  double a4Error = gamma.a4Error;

  double A0 = gamma.A0;

  double Bk = gamma.Bk;
  double Rk1 = gamma.Rk1;
  double Rk2 = gamma.Rk2;
  double Rk3 = gamma.Rk3;

  double Bk2 = gamma.Bk;
  double Rk12 = gamma.Rk1;
  double Rk22 = gamma.Rk2;
  double Rk32 = gamma.Rk3;

  double dmin = -3.14159 / 2.0; 
  double dmax = 3.14159 / 2.0; 

  int points = (dmax - dmin)/step;    
  printf("Points = %i\n", points);


  TGraph* AD= new TGraphErrors(points);

  int num = 0;
  double delta = 0.0; 
  double td = 0.0;

  double rd_1 = 0.0;
  double rd_2= 0.0;

  double A2 = 0.0;
  double A4 = 0.0;

  double a2T = 0.0; 
  double a4T = 0.0; 

  double X2 = 0.0; 

  double ang = 90.0; //90 or 135 ?

  double p2 = ::ROOT::Math::legendre(2,TMath::Cos(TMath::DegToRad()*ang));
  double p4 = ::ROOT::Math::legendre(4,TMath::Cos(TMath::DegToRad()*ang));

  double W_T = 0.0; 
  double W_E = 0.0; 

  double total_error = 0.0;  

  for(int i = 0; i < points; i++){	

    td    = i * step + dmin; 
    delta = TMath::Tan(td);
    rd_1    = (Rk1 + 2*delta*Rk2 + pow(delta,2) * Rk3)/(1+pow(delta,2)); 
    rd_2    = (Rk12 + 2*delta*Rk22 + pow(delta,2) * Rk32)/(1+pow(delta,2)); 
    A2    = Qk*Bk*rd_1; 
    A4    = Qk*Bk2*rd_2; 

 //   W_T   = A0 + A2*p2 + A4*p4; 
//    W_E   = A0*(1+a2E*p2 + a4E*p4); 
      W_T   = A2*p2 + A4*p4;
      W_E   = (a2E/A0)*p2 + (a4E/A0)*p4;
    total_error = W_E * TMath::Sqrt(pow(a2Error,2)+pow(a4Error,2)+pow(0.0708,2)); 

    X2    = TMath::Abs((W_T - W_E))/(2*pow(total_error,2));

    // printf("delta = %f, td = %f, rd = %f, A2 = %f\n",delta,td,rd,A2);

    AD->SetPoint(num,td*180/3.1415,X2);

    //dmin += i * step;	
    num++; 
  }

  return AD; 
}

int Angular_chi(double step = 0.00001){
  //228 values
  double a2E = 0.134;
  double a2Error = 0.084;
  double a4E = 0.313; 
  double a4Error = 0.093;

  double A0 = 0.503;

  double Bk = -11.289;
  double Rk1 = 3.294;
  double Rk2 = -4.253;
  double Rk3 = -5.491;

  double Bk2 = 8.745; 
  double Rk12 = -2.552; 
  double Rk22 = 3.295; 
  double Rk32 = 4.253; 

//  transition gamma = transition( a2E,   a2Error,a4E,a4Error, A0,Bk,    Rk1,   Rk2,   Rk3,   Bk2,   Rk12, Rk22, Rk32);   //228

  //transition gamma = transition(-0.140,0.085,0.124,0.083,0.937,12.736,3.468,-4.477, -5.778,-9.865,-2.686,3.468,4.478);  //494 
transition gamma = transition(0.284,0.092,0.989,0.127, 1.663,15.222,3.753,-4.845,-6.255,-11.791,-8.722,11.260,14.536);//565
//  transition gamma = transition(-0.794,0.177,-1.198,0.137, 1.355,-13.615,-5.95,7.042,8.331,-10.546,-4.609,5.454,6.454); //1001
//  transition gamma = transition(-0.897,0.122,-1.455,0.128,1.980,11.791,5.592,-6.617,-7.829,-9.133,-4.332,5.125,6.064);  //1147



  TGraph * theGraph = ADBestFit(step,gamma);

  theGraph->SetTitle("565 [keV] 19/2\\hbar \\rightarrow 17/2\\hbar ");
  theGraph->GetXaxis()->SetTitle("atan(\\delta) [Deg]");
  theGraph->GetYaxis()->SetTitle("\\Chi ^2");


  theGraph->SetLineColor(2);
  theGraph->SetLineWidth(4);

  auto canvas_1 = new TCanvas("canvas_1", "canvas_1");

  theGraph->Draw("AC");

  canvas_1->SetLogy();
  //	theGraph->GetHistogram()->GetYaxis()->SetRangeUser(0,1000);
  canvas_1->SaveAs("AD565_2.pdf");

  return 1;

}



