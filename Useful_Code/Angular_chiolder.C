#include <stdio.h>
#include <map>
#include <fstream>


#include </home/derosa/packages/GRUTinizer/include/GPeak.h>
#include </home/derosa/packages/GRUTinizer/include/GRootCommands.h>
#include </home/derosa/packages/GRUTinizer/include/TChannel.h>
#include <TGraph.h>

//a2E is the experimentally fitted value of the data, error, A0 is the L=0 solution to eq3.9, Qk is a geometric attenuation factor, Bk is a Clsheb Gordon coeff., Rk1-3 are Racah Coeffients, dmin/max are the range of  


TGraph* ADBestFit(double a2E, double a2Error, double A0, double Qk, 
                  double Bk, double Rk1, double Rk2, double Rk3, int dmin, int dmax, double step){

a2E = 0.134;
a2Error =  .084;
A0 = 0.503;
Qk = 0.03;
Bk = -11.289;
Rk1 = 3.294;
Rk2 = -4.253;
Rk3 = -5.491;
dmin = -20.0; 
dmax = 20.0; 
step = 0.00001;
  int points = (dmax - dmin)/step;    

  TGraph* AD= new TGraphErrors(points);

  int num = 0;
  double delta = 0.0; 
  double td = 0.0;
  double rd = 0.0;
  double A2 = 0.0;
  double a2T = 0.0; 
  double X2 = 0.0; 

//The two arrays that are populated with local variables X2 and td*180/3.1415 at each index num      

  std::vector<double> Chisqr;
  std::vector<double> atanD; 
  for(int i = 0; i < points; i++){

    delta = dmin;

    td    = TMath::ATan(delta); 
    rd    = (Rk1 + 2*delta*Rk2 + pow(delta,2) * Rk3)/(1+pow(delta,2)); 
    A2    = Qk*Bk*rd; 

    a2T   = A2/A0; 

    X2    = TMath::Abs((a2T - a2E))/(2*pow(a2Error,2));

   // printf("delta = %f, td = %f, rd = %f, A2 = %f\n",delta,td,rd,A2);


    AD->SetPoint(num,td*180/3.1415,X2);
    
    dmin += i*step;
    num++; 
    }
    
    return AD; 

  }





