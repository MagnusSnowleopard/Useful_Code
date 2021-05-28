#include <stdio.h>
#include <map>
#include <fstream>

#include <TGraph.h>

//a2E is the experimentally fitted value of the data, error, A0 is the L=0 solution to eq3.9, Qk is a geometric attenuation factor, Bk is a Clsheb Gordon coeff., Rk1-3 are Racah Coeffients, dmin/max are the range of  
struct transition {
   
	double a2E = 0.;
	double a2Error = 0.;
	double A0 = 0.;
	double Bk = 0.;
	double Rk1 = 0.;
	double Rk2 = 0.;
	double Rk3 = 0.;

  transition( double a2E_f, double a2Error_f, double A0_f, double Bk_f, double Rk1_f,	double Rk2_f,	double Rk3_f ){
    
    a2E =     a2E_f;
    a2Error = a2Error_f;   
    A0 =      A0_f;  
    Bk =      Bk_f; 
    Rk1 =     Rk1_f; 
    Rk2 =     Rk2_f; 
    Rk3 =     Rk3_f; 
    }

}; 

TGraph* ADBestFit(double step, transition gamma){
    
	double Qk = 0.03;

  
	double a2E = gamma.a2E;
	double a2Error = gamma.a2Error;
	double A0 = gamma.A0;
	double Bk = gamma.Bk;
	double Rk1 = gamma.Rk1;
	double Rk2 = gamma.Rk2;
	double Rk3 = gamma.Rk3;

	double dmin = -3.14159 / 2.0; 
	double dmax = 3.14159 / 2.0; 

	int points = (dmax - dmin)/step;    
	printf("Points = %i\n", points);
	

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
		//delta = dmin;

		td    = i * step + dmin; 
//		td    = TMath::ATan(delta); 
		delta = TMath::Tan(td);
		rd    = (Rk1 + 2*delta*Rk2 + pow(delta,2) * Rk3)/(1+pow(delta,2)); 
		A2    = Qk*Bk*rd; 
		
		a2T   = A2/A0; 		
		X2    = TMath::Abs((a2T - a2E))/(2*pow(a2Error,2));
		
		// printf("delta = %f, td = %f, rd = %f, A2 = %f\n",delta,td,rd,A2);
		
		AD->SetPoint(num,td*180/3.1415,X2);
	
		//dmin += i * step;	
		num++; 
	  }
	  
	  return AD; 
}

int Angular_chi(double step = 0.00001){

	double a2E = 0.134;
	double a2Error = 0.084;
	double A0 = 0.503;
	double Bk = -11.289;
	double Rk1 = 3.294;
	double Rk2 = -4.253;
	double Rk3 = -5.491;
  
//  transition gamma = transition(a2E,a2Error,A0,Bk,Rk1,Rk2,Rk3); //228
//  transition gamma = transition(-0.140,0.085,0.937,12.735,3.468,-4.477,-5.780); //494
//  transition gamma = transition(0.284,0.092,1.663,15.222,3.753,-4.845,-6.255);  //565
    transition gamma = transition(-0.794,0.117,1.355,-13.615,-5.951,7.042,8.332); //1001
//     transition gamma = transition(-0.897,0.122,1.980,11.791,5.592,-6.617,-8.367); //1147
	TGraph * theGraph = ADBestFit(step,gamma);

	theGraph->SetTitle("1001 [keV] 15/2\\hbar \\rightarrow 11/2\\hbar ");
	theGraph->GetXaxis()->SetTitle("atan(\\delta)");
	theGraph->GetYaxis()->SetTitle("\\Chi ^2");


	theGraph->SetLineColor(2);
	theGraph->SetLineWidth(4);

	auto canvas_1 = new TCanvas("canvas_1", "canvas_1");

	theGraph->Draw("AC");

	canvas_1->SetLogy();
//	theGraph->GetHistogram()->GetYaxis()->SetRangeUser(0,1000);
	canvas_1->SaveAs("AD1001.pdf");

	return 1;

}



