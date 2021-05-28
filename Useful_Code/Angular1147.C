#include <stdio.h>
#include <map>

#include </home/derosa/packages/GRUTinizer/include/GPeak.h>
#include </home/derosa/packages/GRUTinizer/include/GRootCommands.h>
#include </home/derosa/packages/GRUTinizer/include/TChannel.h>
#include <TGraph.h>

double Only_Back(int i, int j, double * original){
  double a, b;
  a = original[j];
  b = ( original[j + i] + original[j - i] ) / 2.0;
  if (b < a) a = b;
  return a;
}

//Attempt to have inhouse TSpectrum background subtraction. 
TH1D* BackSub(TH1D* Hist, int Iterations = 30, char *opt_2 = (char*) "N", char * opt = (char*) "Forwards"){
  //This is for the option
  TString sopt(opt);
  sopt.ToLower();
  TString sopt_2(opt_2);
  sopt_2.ToLower();

  //Turn hist into array
  Int_t first = Hist->GetXaxis()->GetFirst();
  Int_t last  = Hist->GetXaxis()->GetNbins();


  Int_t hist_size = last-first;
  double * original = new double[hist_size];
  double * reference = new double[hist_size];

  for(int i =0; i < hist_size; i++){
    original[i] = Hist->GetBinContent(i + first);
    reference[i] = Hist->GetBinContent(i + first);
  }


  double * new_hist = new double[hist_size];
  if(sopt.Contains("forwards")) {
    for(int i=1; i <= Iterations; i++){
      double start = 10 * i;
      if(i == 1) start =i;
      for(int j = start; j < hist_size - i; j++){
        new_hist[j] = Only_Back(i, j, original);
      }
      for(int j = i; j < hist_size - i; j++){
        original[j] = new_hist[j];
      }
    }
  }

  if(sopt.Contains("back")) {
    for(int i=Iterations; i >= 1; i--){
      double start = 10 * i;
      if(i == 1) start = i;
      for(int j = start; j < hist_size - i; j++){
        new_hist[j] = Only_Back(i, j, original);
      }
      for(int j = i; j < hist_size - i; j++){
        original[j] = new_hist[j];
      }
    }
  }

  TH1D * output_hist = new TH1D(Form("%s_Background",Hist->GetName()),Form("%s_Background",Hist->GetName()),hist_size,0,hist_size);
  TH1D * output_hist_sub = new TH1D(Form("%s_Background_Sub",Hist->GetName()),Form("%s_Background_Sub",Hist->GetName()),hist_size,0,hist_size);
  for(int i=0; i < hist_size; i++){
    output_hist->SetBinContent(i, original[i]);
    output_hist_sub->SetBinContent(i, (reference[i] - original[i]));  
  }

  delete[]reference;
  delete[]original;

  if(sopt_2.Contains("sub")) return output_hist_sub;

  return output_hist;

}



/*

#define G1 1
#define G2 2
#define G3 3
#define G4 4
#define G5 9
#define G6 10
#define G7 11
#define G8 12
#define G9 17
#define G10 18
#define G11 19
#define G12 20
#define G13 23
#define G14 25
#define G15 26
#define G16 27
#define G17 28
#define G18 33
#define G19 34
#define G20 35
#define G21 36
#define G22 39
#define G23 40
#define G24 41 
#define G25 42
#define G26 43
#define G27 44
std::map<int,double> gmap;
void initmap(){
gmap[G1]= 86.0;
gmap[G2]= 86.0;
gmap[G3]= 94.0;
gmap[G4]= 94.0;
gmap[G5]= 134.0;
gmap[G6]= 134.0;
gmap[G7]= 124.0;
gmap[G8]= 124.0;
gmap[G9]= 86.0;
gmap[G10]= 86.0;
gmap[G11]= 94.0;
gmap[G12]= 94.0;
gmap[G13]= 90.0;//not clover 
gmap[G14]= 134.0;
gmap[G15]= 134.0;
gmap[G16]= 124.0;
gmap[G17]= 124.0;
gmap[G18]= 94.0;
gmap[G19]= 86.0;
gmap[G20]= 94.0;
gmap[G21]= 86.0;
gmap[G22]= 45.0;//not clover
gmap[G23]= 45.0;//not clover
gmap[G24]= 94.0;
gmap[G25]= 94.0;
gmap[G26]= 86.0;
gmap[G27]= 86.0;
}
*/
// This code displays 12 over layed graphs of each angle and 2 - 8 fits on top of eachother
//
// This code also displays the distribution of points summed by area after fitting as a function of angle. 
//
// This code also subtracts the non-linear background before fitting the peaks to reduce error. 
TGraph* Angular_Dist4(TH2D* interest, TH2D* norm, double ih_win, double il_win, double i_high, double i_low,
    double nh_win, double nl_win, double  n_high, double n_low, int iterations_of_bg){

  //  iinitmap();
  int detector_number = 60; 

  std::map<int,double>  interest_map;
  std::map<int,double>  interest_map_err;

  std::map<int,double>  norm_map;
  std::map<int,double>  norm_map_err;

  double interest_fit_sum[7]; 
  double interest_fit_sum_err[7];

  double norm_fit_sum[7]; 
  double norm_fit_sum_err[7];

  double angle[7]={45.,86.,90.,94.,124.,134.};

  TCanvas* I = new TCanvas("I","I",1500,800); 
  I->Divide(3,2);
  for(int j =0; j<detector_number;j++){


    TH1D* sumxi = interest->ProjectionX(Form("interest_x%02i",j),j+1,j+1); 

  

    int b = 0; 
    //45 degrees
    if(j == 39 || j == 40) {

      I->cd(1);  
    TSpectrum s_i; 
    TH1* bg_i = s_i.Background(sumxi,40,"BackOrder6"); 
    TH1D* singles_i_bg = (TH1D*)sumxi->Clone("singles_i_bg"); 
    singles_i_bg->Add(bg_i,-1); 

    TH1D* interest_x = singles_i_bg; 

    interest_x->GetXaxis()->SetRangeUser(il_win,ih_win); 
      if(j==39){
        interest_x->SetTitle("Interest at 45 Degrees");
        interest_x->Draw();
      }
      else {
        interest_x->SetLineColor(b+1);
        interest_x->Draw("same");
      }
      if(interest_x->Integral()>100){
        GPeak* peak_i = PhotoPeakFit(interest_x,i_low,i_high);
        interest_map[j] = peak_i->GetSum(); 
        interest_map_err[j] = peak_i->GetSumErr(); 
        b++;
      }   
      if(b>0){
        interest_fit_sum[0] +=interest_map[j];  
        interest_fit_sum_err[0] +=pow(interest_map_err[j],2);
      }
       
    }
    
    int b1 = 0;
    //86 degrees
    if(j == 1 || j == 2 || j == 17 || j == 18 || j == 34 || j  == 36 || j  == 43 || j == 44) {

      I->cd(2);      
    TSpectrum s_i; 
    TH1* bg_i = s_i.Background(sumxi,25,"BackOrder4"); 
    TH1D* singles_i_bg = (TH1D*)sumxi->Clone("singles_i_bg"); 
    singles_i_bg->Add(bg_i,-1); 

    TH1D* interest_x = singles_i_bg; 

    interest_x->GetXaxis()->SetRangeUser(il_win,ih_win); 
      if(j==1){
        interest_x->SetTitle("Interest at 86 Degrees");
        interest_x->Draw();
      }
      else {
        interest_x->SetLineColor(b1);
        interest_x->Draw("same");
      }
      if(interest_x->Integral()>100){
        GPeak* peak_i = PhotoPeakFit(interest_x,i_low,i_high);
        interest_map[j] = peak_i->GetSum(); 
        interest_map_err[j] = peak_i->GetSumErr(); 
        b1++; 
      }   
      if(b1>0){
        interest_fit_sum[1] +=interest_map[j];  
        interest_fit_sum_err[1] +=pow(interest_map_err[j],2);
      }
    }
    //90 degrees
    if(j == 23) {

      I->cd(3);      
    TSpectrum s_i; 
    TH1* bg_i = s_i.Background(sumxi,8,"BackOrder2"); 
    TH1D* singles_i_bg = (TH1D*)sumxi->Clone("singles_i_bg"); 
    singles_i_bg->Add(bg_i,-1); 

    TH1D* interest_x = singles_i_bg; 

    interest_x->GetXaxis()->SetRangeUser(il_win,ih_win); 
      if(j==23){
        interest_x->SetTitle("Interest at 90 Degrees");
        interest_x->Draw();
      }
      if(interest_x->Integral()>100){
      //  GPeak* peak_i = PhotoPeakFit(interest_x,i_low,i_high);
        GPeak* peak_i = PhotoPeakFit(interest_x,i_low,i_high);
        interest_fit_sum[2] = peak_i->GetSum();;  
        interest_fit_sum_err[2] = pow(peak_i->GetSumErr(),2); 
      }   

    }
    int b2 = 0; 
    //94 degrees
    if(j == 3 || j == 4 || j == 19 || j == 20 || j == 33 || j  == 35 || j  == 41 || j == 42) {

      I->cd(4);      
    TSpectrum s_i; 
    TH1* bg_i = s_i.Background(sumxi,26,"BackOrder4"); 
    TH1D* singles_i_bg = (TH1D*)sumxi->Clone("singles_i_bg"); 
    singles_i_bg->Add(bg_i,-1); 

    TH1D* interest_x = singles_i_bg; 

    interest_x->GetXaxis()->SetRangeUser(il_win,ih_win); 
      if(j==3){
        interest_x->SetTitle("Interest at 94 Degrees");
        interest_x->Draw();
      }
      else { 
        interest_x->SetLineColor(b2+1);
        interest_x->Draw("same");
      }
      if(interest_x->Integral()>100){
        GPeak* peak_i = PhotoPeakFit(interest_x,i_low,i_high);
        interest_map[j] = peak_i->GetSum(); 
        interest_map_err[j] = peak_i->GetSumErr(); 
        b2++; 
      }   
      if(b2>0){
        interest_fit_sum[3] += interest_map[j];  
        interest_fit_sum_err[3] +=pow(interest_map_err[j],2);
      }
    }
    int b3 = 0;
    //124 degrees
    if(j == 9 || j == 10 || j == 25 || j == 26) {

      I->cd(5);      
    TSpectrum s_i; 
    TH1* bg_i = s_i.Background(sumxi,40,"BackOrder6"); 
    TH1D* singles_i_bg = (TH1D*)sumxi->Clone("singles_i_bg"); 
    singles_i_bg->Add(bg_i,-1); 

    TH1D* interest_x = singles_i_bg; 

    interest_x->GetXaxis()->SetRangeUser(il_win,ih_win); 
      if(j==9){
        interest_x->SetTitle("Interest at 124 Degrees");
        interest_x->Draw();
      }
      else {   
        interest_x->SetLineColor(b3+1);
        interest_x->Draw("same");
      }
      if(interest_x->Integral()>100){
     //   GPeak* peak_i = PhotoPeakFit(interest_x,i_low,i_high);
        GPeak* peak_i = PhotoPeakFit(interest_x,1143,i_high);
        interest_map[j] = peak_i->GetSum(); 
        interest_map_err[j] = peak_i->GetSumErr(); 
        b3++;
      }   
      if(b3>0){
        interest_fit_sum[4] += interest_map[j];
        interest_fit_sum_err[4] +=pow(interest_map_err[j],2);
      }
    }
    int b4 = 0; 
    //134 degrees
    if(j == 11 || j == 12 || j == 27|| j == 28) {

      I->cd(6);      
    TSpectrum s_i; 
    TH1* bg_i = s_i.Background(sumxi,iterations_of_bg,"BackOrder4"); 
    TH1D* singles_i_bg = (TH1D*)sumxi->Clone("singles_i_bg"); 
    singles_i_bg->Add(bg_i,-1); 

    TH1D* interest_x = singles_i_bg; 

    interest_x->GetXaxis()->SetRangeUser(il_win,ih_win); 
      if(j==11){
        interest_x->SetTitle("Interest at 134 Degrees");
        interest_x->Draw();
      }
      else {  
        interest_x->SetLineColor(b4+1);
        interest_x->Draw("same");
      }
      if(interest_x->Integral()>100){
//        GPeak* peak_i = PhotoPeakFit(interest_x,i_low,i_high);
        GPeak* peak_i = PhotoPeakFit(interest_x,1144,i_high);

        interest_map[j] = peak_i->GetSum(); 
        interest_map_err[j] = peak_i->GetSumErr(); 
        b4++; 
      }   
      if(b4>0){
        interest_fit_sum[5] += interest_map[j];  
        interest_fit_sum_err[5] +=pow(interest_map_err[j],2);
      }
    }
    I->Draw(0); 

  }
  TCanvas* N = new TCanvas("N","N",1500,800); 
  N->Divide(3,2);
  for(int j =0; j<detector_number;j++){

    TH1D* sumxn = norm->ProjectionX(Form("norm_x%02i",j),j+1,j+1); 
    
    TSpectrum s_n; 
    TH1* bg_n = s_n.Background(sumxn,iterations_of_bg,"BackOrder2"); 
    TH1D* singles_n_bg = (TH1D*)sumxn->Clone("singles_n_bg"); 
    singles_n_bg->Add(bg_n,-1); 

    TH1D* norm_x = singles_n_bg; 
 
    norm_x->GetXaxis()->SetRangeUser(nl_win,nh_win); 


    int b = 0; 
    //45 degrees
    if(j == 39 || j == 40) {

      N->cd(1);  
      if(j==39){
        norm_x->SetTitle("Normalization at 45 Degrees");
        norm_x->Draw();
      }
      else {
        norm_x->SetLineColor(b+1);
        norm_x->Draw("same");
      }
      if(norm_x->Integral()>100){
        GPeak* peak_n = PhotoPeakFit(norm_x,n_low,n_high);
        norm_map[j] = peak_n->GetSum(); 
        norm_map_err[j] = peak_n->GetSumErr(); 
        b++; 
      }   
      if(b>0){
        norm_fit_sum[0] += norm_map[j];  
        norm_fit_sum_err[0] +=pow(norm_map_err[j],2);
      } 
    }
    int b1 = 0;
    //86 degrees
    if(j == 1 || j == 2 || j == 17 || j == 18 || j == 34 || j  == 36 || j  == 43 || j == 44) {

      N->cd(2);      
      if(j==1){
        norm_x->SetTitle("Normalization at 86 Degrees");
        norm_x->Draw();
      }
      else {
        norm_x->SetLineColor(b1);
        norm_x->Draw("same");
      }
      if(norm_x->Integral()>100){
        GPeak* peak_n = PhotoPeakFit(norm_x,n_low,n_high);
        norm_map[j] = peak_n->GetSum(); 
        norm_map_err[j] = peak_n->GetSumErr(); 
        b1++;
      }   
      if(b1>0){
        norm_fit_sum[1] += norm_map[j];  
        norm_fit_sum_err[1] +=pow(norm_map_err[j],2);
      }
    }
    //90 degrees
    if(j == 23) {

      N->cd(3);      
      if(j==23){
        norm_x->SetTitle("Normalization at 90 Degrees");
        norm_x->Draw();
      }
      if(norm_x->Integral()>100){
        GPeak* peak_n = PhotoPeakFit(norm_x,n_low,n_high);
        norm_fit_sum[2] = peak_n->GetSum(); 
        norm_fit_sum_err[2] =pow(peak_n->GetSumErr(),2);
      }   

    }
    int b2 = 0; 
    //94 degrees
    if(j == 3 || j == 4 || j == 19 || j == 20 || j == 33 || j  == 35 || j  == 41 || j == 42) {

      N->cd(4);      
      if(j==3){
        norm_x->SetTitle("Normalization at 94 Degrees");
        norm_x->Draw();
      }
      else {   
        norm_x->SetLineColor(b2+1);
        norm_x->Draw("same");
      }
      if(norm_x->Integral()>100){
        GPeak* peak_n = PhotoPeakFit(norm_x,n_low,n_high);
        norm_map[j] = peak_n->GetSum(); 
        norm_map_err[j] = peak_n->GetSumErr(); 
        b2++; 
      }   
      if(b2>0){
        norm_fit_sum[3] += norm_map[j];  
        norm_fit_sum_err[3] +=pow(norm_map_err[j],2);
      }
    }
    int b3 = 0;
    //124 degrees
    if(j == 9 || j == 10 || j == 25 || j == 26) {

      N->cd(5);      
      if(j==9){
        norm_x->SetTitle("Normalization at 124 Degrees");
        norm_x->Draw();
      }
      else {
        norm_x->SetLineColor(b3+1);
        norm_x->Draw("same");
      }
      if(norm_x->Integral()>100){
        GPeak* peak_n = PhotoPeakFit(norm_x,n_low,n_high);
        norm_map[j] = peak_n->GetSum(); 
        norm_map_err[j] = peak_n->GetSumErr(); 
        b3++;
      }   
      if(b3>0){
        norm_fit_sum[4] += norm_map[j];  
        norm_fit_sum_err[4] +=pow(norm_map_err[j],2);
      }
    }
    int b4 = 0; 
    //134 degrees
    if(j == 11 || j == 12 || j == 27|| j == 28) {

      N->cd(6);      
      if(j==11){
        norm_x->SetTitle("Normalization at 134 Degrees");
        norm_x->Draw();
      }
      else {
        norm_x->SetLineColor(b4+1);
        norm_x->Draw("same");
      }
      if(norm_x->Integral()>100){
        GPeak* peak_n = PhotoPeakFit(norm_x,n_low,n_high);
        norm_map[j] = peak_n->GetSum(); 
        norm_map_err[j] = peak_n->GetSumErr(); 
        b4++;
      }   
      if(b4>0){
        norm_fit_sum[5] += norm_map[j];  
        norm_fit_sum_err[5] +=pow(norm_map_err[j],2);
      }
    }
    N->Draw(0); 

  }
   
  //set angular_map.first to an angle corresponding to the detector number 
  TGraph* angular_hist = new TGraphErrors(6);
  int pnt_num=0;
/*  for(int jj =0; jj<detector_number; jj++){
     printf(" Values of each detector for interest,normalization, and the ratio :\n"); 
     printf(" %.02f \t %.02f \t %.02f\n", interest_map[jj], norm_map[jj], (interest_map[jj]/norm_map[jj]));
  }*/
  for(int ii =0; ii< 6;ii++){ //ii< detector_number

    double ang = angle[ii]; //gmap[ii]

    if(ang==0)continue; 

    // printf(" Values at each angle summed after fitting the area, and ratio \n"); 
    printf(" %.02f \t %.02f \t %.02f  \t %.02f \t %.02f \t %.02f\n", angle[ii] ,interest_fit_sum[ii],pow(interest_fit_sum_err[ii],0.5), norm_fit_sum[ii],pow(norm_fit_sum_err[ii],0.5), (interest_fit_sum[ii]/norm_fit_sum[ii]));

    double interest_error = sqrt(pow(1/sqrt(interest_fit_sum[ii]),2)+pow(.05,2));
    double norm_error     = sqrt(pow(1/sqrt(norm_fit_sum[ii]),2)+pow(.05,2));
    double yerror         = sqrt(pow(interest_error,2)+pow(norm_error,2));

    printf("yerror for %i = %f\n",ii, yerror); 
    printf("\n"); 
    angular_hist->SetPoint(pnt_num,ang,(interest_fit_sum[ii]/norm_fit_sum[ii]));
    ((TGraphErrors*) angular_hist)->SetPointError(pnt_num,2*3.1415/180,yerror);
    pnt_num++;
  }

  return angular_hist;  
}












