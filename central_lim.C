/***
This ROOT macro serves as a demonstration of the Central Limit Theorem.
In probability theory, the Central Limit Theorem states that, if you have a
random variable X which follows some arbitrary probability distribution P(X),
if you get experimental values of X and then average them, and then repeat
this process, the resulting distribution of average X's will converge to the 
normal distribution.

This macro lets one to define an arbitrary probability density function (PDF), 
so long as it is never negative. Normalization is handled by the
code. Then, a random variable with this probability distribution is simulated
via discrete inverse transform sampling. The steps for this inverse transform
sampling are as follows:

1. Construct a discrete approximation to the cumulative distribution function
(CDF) by calculating Riemann sums over the PDF. The CDF is just the integral
of the PDF and will range from 0 to 1.

2. Generate a random value between 0 and 1 according to a uniform distribution.

3. Match the closest y-value of the CDF to this random value. Return the
corresponding x-value of the CDF. This is the result of one experiment.

Now that we can simulate a random variable using a probability distribution
of our own definition, we can run the experiment a set number of times and
then average the value, then put that average value into a histogram. We then
repeat many times until our histogram inevitably resembles a normal
distribution.
    
 ***/


/** Probability density function
   
    Cannot be negative, as probability cannot have a negative value.
    Several examples are given in comment blocks. You may also define
    your own, but an improper PDF will throw an error.
 **/
Double_t pdf(Double_t x) {
  //return TMath::Exp((-1*(x-5)*(x-5))/8); //A normal distribution with mean 5, stdev 2
  //return sin(x)+10;  //A bimodal distribution 
  //return x*x;        //A parabolic distribution. 
  //return 1.0;        //A uniform distribution
  //{if (x < 1) return 0; if (x < 5) return 1; return 0;}  //A step distribution
  return TMath::Exp(x)*TMath::Exp((-1*(x-3)*(x-3))/10);  //A lopsided distribution where the most probable value is not the same as the expectation value.
}

Double_t average_value(TF1 *func) {
  double x1 = func->GetXmin();
  double x2 = func->GetXmax();
  return (func->Integral(x1,x2))/(x2-x1);
}

/**
   Main code of the macro.
 **/
void central_lim() { 
  Double_t xmin = 0.;       //Minimum x value of the PDF 
  Double_t xmax = 10.;      //Maximum x value of the PDF

  int numpoints = 1000;      //Points to use in our discrete approximations
  int num_iterations = 100; //Times to run the experiment before averaging 
  int num_means = 10000;     //Number of averages to take


  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  //Create a TF1 from our PDF, normalize it, and ensure it's never negative.
  TF1 *Pdf_t = new TF1("prob_density","[0]*pdf(x)",xmin,xmax);
  Pdf_t->SetParameter(0,1.0);
  double norm_const = 1/(Pdf_t->Integral(xmin,xmax));
  Pdf_t->SetParameter(0,norm_const);
  if (Pdf_t->GetMinimum() < 0.) {
    cout << "Error: Probability density function must always be positive." << endl;
    throw("");
  }
  //Draw our PDF
  Pdf_t->SetTitle("Probability Density Function");
  Pdf_t->Draw();

  
  //Here we create a discrete approximation to the CDF.
  //Calculate points of the CDF by 'integrating' PDF with Riemann sums
  double cdf_x[numpoints];
  double cdf_y[numpoints];

  double area = 0.;
  double x = xmin;
  const double dx = (xmax-xmin)/numpoints;
  for (int i = 0; i < numpoints; i++) {
    area += (Pdf_t->Eval(x))*dx;
    cdf_x[i] = x;
    cdf_y[i] = area;
    x += dx;
  }

  //Plot the CDF
  TGraph *cdf_graph = new TGraph(numpoints,cdf_x,cdf_y);
  c1->cd(2);
  cdf_graph->SetTitle("Cumulative Distribution Function");
  cdf_graph->Draw();

  
  //Here we simulate our random variable using inverse transform sampling.
  TGraph *inv_cdf_graph = new TGraph(numpoints,cdf_y,cdf_x);
  TH1F *trials = new TH1F("h1","trials",numpoints,xmin,xmax);
    
  TRandom3 *r = new TRandom3(1234); //A uniform random variable from 0 to 1
  double rand_val;
  for (int i = 0; i < num_means; i++) {
    TH1F *trial = new TH1F("trial","trial",numpoints/10,xmin,xmax);
    //Run the experiment 'num_iterations' times, then take the average
    for (int j = 0; j < num_iterations; j++) {
      rand_val = r->Uniform(0,1);
      //We find the y-value of our CDF which is closest to rand_var.
      for (int k = 0; k < numpoints-1; k++) {
	if ((cdf_y[k] <= rand_val) and (rand_val <= cdf_y[k+1])) {
	  //The corresponding x-value to this y-value is the experiment result.
	  trial->Fill(cdf_x[k]);
	  break;
	}
      }
    }
    //Now we take the average value of our experiments,
    //and we fill our final histogram with that value.
    trials->Fill(trial->GetMean());
    trial->Delete();
  }
  
  //Drawing the final resulting histogram
  c1->cd(3);
  trials->SetTitle("Average Results");
  trials->Draw();


  //Here we use the PDF to predict the mean and standard deviation of the
  //average results according to the central limit theorem.
  //--The mean should be equal to that of the PDF.
  //--The standard deviation should be equal to that of the PDF divided by the
  //square root of the number of trials.
  TF1 *x_func = new TF1("x*prob_density","[0]*pdf(x)*x",xmin,xmax);
  x_func->SetParameter(0,norm_const);
  TF1 *x2_func =  new TF1("x^2*prob_density","[0]*pdf(x)*x*x",xmin,xmax);
  x2_func->SetParameter(0,norm_const);

  double x_av = average_value(x_func);
  double x2_av = average_value(x2_func);

  double stdev = TMath::Sqrt(x2_av - (x_av)*(x_av));
  stdev /= TMath::Sqrt(num_iterations);

  cout << endl;
  cout << "STATISTICS" << endl;
  cout << "\tPredicted" << "\t" << "Experimental" << endl;
  cout << "Mean: \t" << Pdf_t->Mean(xmin,xmax) << "\t\t" << trials->GetMean() << endl;;
  cout << "Stdev: \t" << stdev << "\t" << trials->GetStdDev() << endl << endl;

  
  c1->Draw();
  c1->SaveAs("plot.png");
}
