void GraphErrorTable(Int_t iMark, Int_t iColor) {

  vector<float> x;
  vector<float> y;
  vector<float> ex;
  vector<float> ey;

  Int_t npt = 0;
  float a,b,ea,eb;

  // read data file
  ifstream in;
  in.open("dat.txt");

  while ( kTRUE ) {

    in >> a >> ea >> b >> eb;
    x.push_back(a);
    ex.push_back(ea);
    y.push_back(b);
    ey.push_back(eb);

    if ( ! in.good() ) break;

      //cout << "x = " << x[npt] << " y = " << y[npt] << endl;

    npt++;

  }

  const Int_t ndim = npt;
  Int_t index = 0;

  Double_t xx[ndim],yy[ndim],exx[ndim],eyy[ndim];

  for (Int_t i=0;i<npt;i++) {

    if (std::abs(y[i]) > 0.5 ) { continue; }
    if (std::abs(x[i]) > 0.3 ) { continue; }
    if ( x[i] < 0.08 || x[i] > 0.09 ) { continue; }
    if (ey[i] > 0.3 ) { y[i] = 0.; ey[i] = 0.; }

    xx[index] = x[i];
    yy[index] = y[i];
    exx[index] = ex[i];
    eyy[index] = ey[i];
    //cout << "# " << i+1 << " x = " << xx[i] << " ex = " << exx[i] << " y = " << yy[i] << " ey = " << eyy[i] << endl;
    index++;
  }   
  
  in.close();

  npt = index-1;
  printf(" found %d points\n", npt);
 
  TGraphErrors *gr1 = new TGraphErrors (npt, xx, yy, exx, eyy); 

  TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 500, 500); 

  c1->cd(); 

  gr1->SetMarkerStyle(iMark);
  gr1->SetMarkerColor(iColor);
  gr1->Draw("AP");

  c1->Print("plot.pdf"); 

}
