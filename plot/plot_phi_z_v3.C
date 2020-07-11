//

TCanvas* c_plot_hits_phiz(0);

//-----------------------------------------------------------------------------
//double 

//-----------------------------------------------------------------------------
void plot_phi_z_v3(const char* Fn) {

  float  x[1000], y[1000], rho[1000], z[1000], phi[1000], phi_fh[1000];
  int    flags[1000], used[1000], done(0), hitIndex, strawId;

  int    parameters_read(0);

  float   X0, Y0, R, chi2, phi_0, dfdz, chi2N;
  float dphi, xdphi[1000], zl, dz, cdfdz, sdfdz, cchi2;

  int np(0), ngood(0), index;

  char c[1000], name[100], s[100];

  FILE* f  = fopen(Fn,"r");

  if (f == 0) {
    Error("plot_hits",Form("missing file %s\n",Fn));
    return -2;
  }

  while ( ((c[0]=getc(f)) != EOF) && !done) {

					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);

      if (parameters_read == 0) {
	parameters_read = 1;
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );

	fscanf(f,"%f",&phi_0);

	fscanf(f,"%s",name);
	fscanf(f,"%s",name);

	fscanf(f,"%f",&dfdz);

	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );

// 	fscanf(f,"%s =",name);
// 	fscanf(f,"%f",&chi2N);

// 	printf("X0 = %12.5f Y0 = %12.5f R = %12.5f  chi2 = %12.5e \n",
// 	       X0, Y0, R, chi2);

      }
      else {
					// read points
	fscanf(f,"%i"  ,&index);
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );

	fscanf(f,"%f"  ,&z[np] );
	//	fscanf(f,"%f"  ,&dz);
	//	fscanf(f,"%i"  ,&dphi);
	fscanf(f,"%f"  ,&phi[np]);

	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%s"  ,name  );

 	// fscanf(f,"%f"  ,&cdfdz);
 	// fscanf(f,"%f"  ,&dphi);
 	// fscanf(f,"%f"  ,&xdphi[np]);
 	// fscanf(f,"%f"  ,&cdfdz);
	//	fscanf(f,"%f"  ,&phi_fh[np]);
	
// 	printf("name = %s np=%3i flags=%08x %4i x[np], y[np], z[np] = %10.3f %10.3f %10.3f \n",
// 	       name,np,flags[np], flags[np],x[np],y[np],z[np]);
	if (flags[np] < 256) ngood++;
	np++;
      }

    }
					// skip line
    fgets(c,100,f);
  }

  fclose(f);


  dfdz = 1./dfdz;

  //  printf(">> np = %i, ngood = %3i\n",np,ngood);

  //  TGraph* gr_phiz = new TGraph(np,z,phi);

  if (c_plot_hits_phiz != 0) delete c_plot_hits_phiz;

  c_plot_hits_phiz = new TCanvas(Form("c_plot_hits_phiz_%s",name),Form("c_phiz %s",name),1200,600);
//-----------------------------------------------------------------------------
// plot PHI-Z picture
//-----------------------------------------------------------------------------
  c_plot_hits_phiz->cd(1)->SetGrid();

  TH2F* h2_phiz = new TH2F("h2_phiz","phiZ VIEW; z [mm]; helix-#phi [rad]",2200,-1600,2800,80,-20,20);
  h2_phiz->SetStats(0);
  h2_phiz->Draw();

  float yfit[200];
  int color;

  for (int i=0; i<np; i++) {
    TMarker* m = new TMarker(z[i],phi[i],2);
    //    if (flags[i] >= 256) color = kBlack;
    //    if ( (used[i] < 1)  || (xdphi[i] > 10.0) )   color = kBlack;
    //    else
    color = kRed;
    m->SetMarkerColor(color);
    m->SetMarkerSize(0.7);
    m->Draw();
  }

  TF1 *yf = new TF1("yf","[0]+x*[1]",-1600., 2800.);
  yf->SetParameters(phi_0, dfdz);
  yf->Draw("same");
}
