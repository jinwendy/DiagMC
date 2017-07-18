//Structure of data profile
struct profile {

	double min;       //Lower boundary of the parameter
	double max;       //Upper boundary of the parameter
	int nbin;         //Number of bin subdividing the parameter
	int underflow;    //Weight entered above the upper boundary
	int overflow;     //Weight entered below the lower boundary
	//struct cmplx *c;        //Array of accumulated values
	//struct cmplx *c2;       //Array of accumulated squared values
	complex <double> *c;
	complex <double> *c2;
	int *n;          //Array of the number of entries per bin 
	int nentries;    //Total number of entries
};


//multi_dim profile
void prof_fill_multi(struct profile_multi *h,double x[],complex <double> w);
void prof_print_multi(FILE *stream,struct profile_multi *h);
void prof_print_multi_qw(FILE *stream,struct profile_multi *h);
struct profile_multi prof_four_multi(struct profile_multi *h, int isign);

void prof_bincoord(struct profile_multi *h,int bin, double x[]);

complex <double> read_bin(struct profile_multi *h, double x);
complex <double> read_binx(struct profile_multi *h, double x[]);

//==============================================================
struct profile prof_create(int n, double min, double max) {
	/*This function creates and initializes a profile structure
	with n bins, and for a parameter ranging from min to max. */
	struct profile rtn;
	if (n<0) {
		fprintf(stderr, "In prof_create, the number of bins must be >0.\n");
		exit(0);
	}
	if (max <= min) {
		fprintf(stderr, "In prof_create, we need to have max>min\n");
		exit(0);
	}
	rtn.max = max;
	rtn.min = min;
	rtn.nbin = n;
	rtn.underflow = 0;
	rtn.overflow = 0;
	rtn.nentries = 0;
	rtn.c = new complex <double>[n];
	rtn.c2 = new complex <double>[n];
	rtn.n = new int[n];

	complex <double> zero(0.0, 0.0);
	for (int i = 0; i<n; i++) {
		rtn.c[i] = zero;
		rtn.c2[i] = zero;
		rtn.n[i] = 0;
	}

	return rtn;
}
//==============================================================
void prof_fill(struct profile *h, double x, complex <double> w){
	/*Updates the averages in profile h for a new entry for
	parameter x with the complex value w*/
	h->nentries++; //Increment the total number of entries
	//If x out of range increment the corresponding counter
	if (x < h->min) { h->underflow++; /*= h->underflow + 1;*/ return; }
	if (x > h->max) { h->overflow++; /* = h->overflow + 1;*/ return; }

	//Identifies which bin to update, kth
	int k = int(floor(h->nbin*(x - h->min) / (h->max - h->min)));//find the kth bin, which this w falls into
	double n_onebin = h->n[k];//read out the number of datas in this particular bin

	//new Mean value
	h->c[k] = (h->c[k] * n_onebin + w) / (n_onebin + 1);
	//new Mean of the square 
	double re = (real(h->c2[k]) * n_onebin + real(w)*real(w)) / (n_onebin + 1);
	double im = (imag(h->c2[k]) * n_onebin + imag(w)*imag(w)) / (n_onebin + 1);
	complex <double> temp(re, im);
	h->c2[k] = temp;

	//Count the number of entries in that bin
	h->n[k] = h->n[k] + 1;

	return;
}
//==============================================================
void prof_print(FILE *stream, struct profile h) {
	/*Print data in ptofile h to channel stream*/
	//  fprintf(stderr,"MIN=%f\nMAX=%f\n",h.min,h.max);
	//fprintf(stderr,"NBIN=%d\n",h.nbin);
	for (int i = 0; i<h.nbin; i++){
		double sr = real(h.c2[i]) - real(h.c[i])*real(h.c[i]);//discrepancy,real
		if (sr>0) sr = sqrt(sr); else sr = 0;
		double si = imag(h.c2[i]) - imag(h.c[i])*imag(h.c[i]);//discrepancy,imaginary part
		if (si > 0) si = sqrt(si); else si = 0;
		fprintf(stream, "%3d %d %g %g %g %g %g %g\n", i + 1, h.n[i],
			h.min + (i + 0.5)*(h.max - h.min) / h.nbin,
			0.5*(h.max - h.min) / h.nbin,
			real(h.c[i]), imag(h.c[i]), sr, si);
	}
	//  fprintf(stderr,"UNDERFLOW=%d\nOVERFLOW=%d\n",h.underflow,h.overflow);
	// fprintf(stderr,"ENTRIES:%d\n",h.nentries);
	return;
}
//==============================================================

//==============================================================
void prof_fill_multi(struct profile_multi *h,double x[],complex <double> w){
  /*Updates the averages in profile h for a new entry for 
    parameter x with the value w*/
  double n_onebin = 0.0;
  h->nentries++; //Increment the total number of entries
  double re = 0.0;
  double im = 0.0;
  
  //If x out of range increment the corresponding counter
  for(int i=0;i<h->dim;i++) 
    if(x[i]<h->min[i] || x[i]>h->max[i]) 
      {h->outrange=h->outrange+1;return;}
  
  //Identifies the bin to update
  int k=0;
  int pwr=1;
  for(int i=h->dim-1;i>=0;i--){//starts from the most right dimension eg: for (r,t), starts from t
    int j=int(floor(h->nbin[i]*(x[i]-h->min[i])/(h->max[i]-h->min[i]))); //the bin order in each dimension (interval = max-min/nbin)
    k=k+pwr*j;
    pwr=pwr*h->nbin[i];//for next dimension
  }//end of loop
  
  n_onebin = h->n[k];//read out the number of datas in this particular bin

  //new Mean value
  //h->c[k] = h->c[k]+ w;
  h->c[k] = (h->c[k] * n_onebin+ w) / (n_onebin+ 1.0);
  
  //new Mean of the square 
  re = (real(h->c2[k]) * n_onebin + real(w)*real(w)) / (n_onebin + 1.0);//update the values in the corresponding bin
  im = (imag(h->c2[k]) * n_onebin + imag(w)*imag(w)) / (n_onebin + 1.0);
  //re = (real(h->c2[k]) + real(w)*real(w));//update the values in the corresponding bin
  //im = (imag(h->c2[k]) + imag(w)*imag(w));
  
  complex <double> temp(re, im);
  h->c2[k] = temp;
  
  //Count the number of entries in that bin
  h->n[k]=h->n[k]+1;
  
 // printf("real: %f", real (h->c[k]));
  
  return;
}
//==============================================================
void prof_print_multi(FILE *stream,struct profile_multi *h) {
  /*Print data in ptofile h to channel stream*/

  //Loop over all the bins
  for(int i=0;i<h->ntot;i++){
    //Standard deviation over the real part
    double sr=real(h->c2[i])-real(h->c[i])*real(h->c[i]); 
    if(sr>0) sr=sqrt(sr); else sr=0;
    //Standard deviation over the imaginary part
    double si=imag(h->c2[i])-imag(h->c[i])*imag(h->c[i]);
    if(si>0) si=sqrt(si); else si=0;


    //Print the bin number and the number of entries in that bin
    fprintf(stream, "%d %d ",i+1,h->n[i]);
    //Print the bin coordinate
    double bincrd[h->dim];prof_bincoord(h,i,bincrd);
    for(int d=0;d<h->dim;d++) fprintf(stream,"%g ",bincrd[d]);
    
    //Print the bin value and standard deviation
    fprintf(stream,"%g %g %g %g \n",real(h->c[i]),imag(h->c[i]),sr,si);
  }
  
  return;
}
//==============================================================
void prof_print_multi_qw(FILE *stream,struct profile_multi *h) {
  /*Print data in ptofile h to channel stream, and the data is in (q,omega) space, got after the FT.*/

  //Loop over all the bins
  for(int i=0;i<h->ntot;i++){
    //Standard deviation over the real part
    double sr=real(h->c2[i])-real(h->c[i])*real(h->c[i]); 
    if(sr>0) sr=sqrt(sr); else sr=0;
    //Standard deviation over the imaginary part
    double si=imag(h->c2[i])-imag(h->c[i])*imag(h->c[i]);
    if(si>0) si=sqrt(si); else si=0;

    //Print the bin number and the number of entries in that bin
    fprintf(stream, "%d %d ",i+1,h->n[i]);
		
    //Print the bin coordinate
    double bincrd[h->dim];prof_bincoord(h,i,bincrd);
    for(int d=0;d<h->dim;d++) fprintf(stream,"%f ",bincrd[d]);
    
    //Print the bin value and standard deviation
    fprintf(stream,"%g %g %g %g \n",real(h->c[i]),imag(h->c[i]),sr,si);
  }
  
  return;
}


//===============================================================
complex <double> read_bin(struct profile_multi *h, double x){//for 1D profile only!
	/*This function takes a double value, (time) and a data profile as arguments. And returns the data value corresponds to that time intervale x.*/
	complex <double> rtn;
	if (x>h->max[0] ||x<h->min[0]) printf("outside range\n");
	//Identifies which bin to update, kth
	int k = int(floor(h->ntot*(x - h->min[0]) / (h->max[0] - h->min[0])));//find the kth bin, which this x falls into

	rtn = h->c[k]; //read out the value in kth bin
	return rtn;
}
//===============================================================
complex <double> read_binx(struct profile_multi *h, double x[]){//for 2D and above profile
	/*This function takes a double array (containing the position and time) and a data profile as arguments. 
	And returns the data value corresponds to that time intervale t and position x.*/
	complex <double> rtn;
	
	//int n = sizeof(&x) / sizeof(int); //length of the array x
	//printf("length of array x:%d\n",n);
	
	//if (n != h->dim) errorquit("not a instant interaction dimension in read_binx");//profile h and double x have the same dimensions
	
	//Identifies which bin to update, kth in each direction
	int dim = h->dim;
	//printf("dim :%d\n",dim);
	int k[] = {0,0,0,0};
	int kt = 1;
	for (int i=0; i<dim; i++){
		k[i]=int(floor(h->nbin[i]*(x[i] - h->min[i]) / (h->max[i] - h->min[i])));//find the kth in position, which this x[i] falls into
		//printf(" k[i]:%d\n",k[i]);
		kt = kt*k[i];//bin number in the profile
	}
	
	rtn = h->c[kt]; //read out the value in kth bin
	
	return rtn;
}
//===============================================================