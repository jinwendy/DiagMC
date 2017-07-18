/*Functions to do FFT*/
 

#define N_rescale 2000 //The uplimit of the rescale factor
#define R_wf  1.0 //The ratio between ~W and F

//Structure for FT
struct ft_comp {
	int nentries;      //Number of entries in the coefficient evaluation
	int dim;           //Number of components
	complex <double> *a;   //Array of complex numbers
};



complex <double> init_complex(complex <double> c);
struct ft_comp decomp_init(int n);
int Fourier_T(struct ft_comp *ft, complex <double> c, double t);

struct profile prof_create(int n, double min, double max);
void prof_fill(struct profile *h, double x, complex <double> w);
void prof_print(FILE *stream, struct profile h);
struct profile prof_four(struct profile s, int isign);
void fourn(double data[], unsigned long nn[], int ndim, int isign);

struct profile_multi minusbeta(struct profile_multi *h);
struct profile_multi betaback(struct profile_multi *h);



void init_profile_multi(struct profile_multi *h);


struct profile_multi bold_green(struct profile_multi *bare, struct profile_multi *self);
struct profile_multi wequation(struct profile_multi *bare, struct profile_multi *polar);
struct profile_multi wtilde(struct profile_multi *bare, struct profile_multi *w);
struct profile_multi chi(struct profile_multi *bare, struct profile_multi *polar);
void newf(struct profile_multi *wtilde, struct profile_multi *unphys_f);
void sum_rule(struct profile_multi ch, struct profile_multi *pi, struct profile_multi *j);
void correct_j(struct profile_multi *jft);


void complement_profile_multi(struct profile_multi *sp);

complex <double> sum_up(struct profile_multi *h);
complex <double> lambda(struct profile_multi *j, struct profile_multi *polar, double lambda);
double rescale(struct profile_multi *j, struct profile_multi *polar);
double Qn2(struct profile_multi *green);

void times_profile(struct profile *h, double d);

void scale_profile_multi(struct profile_multi *h, double d);
void normalize_profile(struct profile_multi *h, int har_or_g, double d);
void new_normalize_profile(struct profile_multi *sp, struct profile_multi *hg, double d);

struct profile bg(struct profile *self);
struct profile bg(struct profile *g, struct profile *self);
struct profile bw(struct profile *polar);
struct profile wtilde(struct profile *inter);


//Functions to do integration 
double tiner(int n);
int reverse_tiner(double t);
struct profile ite_g(struct profile *oldg, struct profile *self);

void print_profile(struct profile_multi *h);

//============================================================
complex <double> init_complex(complex <double> c){
	/*This function initializes a complex number, real part=imag part=0.*/
	std::complex<double> zero0 (0,0);
	c = zero0;
	return c;
}
//==========================================================================
struct ft_comp decomp_init(int n) {
	/*Initializes an instance of decomp using 2n+1 as the number of complex
	coeffcients.*/
	struct ft_comp rtn;
	rtn.nentries = 0;
	rtn.dim = 2 * n + 1;
	rtn.a = new complex <double>[rtn.dim]; //Fourier components (2n+1)
	for (int i = 0; i < rtn.dim; i++){
		rtn.a[i] = init_complex(rtn.a[i]);
		//printf("ith component %d\n", i);
		//print_cmplx(rtn.a[i]);
	}
	return rtn;
}
//==========================================================================
int Fourier_T(struct ft_comp *ft, complex <double> c, double time) {
	/*For the G sectr.This function updates the fourier decompositions *ft for a complex
	datapoint c at parameter t with t in the interval [-1,1].*/
	if (time<-TEM || time>TEM) return DECOMP_NOTOK;
	ft->nentries++; //Increment the number of entries
	int n = (ft->dim - 1) / 2;
	for (int i = -n; i <= n; i++){
		double wm = DECOMP_PI*(i);
		double phi = 1.0*wm*time;
		double cphi = cos(phi);
		double sphi = sin(phi);
		complex <double> temp(cos(phi), sin(phi));
		ft->a[i + n] = ft->a[i + n] + c*temp;
	}
	return DECOMP_OK;
}
//==========================================================================


/*Fast Fourier Transformation_new.*/
//==============================================================
struct profile prof_four(struct profile s, int isign) {
	/*Create a profile with the fourier transform of in when isign=1
	and the inverse-fourier transform for i=-1*/

	//Create an array to be used to compute the fourier transform
	double *d;
	d = new double[2 * s.nbin];
	for (int i = 0; i<s.nbin; i++) {
		d[2 * i] = real(s.c[i]);     //Real part
		d[2 * i + 1] = imag(s.c[i]); //Imaginary part
		//printf("%lf  %lf\n", d[2*i],d[2*i+1] );
	}
	
	//Call the FFT algorithm
	unsigned long ndim[1]; ndim[0] = s.nbin; /*Array containing the numbers
											 of sample in the different
											 directions. */
	fourn(d, ndim, 1, isign);//FT or inverse FT
	
	//Store the result in a profile
	struct profile rtn;
	rtn = prof_create(s.nbin, s.min, s.max);
	for (int i = 0; i<s.nbin; i++) {
		double nbin = s.nbin;
		printf("nbin:%f\n",nbin);
		double re = 0.0;
		double im = 0.0;
		if(isign == 1){
			printf("FT\n");	
			re = d[2 * i];
			im = d[2 * i + 1];
		}
		else{//inverse FT, the data is acutal one times the product of the lengths of all dimensions.
			 re = d[2 * i]/nbin;
			 im = d[2 * i + 1]/nbin;
		 }
			
		complex <double> temp(re, im);
		rtn.c[i] = temp;
	}

	return rtn;
}
//==============================================================
void fourn(double data[], unsigned long nn[], int ndim, int isign) {
	/*
	Replaces data by its ndim-dimensional discrete Fourier transform,
	if isign is input as 1.
	nn[1..ndim] is an integer array containing the lengths of each dimension
	(number of complex values), which MUST all be powers of 2.
	data is a real array of length twice the product of these lengths,
	in which the data are stored as in a multidimensional complex array:
	real and imaginary parts of each element are in consecutive locations,
	and the rightmost index of the array increases most rapidly as one
	proceeds along data.
	For a two-dimensional array, this is equivalent to storing the array by
	rows.
	If isign is input as -1, data is replaced by its inverse transform times
	the product of the lengths of all dimensions.
	*/


	int idim;
	unsigned long i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
	unsigned long ibit, k1, k2, n, nprev, nrem, ntot;
	double tempi, tempr;
	//Double precision for trigonometric references
	double theta, wi, wpi, wpr, wr, wtemp;
	//Compute the total number of complexes values
	for (ntot = 1, idim = 0; idim<ndim; idim++)
		ntot *= nn[idim];


	nprev = 1;
	//Main loop over the dimensions
	for (idim = ndim - 1; idim >= 0; idim--) {
		n = nn[idim];
		nrem = ntot / (n*nprev);
		ip1 = nprev << 1;
		ip2 = ip1*n;
		ip3 = ip2*nrem;
		i2rev = 0;
		//This is the bit-reversal section of the routine.
		for (i2 = 0; i2<ip2; i2 += ip1) {
			if (i2 < i2rev) {
				for (i1 = i2; i1<i2 + ip1 - 1; i1 += 2) {
					for (i3 = i1; i3<ip3; i3 += ip2) {
						i3rev = i2rev + i3 - i2;
						SWAP(data[i3], data[i3rev]);
						SWAP(data[i3 + 1], data[i3rev + 1]);
					}
				}
			}
			ibit = ip2 >> 1;
			while (ibit >= ip1 && i2rev + 1 > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		//Here begins the Danielson-Lanczos section of the routine.
		ifp1 = ip1;
		while (ifp1 < ip2) {
			ifp2 = ifp1 << 1;
			//Initialize for the trig. recurrence.
			theta = isign*6.28318530717959 / (ifp2 / ip1);
			wtemp = sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi = sin(theta);
			wr = 1.0;
			wi = 0.0;
			for (i3 = 0; i3<ifp1; i3 += ip1) {
				for (i1 = i3; i1<i3 + ip1 - 1; i1 += 2) {
					for (i2 = i1; i2<ip3; i2 += ifp2) {
						//Danielson-Lanczos formula:
						k1 = i2;
						k2 = k1 + ifp1;
						tempr = (double)wr*data[k2] - (double)wi*data[k2 + 1];
						tempi = (double)wr*data[k2 + 1] + (double)wi*data[k2];
						data[k2] = data[k1] - tempr;
						data[k2 + 1] = data[k1 + 1] - tempi;
						data[k1] += tempr;
						data[k1 + 1] += tempi;
					}
				}
				//Trigonometric recurrence.
				wr = (wtemp = wr)*wpr - wi*wpi + wr;
				wi = wi*wpr + wtemp*wpi + wi;
			}
			ifp1 = ifp2;
		}
		nprev *= n;
	}
}
//==============================================================


/*Multi-dimensional profile*/

//==============================================================
struct profile_multi prof_four_multi(struct profile_multi *h, int isign) {
	/*Create a profile with the fourier transform of in when isign=1
	and the inverse-fourier transform for i=-1*/
	
    //Create the various arrays needed to run fourn
    int ndim=h->dim;
    unsigned long nn[ndim];
    for(int i=0;i<ndim;i++) nn[i]=h->nbin[i];//number of bins in each dimension
    double data[2*h->ntot];//a data array to store the c[i] in the file, length of the data is 2*#of c[i], both imag and real part

    for(int i=0;i<h->ntot;i++){//each complex value occupies two sequential locations
      data[2*i] = real(h->c[i]);//real part
      data[2*i+1] = imag(h->c[i]);//followed by imaginary part
    }

	
	//Call the FFT algorithm
	fourn(data, nn, ndim, isign);/* if isign=1 :replaces data by its ndim-dimensional discrete FT
	                                if isign=-1:replaces data by its inverse FT, times the product of
		                                        the length of all dimensions.*/
	//Store the result in a profile
	struct profile_multi rtn;
	rtn = prof_create_multi(h->dim, h->nbin, h->min, h->max);
	//for(int i=0;i<ndim;i++) printf("min: %g\n",h->min[i]);
	
	for (int i = 0; i<h->ntot; i++) {
		
		double nbin = h->ntot;//total number of bins
		
		double re = 0.0;
		double im = 0.0;
		if(isign == 1){
			re = data[2 * i];
			im = data[2 * i + 1];
		}
		else{//inverse FT, the data is acutal one times the product of the lengths of all dimensions.
			
			//if (i==0) printf("length of all dimensions: %f\n",nbin);
	
			 re = data[2 * i]/nbin;
			 im = data[2 * i + 1]/nbin;
			 //printf("Data %f Re :%f\n",data[2 * i],re);
		 }
		
		rtn.c[i] = assign_complx(re,im);
		rtn.c2[i]= assign_complx(0.0,0.0);
	    rtn.n[i]=1;
	    rtn.nentries++;
	}

	return rtn;
}
//==============================================================

//==============================================================
void prof_bincoord(struct profile_multi *h,int bin, double x[]){
  /*
    Array x is made to contain the coordinates if the n-th 
    bin in profile *h
  */

  //Products of the successive numbers of bins in the different dimensions
  int *pwr; pwr=new int[h->dim];
  pwr[h->dim-1]=1;
  for(int i=h->dim-2;i>=0;i--){
    if(i>=0)    pwr[i]=pwr[i+1]*h->nbin[i+1];
  }
  
  //Integer coordinates of the bin
  int *k; k=new int[h->dim];
  int icp=bin;
  for(int d=0;d<h->dim;d++) { //Reversed loop
    k[d]=icp/pwr[d]; //Dimension d coordinate
    icp=icp%pwr[d];  //Remainder for other coordinates
    x[d]=h->min[d]+(k[d]+0.5)*(h->max[d]-h->min[d])/h->nbin[d];
  }

  return;
}
//==============================================================
//====================================================
void init_profile_multi(struct profile_multi *h){
	/*This function initialize a profile.*/
	for (int i = 0; i < h->ntot; i++){
		//double n = h->n[i];
		complex <double> zero(0.0,0.0);
		h->n[i] = 0;
		h->c[i] = zero;
		h->c2[i] = zero;
	}
	return;
}


/*Functions to manipulate the profile, save information for FFT*/
//===============================================================
/*Dyson equation*/
struct profile_multi bold_green(struct profile_multi *bare, struct profile_multi *self){
	/*This function calculate the frequence dependent bold green's function by Dyson equation.*/
	
	//create a 1D profile
	int nft[1] = {2*N_t};
	double min[1] = {-BETA};
	double max[1] = {BETA};
	
	struct profile_multi rtn = prof_create_multi(1, nft, min, max);
	
	for (int i = 0; i < bare->ntot; i++){/*ith Fourier component of bold green's function by Dyson*/
		rtn.c[i] = bare->c[i] / (1.0 - bare->c[i] * self->c[i]);
		rtn.c2[i] = self->c2[i];
		//print_cmplx(rtn.c2[i]);
	}
	return rtn;
}
//===============================================================
/*W equation (7 ):first part in PRB 87,024407(2013)*/
struct profile_multi wequation(struct profile_multi *bare, struct profile_multi *polar){
	/*This function calculate the frequence dependent W equation.*/
	
	//create a 2D profile
	int nft[2] = {N, N_t};
	double min[2] = {XL, Tmin};
	double max[2] = {XU, BETA};
	
	struct profile_multi rtn = prof_create_multi(2, nft, min, max);
	
	for (int i = 0; i < bare->ntot; i++){/*ith Fourier component of bold green's function by Dyson*/
		rtn.c[i] = bare->c[i] / (4.0 - bare->c[i] * polar->c[i]);
		rtn.c2[i] = polar->c2[i];
	}
	return rtn;
}
//===============================================================
/*W equation (7 ):second part in PRB 87,024407(2013)*/
struct profile_multi wtilde(struct profile_multi *bare, struct profile_multi *w){
	/*This function calculate the frequence dependent ~W.*/
	
	//create a 2D profile
	int nft[2] = {N, N_t};
	double min[2] = {XL, Tmin};
	double max[2] = {XU, BETA};
	
	struct profile_multi rtn = prof_create_multi(2, nft, min, max);
	
	for (int i = 0; i < bare->ntot; i++){/*ith Fourier component of bold green's function by Dyson*/
		rtn.c[i] = w->c[i] - (bare->c[i]/4.0);
		rtn.c2[i] = w->c2[i];
	}
	return rtn;
}
//===============================================================
/*Chi equation (8) in PRB 87,024407(2013)*/
struct profile_multi chi(struct profile_multi *bare, struct profile_multi *polar){
	/*This function calculate the frequence dependent chi.*/
	if (bare->dim != polar->dim) errorquit("J and PI do not have same dimensions");
	
	int d = bare->dim;
	
	//create a 2D profile
	int nft[2] = {N, N_t};
	double min[2] = {XL, Tmin};
	double max[2] = {XU, BETA};
	/*int nft[d] = new int[d];	
	double min[d] = new double[d];
	double max[d] = new double[d];
	
	nft[d] = bare->nbin;
	min[d] = bare->min;
	max[d] = bare->max;*/
	
	struct profile_multi rtn = prof_create_multi(2, nft, min, max);
	
	for (int i = 0; i < bare->ntot; i++){/*ith Fourier component of bold green's function by Dyson*/
		rtn.c[i] = polar->c[i] / (4.0 - bare->c[i] * polar->c[i]);
		rtn.c2[i] = polar->c2[i];
	}
	return rtn;
}
//===============================================================
void newf(struct profile_multi *wtilde, struct profile_multi *unphys_f){
	/*This function update the F-profile according to the one of W.*/
	if (unphys_f->ntot != wtilde->ntot ) errorquit("F and ~W do not have same dimensions");//If they have different dimensions
	double factor = (double) R_wf;
	
	complex <double> sum(0.0,0.0);
	complex <double> average(0.0,0.0);
	
	for (int i = 0; i < wtilde->ntot; i++){
		sum = sum + wtilde->c[i];
	}
	average = sum/double(wtilde->ntot);
	
	for (int i = 0; i < unphys_f->ntot; i++){//update the F-profile
		if (abs(wtilde->c[i])< 1E-3) printf("newf\n");//unphys_f->c[i]= average;
		else{
			//printf("there : %d\n",i);
			//print_cmplx(wtilde->c[i]);
			unphys_f->c[i] = wtilde->c[i];
			//print_cmplx(unphys_f->c[i]);
			unphys_f->c2[i] = wtilde->c2[i];
			}
	}
	return;
}
//===============================================================
void sum_rule(struct profile_multi ch, struct profile_multi *pi, struct profile_multi *j){
	/*After solving the Dyson equation, we check the value of Chi(r=0,t=0), adjust the scaling factor
	 for Pi, and go back to solving the Dyson equation again until the sum rule is satisfied with three-digit accuracy.*/
	/*Chi  is the profile of Chi in (r,t) space*/
	complex <double> ch0 = ch.c[0];//read Chi(r=0,t=0) for the profile, which is the first element
	
	double factor = 1.5;
	int number = 0;
	
	//print_cmplx(ch.c[0]);
	
	struct profile_multi chi_qw = chi(j, pi);
	
	while (number < 20) {
		printf("=====\n");
		printf("number: %d\n",number);
		scale_profile_multi(pi, factor);
		chi_qw = chi(j,pi);
		ch = prof_four_multi(&chi_qw, -1);//Inverse FT
		print_cmplx(ch.c[0]);
		
		printf("diff:%g\n",fabs(imag(ch.c[0])-SUM));
		number++;
	}
	
	//print_cmplx(ch.c[0]);
return;	
	
}
//===============================================================
void correct_j(struct profile_multi *jft){
	/*This function correct the expression of J(q,omega), which is independent of Omega.*/
	double *data;
	data = new double[2*N];
	for (int i = 0; i< N; i++) {
		data[2 * i] = (double)i*J0;//Real part
		data[2 * i + 1] = 0.0; //Imaginary part
		//printf("%d: %lf  %d:%lf\n", 2*i, data[2*i], 2*i+1,data[2*i+1] );
	}
	unsigned long num[1]; num[0] = N;
	complex <double> temp(0.0, 0.0);
	
	fourn(data, num, 1, 1);
	
	for (int j = 0; j < jft->ntot; j++){
	 if (j < N_t) temp = assign_complx(data[0],data[1]);
	 else temp = assign_complx(data[2],data[3]);
	 //printf("here\n");
	 jft->c[j] = temp;
		
	 } 
return;		 
}

complex <double> sum_up(struct profile_multi *h){
	/*This function caculate the sum up of a profile.*/
	
	complex <double> result(0.0,0.0);
	
	for (int i = 0; i < h->ntot; i++){
		result = result + h->c[i];
	}
	
	result = result/ BETA / ((double)N);
	
	printf("Real : %f\n",real(result));
	
    return result;
}
//===============================================================
complex <double> lambda(struct profile_multi *j, struct profile_multi *polar, double lambda){
	/*This function multiply polarization by a factor lambda, and returns the sum of Chi, which is a complex number.*/
	
	complex <double> result(0.0,0.0);
	
	for (int i = 0; i < polar->ntot; i++){/*ith Fourier component of bold green's function by Dyson*/
		result =  result + lambda * polar->c[i] / (4.0 - lambda *j->c[i] * polar->c[i]);
	}
	
	result = result/ BETA / ((double)N);
	
	//printf("Real : %f\n",real(result));
	
    return result;
}
//===============================================================
double rescale(struct profile_multi *j, struct profile_multi *polar){
	/*This function find the factor which after rescaling polarization, satisfies the sum-rule.*/ 
	
	double lamb = 0;
	//complex <double> previous(0.0,0.0);
	//complex <double> now(0.0,0.0);
	
	double previous = 0;
	double now = 0;
	
	double precise = 0;
	
	previous = real(lambda(j,polar,1.0));
	
	for (int i = 2; i< N_rescale; i++){//looking for lambda roughly
		lamb = 0.5*(double)i;
		now = real(lambda(j,polar,lamb));
		//printf("lambda: %g Real-sum: %g\n",lamb,real(lambda(j,polar,lamb)));
		//fprintf(f_scale,"%g %g\n",lamb,real(lambda(&prf_jj,&prf_pp,lamb)));
		if (previous < 0.25 && now > 0.25){//the solution between these two numbers i-1 and i
			precise = 0.5*(double)(i-1);
			break;
		}
		previous = now;
	}
	//printf("I am here\n");
	
	printf("i:%f\n",precise);
	if (precise == 0) {
		lamb = 1.0;
		return 1;}
	
	for (int i = 0; i<200; i++){//looking for lambda more precisely
		lamb = precise + 0.005*((double)i);
		now = real(lambda(j,polar,lamb));
		
		if (fabs(now - 0.25)< 0.002) {
			printf("Lambda : %f\n", lamb);
			break;
		}
		//printf("lambda: %g Real-sum: %g\n",lamb,real(lambda(j,polar,lamb)));
		//fprintf(f_scale,"%g %g\n",lamb,real(lambda(&prf_jj,&prf_pp,lamb)));
	}
	
	return lamb;
}
//===============================================================

double time_bin(struct profile_multi *h, int bin){
	/*This funtion calculates the time coordinate of the corresponding bin */
	double rtn = 0;
	int n_time = bin % N_t;
	rtn = h->min[1]+(h->max[1]-h->min[1])/N_t*(bin+0.5); 
	
	return rtn;
}

//====================================================
void scale_profile_multi(struct profile_multi *h, double d){
	for (int i = 0; i < h->ntot; i++){
		//double n = h->n[i];
		h->c[i] = h->c[i] * d;
		h->c2[i] = h->c2[i] * d*d;
	}
	return;
}
//====================================================
void normalize_profile(struct profile_multi *h, int har_or_g, double d){
	/*Normalize profile h by normal har_or_g, number of hartree or gg diagrams.*/
	double n1 = 0.0;
	double n2 = double(har_or_g);
	
	for (int i = 0; i < h->ntot; i++){
		n1 = double(h->n[i]);
		
		//h->c[i] = h->c[i] * d / n2;
		//h->c2[i] = h->c2[i] * d * d /n2/n2;
		
		
		h->c[i] = h->c[i] * d * n1/ n2;
		h->c2[i] = h->c2[i] * d * d *n1 *n1/n2/n2;
	}
	return;
}
//====================================================
//====================================================
void new_normalize_profile(struct profile_multi *sp, struct profile_multi *hg, double d){
	/*Normalize profile sp(self-energy or polarization) by profile of (har_or_g),the profile only contains 
	times it appears in correponding bin.*/
	double n1 = 0.0;
	
	int bin = 0;
	/*for (int i = 0; i < hg->ntot; i++){
		n1 = n1 + double (hg->n[bin]);
	}
	printf("n1:%f\n",n1);
	n1 = n1/(hg->ntot);*/
	/*do{
		bin = rand() % hg->ntot; 
	} while (bin == 0);
	
	n1 = hg->n[bin];*/

	for (int i = 0; i < sp->ntot; i++){
		//read out the number of appearance in each bin
		if (i < N_t) n1=hg->n[i];
		else n1=hg->n[i-N_t];
		
		//sp->c[i] = sp->c[i] * d / n1;
		//sp->c2[i] = sp->c2[i] * d * d /n1/n1;
		
		if (sp->n[i] != 0){
			sp->c[i] = sp->c[i] * d / n1;
			sp->c2[i] = sp->c2[i] * d * d /n1/n1;
		}
		else{
			//printf("empty bin\n");
		}
		
		printf("i: %d nn:%f nsp: %d\n", i, n1, sp->n[i]);
			
		//h->c[i] = h->c[i] * d * n2/ n1;
		//h->c2[i] = h->c2[i] * d * d *n2 *n2/n1/n1;
	}
	return;
}
//====================================================
double Qn2(struct profile_multi *green){
	/*This function caculate the sum up of a profile.*/
	
	//complex <double> result(0.0,0.0);
	double result = 0;
	double tnow = 0;
	double tn = 0;
	double dt = (BETA - Tmin)/N_t;
	
	for (int i = 0; i < N_t; i++){
		tnow = Tmin + (i+1/2)*dt;
		tn = BETA - tnow;
		result = result + abs(-green_f(green,tnow)*green_f(green,tn))*dt;
		//printf("Qn2 i : %d\n",i);
	}

	printf("Real of Qn2 : %f\n",result);
	
    return result;
}
//====================================================
void complement_profile_multi(struct profile_multi *sp ){
	/*This function complement the profile (self-energy of polarization), if there is one or 
	more than one bin without data, calculate the data by two neighboring bins.*/
	
	/*Here if n[i]=0, then there is no data wrote on that bin, we only need to deal with this kind of bin.*/
	int n_sites = N;
	int bin = 0;//identify the bin number in the data file 
	
	int b1 = 0;
	int b2 = 0;

	complex <double> slope(0.0,0.0);
	double dt = (BETA- Tmin)/N_t;
	
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N_t; i++){
			bin = i + j * N_t;
			b1 = b2 = bin;
			if (sp->n[bin] == 0){//only manipulate the bin without data
				printf("empty_bin :%d\n",bin);
				if (i < N_t/2) {
					do {
						b1++;
					}while (sp->n[b1] == 0);//find the next two bins have data
					b2 = b1+1;
					do {
						b2++;
					}while (sp->n[b2] == 0);	
				}
				else{
					do {
						b2--;
					}while (sp->n[b2] == 0);
					b1 = b2 -1 ;
					do {
						b1--;
					}while (sp->n[b1] == 0);	
					
				}
				if (b2 < b1) errorquit("b2<b1");
				/*else if (i == N_t -1){
					do {
						b2--;
					}while (sp->n[b2] == 0);
					b1 = b2 -1 ;
					do {
						b1++;
					}while (sp->n[b1] == 0);	
					
				} 
				else{
					do {
						b1--;
					}while (sp->n[b1] == 0);
					do {
						b2++;
					}while (sp->n[b1] == 0);	
				}*/
				slope = (sp->c[b2]-sp->c[b1])/(dt*(b2-b1));
				sp->c[bin] = sp->c[b1]+slope * dt *double(bin-b1);
			}	
		}	
	}
	
	return;
}


/*Simple expansion, Improper self-energy and polarization operator*/
struct profile bg(struct profile *self){
	/*This function calculate the bold green's function of improper self energy.*/
	struct profile rtn = prof_create(N_t, -BETA, BETA);
	complex <double> g(0.0, 1.0);
	double tiner = 0;

	for (int i = 0; i < self->nbin; i++){
		tiner = self->min + i*(self->max - self->min) / self->nbin;//boundary of each bin in time axis
		g = green(tiner);//bare green's function

		//calculate bold green's function
		rtn.c[i] = g*(1.0 + self->c[i]);
		rtn.c2[i] = (rtn.c[i])*(rtn.c2[i]);

	}
	return rtn;
}
//===============================================================
struct profile bg(struct profile *g, struct profile *self){
	/*This function calculate the bold green's function of improper self energy.*/
	struct profile rtn = prof_create(N_t, -BETA, BETA);

	for (int i = 0; i < self->nbin; i++){
		//calculate bold green's function
		rtn.c[i] = g->c[i] * (1.0 + self->c[i]);
		rtn.c2[i] = (rtn.c[i])*(rtn.c2[i]);
	}
	return rtn;
}
//===============================================================
struct profile bw(struct profile *polar){
	/*This function calculate the effective interaction for reducible diagrams.*/
	struct profile rtn = prof_create(N_t, -TEM, TEM);

	double j = J / 4;
	for (int i = 0; i < polar->nbin; i++){
		//calculate effective interaction
		rtn.c[i] = j*(1.0 + polar->c[i]);
		rtn.c2[i] = (rtn.c[i])*(rtn.c2[i]);
	}

	return rtn;
}
//===============================================================
struct profile wtilde(struct profile *inter){
	/*This function calculate the effective interaction for reducible diagrams.*/
	struct profile rtn = prof_create(N_t, -TEM, TEM);

	double j = J / 4;
	for (int i = 0; i < inter->nbin; i++){
		//calculate effective interaction
		rtn.c[i] = inter->c[i] - j;
		rtn.c2[i] = (rtn.c[i])*(rtn.c2[i]);
	}

	return rtn;
}
//===============================================================

//Functions to do integration 
//===============================================================
double tiner(int n){
	/*Read out the time according to the bin number n, the time is the middle of the bin interval*/
	double rtn = 0;
	rtn = Tmin + (n+0.5)*(BETA-Tmin)/N_t;
	
	return rtn;
	
};
//===============================================================
int reverse_tiner(double t){
	/* Read out the bin number according to the time, and the time is the middle of the bin interval*/
	int k = 0;
	k = int(floor(N_t*(t -Tmin) / (BETA-Tmin)));//find the kth bin, which this w falls into
	return k;
};
//===============================================================
struct profile ite_g(struct profile *oldg, struct profile *self){
	/*This function calculate the new g profile according to the iteration intergral form of propagator.*/
	struct profile rtn = prof_create(N_t, -TEM, TEM);
	int m = 0;
	int n = 0;
	int k = 0;
	double delta = 2*BETA/N_t;
	
	double t = 0;
	double tau = 0;
	double t1 = 0;
	double t2 = 0;
	
	complex <double> zero(0.0, 0.0);
	complex <double> cinter1 = zero;//initialize the complex number
	complex <double> cinter = zero;//initialize the complex number
	
	for(int i = 0; i < N_t; i++){
		
		t = tiner(i);//calculate t
		
		for(int j = 0; j < N_t; j++){
			tau = tiner(j);//calcuate tau
			//Determine the two limits
			if(tau > 0){//if tau>0, 
				m = 0;
				n = reverse_tiner(BETA - tau);
				
			}
			else{//if tau<0
				m = reverse_tiner(-BETA - tau);
				n = N_t - 1;
			}
			
			//printf("=======\n");
			for(int l = m; l < n + 1; l++){
				t1 = tiner (l);
				t2 = t - t1 - tau;
				
				if (t2 > BETA) t2 = t2 - 2 * BETA;
				if (t2 < -BETA) t2 = t2 + 2 * BETA;
				
				k = reverse_tiner(t2);
				
				cinter1 = cinter1 + oldg->c[l]*oldg->c[k]*delta;
				
			}//end of loop l
			
			cinter = cinter + self->c[j]*cinter1*delta;
			cinter1 = zero;
			
		}//end of loop j
		rtn.c[i] = oldg->c[i] + cinter;
		cinter = zero;
	}//end of loop i
	
	return rtn;
	
}
//=========================================================

//make complete profile from (0,\beta) to (-\beta,\beta)
//==============================================================
struct profile_multi minusbeta(struct profile_multi *h){
	/*Construct the diagram	to have 2times number of bins in profile h, fill in the first part with (-beta, 0).*/
	
	double newmin[2] = {XL,-BETA};
	int newnbin[1]= {2*N_t};
	
	struct profile_multi rtn;
	rtn = prof_create_multi(h->dim, newnbin, newmin, h->max);
	
	//printf("Number of bins in h:%d\n",h->ntot);
	
	//Fill profile for tau from -beta to 0
	for (int i = 0; i<h->ntot; i++) rtn.c[i]=-1.0*h->c[i];
	
    //Fill profile for tau from 0 to beta
	for (int i = h->ntot; i<2*h->ntot; i++){
		rtn.c[i]=h->c[i-N_t];
		//printf("i:%d\n",i);
	}
	
	return rtn;
}
//==============================================================
struct profile_multi betaback(struct profile_multi *h){
	/*Construct the diagram	*/
	
	double newmin[2] = {XL,Tmin};
	int newnbin[1]= {N_t};
	
	struct profile_multi rtn;
	rtn = prof_create_multi(h->dim, newnbin, newmin, h->max);
	//printf("Number of bins in h-back:%d\n",h->ntot);
	
    //Fill profile for tau from 0 to beta
	for (int i = 0; i<rtn.ntot; i++){
		rtn.c[i]=h->c[i+N_t];
	}
	
	return rtn;
}
//==============================================================

void print_profile(struct profile_multi *h){
	/*This function prints out all the element in a profile*/
	
	printf("========print a profile=======\n");
	
	for (int i = 0; i < h->ntot; i++){
		printf("i : %d", i);
		printf(" ");
		print_cmplx(h->c[i]);
	}
	return;
}

