//----------------------------------------------------------------------------------------------------
/*
This project , do simulation by use of linked list. save information of lines/ vertices to each node of the two list.
*/
//--------------------------------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE
#include <stdlib.h>   /*rand, srand*/
#include <stdio.h>    /*printf*/
#include <string.h>
#include <math.h>
#include <time.h>     /*time*/
//#include <random>
#include <complex>    /* Standard Library of Complex Numbers */
#include <iostream>
using namespace std;

//#include"linkedlist.h"
#include"newh.h"
#include"linkedlist.cpp"
#include"momenta.c"
#include"create.cc"
#include"test.c"
#include"profile.cpp"
#include"calculation.cpp"
#include"fft.cpp"
#include"updates1.cpp"
#include"1storder.cpp"



#define N_update 10000000 //number of update in update loop
#define N_simul 1       //number of smimulation in simulation loop
#define N_iteration 1   //iteration order
#define PI 3.14159265358979  //Pi

#define N_in      48   //types of updates be included

#define G 1            //indicates self-energy sector
#define P 2            //indicates polarization sector

#define PLUS 1  
#define MINUS -1

int imp_interation(struct diag *t);
void test_2d_fft();
void test_randomt();
//---------------------------------------------------------------------------
int main()
{   

	//srand(time(NULL));//initialize the random number generator, by setting the seed using function time(<time.h>)
    //printf("%f\n",1.0/0.0);
	struct diag tad; //declare a tadpole
	idiag(&tad);
	imp_interation(&tad);
	
	//test_randomt();
	return 0;
}

//=============================================================================
void test_randomt(){
	FILE *f_s;//file pointers
	f_s = fopen("tests.txt", "w+");//create and open a file, save the raw data of  bare green's function in t space
	//f_w = fopen("testw.txt", "w+");//save
	double td = 0.0;
	double t1 = 0.0;
	double t2 = 0.0;
	double ss[1] = {0.0};
	complex <double> self(1.0,0);
	
	//1D profile: t-varible
	int nft[1] = {N_t};
	double min[1] = {Tmin};
	double max[1] = {BETA};
	
	//struct profile_multi prf_g = prof_create_multi(1, nft, min, max);//profile which saves the data of green's function
	struct profile_multi prf_s = prof_create_multi(1, nft, min, max);//profile which saves the data of self energy
	
	for (int i = 1; i<100000; i++){
		t1 = rand_time();
		t2 = rand_time();
		td = t1-t2;
		if (td < 0 ) td = td + BETA;
		ss[0]= td;
		prof_fill_multi(&prf_s, ss, self);
	}
	prof_print_multi(f_s, &prf_s);
	fclose(f_s);
	
	
}
//=============================================================================
int imp_interation(struct diag *t){
	/*This function runs the simulation, only propose to apply first 8 updates.*/ /*No FFT*/
	int update = NOTDEFINE;
	int update_applied = NOTDEFINE; //indientifier to indicate whether this update been applied
	long int n_update = 0;               //counter of number updates applied to the diagram
	int k = 0;                      //lables number of attempted upgrades so far
	//counter for how many update of each type we actually perform
	int c, d, ch, dh, mp, mi, com, dummy, in, rem, dr, udr, rec, mt;
	c = d = ch = dh = mp = mi = com = dummy = in = rem = dr = udr = rec = mt = 0;
	int data_set = 0; //count number of data set(lines) we print out in the data file 
	int irredu = 0;   //

	double normal_factor_s1 = 0; //normalization factor for self-energy part
	double normal_factor_s2 = 0; //normalization factor for polarization operator part

	int iteration = 0; //number of iteration

	int m1, m2, m3, m4, m5, m6; //number of 1st, 2nd , 3rd and 4th order diagram
	m1 = m2 = m3 = m4 = m5 = m6 = 0;

	//1D profile: t-varible
	int nft[1] = {N_t};
	double min[1] = {Tmin};
	double max[1] = {BETA};
	
	struct profile_multi prf_g = prof_create_multi(1, nft, min, max);//profile which saves the data of green's function
	struct profile_multi prf_s = prof_create_multi(1, nft, min, max);//profile which saves the data of self energy
	
	struct profile_multi prf_hartree = prof_create_multi(1, nft, min, max);//profile which saves the data of green's function
	struct profile_multi prf_ggdiag = prof_create_multi(1, nft, min, max);//profile which saves the data of self energy
	
	struct profile_multi prf_gg = prof_create_multi(1, nft, min, max);//data of the fourier transform of g
	struct profile_multi prf_ss = prof_create_multi(1, nft, min, max);//data of the fourier transform of s
	struct profile_multi prf_ge = prof_create_multi(1, nft, min, max);//profile of effective green's func after last iteration in w space
	struct profile_multi prf_get = prof_create_multi(1, nft, min, max);//profile of effective green's func after last iteration in t space
	
	//2D profile: t and x-varibles
	int nft_i[2] = {N,N_t};
	double min_i[2] = {XL,Tmin};
	double max_i[2] = {XU,BETA};
	
	//struct profile_multi prf_i = prof_create_multi(2, nft_i, min_i, max_i);//save the data of bare interaction in 
	struct profile_multi prf_p0 = prof_create_multi(2, nft_i, min_i, max_i);//polarization before iteraction, GG diagram
	struct profile_multi prf_p = prof_create_multi(2, nft_i, min_i, max_i);//save the data of polarization
	struct profile_multi prf_j = prof_create_multi(2, nft_i, min_i, max_i);//save the data of polarization
	
	//struct profile_multi prf_ii = prof_create_multi(2, nft_i, min_i, max_i);//data of the fourier transform of i
	struct profile_multi prf_pp = prof_create_multi(2, nft_i, min_i, max_i);//data of the fourier transform of p
	struct profile_multi prf_jj = prof_create_multi(2, nft_i, min_i, max_i);//data of the fourier transform of j
	
	struct profile_multi prf_we = prof_create_multi(2, nft_i, min_i, max_i);//data of the W-equation
	struct profile_multi prf_wet = prof_create_multi(2, nft_i, min_i, max_i);//data of the W-equation
	
	struct profile_multi prf_wt = prof_create_multi(2, nft_i, min_i, max_i);//data of  ~W in (q,omega) space

	struct profile_multi prf_chi = prof_create_multi(2, nft_i, min_i, max_i);//data of chi in (q,omega) space
	struct profile_multi prf_chichi = prof_create_multi(2, nft_i, min_i, max_i);//data of chi in (r,t) space
	
	
	FILE *f_green, *f_self, *f_polar, *f_selfw, *f_greenw,*f_polarw, *f_greene, *f_greenet, *f_wequation, *f_j, *f_jw, *f_wt, *f_chi,
		*f_wtwt, *f_chichi, *f_f0, *f_p0;//file pointers
		
	f_green = fopen("green.txt", "w+");//create and open a file, save the raw data of  bare green's function in (r,t) space
	f_self = fopen("self.txt", "w+");//save the raw data of self energy in (r,t) space
	f_polar = fopen("polar.txt", "w+");//create and open a file, save the polrization  in (r,t) space
	f_greenw = fopen("green-ft.txt", "w+");//create and open a file, save the green function in (q,w) space
	f_selfw = fopen("self-ft.txt", "w+");//create and open a file, save the self-energy in (q,w) space
	f_polarw = fopen("polar-ft.txt", "w+");//create and open a file, save the polarization in (q,w) space
	f_greene = fopen("green-e.txt", "w+");//create and open a file, save the effective green function in (q,w) space
	f_greenet = fopen("green-et.txt", "w+");//create and open a file, save the effective green function in (r,t) space
	f_wequation = fopen("wequation.txt", "w+");//create and open a file, save the w equation (7) in PRB 87,024407(2013) (q,w)
	f_j = fopen("j.txt", "w+");//create and open a file, save data of of J-profile in (r,t) space
	f_jw = fopen("j-ft.txt", "w+");//create and open a file, save fourier transformation of J-profile in (q,w) space
	f_wt = fopen("wt.txt", "w+");//create and open a file, save ~W in (q,omega) space
	f_chi = fopen("chi-ft.txt", "w+");//create and open a file, save Chi in (q,omega) space
	f_wtwt = fopen("wt-ft.txt", "w+");//create and open a file, save ~W in (r,t) space
	f_chichi = fopen("chi.txt", "w+");//create and open a file, save Chi in (r,t) space
	f_f0 = fopen("f.txt", "w+");//create and open a file, save F coupling in the unphysical sector in (r,t) space
	f_p0 = fopen("p0.txt", "w+");//create and open a file, save Pi_0 polarization in (r,t) space, first order diagram
	

    FILE *f_chi00;
	f_chi00= fopen("chi(0-0).txt", "w+");

	/*FILE *number_order;
	number_order = fopen("order.txt", "w+");//save counts of each diagram_order, from 1 to 4*/
	
	FILE *f_scale;//file pointer for rescale Pi to satisfy the sum rule
	f_scale= fopen("rescale.txt", "w+");
	
    FILE *f_iter;
	f_iter= fopen("iteration.txt", "w+");
	

	double time0 = 0.0;
	complex <double> g(0.0, 0.0);
	complex <double> self(0.0, 0.0);
	complex <double> polar(0.0, 0.0);
	complex <double> zero(0.0, 0.0);

	int nm = 0;

	int n_hartree = 0;
	int n_gg = 0;
	int n_fake_hartree = 0;
	int n_fake_gg = 0;
	
	int n_gsector = 0;
	int n_psector = 0; 
	
	int hartree_p = 0;
	int hartree_n = 0;

    int irr_c, irr_d, irr_ch, irr_dh, irr_com , irr_dum, irr_rec, irr_mt;
	irr_c = irr_d = irr_ch = irr_dh = irr_com = irr_dum = irr_rec = irr_mt = 0;
	
	//array for 4 profiles
	double gg[1] = {0.0};
	double ss[1] = {0.0};
	double ii[2] = {0.0,0.0};
	double pp[2] = {0.0,0.0};
	double jj[2] = {0.0,0.0};
	
	//Construct the profile for J(r)
	complex <double> temp(0.0,0.0);
	double tnow = 0.0;
	
	complex <double> complx_zero(0.0, 0.0);
	
	struct node* old_measure;
	struct node* measure;
	
	int rand_meas  = 0;
	
	
	idiag(t);     //initialize the diagram
    //startpoint(t);//create a simplest diagram to get start
	//create_tadpole_ladder_phy_s1(t);
    firstdiag(t);
	test_diag(t); //check the created diagram 
	assign_parity(t);//store the parity of the diagram
	if (hartree_or_g(t) == OK) printf ("GG\n");
	else printf("No GG\n");
	//printf("I am here\n");
	
	//construct the J-profile (r)
	for (int j = 0; j < nft_i[0]; j++){
	
		jj[0] = (double)j;//the x-component of J-profile (bare interaction)
	
		if (jj[0]== 1.0){// only when the position difference is 1.0, the thing is non-zero
		
			for (int i=0; i < prf_j.nbin[1]; i++){//t - coordinates, delta function in t (when t<0.03)
			    //tnow = i* BETA / N_t;
				tnow =prf_j.min[1]+(prf_j.max[1]-prf_j.min[1])/prf_j.nbin[1]*(i+0.5);//t_i=i*Delta_t
		
				//for green's function
				if (fabs(tnow)<0.03) temp = assign_complx(1.0*J0,0.0); 
				else temp = assign_complx(0.0,0.0);
		
				jj[1] = tnow;//the t-component of J-profile
	
				prof_fill_multi(&prf_j,jj,temp);//write information into profile
			}
		
		}
	
	}
	//interacton J (q,omega)
	prf_jj = prof_four_multi(&prf_j, 1); 
	//correct_j(&prf_jj);
	
	//construct the G-profile for the zero-th iteration, which is just the G0
	for (int i = 0;i< t->prf_g.ntot; i++ ){
		tnow =t->prf_g.min[0]+(t->prf_g.max[0]-t->prf_g.min[0])/t->prf_g.nbin[0]*(i+0.5); 
		//printf("tnow:%g\n",tnow);
		t->prf_g.c[i] = green(tnow);
		//print_cmplx(t->prf_g.c[i]);
	}
	
	//construct the Pi-profile for the zero-th iteration, from GG diagram
	for (int i = 0 ;i< N_t; i++ ){
		prf_p0.c[i] = assign_complx(0.0,2.0/PI/TEM);
	}
	/*prf_pp = prof_four_multi(&prf_p0, 1); 
	//calculate the w equation
	prf_we = wequation(&prf_jj, &prf_pp);//Dyson equation 2, W, effective interaction in (q,omega)space
	prf_wet= prof_four_multi(&prf_we, -1);// Effective W in (r,t) space

	//calculate the ~W value for next iteration
	prf_wt = wtilde(&prf_jj, &prf_we);//in (q,omega)space
	t->prf_w = prof_four_multi(&prf_wt, -1); //Inverse FT, Wtilde, the screening interaction in (r,t) space*/
		
	//construct the ~W-profile for the zero-th iteration, which is just a constant W0
	for (int i = 0;i< (t->prf_w.ntot); i++ ){
		if (i < N_t) (t->prf_w).c[i] = assign_complx(W0,0.0);
		else (t->prf_w).c[i] = assign_complx(W0,0.0);
		//printf("%d\n",i);
	}
	
	//construct the F-profile for the zero-th iteration, which is just a constant F0
	for (int i = 0;i< t->prf_f.ntot; i++ ){
		t->prf_f.c[i] = assign_complx(F0,0.0);
		//printf("i: %d f:%g\n",i,real(t->prf_f.c[i]));
	}
	
	//loop to do interation, the identificator for interation order is integer "interation"
	for (int iteration = 0; iteration < N_iteration; iteration++){
		
		printf("======================\n");
		printf("ITERATION ORDER: %d\n",iteration);
		
		//In each interation, propose random updates by "N_updates" times
		for (unsigned long long j = 0; j < N_update; j++){//update loop in each simulation
			//update = pro_updates(); // randomly generate integer form 0 to 13, according to probability distribution
			update = 1 + rand() % N_in; // randomly generate integer form 0 to 10, each integer stands for a update
			//printf("update %d\n",update);
			switch (update){
			case 1 ... 9:
			    //printf("%d\n", 1);
				update_applied = Create(t);
				if (OK == update_applied) {
					//printf("%d\n", 1);
					c++;
					if (NOTOK == check_irreducibility(t)) irr_c++;
				}
				break;
			case 10 ... 18:
			    //printf("%d\n", 2);
			    update_applied = Delete(t);
				if (OK == update_applied) {
					//printf("%d\n", 2);
					d++;
					if (NOTOK == check_irreducibility(t)) irr_d++;
				}
				break;
			case 19 ... 27:
			    //printf("%d\n", 3);
				update_applied = Create_H(t);
				//printf("%d\n", update_applied);
				if (OK == update_applied) {
					//printf("%d\n", 3);
					ch++;
					change_parity(t);
					if (NOTOK == check_irreducibility(t)) irr_ch++;
				}
				
				break;
			case 28 ... 36:
			    //printf("%d\n", 4);
			    update_applied = Delete_H(t);
				if (OK == update_applied) {
					//printf("%d\n", 4);
					dh++;
					change_parity(t);
					if (NOTOK == check_irreducibility(t)) irr_dh++;
				}
				
				break;
			case 37 ... 41:
			    //printf("%d\n", 5);
				update_applied = Move_P(t);
				if (OK == update_applied) {
					//printf("%d\n", 5);
					mp++;
					//if (NOTOK == check_irreducibility(t)) printf("reduc-p\n");
				}
				break;
			case 42:
			    update_applied = Move_I(t);
				if (OK == update_applied) {
					mi++;
					//printf("%d\n", 6);
				}
				break;
			case 43 ... 44:
			    //printf("%d\n", 7);
				update_applied = Commute(t);
				if (OK == update_applied) {
					//printf("%d\n", 7);
					com++;
					change_parity(t);//change number of fermion loops
					if (NOTOK == check_irreducibility(t)) irr_com++;
				}
				
				break;
			case 45:
				update_applied = Dummy(t);
				if (OK == update_applied) {
					//printf("%d\n", 8);
					dummy++;
					if (NOTOK == check_irreducibility(t)) irr_dum++;
				}
				break;
			case 97:
				//printf("%d\n", 9);
			    update_applied = Insert(t);
				if (OK == update_applied) in++;
				break;
			case 98:
			    //printf("%d\n", 12);
			    update_applied = Remove(t);
				if (OK == update_applied) rem++;
				break;
			case 49:
				//printf("%d\n", 11);
			    update_applied = Dress(t);
				if (OK == update_applied) dr++;
				break;
			case 50:
				//printf("%d\n", 12);
			    update_applied = Undress(t);
				if (OK == update_applied) udr++;
				break;
			case 46 ... 47:
			    //printf("%d\n", 13);
			    update_applied = Recolor(t);
			    if (OK == update_applied) {
					//printf("%d\n", 13);
					rec++;
					if (NOTOK == check_irreducibility(t)) irr_rec++;
				}
				break;
			case 48:
			    //printf("%d\n", 14);
				update_applied = Move_T(t);
				if (OK == update_applied) {
					//printf("%d\n", 14);
					mt++;
					if (NOTOK == check_irreducibility(t)) irr_mt++;
				}
				break;
			}      //  propose a random update
			//if (OK == update_applied) test_detached(t);

			//print_information(t);
			if (t->n_wo == 0){
				test_diag(t); //test diagram after each update
				//diag_print(t);		
			}

			if (OK == update_applied) {//if the update been applied succesfully
				k++;    //increase k by 1,attempted upgrades so far
			}

			//print out a number every 100000 attempts
			n_update++;
			if (t->n_wo == 0){//if the diagram is in physical sector
				if (diag_order(t) == 1) m1++;
				else if (diag_order(t) == 2) m2++;
				else if (diag_order(t) == 3) m3++;
				else if (diag_order(t) == 4) m4++;
				else if (diag_order(t) == 5) m5++;
				else if (diag_order(t) == 0) printf("zero-th order\n");
				else m6++;
			}

			if (n_update % 1000000 == 0) {
				printf("%ld/%ld \n", n_update, (long int)(iteration+1)*N_update);
				//printf("success update:%d ", k);
				printf("Physical: %d %d %d %d %d %d\n", m1, m2, m3, m4, m5, m6);
			}
			
            //MEASURE DIAGRAM;update data profile, which stores information of Sigma, Pi
			if (t->n_wo == 0 && diag_order(t) < (Max_meas+1)){
				time0 = time_diff(t->meas_line);
				//printf("timediff:%f\n",time0);
				double presign = 1.0;
				if (time0 < 0) {
					time0 = BETA + time0; // t=t+Beta
					presign = -1.0;   // flip sign of G and Sigma
				}
				
				if (sector(t) == G){//s=1, G-sector
					//Measure the self-energy 
					n_gsector++;
					//printf
					//print_cmplx(self_energy(t));
				
					if (diag_order(t) == 1){
					    //number of Hartree diagrams
					    if (hartree_or_g(t) == OK) {
							n_hartree++;	
							if (fake_hartree(t) == OK) printf("fake hartree");//test of hartree diagram
							hartree(t);
							//ss[0] = time0;
						    ss[0] = rand_time();
					
							prof_fill_multi(&prf_hartree, ss, zero);
					
							self = presign * self_energy(t);
							//self = self_energy(t)
							//print_cmplx(self);
		
							//save self energy to profile
							//prof_fill_multi(&prf_s, ss, self);
							//printf("loops:%d\n",diag_order(t));
							
							if (coupling(t->meas_line->data.l.head->data.v.pend)>0) hartree_p++;
							else hartree_n++;
						}
						else {
							ss[0] = time0;
							self = presign * self_energy(t);
							prof_fill_multi(&prf_s, ss, self);
						}
					}
					else{
						ss[0] = time0;
						self = presign * self_energy(t);
						prof_fill_multi(&prf_s, ss, self);
				
					}
				}//measuring G-sector
				else{//s=2, P-sector
					
					//MEASURE Polarization operator
					n_psector++;
					//printf("P sector\n");
					if((t->meas_line)->data.l.phys.screen == J) printf("GG with J\n");
					struct node* p; //find the measuring pline
					p = t->meas_line;
			
					pp[1] = time0;
					polar = contribution(t);
					
					//print_cmplx(polar);
				
					if (position_diff(p) == 0) {
			
						//number of GG diagrams
						if (diag_order(t) == 1){
							if (hartree_or_g(t) == OK /*&& time0 < BETA/N_t*/) {
								//if (pi_sign(t) == -1) printf("wrong sign\n");
								ggdiag(t);
								//if (fabs(time0) < (BETA-Tmin)/N_t) n_gg++;//account for gg_diagram within a time interval
								n_gg++;
							    ss[0] = time0;
								prof_fill_multi(&prf_ggdiag, ss, zero);
								
								pp[0] = 0.0;
								prof_fill_multi(&(prf_p), pp, polar);
							}
						
						}
						else{
							if (fabs(p->data.l.phys.p.x)+fabs(p->data.l.phys.p.y)+fabs(p->data.l.phys.p.z) < 1E-4) printf("00\n");
							else prof_fill_multi(&(prf_p), pp, polar);
						}
			
					}
					if (abs(position_diff(p)) == 1) {
						//printf("1\n");
						pp[0] = 1.0;
						if (fabs(p->data.l.phys.p.x)+fabs(p->data.l.phys.p.y)+fabs(p->data.l.phys.p.z) < 1E-4);//printf("00\n");
						else prof_fill_multi(&(prf_p), pp, contribution(t));
						/*if (count(t) == 1) {
							if ((p->data.l.head->data.v.ein)==(p->data.l.head->data.v.eout)) {
								printf("diagram order is not one\n");
								diag_print(t);
							}
						}
						else prof_fill_multi(&prf_p, pp, polar);*/
					}
				
				}//measuring P-sector
				
		

				data_set++;
				
				
			}
			
			//initialize identifier, prepare for next loop
			update_applied = NOTDEFINE; //re-initialize it in every loop
		}
		/*printf("POLARIZATION\n");
		print_profile(&prf_p);
		printf("SELF-ENERGY\n");
		print_profile(&prf_s);*/
       
		//calculate normalization factor after first set of N_updates updates
		if (iteration == 0) {
			normal_factor_s1 = N * 2 * fabs(J0/4) * abs(green(BETA));//the normalization factor for self-energy
			//printf("Qn(1),%f\n",normal_factor_s1);
			normal_factor_s2 = 2*BETA;
			//printf("Qn(2),%f\n",normal_factor_s2);
		}
		else {//not first iteration
			normal_factor_s1 =  N * 2 * fabs(J0/4) * abs(green_f(&(t->prf_g), 0.0));
			//printf("normal_factor_s1 %f\n",normal_factor_s1);
			print_profile(&(t->prf_g));
			normal_factor_s2 = 2*Qn2(&(t->prf_g));
			//printf("normal_factor_s1 %f\n",normal_factor_s2);
		}
		
		
		printf("n_hartree:%d\n",n_hartree);
		printf("n_gg:%d\n",n_gg);
		//Update profiles only when n_hartree != 0 and n_gg != 0
		if (n_hartree != 0 && n_gg != 0){
			//normalize the self_energy in space-time domain
			printf("G-sector:%d\n",n_gsector);
			printf("n_hartree:%d\n",n_hartree);
			normalize_profile(&prf_s, n_hartree, normal_factor_s1);
			//new_normalize_profile(&prf_s, &prf_hartree, normal_factor_s1);
		
			//normalize the polarization operator
			printf("P-sector:%d\n",n_psector);
			printf("n_gg:%d\n",n_gg);
			printf("======Normalize Polarization======\n");
			//if(iteration !=0){
			normalize_profile(&prf_p, n_gg, normal_factor_s2);	
			//new_normalize_profile(&prf_p, &prf_ggdiag, normal_factor_s2);
			
			printf("P(r=0,t=0) order:%f\n",real(prf_p.c[0]));	
			
			//complement_profile_multi(&prf_s);
			//complement_profile_multi(&prf_p);	
		
			//prof_print_multi(f_f0, &t->prf_f);
			
			
			/*printf("POLARIZATION\n");
			print_profile(&prf_p);
			printf("SELF-ENERGY\n");
			print_profile(&prf_s);*/
			
		
			//Fourier transformation of (FFT) green function and self-energy
			prf_g = t->prf_g;
			
			struct profile_multi prf_e = minusbeta(&prf_g);//extand G to tau<0
			struct profile_multi prf_snew = minusbeta(&prf_s);//extand \Sigma to tau<0
			
			/*printf("E\n");
			print_profile(&prf_e);
			printf("SE\n");
			print_profile(&prf_snew);*/
			
			prf_e = prof_four_multi(&(prf_e),1);
			prf_snew = prof_four_multi(&(prf_snew),1);
			
	        //prf_gg = prof_four_multi(&(t->prf_g),1);//bare propagator
			//prf_ss = prof_four_multi(&prf_s, 1); //self energy
			
			printf("======FFT of P======\n");
			prf_pp = prof_four_multi(&prf_p, 1); //ploarization
			printf("PP(q=0,omega=0) order:%g\n",real(prf_pp.c[0]));

			//calculate bold green's function by Dyson equation
		    //prf_ge = bold_green(&prf_gg, &prf_ss);  //Dyson equation, bold green's function in (q,omega)space
		    struct profile_multi prf_gege = bold_green(&prf_e, &prf_snew);
			//printf("GE-FT\n");
			//print_profile(&prf_gege);
			
			//update the bold propagator profile
			//t->prf_g = prof_four_multi(&prf_ge, -1);//back to (r,t)space
			
			struct profile_multi prf_over = prof_four_multi(&prf_gege, -1);//back to (r,t)space
			//printf("GE\n");
			//print_profile(&prf_gege);
			
			printf("bins in GeGe in tau:%d\n", prf_gege.ntot);
			t->prf_g = betaback(&prf_over);	//set t back to interval (0,BETA) again
			//printf("GE\n");
			//print_profile(&t->prf_g);
			printf("bins in bold G in tau:%d\n", t->prf_g.ntot);
			
			printf("======\n");
		
			//Rescale polarization to satisfy the sum-rule
			double lamb = rescale(&prf_jj,&prf_pp);
			//double lamb = 2.0-1*iteration;
		    //scale_profile_multi(&prf_pp,lamb);
			//prf_p = prof_four_multi(&prf_pp, -1); //ploarization
		
			//Susceptibility
			prf_chi = chi(&prf_jj, &prf_pp);
			prf_chichi = prof_four_multi(&prf_chi, -1);//Inverse FT
		
			//calculate the w equation
			prf_we = wequation(&prf_jj, &prf_pp);//Dyson equation 2, W, effective interaction in (q,omega)space
			prf_wet= prof_four_multi(&prf_we, -1);// Effective W in (r,t) space
		
			//calculate the ~W value for next iteration
			prf_wt = wtilde(&prf_jj, &prf_we);//in (q,omega)space
			printf("here\n");
			t->prf_w = prof_four_multi(&prf_wt, -1); //Inverse FT, Wtilde, the screening interaction in (r,t) space
		    
			double gx[1] = {0.0};
			
			double xx[2] = {0.0,0.0};
			double yy[2] = {1.0,0.0};
			
			printf("Pi(tau) order:%g\n",real(read_binx(&prf_p, xx)));
			printf("Pi(omega) order:%g\n",real(read_binx(&prf_pp, xx)));
			printf("~W (tau):%g\n",real(read_binx(&t->prf_w, xx)));
			printf("~W (omega):%g\n",real(read_binx(&prf_wt, xx)));
			printf("~W (r=1,t = 0):%g\n",real((t->prf_w).c[32]));
			printf("~W (r=0,t = 0):%g\n",real((t->prf_w).c[0]));
			printf("Chi(t=0,r=1) order:%g\n",real((prf_chichi).c[N_t]));
			printf("G(t=0) order:%g\n",real((t->prf_g).c[0]));
			
		
		
			fprintf(f_chi00,"%d %f\n",iteration,real(read_binx(&prf_chi, xx)));
			fprintf(f_iter,"%d %d %d %.4f %.4f %.4f %.4f %.4f %.4f %.6f %.4f\n",iteration,n_hartree,n_gg,lamb,real(read_binx(&prf_chi, xx)),real(read_binx(&prf_chi, yy)),
			real(read_binx(&prf_chi, xx))+real(read_binx(&prf_chi, yy)),real((t->prf_g).c[0]),
			real((prf_chichi).c[N_t]),real((t->prf_w).c[N_t]),real((prf_wet).c[N_t]));
		
			//Update F for the next iteration
			//newf(&(prf_wet), &(t->prf_f));//update F
		
		}
		//printf("~~~~~ after normalization\n");
		//printf("real p:%f\n",real(prf_p.c[20]));
		//clear self-energy and polarization profiles for the next iteration
		if (iteration < (N_iteration - 1)){
			init_profile_multi(&prf_s);
			printf("initialize p\n");
			init_profile_multi(&prf_p);
			init_profile_multi(&prf_hartree);
			init_profile_multi(&prf_ggdiag);
			n_hartree = 0;
			n_gg = 0;
			n_gsector = 0;
			n_psector = 0;
		}
		/*for (int i=0; i<prf_p.ntot; i++){
			printf("i:%d\n",i);
			printf("real p:%f\n",real(prf_p.c[i]));
		}*/
		printf("applied updates in each iteration:%d\n", k);
	}
	
	//printf("%f\n",real(t->prf_w.c[10]));
	prof_print_multi(f_wt, &(t->prf_w));
	
	//print out profile information
	prof_print_multi(f_green, &prf_g);
	prof_print_multi(f_self, &prf_s);
	prof_print_multi(f_polar, &prf_p);
	prof_print_multi(f_j, &prf_j);

	prof_print_multi(f_greenw, &prf_gg);
	prof_print_multi(f_selfw, &prf_ss);
	prof_print_multi(f_polarw, &prf_pp);
	prof_print_multi(f_jw, &prf_jj);

	prof_print_multi(f_greene, &prf_ge);
	prof_print_multi(f_greenet, &(t->prf_g));
	prof_print_multi(f_wequation, &prf_we);

	
	prof_print_multi(f_chi, &prf_chi);

	prof_print_multi(f_wtwt, &prf_wt);
	prof_print_multi(f_chichi, &prf_chichi);

	prof_print_multi(f_f0, &(t->prf_f));
	prof_print_multi(f_p0, &prf_p0);
	
	

	//close file,save data
	fclose(f_green);
	fclose(f_self);
	fclose(f_polar);
	fclose(f_polarw);
	fclose(f_selfw);
	fclose(f_greenw);
	fclose(f_greene);
	fclose(f_greenet);
	fclose(f_wequation);
	fclose(f_j);
	fclose(f_jw);
	fclose(f_wt);
	fclose(f_chi);
	fclose(f_wtwt);
	fclose(f_chichi);
	
	fclose(f_scale);
	fclose(f_f0);
	fclose(f_p0);
	fclose(f_chi00);
	fclose(f_iter);
	
	//prof_print_multi(f_wt, &(t->prf_w));
	

	printf("applied updates:%d\n", k);
	printf("create:  %d\n", c);
	printf("delete:  %d\n", d);
	printf("create_H:%d\n", ch);
	printf("delete_H:%d\n", dh);
	printf("move-p:  %d\n", mp);
	printf("move-i:  %d\n", mi);
	printf("dummy:   %d\n", dummy);
	printf("commute: %d\n", com);
	printf("insert:  %d\nREMOVE: %d\nDRESS: %d\n"
		"UNDRESS: %d \nRECOLOR: %d \nMOVE-T: %d\n", in, rem, dr, udr, rec, mt);

	printf("Order of diagram:%d\n", diag_order(t));
	printf("first order diagram:%d\n", m1);

	printf("number of Hartree diagram with J: %d\n", n_hartree);
	printf("number of GG diagram with ~W:%d\n", n_gg);
	printf("total number of diagram :%d\n", data_set);
	
	/*diag_print(t);
	test_diag(t);
	printf("delte:%d\n",Delete(t));
	printf("delte-h:%d\n",Delete_H(t));
	printf("dummy:%d\n",Dummy(t));
	//printf("delte:%d\n",Delete(t));
	printf("Irreducible:%d\n",check_irreducibility(t));*/
	
	/*printf("========CHECK IRREDUCIBILITY=====\n");
	printf("create:  %d\n", irr_c);
	printf("delete:  %d\n", irr_d);
	printf("create_H:%d\n", irr_ch);
	printf("delete_H:%d\n", irr_dh);
	printf("dummy:   %d\n", irr_dum);
	printf("commute: %d\n", irr_com);
	printf("RECOLOR: %d\nMOVE-T: %d\n", irr_rec, irr_mt);*/
	

	return 1;
}
//=============================================================================
void test_2d_fft(){
	printf("main N_t:%d\n",N_t);
	FILE *f_g, *f_i, *f_ftg, *f_fti,*f_gre;//file pointers
	f_g = fopen("testg.txt", "w+");//create and open a file, save the raw data of  bare green's function in t space
	f_i = fopen("testi.txt", "w+");//save
	f_ftg= fopen("gg.txt", "w+");
	f_fti= fopen("ii.txt", "w+");
	f_gre= fopen("ggre.txt", "w+");
	
	
	
	int nft[1] = {N_t};
	double min[1] = {Tmin};
	double max[1] = {BETA};

	
	struct profile_multi prf_g = prof_create_multi(1, nft, min, max);//profile which saves the data of green's function
	
	int nft_i[2] = {2,N_t};
	double min_i[2] = {-0.5,Tmin};
	double max_i[2] = {1.5,BETA};
	
	struct profile_multi prf_i = prof_create_multi(2, nft_i, min_i, max_i);//save the data of bare interaction
	
	struct profile_multi prf_gg = prof_create_multi(1, nft, min, max);//data of the fourier transform of g
	struct profile_multi prf_it = prof_create_multi(2, nft_i, min_i, max_i);//data of the fourier transform of i
	
	
	double t = 0.0;
	double gg[1] = {0.0};
	double ii[2] = {0.0,0.0};
	complex <double> temp(0.0,0.0);
	
	for (int i=0; i < N_t; i++){//time coordinates
	    t = prf_g.min[0]+(prf_g.max[0]-prf_g.min[0])/prf_g.nbin[0]*(i+0.5);//t_i=i*Delta_t
		//printf("%f\n",2.0/N_t);
		//printf("%d %f %f\n",i,t,real(green(t)));
		printf("I:%d tnow:%f \n",i,t);
		//for green's function
		if (fabs(t)<0.03) temp=assign_complx(1.0,0.0); 
		else temp=assign_complx(0.0,0.0);
		
		gg[0] = t;
		ii[1] = t;
		
		prof_fill_multi(&prf_g,gg,green(t));
		
		
		ii[0]=1.0;
		prof_fill_multi(&prf_i,ii,temp);
		
		ii[0]= 0.0;
		temp=assign_complx(0.0,0.0);
		prof_fill_multi(&prf_i,ii,temp);
	}
	

	
	prof_print_multi(f_g, &prf_g);
	prof_print_multi(f_i, &prf_i);
	
	struct profile_multi four_g = prof_four_multi(&prf_g, 1); //bare propagator
	struct profile_multi four_i = prof_four_multi(&prf_i, 1); //bare propagator
	
	prof_print_multi(f_ftg, &four_g);
	prof_print_multi(f_fti, &four_i);
	
	struct profile_multi four_gre = prof_four_multi(&four_g, -1); //bare propagator
	prof_print_multi(f_gre, &four_gre);
	
	fclose(f_g);
	fclose(f_i);
	fclose(f_ftg);
	fclose(f_fti);
	fclose(f_gre);

	
	return;
}

/*The integration for lowest-order GG diagram, to do normalization, in Eq.(15)*/
double gg_value(struct profile *h){
	/*This function calculate the integration in the Eq.(15)*/
	double rtn = 0;
	int number = h->nbin / 2; //a integer number which is half of the number of bins in profile
	complex <double> g1(1.0, 0.0);
	complex <double> g2(1.0, 0.0);
	//printf("number:%d\n", number);
	for (int i = 0; i < number; i++){
		g1 = h->c[number + i + 1];
		g2 = h->c[number - i];

		rtn = rtn + abs(g1*g2);
		//printf()
		printf("g1*g2 %f\n", abs(g1*g2));
	}

	return rtn;
}

/*Functions to do Normalization*/
//====================================================
void times_profile(struct profile *h, double d){
	for (int i = 0; i < h->nbin; i++){
		double n = h->n[i];
		h->c[i] = h->c[i] * d*n;
		h->c2[i] = h->c2[i] * d*d*n*n;
	}
	return;
}
//====================================================




//====================================================
void diag_print(struct diag *t) {
	char measure = 'i';
	if (t->meas_line->data.l.type == ELECTRON) measure = 'E';
	else measure = 'P';
	printf("===============================\n");
	printf("Number of electron lines: %d\n", t->n_el);
	printf("Number of phonon lines: %d\n", t->n_pl);
	printf("Number of measuring lines: %d, type: %c\n", t->n_ml, measure);
	printf("Number of vertices: %d\n", t->n_ver);
	printf("Number of worms: %d\n", t->n_wo);
	printf("ID of measuring line: %d\n", t->meas_line->data.id); //ID of measuring line
	
	

	if (t->s_worm != NULL && t->t_worm != NULL){//print the information of S, T worms
		printf("ID of S worm: %d\n", t->s_worm->data.id);
		printf("ID of T worm: %d\n", t->t_worm->data.id);
	}
	else{
		printf("No worms, we are in a physical sector.\n");
	}
	//Start by printing the vertex
	reset_current(&t->vertex);
	if (t->vertex.current == NULL) {
		printf("Duh! Novertex in the diagram!\n");
	}
	else
		do{
			char ver_type = 'i';
			if (t->vertex.current->data.v.type == NO) ver_type = 'N';
			if (t->vertex.current->data.v.type == S) ver_type = 'S';
			if (t->vertex.current->data.v.type == T) ver_type = 'T';
			printf("Vertex %d In:%d Out:%d Inter:%d Position:%d Time:%f type:%c\n",
				t->vertex.current->data.id,
				t->vertex.current->data.v.ein->data.id,
				t->vertex.current->data.v.eout->data.id,
				t->vertex.current->data.v.pend->data.id,
				t->vertex.current->data.v.phys.r,
				t->vertex.current->data.v.phys.t,
				ver_type);
		} while (next_one(&t->vertex) == OK);

		//Then print the lines
		reset_current(&t->line);
		if (t->line.current == NULL) {
			printf("Duh! No lines in the diagram!\n");
		}
		else
			do{
				printf("Line %d Head:%d Tail:%d",
					t->line.current->data.id,
					t->line.current->data.l.head->data.id,
					t->line.current->data.l.tail->data.id);
				if (t->line.current->data.l.type == ELECTRON){
					if (t->line.current->data.l.phys.s == UP) printf(" spin:%s", "UP");
					else printf(" spin:%s", "DOWN");
				}
				else{
					if (t->line.current->data.l.phys.screen == J) printf(" screen:%s", "J");
					if (t->line.current->data.l.phys.screen == W) printf(" screen:%s", "W");
					//if (t->line.current->data.l.phys.screen == F) printf(" screen:%s", "F");
				}
				printf(" ");
				print_momentum(t->line.current->data.l.phys.p);
			} while (next_one(&t->line) == OK);

			printf("===============================\n");

			return;
}
//==========================================================

//===========================================================
struct node *rand_eline(struct list *l){
	/* This function choose one propagator from the line list randomly, and returns the pointer of this chosen line.*/
	int order = 0;        //order of which node in the list will be chosen

	do{
		order = rand() % (l->n_nodes) + 1;  //choose a node in the list at random, random number between 1 and node_number
		nth_node(l, order);                 //set the current pointer to the chosen node 
	} while (l->current->data.l.type != ELECTRON); //as long as the chosen line is not electron, we choose again.

	return l->current;  // returns the pointer of chosen line
}
//===========================================================
struct node *rand_pline(struct list *l){
	/* This function choose one interaction line from the line list randomly, and returns the pointer of this chosen line.*/
	int order = 0;        //order of which node in the list will be chosen

	do{
		order = rand() % (l->n_nodes) + 1;  //choose a node in the list at random, random number between 1 and node_number
		nth_node(l, order);                 //set the current pointer to the chosen node 
	} while (l->current->data.l.type != PHONON); //as long as the chosen line is not phonon, we choose again.

	return l->current;  // returns the pointer of chosen line
}
//===========================================================
struct node *rand_pline_meas(struct list *l){
	/* This function choose one interaction line from the line list randomly, and returns the pointer of this chosen line.*/
	int order = 0;        //order of which node in the list will be chosen

	do{
		order = rand() % (l->n_nodes) + 1;  //choose a node in the list at random, random number between 1 and node_number
		nth_node(l, order);                 //set the current pointer to the chosen node 
	} while ((l->current->data.l.type == ELECTRON) || (l->current->data.l.type == PHONON &&l->current->data.l.phys.screen == J)); //as long as the chosen line is not phonon, we choose again.

	return l->current;  // returns the pointer of chosen line
}
//===========================================================
struct node *rand_ver(struct list *l){
	/* This function choose one vertex from the vertex list randomly, and returns the pointer of this chosen line.*/
	int order = 0;    //order of which node in the list will be chosen

	order = rand() % (l->n_nodes) + 1; //choose a node in the list at random, random number between 1 and node_number
	nth_node(l, order);                //set the current pointer to the chosen node

	return l->current;  // returns the pointer of chosen vertex
}
//===========================================================
struct node *rand_meas_line(struct list *l,int s){
	/* This function choose one line (either eline/pline) from the line list randomly, and make the new chosen line to be measuring line
	then returns the pointer of this chosen line.*/
	
	struct node *rtn;//the pointer to be returned
	
	int order = 0;        //order of which node in the list will be chosen

	if(s == 1) rtn = rand_eline(l);//if s=1 chose a eline
	if(s == 2) {//if s=-1, chose a ~W pline
		do{
			order = rand() % (l->n_nodes) + 1;  //choose a node in the list at random, random number between 1 and node_number
			nth_node(l, order);                 //set the current pointer to the chosen node 
		} while ((l->current->data.l.type == ELECTRON) || (l->current->data.l.type == PHONON &&l->current->data.l.phys.screen == J)); //as long as the chosen line is not phonon, we choose again.
		rtn = l->current;
	}
	
	rtn->data.l.measure = MEASURE;//make the new line to be a meauring line

	return rtn;  // returns the pointer of chosen line
}
//===========================================================
int rand_worm_tobe(struct node *ver){
	/*this function assign the special type(S or T) to a vertex at random. */
	int ran = rand() % 2; //generate random number 0 or 1
	if (0 == ran) ver->data.v.type = S;  //case 0, vertex is a S worm 
	else ver->data.v.type = T;           //case 1, vertex is a T worm

	return OK;
}
//===========================================================
struct node *select_worm(struct diag *t){
	/*This function radomly choose a worm to be delete/modify, returns a pointer to this worm.*/
	int ran = rand() % 2;
	if (ran == 0) return t->s_worm;   //case 0, returns the pointer of S worm
	else return t->t_worm;            //case 1, returns the pointer of T worm
}
//===========================================================
struct node *pair_ver(struct node *ver){
	/*This function read out the pair vertex of the argument, which is connect to the argument vertex by a interaction line.*/
	struct node *p = ver->data.v.pend;
	if (ver->data.v.p_in_out == IN) return p->data.l.tail;
	else {
		return p->data.l.head;
	}
}
//===========================================================
int tildeW_present(struct list *l){
	/* This function tells whether there is at least one ~W line in the line list.*/
	int rtn = NOTOK;        //order of which node in the list will be chosen
	reset_current(l); //set pointer to the head of the list
	for (int i = 0; i < l->n_nodes; i++){
		if (l->current->data.l.type == PHONON &&l->current->data.l.phys.screen == W) {
			rtn = OK;
			break;
		}
		else next_one(l);		
	}
	
	return rtn;
}
//===========================================================

//===========================================================
int check_irreducibility(struct diag *t){
	/*An easy way to ensure irreducibility is to check that no two lines (same type, either electron or phonon)
	in the diagram have the same momentum.*/

	struct momentum p_temp, delta;
	enum line_type line_type_temp;
	enum spin spin_temp;
	enum pline screen_temp;
	int worm_end_temp = -1;

	reset_current(&t->line);//move to the first element(node) in the line list
	
	for (int i = 1; i < t->line.n_nodes; i++){//only need to account for the n-1 nodes
		
		/*Read out the line type and momentum on the temporary line */
		line_type_temp = t->line.current->data.l.type;

		if (line_type_temp == ELECTRON) spin_temp = t->line.current->data.l.phys.s; //spin of the temp electron line
		else {
			screen_temp = t->line.current->data.l.phys.screen; //screen of the temp interaction line
			worm_end_temp = tell_end_worm(t->line.current);//tell whether this pline is F type
		}
		
		p_temp = t->line.current->data.l.phys.p;//momentum on the temp line
		
		//compare to the rest of the lines
		for (int j = i; j < t->line.n_nodes; j++){//loop thrugh all the lines below temp lines
			next_one(&t->line);
			if (t->line.current->data.l.type == line_type_temp){//begin if 1
				
				delta = momentum_subtract(p_temp, t->line.current->data.l.phys.p);//momntum difference
				
				/*if (line_type_temp == ELECTRON && t->line.current->data.l.phys.s != spin_temp);
				else{
					if ((fabs(delta.x) + fabs(delta.y) + fabs(delta.z)) < 1E-4) return NOTOK;
				}*/
				
				if (line_type_temp == ELECTRON && t->line.current->data.l.phys.s == spin_temp){//ELECTRON WITH SAME SPIN
					if ((fabs(delta.x) + fabs(delta.y) + fabs(delta.z)) < 1E-5) return NOTOK;	
				}//END OF ELECTRON
				
				else if (line_type_temp == PHONON){//PHONON
					
					if ((worm_end_temp == OK || tell_end_worm(t->line.current) == OK) && t->line.current->data.l.phys.screen == screen_temp){
							if ((fabs(delta.x) + fabs(delta.y) + fabs(delta.z)) < 1E-5) return NOTOK;
						}
					else {
						if ((fabs(delta.x) + fabs(delta.y) + fabs(delta.z)) < 1E-5) return NOTOK;
					}
					
				}//END OF PHONON
				
			}//end if 1
			
		}//end of loop j
		
		nth_node(&t->line, i + 1);
	}//end of loop i
	return OK;
}
//===========================================================
int check_irreducibility_old(struct diag *t){
	/*This function check the irreducibility of a diagram. If there is no line has the same momentum as the
	measuring line, it is irreducible.*/

	//read out the momentum on measuring line
	struct momentum p;
	p = t->meas_line->data.l.phys.p;

	int result = OK;
	int spin, screen;
	spin = screen = 0;

	//read out the type of measuring line, pline or eline, self-energy or polarization sector
	int line_type = t->meas_line->data.l.type;

	if (line_type == PHONON){
		//if the measuring line is pline, read the screen_component
		screen = t->meas_line->data.l.phys.screen;
	}
	else{
		//if the measuring line is eline, read the spin_component, 
		spin = t->meas_line->data.l.phys.s;
	}

	reset_current(&t->line);  //set the current node to be the first node on the line list
	//check whether there is a second line has the same momentum as the measuring line

	if (line_type == PHONON){
		for (int i = 0; i < t->line.n_nodes; i++){
			if (t->line.current->data.l.type == PHONON && t->line.current->data.l.measure == REGULAR &&
				t->line.current->data.l.phys.screen == screen){
				if (fabs(p.x - t->line.current->data.l.phys.p.x) + fabs(p.y - t->line.current->data.l.phys.p.y)
					+ fabs(p.z - t->line.current->data.l.phys.p.z) < 1E-3){
					result = NOTOK;
					break;
				}
			}
			next_one(&t->line);
		}
	}
	else{
		for (int i = 0; i < t->line.n_nodes; i++){
			if (t->line.current->data.l.type == ELECTRON && t->line.current->data.l.measure == REGULAR &&
				t->line.current->data.l.phys.screen == spin){
				if (fabs(p.x - t->line.current->data.l.phys.p.x) + fabs(p.y - t->line.current->data.l.phys.p.y)
					+ fabs(p.z - t->line.current->data.l.phys.p.z) < 1E-3){
					result = NOTOK;
					break;
				}
			}
			next_one(&t->line);
		}

	}
	return result;
}
//===========================================================
/*Change interaction line from F to ~W/J*/
//===========================================================
int change_pline_type(struct node *pline){
	//This function make two ends vertices to be normal

	struct node *vhead, *vtail;
	vhead = pline->data.l.head;
	vtail = pline->data.l.tail;

	//make two ends vertices to be normal
	vhead->data.v.type = NO;
	vtail->data.v.type = NO;

	return OK;
}
//===========================================================
int position_diff(struct node *pline){
	/*This function reads out the position difference at two ends.*/
	int rtn = 2;

	struct node *vhead, *vtail;
	vhead = pline->data.l.head;
	vtail = pline->data.l.tail;

	if ((vhead->data.v.phys.r - vtail->data.v.phys.r) == 1) rtn = 1;
	else if ((vhead->data.v.phys.r - vtail->data.v.phys.r) == -1) rtn = -1;
	else if ((vhead->data.v.phys.r - vtail->data.v.phys.r) == 0) rtn = 0;

	return rtn;
}
//============================================================
//===========================================================
int tell_end_worm (struct node *pline){
	/*This function test whether there is worm on one of pline's end vertex. if there is a worm end, returns OK. otherwise,
	return NOTOK.*/
	
	//if (t->n_wo == 0) return NOTOK;//NO WORMS
	
	struct node *vhead, *vtail;
	vhead = pline->data.l.head;
	vtail = pline->data.l.tail;
	
	if (vhead->data.v.type == NO && vtail->data.v.type == NO) return NOTOK;//no worms at to ends
	
	else return OK;
}
//===========================================================

//Rules of dummy lines:
//1. the dummy interaction line can not be removed or created in any update
//2. if the dummy propagator originating from vertex A is modified by adding/removing an intermediate vertex C, 
//   then A always remains the originating vertex of the new dummy line.
//Note: for interaction dummy line, we can only have tilde_W type, J type pline cannot be a dummy line.

//==================Diagram Order============================
int diag_order(struct diag *t){
	/*This function counts the order of a diagram, namely the number of plines in diagram.*/
	struct list *l;
	l = &t->line;
	reset_current(l);    //set current pointer to the head of line list
	int rtn = 0;
	for (int i = 0; i < t->line.n_nodes; i++){ //loop through the whole line list
		//count the number of plines //when it is not a measuring line
		if (l->current->data.l.type == PHONON /*&& l->current->data.l.measure == REGULAR*/) rtn++;
		next_one(l);
	}
	return rtn;
}
//===========================================================

//==================Number of fermion loops==================
int count(struct diag *t){
	/*This function count the number of loops in the diagram.*/

	int first_lin = NOTDEFINE;
	int loop_order = 0;

	struct list *l;
	l = &t->line;
	reset_current(l);    //set current pointer to the head of line list
	
	for (int i = 0; i < t->line.n_nodes; i++){ //loop through the whole line list

		if (l->current->data.l.type == ELECTRON && l->current->data.print == 0){
			first_lin = l->current->data.id;//set the first eline index number of this loop
			loop_order++;
			l->current->data.print = 1;
			int spin = l->current->data.l.phys.s;
			struct node *v1, *e1;
			v1 = l->current->data.l.head;
			do{
				e1 = v1->data.v.eout;
				e1->data.print = 1;
				//break;?
				v1 = e1->data.l.head;
			} while (v1->data.v.eout->data.id != first_lin);

		}

		next_one(l);
	}

	//then set the identifier .print back to initialization, so we can call this function again and again
	reset_current(l);    //set current pointer to the head of line list
	for (int i = 0; i < t->line.n_nodes; i++){
		if (l->current->data.print != 0) l->current->data.print = 0;
		next_one(l);
	}

    if (t->meas_line->data.l.type == ELECTRON) loop_order--;   //number of fermion loops = virtual loops - 1
	return loop_order;
}
//===========================================================
void assign_parity(struct diag *t){
	/*This function assign the parity to a diagram.*/
	if (count(t) % 2 == 0) t->parity = 0;//modulus operator(%) which returns remainder,# of fermion loops is even number
	else t->parity = 1;
}
//===========================================================
void change_parity(struct diag *t){
	/*Change the parity of a diagram, between 0 to 1.*/
	if (t->parity == 0) t->parity = 1; 
	else t->parity = 0;
}
//===========================================================

//===========================================================
//Function to tell the system is in sector G or P
//===========================================================
int sector(struct diag *t){
	/*This function tells a diagram is in a G-sector or P-sector.*/
	if (t->n_wo != 0) {
		errorquit("not in Physical sector, cannot decide G or P");
		return 1;//make sure diagram in physical 
	}
	int rtn = 0;
	if (t->meas_line->data.l.type == ELECTRON) rtn = G; //measuring line is electron, in G-sector
	else rtn = P;//measuring line is phonon

	return rtn;
}
//===========================================================

double time_diff(struct node *line){
	/*Calculate the time difference of a line.*/
	struct node *v1, *v2;     //two end vertices of e1
	v1 = line->data.l.head;
	v2 = line->data.l.tail;
	double time0 = v1->data.v.phys.t - v2->data.v.phys.t;

	return time0;
}
//===========================================================


//===========================================================
//Functions to manipulate a complex number
//===========================================================
void print_cmplx(complex <double> c){
	/*Print out the real and imgaginary parts of a complex number.*/
	printf("Re: %f  Im: %f\n", real(c), imag(c));
}
//===========================================================

int print_information(struct diag *t){
	/*This function print the information of a diagram.*/

	int loop_order, line_id;
	loop_order = 0;
	struct list *l;
	l = &t->line;     //l to be the line list in the diagram
	reset_current(l); //set current pointer to the first node in the list
	int n = t->n_el;
	for (int i = 0; i < n; i++){

		if (l->current->data.l.type == PHONON){
			if (l->current == l->last) break;
			next_one(l);
		}
		if (l->current->data.l.type == ELECTRON && l->current->data.print != 1){
			//printf("current line node:%d\n", l->current->data.id);
			//printf("did I\n");
			loop_order++;        //the order of loop, nth loop in the diagram
			int n_el_inloop = 1;   //number of elines in this loop
			printf("Loop #:%d \n", loop_order);

			//print the infromation of first eline in this loop
			l->current->data.print = 1;
			line_id = l->current->data.id;
			printf("EL:%d           P[%g,%g,%g]", line_id, l->current->data.l.phys.p.x, l->current->data.l.phys.p.y, l->current->data.l.phys.p.z);

			//spin component of the eline
			if (l->current->data.l.phys.s == UP) printf(" Spin UP\n");
			else printf(" Spin DOWN\n");

			//print the vertex which first eline comes into
			struct node *v1, *p1;
			v1 = l->current->data.l.head; //which vertex first eline coming into 
			p1 = v1->data.v.pend;         //pline connects to above vertex
			printf("V:%d  PL:%d ", v1->data.id, p1->data.id);
			if (v1->data.v.p_in_out == IN) printf(" IN");
			else printf(" OUT");
			printf("  P[%g,%g,%g]", p1->data.l.phys.p.x, p1->data.l.phys.p.y, p1->data.l.phys.p.z);
			//time of this vertex
			printf("  t:%g ", v1->data.v.phys.t);
			//space position of this vertex
			printf("  r:%d \n", v1->data.v.phys.r);
			v1->data.print = 1;
			p1->data.print = 1;
			do{
				//print the eline coming out off the above vertex
				struct node *e1;
				e1 = v1->data.v.eout;
				n_el_inloop++;
				e1->data.print = 1;
				//printf the information of this eline
				printf("Loop #:%d \n", loop_order);
				printf("EL:%d           P[%g,%g,%g]", e1->data.id, e1->data.l.phys.p.x, e1->data.l.phys.p.y, e1->data.l.phys.p.z);
				//spin component of the eline
				if (e1->data.l.phys.s == UP) printf(" Spin UP\n");
				else printf(" Spin DOWN\n");
				v1 = e1->data.l.head; //which vertex first eline coming into 
				p1 = v1->data.v.pend;         //pline connects to above vertex
				printf("V:%d  PL:%d ", v1->data.id, p1->data.id);
				if (v1->data.v.p_in_out == IN) printf(" IN");
				else printf(" OUT");
				printf("  P[%g,%g,%g]", p1->data.l.phys.p.x, p1->data.l.phys.p.y, p1->data.l.phys.p.z);
				//time of this vertex
				printf("  t:%g ", v1->data.v.phys.t);
				//space position of this vertex
				printf("  r:%d \n", v1->data.v.phys.r);
				v1->data.print = 1;
				p1->data.print = 1;

				printf("loop closed\n");
				printf("N propagators in this loop: %d\n", n_el_inloop);
				printf("----------------\n");
			} while (v1->data.v.eout->data.id != line_id);
			if (l->current->data.id == l->last->data.id) break;
			next_one(l);
		}


	}


	printf("--------------------------------------\n");
	//then set the identifier .print back to initialization, so we can call this function again and again
	reset_current(l);
	int n_lin_node = l->n_nodes;
	for (int i = 0; i < n; i++){
		l->current->data.print = 0;
		next_one(l);
	}

	for (int i = 0; i < t->n_ver; i++){
		t->vertex.current->data.print = 0;
		next_one(&t->vertex);
	}

	printf("Number of electron lines: %d\n", t->n_el);
	printf("Number of phonon lines: %d\n", t->n_pl);
	printf("Number of measuring line: %d\n", t->n_ml);
	printf("Number of vertices: %d\n", t->n_ver);
	printf("Number of worms: %d\n", t->n_wo);
	printf("ID of measuring eline: %d\n", t->meas_line->data.id);

	if (0 == t->n_wo) printf("No worms in diagram\n");
	if (2 == t->n_wo){
		printf("ID of I-worm:%d\n", t->s_worm->data.id);
		printf("ID of M-worm:%d\n", t->t_worm->data.id);
	}
	printf("------------------------------------\n");

	return OK;
}
//==========================================================
