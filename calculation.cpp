//===========================================================
//Functions calculate the Green's function and self-energy
//=========================================================== 

#define PLUS 1  
#define MINUS -1

#define G 1           //indicates self-energy sector
#define P 2           //indicates polarization sector



complex <double> green(double time0);
complex <double> green_f(struct profile_multi *h, double time0);
complex <double> green_0(struct node *eline);
complex <double> green_bold(struct profile_multi *h, struct node *eline);
complex <double> w_tilde(struct profile_multi *h, struct node *pline);
complex <double> f_worm(struct profile_multi *h, struct node *pline);
double interaction(struct profile_multi *h, struct node *pline);

double coupling(struct node *pline);
double coupling_M(struct node *pline);

complex <double> digram_value_D(struct diag *t);
complex <double> self_energy(struct diag *t);

int pi_sign(struct diag *t);

double phase_phi(struct diag *t);
complex <double> contribution(struct diag *t);

//===========================================================
complex <double> green(double time0){
	/*This funtion calculate the green's function of a propagator, returns a complex value*/

	complex <double> result;       //declare a complex number to be return
	complex <double> de(-1.0, 1.0);    //define a constant in the denominator of the green's function, class in C++

	complex <double> mu(0.0, -PI*TEM / 2); //chemical potential

	//calculate the green's function
	if (time0 >= 0){
		result = exp(mu*time0) / de;
	}
	else{
		time0 = BETA + time0;
		result = -1.0 * exp(mu*time0) / de;
	}

	return result;
}
//===========================================================
complex <double> green_f(struct profile_multi *h, double time0){
	/*Read out the bold propagator according to time*/
	complex <double> rtn(1.0, 0.0);

	if (h == NULL) rtn = time0;
	else rtn = read_bin(h, time0);
	return rtn;
}
//===========================================================
complex <double> green_0(struct node *eline){
	/*This funtion calculate the free green's function of a propagator, taking a eline as argument, and returns a complex value*/
	complex <double> result(1.0, 0.0);  //declare a complex number to be return
	//if (e1->data.l.measure == MEASURE) return result;

	struct node *v1, *v2;     //two end vertices of e1
	v1 = eline->data.l.head;
	v2 = eline->data.l.tail;
	double time0 = v1->data.v.phys.t - v2->data.v.phys.t; //time interval on e1

	result = green(time0); //calculate the green's function, taking time as an arguement

	return result;
}
//===========================================================
complex <double> green_bold(struct profile_multi *h, struct node *eline){
	/*This funtion calculate the bold green's function of a propagator, taking a eline as argument, and returns a complex value
	obtained by reading the bold green function profile*/
	complex <double> result(1.0, 0.0);  //declare a complex number to be return
	//if (e1->data.l.measure == MEASURE) return result;

	double tnow = time_diff (eline);//time difference at two ends of an electron
	
	if (tnow < 0){
		tnow = tnow + BETA;
		result = -1.0 * green_f(h, tnow);
	}
	else result = green_f(h, tnow); //calculate the green's function, taking time as an arguement

	return result;
}
//===========================================================
complex <double> w_tilde(struct profile_multi *h, struct node *pline){
	/*This funtion calculate the free green's function of a propagator, taking a eline as argument, and returns a complex value*/
	complex <double> result(1.0, 0.0);  //declare a complex number to be return
	//if (e1->data.l.measure == MEASURE) return result;

	double x[2] = {0,0};
	x[0] = double(position_diff(pline));
	x[1] = time_diff (pline);
	
	if (x[1]<0){//if time_diff <0, set t =t +BETA
		x[1] = x[1]+ BETA;
	}

	result = read_binx(h, x); //calculate the green's function, taking time as an arguement

	return result;
}
//===========================================================
complex <double> f_worm(struct profile_multi *h, struct node *pline){
	/*This funtion calculate read out the F(r,t), taking a pline as argument, and returns a complex value*/
	complex <double> result(1.0, 0.0);  //declare a complex number to be return
	
   /*	struct node *vend1, *vend2;
	
	vend1 = pline->data.l.head;
	vend2 = pline->data.l.tail;
	
	if ((vend1->data.v.type == NO) && (vend2->data.v.type == NO)) {//if the pline has no worm on its end
		errorquit("No worm vertex");
	}*/

	double x[2] = {0,0};
	x[0] = double(position_diff(pline));
	x[1] = time_diff (pline);
	
	if (x[1]< 0){//if time_diff <0, set t =t +BETA
		x[1] = x[1]+ BETA;
	}

	result = read_binx(h, x); //calculate the green's function, taking time as an arguement

	return result;
}
//===========================================================
double coupling_M(struct node *pline){
	/*This function returns the coupling Matrix M in Eq.(6).*/
	
	struct node *v1, *v2, *e1_in, *e1_out, *e2_in, *e2_out;
	double rtn = 0;

	v1 = pline->data.l.head;
	v2 = pline->data.l.tail;

	e1_in = v1->data.v.ein;
	e1_out = v1->data.v.eout;
	e2_in = v2->data.v.ein;
	e2_out = v2->data.v.eout;

	int s1_in, s1_out, s2_in, s2_out; //spin components on four elines
	s1_in = e1_in->data.l.phys.s;
	s1_out = e1_out->data.l.phys.s;
	s2_in = e2_in->data.l.phys.s;
	s2_out = e2_out->data.l.phys.s;

	//diagonal coupling
	if (s1_in == s1_out && s2_in == s2_out){
		if (s1_in == s2_in) rtn = 1.0;
		else rtn = -1.0;
	}

	//spin flip term
	if (s1_in == s2_out && s1_out == s2_in && s1_in != s2_in) rtn = (2.0);
	
	if (rtn == 0) printf("wrong J type, spin not conserved\n");
	//printf("J %f\n", rtn);
	return rtn;

}
//===========================================================
double coupling(struct node *pline){
	/*This function returns the magnitude of interaction coupling, depending on the type of the interaction line.*/
	double result = 0;
    struct node *vend1, *vend2;
	
	vend1 = pline->data.l.head;
	vend2 = pline->data.l.tail;
	
	if (vend1->data.v.type != NO || vend2->data.v.type != NO) return F0; // if pline has a worm vertex on its end
	
	int screen = pline->data.l.phys.screen;
	
	if (screen == J) result = coupling_M(pline)*J0/4.0;
	else result = W0;//read from the ~W profile

	return result;
}
//===========================================================
double interaction(struct profile_multi *h, struct node *pline){
	/*This function returns the magnitude of interaction coupling, depending on the type of the interaction line.*/
	double rtn;
    struct node *vend1, *vend2;
	
	vend1 = pline->data.l.head;
	vend2 = pline->data.l.tail;
	
	if (vend1->data.v.type != NO || vend2->data.v.type != NO) {//if the pline has a worm on its end
		//rtn = f_worm(h,pline);
		//return rtn; // if pline has a worm vertex on its end
		return F0;
	}
	
	int screen = pline->data.l.phys.screen;
	
	if (screen == J) rtn = coupling_M(pline)*J0/4.0;
	//if (screen == J) rtn = assign_complx(J0/4,0.0);//coupling_M(pline);
	else rtn = coupling_M(pline)*real(w_tilde(h,pline));//read from the ~W profile

	return rtn;
}
//===========================================================

//===========================================================
complex <double> digram_value_D(struct diag *t){
	/*This function calculate the value of the digram, products of Green's function of each propagator and interaction coupling.*/
	if (t->n_wo != 0) return 1;//make sure diagram in physical sector

	complex <double> rtn(1.0, 0.0);
	//double couple = 0;
	struct list *l;
	l = &t->line;
	reset_current(l);    //set current pointer to the head of line list

	for (int i = 0; i < t->line.n_nodes; i++){ //loop through the whole line list
		//printf("i:%d", i);
		//print_cmplx(rtn);
		if (l->current->data.l.measure == MEASURE) rtn = rtn;  //if this node is a measuring line
		else {                                                 //if this line is a regular line
			if (l->current->data.l.type == ELECTRON) {//for a propagator, calculate green's function
				//rtn = green_0(l->current)*rtn;
			      rtn = green_bold(&(t->prf_g),l->current)*rtn;
			}
			else {//for an interaction line, calculate the coupling
				//couple = coupling(l->current);
				rtn = interaction(&(t->prf_w),l->current)*rtn;
			}
		}

		next_one(l);
	}
	//printf("D\n");
    //print_cmplx(rtn);
	return rtn;
}
//===========================================================
double phase_phi(struct diag *t){
	/*Calculate the phase factor.*/
	double rtn = 0;
	//int loop = t->parity;
	int loop = count(t);
	int order = diag_order(t);
	rtn = arg(digram_value_D(t)) + PI*(loop + order);
	return rtn;
}
//===========================================================
int pi_sign(struct diag *t){
	/*This function reads out the spin components of the polarization.*/
	if (t->meas_line->data.l.type != PHONON){
		errorquit("not a P-sector, pi_sign is void");
		return 1;
	}
	//read out the two vertices which connect the measuring pline
	struct node *v1, *v2;
	v1 = t->meas_line->data.l.head; 
	v2 = t->meas_line->data.l.tail;

	//read out the two propagators (outgoing ones)
	struct node *e1, *e2;
	e1 = v1->data.v.eout;
	e2 = v2->data.v.eout;

	int s1, s2; //spin of e1 and e2
	s1 = e1->data.l.phys.s;
	s2 = e2->data.l.phys.s;

	int result = s1*s2;
	
	//if (result != 1 && result !=-1) printf("Nah\n");
	
	if (result == PLUS) return PLUS;
	else return MINUS;
}
//===========================================================
complex <double> self_energy(struct diag *t){
	/*diagrammatic contribution from the configuration space.*/
	complex <double> rtn(1.0, 0.0);
	complex <double> i_unit(0.0, 1.0); //chemical potential
	rtn = exp(i_unit*phase_phi(t));

	return rtn;
}
//===========================================================
//for G/P sector, Function to get digrammatic contribution from the configuration space.*/
complex <double> contribution(struct diag *t){
	/*diagrammatic contribution from the configuration space.*/
	complex <double> rtn(1.0, 0.0);
	complex <double> i_unit(0.0, 1.0); //Unit imaginary 
	rtn = exp(i_unit*phase_phi(t));

	if (sector(t) == P) {
		if (MINUS == pi_sign(t)) rtn = -1.0*rtn;//if two propagator have different signs
	}
	return rtn;
}
//===========================================================