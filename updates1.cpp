/*14 updates*/
//#include"create.cc"


/*The probability of each update.*/
#define u1 0.125       //Create, probability
#define u2 0.125       //Delete, probability
#define u3 0.125       //Create-H, probability
#define u4 0.125       //Delete-H, probability
#define u5 0.05        //Move-P, probability
#define u6 0.05        //Move-I, probability
#define u7 0.05        //Commute, probability
#define u8 0.05        //Dummy, probability
#define u9 0.05        //Insert, probability
#define u10 0.05       //Remove, probability
#define u11 0.05       //Dress, probability
#define u12 0.05       //Undress, probability
#define u13 0.05       //Recolor, probability
#define u14 0.05       //Move-T, probability



int Create(struct diag *t);
int Delete(struct diag *t);
int Create_H(struct diag *t);
int Delete_H(struct diag *t);
int Move_P(struct diag *t);
int Move_I(struct diag *t);
int Commute(struct diag *t);
int Dummy(struct diag *t);
int Dummy_e(struct diag *t);
int Insert(struct diag *t);
int Remove(struct diag *t);
int Dress(struct diag *t);
int Undress(struct diag *t);
int Recolor(struct diag *t);
int Move_T(struct diag *t);


//UPDATES
//===========================================================
/*The pair of complementary updates "Create-Delete" switches between physical and unphysical sectots by
inserting/ reamoving a pair of special verices connected by the particle propagator.*/
//===========================================================
int Create(struct diag *t){
	/*This update, the propagator and the type of special vertex to be placed at its tail are selected at random. One has to verify that flipping the
	propagator spin is consistent with the worm rules or reject the proposal. The minssing momentum p_w at the worm vertex is selected at random.*/

	if (t->n_wo != 0) return 1; //make sure we have no worms, in a physical sector
	
	if (check_irreducibility(t) == NOTOK) return 1;

	struct node *eline, *a, *b;
	eline = rand_eline(&t->line); //randomly select a propagator
	a = eline->data.l.tail;       //tail vertex pointer of the selected propagator
	b = eline->data.l.head;       //head vertex pointer of the selected propagator
	enum spin spin = eline->data.l.phys.s;//save the original spin component 

	//printf("selected eline id:%d\n", eline->data.id);

	//if selected eline form a bubble loop, reject update
	if (a == b) return 1;

	//find two plines connects to vertice A and B
	struct node *p1, *p2;
	p1 = a->data.v.pend;
	p2 = b->data.v.pend;

	//two worms cannot be connected by a interaction line
	if (p1 == p2) return 1;

	//find two vertices pairing with vertex A and B
	struct node *c, *d;
	c = pair_ver(a);
	d = pair_ver(b);

	//Calculate the probability for the update to take place
	int n = diag_order(t);  //order of the diagram 

	//Calculate the probability denominator, before update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline)*interaction(&(t->prf_w),p1)*interaction(&(t->prf_w),p2);
	
	//print_cmplx(green_bold(&(t->prf_g),eline));
    //print_cmplx(interaction(&(t->prf_w),p1));
	
	//printf("A: %d B: %d  P1 :%f\n",a->data.v.type, b->data.v.type,real(unphys_F(&(t->prf_f),p1)));
	//Calculate the probability numerator, after update
	complex <double> ratio_nu = green_bold(&(t->prf_g),eline)*f_worm(&(t->prf_f),p1)*f_worm(&(t->prf_f),p2);
	//complex <double> ratio_nu = green_bold(&(t->prf_g),eline)*F0*F0;

	//probability
	double ratio = abs(ratio_nu / ratio_de)*(2 * n)*u2 / u1;
	//printf("create:%f\n", ratio);

	if (rand_01() > ratio) { //if probability is too small, do nothing
		return 1;
	}

	//assign speciality to the tail vertex 
	rand_worm_tobe(a);
	int worm_type = a->data.v.type; //save special vertex type.

	//assign speciality to the head vertex 
	if (worm_type == S) b->data.v.type = T;
	else b->data.v.type = S;

	//check spin component on selected eline, whether it is consistent with the worm type on its end
	if ((eline->data.l.phys.s == UP && a->data.v.type == T) || (eline->data.l.phys.s == DOWN && a->data.v.type == S)) {
		a->data.v.type = NO; //set vertex A and B back to normal
		b->data.v.type = NO;
		return 1;
	}

	//flip spin, now a and b are worms
	if (worm_type == S) eline->data.l.phys.s = DOWN;
	else eline->data.l.phys.s = UP;

	//assign random excess momentum
	struct momentum delta = momentum_init(0.0, 0.0, 0.0); //initialize momentum excess on I worm
	delta = assign_ran_momentum();

	//change the momentum on the selected eline, and create two worms 
	if (worm_type == S) {
		eline->data.l.phys.p = momentum_subtract(eline->data.l.phys.p, delta);
		t->s_worm = a; //make the tail end a S worm
		t->t_worm = b; //make the head end a T worm
	}
	else {
		eline->data.l.phys.p = momentum_add(eline->data.l.phys.p, delta);
		t->s_worm = b; //make the head end a S worm
		t->t_worm = a; //make the tail end a T worm
	}

	//printf("S,T WORM : %d %d\n", t->s_worm->data.id, t->t_worm->data.id);
	t->n_wo = t->n_wo + 2;

	//the excess momentum = difference with incoming and outgoing momenta in I worm
	t->delta = momentum_excess_ver(t->s_worm);
	//print_momentum(t->delta);

	//change the screen type of two plines, no need
	//p1->data.l.phys.p_f = F;
	//p2->data.l.phys.p_f = F;

	return OK;
}
//===========================================================
int Delete(struct diag *t){
	/*This delete update, one selects a wrom at random, checks that the outgoing propagator arrives at the other worm, and proposes to remove
	the worms from the diagram.*/
	if (t->n_wo != 2) return 1; // make sure there are two worms in the diagram
	if (check_irreducibility(t) == NOTOK) return 1;
	
	int n = diag_order(t);//diagram order

	struct node *a, *b, *eline1;
	a = select_worm(t);//select a worm at random, worm 1 A
	int worm_type = a->data.v.type;//save the worm type of vertex A to a int
	//printf("vertex A :%d\n", a->data.id);

	eline1 = a->data.v.eout; //eline1 would be the outgoing propagator of worm 1
	enum spin s1 = eline1->data.l.phys.s; //stores the original spin on the selected eline

	b = eline1->data.l.head;//head vertex of eline1
	if (b->data.v.type == NO || a == b) return 1; //if vertex B is not a worm, reject update

	//if we have different spin configuration as Fig.5, reject this proposal
	if ((worm_type == S && s1 == UP) || (worm_type == T && s1 == DOWN)) return 1;

	//find two plines connects to vertice A and B
	struct node *p1, *p2;
	p1 = a->data.v.pend;
	p2 = b->data.v.pend;
	if (p1 == p2) return 1;//two worms cannot be connected by an interaction line

	//calculate the denominator of probability, before the update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*f_worm(&(t->prf_f),p1)*f_worm(&(t->prf_f),p2);
	//complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*F0*F0;
	//printf("delete-de: %g  %g ", real(ratio_de),imag(ratio_de));

	/*Now, make changes to diagram*/

	//1.flip the spin, need for calculation the probability
	if (worm_type == S) {
		eline1->data.l.phys.s = UP;
	}
	else {
		eline1->data.l.phys.s = DOWN;
	}

	//2.change the type of two plines, including 3.elimination of worms
	change_pline_type(p1);
	change_pline_type(p2);

	complex <double> ratio_nu = green_bold(&(t->prf_g),eline1)*interaction(&(t->prf_w),p1)*interaction(&(t->prf_w),p2);
	//printf("delete-nu: %g  %g ", real(ratio_de),imag(ratio_nu));

	//Calculate the probability
	double ratio = abs(ratio_nu / ratio_de) / (2 * n)*u1 / u2;
	//printf("delete:%g\n", ratio);

	if (ratio >= 1) ratio = 1;
	else{
		if (rand_01() > ratio) {//if probability is too small, remove all the changes

			//1.flip spin back to orignal
			eline1->data.l.phys.s = s1;

			//2.interaction line type
			//p1->data.l.phys.p_f = F;
			//p2->data.l.phys.p_f = F;

			//2.reset two worms
			if (worm_type == S){
				a->data.v.type = S;
				b->data.v.type = T;
			}
			else{
				a->data.v.type = T;
				b->data.v.type = S;
			}

			return 1;
		}
	}

	//if probability is large enough, do the update

    //momentum before change
	struct momentum m1 = eline1->data.l.phys.p;
	//change momenta
	if (worm_type == S) eline1->data.l.phys.p = momentum_add(eline1->data.l.phys.p, t->delta);
	else eline1->data.l.phys.p = momentum_subtract(eline1->data.l.phys.p, t->delta);
	
	if (check_irreducibility(t) == NOTOK) {
		eline1->data.l.phys.s = s1;

		//2.interaction line type
		//p1->data.l.phys.p_f = F;
		//p2->data.l.phys.p_f = F;

		//2.reset two worms
		if (worm_type == S){
			a->data.v.type = S;
			b->data.v.type = T;
		}
		else{
			a->data.v.type = T;
			b->data.v.type = S;
		}
		
		//recover momentum
		eline1->data.l.phys.p = m1;

		return 1;
		
	}


	//change the type of the two interaction lines
	t->s_worm = NULL;            //remove S worm
	t->n_wo--;                   //decrease number of worms by 1
	t->t_worm = NULL;            //remove T worm
	t->n_wo--;                   //decrease number of worms by 1
	t->delta = momentum_init(0.0, 0.0, 0.0); //no momentum excess

	return OK;

}
//===========================================================
/*The pair "Create-H/Delete-H" swithches between physical and unphysical with an additional ingredient:
it increase/decrease diagram order by attaching a Hartree-type bubble to the exsiting graph.
One of the worms is placed on the bubble vertex.*/
//===========================================================
int Create_H(struct diag *t){
	/*In create_H, the propagator and the type of special vertex to be placed at its head is selected at random.This digram remains irreducible, because
	ont of the worms is placed on the bubble vertex.*/
	/*In this update, a propagator line (going from vertex A to vertex B) and whether to place S or T on vertex B is decided at random.*/

	if (t->n_wo != 0) return 1;   //make sure we have no worms, in a physical sector
	if (diag_order(t) >= Max_order) return 1; //when diagram has up-limit order, reject the updates which will increase order 
	if (check_irreducibility(t) == NOTOK) return 1;

	int n = diag_order(t);        //order of the diagram 

	struct node *eline, *pline, *a, *b;
	eline = rand_eline(&t->line); //randomly select a propagator
	a = eline->data.l.tail;       //tail vertex pointer of the selected propagator
	b = eline->data.l.head;       //head vertex pointer of the selected propagator
	pline = b->data.v.pend;       //the interaction line connecting vertex B

	//printf("pline %d\n",pline->data.id);

	int spin, screen;
	spin = eline->data.l.phys.s;   //original spin on the selected eline
	screen = pline->data.l.phys.screen; //original screen type on the pline


	//calculate the denominator of the acceptance ratio, before the update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline)*interaction(&(t->prf_w),pline);
	//printf("pline-type: %d\n", pline->data.l.phys.screen);
	//printf("pline-value: %f\n",interaction(&(t->prf_w),pline));

	//change diagram
	//assign speciality(type of worm: S/T) at random to the head vertex v_head of selected eline
	rand_worm_tobe(b);
	int worm_type = b->data.v.type;//save the specility of Vertex B to an int
	//printf("vertex B worm,%d\n", b->data.v.type);

	//check the worm rules, if the proposal is inconsistent with worm rules, reject update
	//if vertex B is T, the worm on the bubble should be S, then the incoming line spin must to be UP, 
	if (worm_type == T){
		if (spin != UP) {
			b->data.v.type = NO;
			return 1;
		}
	}
	else{
		if (spin != DOWN) {
			b->data.v.type = NO;
			return 1;
		}
	}

	//change diagram

	//create a new eline and a new vertex C to insert into the selected eline
	struct node *e_insert, *c;
	c = create_ver(t);
	c->data.v.phys.t = rand_time();  //time of vertex C is generated form probability distribution
	c->data.v.phys.r = a->data.v.phys.r;//position of vertex C

	//creat a new eline, spin component is decided according to worm rule
	if (worm_type == T){
		e_insert = create_lin(t, ELECTRON, REGULAR, DOWN);
	}
	else{
		e_insert = create_lin(t, ELECTRON, REGULAR, UP);
	}

	//create 1 eline + 1 pline + 1 vertex, Hartree bubble
	struct node *eline1, *p1, *v1;

	eline1 = create_e_rand(t);//create a eline, spin is decided randomly   
	p1 = create_p_rand(t);//create a pline, J/W is decided randomly 
	//printf("new-p-type: %d\n", p1->data.l.phys.screen);
	v1 = create_ver(t);   //create a vertex

	//assign speciality to the new vertex V1
	if (worm_type == T) v1->data.v.type = S;
	else v1->data.v.type = T;

	//assign time and position to v1,according to the type of interaction p1
	if (p1->data.l.phys.screen == J) {
		v1->data.v.phys.t = c->data.v.phys.t;

		//in two spin system, only one neighbour, two sites 0 or 1
		if (c->data.v.phys.r == 0) v1->data.v.phys.r = 1;
		else v1->data.v.phys.r = 0;

		//int r = rand() % 2;   //r +/- 1, two neighbor sites
		//if (r == 0) v1->data.v.phys.r = 1 + c->data.v.phys.r; //r + 1
		//else v1->data.v.phys.r = -1 + c->data.v.phys.r; //r - 1
	}
	else {
		v1->data.v.phys.t = rand_time();
		v1->data.v.phys.r = rand() % N;  //random position, interaction is ~W type, here in two spin system N is 2
	}

	connect_lin(eline1, v1, v1); //connect bubble lines and vertex, create a bubble

	connect_lin(p1, v1, c); //connect two new vertices, the momentum direction in the pline is decided randomly

	//connect new eline and new vertices along the selected eline
	connect_lin(eline, c, a);
	connect_lin(e_insert, b, c);

	//calculate the numerator acceptance ratio,after the update
	complex <double> ratio_nu = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),eline)*green_bold(&(t->prf_g),e_insert)//*F0*F0;
		*f_worm(&(t->prf_f),p1)*f_worm(&(t->prf_f),pline);
	//printf("DE: %g %g NU: %g %g \n",real(ratio_de),imag(ratio_de),real(ratio_nu),imag(ratio_nu));

	//calculate the probability distribution in Eq (19) and (20)
	double factor;
	if (p1->data.l.phys.screen == J) factor = Z;
	else factor = 1.0 * N / TEM;
	//printf("factor:%f\n", factor);
	//printf("d/d%f\n",abs(ratio_nu / ratio_de));

	//calculate the acceptance ratio
	double ratio;
	ratio = abs(ratio_nu / ratio_de)*pow(2, 3)*n / TEM * factor*u4 / u3;
	//printf("create-h:%f\n", ratio);

	if (ratio >= 1) ratio = 1;
	else{
		if (rand_01() > ratio) { //if probability is too small, remove all modification

			delete_lin(t, e_insert);
			delete_ver(t, c);
			delete_lin(t, eline1);
			delete_lin(t, p1);
			delete_ver(t, v1);

			connect_lin(eline, b, a);

			b->data.v.type = NO;

			return 1;

		}
	}

	//send momentum to every line
	eline1->data.l.phys.p = assign_ran_momentum();  //send random momenta to eline1
	p1->data.l.phys.p = assign_ran_momentum();  //send random momenta to p1
	//print_momentum(p1->data.l.phys.p);

	//calculate momentum on e_insert propagator 
	e_insert->data.l.phys.p = momentum_subtract(eline->data.l.phys.p, p1->data.l.phys.p);

	if (worm_type == S) {
		t->s_worm = b;
		t->t_worm = v1;
	}
	else{
		t->s_worm = v1;
		t->t_worm = b;
	}
	t->n_wo = 2;

	//set the interaction type to be F
	//pline->data.l.phys.screen = F;
	//p1->data.l.phys.screen = F;

	//the excess momentum = difference with incoming and outgoing momenta in S worm
	t->delta = momentum_excess_ver(t->s_worm);
	//printf("id of pline %d p1 %d\n", pline->data.id, p1->data.id);
	//printf("ID s:%d\n", t->s_worm->data.id);
	return OK;
}
//===========================================================
int Delete_H(struct diag *t){
	/*In this update, a random choice id made as to what type of special vertex must be on the Hartree bubble provided that
	the overall topology of lines is identical to that on the tight-hand side of Fig.6. The proposal is to remove worms and
	the bubble from the diagram. It is rejected if either the propagator line originating from vertex C or the interaction
	line attached to C is the dummy one.*/

	if (t->n_wo != 2) return 1; // make sure there are two worms in the diagram
	if (check_irreducibility(t) == NOTOK) return 1;

	//check topology, 

	//select a worm at random
	struct node *v1, *eline;;
	v1 = select_worm(t);
	eline = v1->data.v.ein; //eline on the H-bubble

	if (v1->data.v.eout != eline) return 1;//check if the selected worm is on bubble, if not, reject
	int worm_type = v1->data.v.type;//the type of worm(S/T) on vertex V1

	//find vertex c, C is on the other end of the pline connecting the bubble, it is the pair vertex of the worm on bubble
	struct node *c;
	c = pair_ver(v1);

	if (c->data.v.type != NO) return 1; //worms cannot be connecting by a interaction line

	//find  the elines and pline connect to vertex C
	struct node *eline1, *e2, *p1;
	eline1 = c->data.v.eout;  //eline coming out from C
	e2 = c->data.v.ein;   //eline coming into C
	p1 = c->data.v.pend;  //pline connects to C

	//printf("p1:%d\n", p1->data.l.phys.screen);
	if (eline1 == e2) return 1;

	//if eline comes out c is measuring eline or the interaction line connecting c is measuring pline, updated is rejected
	if (eline1->data.l.measure == MEASURE || p1->data.l.measure == MEASURE || eline->data.l.measure == MEASURE) return 1;

	//find vertex A and B
	struct node *a, *b;
	b = eline1->data.l.head;
	a = e2->data.l.tail;

	if (b->data.v.type == NO) return 1; //if this vertex is not a worm, reject update

	if (c == b || c == a) return 1;//if vertex C is same as vertex A or vertex B,reject
	//if (a == b) return 1;

	//find vertex D, the pair vertex of B
	struct node *d;
	d = pair_ver(b);
	//if (d == a) return 1;
	
	//pline connecting vertex B
	struct node *pline;
	pline = b->data.v.pend;

	//calculate the denominator of probability, before the update
	int n = diag_order(t);

	complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2)*green_bold(&(t->prf_g),eline)//*F0*F0;//
	*f_worm(&(t->prf_f),p1)*f_worm(&(t->prf_f),pline);

	//calculate the numerator of probability, after the update

	//change the the two worms to be normal
	change_pline_type(pline);//change the pline of B to normal, J or W~

	connect_lin(e2, b, a);
	//coupling pline depends the spin component of incoming line 
	complex <double> ratio_nu = green(b->data.v.phys.t - a->data.v.phys.t)*interaction(&(t->prf_w),pline);
	//printf("type of pline after change:%d\n", pline->data.l.phys.screen);

	//Calculate the probability distribution
	change_pline_type(p1);


	double factor;
	if (p1->data.l.phys.screen == J) factor = 1.0 / Z;
	else factor = 1.0 / N * TEM;
	//printf("factor %f\n", factor);

	//calculate probability
	double ratio;
	ratio = abs(ratio_nu / ratio_de)*TEM / (n - 1) / 8 * factor * u3 / u4;
	//printf("Id of A:%d, B:%d C:%d v1:%d d:%d\n", a->data.id, b->data.id, c->data.id, v1->data.id, d->data.id);
	//printf("delete-h:%f\n", ratio);

	if (ratio >= 1) ratio = 1;
	else{
		if (rand_01() > ratio) { //if probability is too small,remove all the change
			//pline->data.l.phys.screen = F;
			//p1->data.l.phys.screen = F;

			if (worm_type == S) {
				v1->data.v.type = S;
				b->data.v.type = T;
			}
			else {
				v1->data.v.type = T;
				b->data.v.type = S;
			}
			connect_lin(e2, c, a);
			connect_lin(eline1, b, c);
			return 1;
		}

	}
	//printf("probability:%f\n", ratio);

	//delete the bubble, a eline and a vertex and p1 
	delete_lin(t, eline);   //delete eline on bubble
	delete_ver(t, v1);      //delete v1
	delete_lin(t, p1);      //delete p1

	//delete vertex C and the eline coming from C to vertex b
	delete_lin(t, eline1);
	delete_ver(t, c);

	//change the diagram 
	t->s_worm = NULL;            //remove I worm
	t->n_wo--;                   //decrease number of worms by 1
	t->t_worm = NULL;            // remove M worm
	t->n_wo--;                   //decrease number of worms by 1
	t->delta = momentum_init(0.0, 0.0, 0.0); //no momentum excess

	//printf("measuring line after:%d\n", t->meas_line->data.id);

	// if (pline->data.l.phys.screen == J);
	//printf("value:%f\n", interaction(&(t->prf_w),pline));

	return OK;
}
//===========================================================
int Move_P(struct diag *t){
	/*Move worm along propagator line.This self-complementary update, one selects at random S or T and proposes to shift the
	selected worm along the incomging or outgoing propagator line, decide again randomly between the two choices. 4 ways*/

	if (t->n_wo != 2) return 1; //make sure there are two worms in the diagram
	if (t->s_worm == NULL || t->t_worm == NULL) return 1;

	struct node *a, *b, *eout, *ein, *c, *d, *eline, *pline;
	int move_direction;

	a = select_worm(t);//select a worm at random (S or T)
	eout = a->data.v.eout; //propagator coming out from Vertex A
	ein = a->data.v.ein;   //propagator coming into Vertex A
	enum ver_worm w_type = a->data.v.type; //stores the original worm type of Vertex A
	if (eout == ein) return 1;//if this worm on bubble,reject update
	c = pair_ver(a);

	//select move direction at random (left or right), along ein or eout
	move_direction = rand() % 2; //generate 0 or 1.  0 move worm the right,along eout;1 move worm to the left

	//find eline (along which propergator the worm moves) and vertex B
	if (move_direction == 0){//move to the right, along outgoing propergator line
		eline = eout;
		b = eout->data.l.head; //vertex B is the head vertex of eout
	}
	else {//move to the left
		eline = ein;
		b = ein->data.l.tail; //vertex B is the tail vertex of ein, 
	}
	enum spin s = eline->data.l.phys.s;//stores the original spin on eline

	pline = b->data.v.pend;//pline connecting vertex B
	d = pair_ver(b);//find the pair vertex of vertex B
	if (c == d) return 1;  //if vertices C and D are the same, reject
	if (c == b) return 1;

	//check if vertex B or vertex D is the other worm, if it is, reject update
	if (w_type == S){
		if (b == t->t_worm || d == t->t_worm) return 1;
	}
	else{
		if (b == t->s_worm || d == t->s_worm) return 1;
	}

	//check spin, worm rule

	//check spin, if moving to the right, along eout
	if (move_direction == 0){//along eout
		if (w_type == S){
			// if spin of eline coming into B is UP, reject
			if (eline->data.l.phys.s == UP) return 1;
		}
		if (w_type == T){
			// if spin of eline coming into B is DOWN, reject
			if (eline->data.l.phys.s == DOWN) return 1;
		}
	}

	//check spin, if moving to the left
	else{//along ein, along ein
		if (w_type == S){
			// if spin of eline coming into B is UP, reject
			if (eline->data.l.phys.s == DOWN) return 1;
		}
		if (w_type == T){
			// if spin of eline coming into B is DOWN, reject
			if (eline->data.l.phys.s == UP) return 1;
		}
	}

	//pline connects to Vertex A
	struct node *pline2;
	pline2 = a->data.v.pend;
	
	//caluculate probability denominator, before the update
	//complex <double> ratio_de = green_bold(&(t->prf_g),eline)*F0*interaction(&(t->prf_w),pline);
    complex <double> ratio_de = green_bold(&(t->prf_g),eline)*f_worm(&(t->prf_w),pline2)*interaction(&(t->prf_w),pline);
	//calculate probability numerator, after the update

	//change the spin component on eline
	if (move_direction == 0){//move to the right
		if (w_type == S) eline->data.l.phys.s = UP;
		else eline->data.l.phys.s = DOWN;
	}
	else{
		if (w_type == S) eline->data.l.phys.s = DOWN;
		else eline->data.l.phys.s = UP;
	}

	
	//change the Vertex A and B type
	if (w_type == S) b->data.v.type = S;
	else b->data.v.type = T;

	a->data.v.type = NO;

	//change the pline2 type from F back to W 
	change_pline_type(pline2);

	complex <double> ratio_nu = green_bold(&(t->prf_g),eline)*f_worm(&(t->prf_w),pline)*interaction(&(t->prf_w),pline2);

	//calculate the probability
	double ratio = abs(ratio_nu / ratio_de);
	//printf("ratio of move-p %f\n", ratio);

	if (rand_01() > ratio) { //if probability is too small,remove all the change
		//set the eline spin back to normal
		eline->data.l.phys.s = s;

		//set the Vertex A back to be the original worm
		a->data.v.type = w_type;

		//set the Vertex B back to be normal 
		b->data.v.type = NO;

		//set pline2 (connecting A) to F
		//pline2->data.l.phys.screen = F;

		return 1;
	}

	//read out the momentum of eout and e2
	struct momentum p1, p2;
	p1 = eout->data.l.phys.p;
	p2 = ein->data.l.phys.p;

	//change momentum on vertex B
	if (w_type == S){
		if (move_direction == 0) eout->data.l.phys.p = momentum_add(p1, t->delta);
		else ein->data.l.phys.p = momentum_subtract(p2, t->delta);
	}
	else{
		if (move_direction == 0) eout->data.l.phys.p = momentum_subtract(p1, t->delta);
		else ein->data.l.phys.p = momentum_add(p2, t->delta);
	}
	
    //enum pline scr_pline = pline->data.l.phys.screen;
	
	//make the pline which connecting the new worm to be F type
	//pline->data.l.phys.screen = F;


   // printf("check-p:%d\n",check_irreducibility(t));
	
	//check irreducibility, activated
	if (check_irreducibility(t) == NOTOK) {
		//printf("irreducible move_p\n");
		
		//set the eline spin back to normal
		eline->data.l.phys.s = s;

		//set the Vertex A back to be the original worm
		a->data.v.type = w_type;

		//set the Vertex B back to be normal 
		b->data.v.type = NO;

		//set pline2 (connecting A) to F
		//pline2->data.l.phys.screen = F;
		
	    //recover the momentum
		eout->data.l.phys.p = p1;
		ein->data.l.phys.p = p2;
		
		//pline->data.l.phys.screen = scr_pline;
		
		return 1;
	}

	
	//change the s_worm and t_worm identifier in diagram
	if (w_type == S) t->s_worm = b;
	else t->t_worm = b;

	return OK;
}
//===========================================================
int Move_I(struct diag *t){
	/*Move worm along interaction line.This is self-complementary update, one selects at random S or T and
	propose to shift the selected worm along the interaction line to vertex B.*/

	if (t->n_wo != 2) return 1; // make sure there are two worms in the diagram

	struct node *a, *b, *p1;

	//select a worm at random (S or T)
	a = select_worm(t);

	//save the worm_type on Vertex A
	int worm_type = a->data.v.type;

	p1 = a->data.v.pend; //move worm along this interaction line
	b = pair_ver(a);     //Vertex B, the proposal new worm

	//for the probability, always 1 in move_I

	if (b->data.v.type != NO) return 1;

	//read out the momenta on p1
	struct momentum q;
	q = p1->data.l.phys.p;

	//change momentum on interaction line
	if (worm_type == S){
		if (a->data.v.p_in_out == IN) p1->data.l.phys.p = momentum_subtract(q, t->delta);
		else p1->data.l.phys.p = momentum_add(q, t->delta);
	}
	else {
		if (a->data.v.p_in_out == IN) p1->data.l.phys.p = momentum_add(q, t->delta);
		else p1->data.l.phys.p = momentum_subtract(q, t->delta);
	}

	//printf("check-i:%d\n",check_irreducibility(t));
	//check irreducibility, ACTIVATED!
	if (check_irreducibility(t) == NOTOK) {
		//printf("irreducible move_I %d\n", NOTOK);
		//printf("========\n");
		
		p1->data.l.phys.p = q;//recover the momentum
		
		return 1;
	}
	
	//Apply the update
	if (worm_type == S){
		b->data.v.type = S;  //set speciality to vertex B
		t->s_worm = b;
	}
	else {
		b->data.v.type = T;
		t->t_worm = b;
	}

	a->data.v.type = NO;    //set selected vertex back to normal

	//if (t->t_worm->data.v.ein == t->s_worm->data.v.eout) printf("move-i : two worms on one loop\n");
	//printf("ID s:%d\n", t->s_worm->data.id);

	return OK;
}
//===========================================================
int Commute(struct diag *t){
	/*randomly connecting outgoing propagator to incoming ones. The self-complementary update proposes to swap destination vertices for
	propagator line originating at S or T. This proposal is valid only all vertices have the same space coordinate and both propagator
	have the same spin index.*/
	if (t->n_wo != 2) return 1; // make sure there are two worms in the diagram
	if (check_irreducibility(t) == NOTOK) return 1;
	

	//check the position of two worms, if they not same, reject update
	if (t->s_worm->data.v.phys.r != t->t_worm->data.v.phys.r) return 1;

	struct node *eline1, *e2, *v1, *v2;

	eline1 = t->s_worm->data.v.eout;  //propagator coming out from S worm
	e2 = t->t_worm->data.v.eout;  //propagator coming out from T worm

	//spin indexes of two elines
	int s1, s2;
	s1 = eline1->data.l.phys.s;
	s2 = e2->data.l.phys.s;

	//check the spins on two elines, if they not same, reject update
	if (s1 != s2) return 1;

	//read out momenta on the incomging ends of S/T 
	struct momentum p1, p2;
	p1 = eline1->data.l.phys.p;
	p2 = e2->data.l.phys.p;

	//read out destination end of the two elines, besides worms
	v1 = eline1->data.l.head;
	v2 = e2->data.l.head;

	if (v1 == v2) return 1;
	//if (v1 == t->t_worm) return 1;
	//if (v2 == t->s_worm) return 1;
	//if (v1 == t->s_worm || v1 == t->t_worm) return 1;
	//if (v2 == t->s_worm || v2 == t->t_worm) return 1;

	//check the position of v1 and v2, if they not same, reject update
	if (v1->data.v.phys.r != v2->data.v.phys.r) return 1;

	//caluculate probability denominator, before the update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2);

	//calculate probability numerator, after the update

	//swaps heads of eline1 and e2 
	connect_lin(eline1, v2, t->s_worm);
	connect_lin(e2, v1, t->t_worm);

	complex <double> ratio_nu = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2);

	//probability 
	double ratio = abs(ratio_nu / ratio_de);
	//printf("commute %f\n", ratio);

	if (ratio >= 1) ratio = 1;
	else{//if ratio smaller than 1, probability is ratio, need compare it to random number
		if (rand_01() > ratio) { //if probability is too small,remove all the change

			connect_lin(eline1, v1, t->s_worm);
			connect_lin(e2, v2, t->t_worm);
			return 1;
		}
	}

	//swap momentum on two elines
	eline1->data.l.phys.p = p2;
	e2->data.l.phys.p = p1;

	//change momentum excess delta
	t->delta = momentum_add(t->delta, momentum_subtract(p1, p2));

	return OK;
}
//===========================================================
int Dummy(struct diag *t){
	/*We select one of the vertices at random, say vertex A, and make a random decision whether the new dummy line should be the interaction or the
	propagator line originationg from A.*/
	if (check_irreducibility(t) == NOTOK) return 1;
	

	struct node *m_line;        //read out the measuring line pointer
	m_line = t->meas_line;

	struct node *v1, *eline1, *p1;
	v1 = rand_ver(&t->vertex);  //choose a vertex at random

	//read out the eline and pline connecting to this vertex 
	eline1 = v1->data.v.eout;
	p1 = v1->data.v.pend;

	//generate a const 0 or 1. //0 pline to be new dummy interaction line, 1 eline coming out from v1 is the new dummy propagator 
	int pline_or_eline = rand() % 2;

	if (pline_or_eline == 0) {//p1 to be the dummy line
		if (p1->data.l.phys.screen == J) return 1; //if pline is J type, it could not be a dummy line, exit
		//if pline is connecting to a bubble, it cannot be a measuring line
		if (pair_ver(v1)->data.v.ein == pair_ver(v1)->data.v.eout || v1->data.v.ein == v1->data.v.eout) return 1;

		m_line->data.l.measure = REGULAR;        //set the old measuring line back to regular
		p1->data.l.measure = MEASURE;            //set the selected pline to be the new measuring line

		t->meas_line = p1;                       //change the property of the diagram
	}
	else{//eline1 to be the dummy line
		m_line->data.l.measure = REGULAR;
		eline1->data.l.measure = MEASURE;

		t->meas_line = eline1;
	}
	//probability is always 1

	return OK;

}
//===========================================================
int Dummy_e(struct diag *t){
	/*Only for dummy propagators.*/

	struct node *m_line;        //read out the measuring line pointer
	m_line = t->meas_line;

	struct node *v1, *eline1;
	v1 = rand_ver(&t->vertex);  //choose a vertex at random
	eline1 = v1->data.v.eout;//read out the eline come out from this vertex 

	if (eline1->data.l.measure == MEASURE) return 1;//if the chosen propagator happens to be a measuring line, reject

	//change the measuring line
	m_line->data.l.measure = REGULAR;
	eline1->data.l.measure = MEASURE;
	t->meas_line = eline1;

	return OK;
}
//===========================================================
//The parir "insert/remove" is to increase/decrease the diagram oder by inserting/removing 
// a ladder-type structure.
//-----------------------------------------------------------
int Insert(struct diag *t){
	/*In "insert", we make a random choice between S and T to start the construction from the special vertex V1.
	Next, we identify  vertex A as the destination vertex of the propagator originating from vertex V1, and vertex B
	as the originating vertex for the propagator with destination V2, which is the other worm end.*/

	if (t->n_wo != 2) return 1; //this update can be only applied in unphysical sector
	if (diag_order(t) >= Max_order) return 1; //when diagram has up-limit order, reject the updates which will increase order
	if (check_irreducibility(t) == NOTOK) return 1;

	struct node *v1, *v2, *a, *b;
	struct node *eline1, *e2;
	int worm_type;

	v1 = select_worm(t);//select a worm at random (S or T)
	worm_type = v1->data.v.type;
	eline1 = v1->data.v.eout;//eline coming from v1
	a = eline1->data.l.head; //vertex A, which is the head vertex of the eline(eline1) coming from V1

	//select another worm to be V2
	if (worm_type == S) v2 = t->t_worm;
	else v2 = t->s_worm;

	//if vertex A and the orther worm are the same, (two worms connected directly by an eline), reject proposal
	if (a == v2) return 1;

	//find vertex  B
	e2 = v2->data.v.ein;  //eline coming into V2
	b = e2->data.l.tail;  //originating vertex of e2

	//caluculate probability denominator, before the update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2);

	//apply update
	int spin_1, spin_2;
	struct node *new_eline1, *new_e2, *new_v1, *new_v2, *new_p;

	//read out the spin components
	spin_1 = eline1->data.l.phys.s;
	spin_2 = e2->data.l.phys.s;

	//create two new elines and two new vertices
	new_eline1 = create_lin(t, ELECTRON, REGULAR, spin_1);
	new_e2 = create_lin(t, ELECTRON, REGULAR, spin_2);
	new_v1 = create_ver(t);
	new_v2 = create_ver(t);

	//create a pline, screen depends on the position of two vertices
	if (v1->data.v.phys.r == v2->data.v.phys.r) new_p = create_lin(t, PHONON, REGULAR, W);
	else new_p = create_p_rand(t);

	//time components of two new vertices
	double t_new = BETA* rand_01(); //generate the new time

	//time of v1
	new_v1->data.v.phys.t = t_new;

	//time of v2,according to the screen type of new pline
	if (new_p->data.l.phys.screen == J) new_v2->data.v.phys.t = t_new; //J type connection, two vertices have same time
	else new_v2->data.v.phys.t = BETA* rand_01();

	//assign real space position to 2 new vertices
	new_v1->data.v.phys.r = v1->data.v.phys.r;
	new_v2->data.v.phys.r = v2->data.v.phys.r;

	//connect new lines and vertices
	connect_lin(new_p, new_v2, new_v1);

	connect_lin(new_eline1, new_v1, v1);
	connect_lin(new_e2, v2, new_v2);

	connect_lin(eline1, a, new_v1);
	connect_lin(e2, new_v2, b);

	//calculate probability numerator, after the update
	complex <double> ratio_nu = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2)*green_bold(&(t->prf_g),new_eline1)*green_bold(&(t->prf_g),new_eline1)*interaction(&(t->prf_w),new_p);

	//calculate the factor for probability, according to the screen type of new pline
	double factor = 1.0;
	//if screen of new pline is J, factor keeps 1, if is ~W, 1/T
	if (new_p->data.l.phys.screen == W) factor = 1.0 / TEM;

	//probability
	double ratio = abs(ratio_nu / ratio_de)*factor * 2 / TEM *u10 / u9;
	//printf("insert :%f\n", ratio);

	if (rand_01() > ratio) { //if probability is too small, remove all modification
		delete_lin(t, new_eline1);
		delete_lin(t, new_e2);

		delete_ver(t, new_v1);
		delete_ver(t, new_v2);

		delete_lin(t, new_p);

		connect_lin(eline1, a, v1);
		connect_lin(e2, v2, b);

		return 1;
	}

	//read out momentum on original elines
	struct momentum p1, p2;
	p1 = eline1->data.l.phys.p;
	p2 = e2->data.l.phys.p;

	//momentum on new lines
	new_p->data.l.phys.p = assign_ran_momentum();  //send random momenta to new_p
	new_eline1->data.l.phys.p = momentum_add(p1, new_p->data.l.phys.p);
	new_e2->data.l.phys.p = momentum_add(p2, new_p->data.l.phys.p);

	//change momentum excess
	t->delta = momentum_excess_ver(t->s_worm);

	//if one of eline1,e2 is the dummy line, the new dummy line has to have the same originating vertex
	//we only consider eline1, because we did not change the originating vertex of line e2
	if (eline1->data.l.measure == MEASURE){
		eline1->data.l.measure = REGULAR;
		new_eline1->data.l.measure = MEASURE;
		t->meas_line = new_eline1;
	}

	return OK;
}
//===========================================================
int Remove(struct diag *t){
	/*Select one of the special vertice V1 at random and verify that the topology of lines connecting
	it to the other special vertex V2, as well as lin parameters. */

	if (t->n_wo != 2) return 1; //this update can be only applied in unphysical sector
	if (check_irreducibility(t) == NOTOK) return 1;

	struct node *v1, *v2, *a, *b, *c, *d;
	struct node *eline1, *e2, *eline1_move, *e2_move, *p1;

	//select a worm at random (S or T)
	v1 = select_worm(t);

	eline1_move = v1->data.v.eout;//find the eline coming from v1 (eline1_move to removed)
	c = eline1_move->data.l.head; //find vertex C, which is the head vertex of the eline coming from V1
	eline1 = c->data.v.eout;      //find the eline coming from C
	a = eline1->data.l.head;      //find vertex A,  which is the head vertex of the eline coming from C
	p1 = c->data.v.pend;      //find p1, interaction line connecting vertex C an vertex D

	if (c == v1) return 1;
	if (c == a) return 1;

	//if vertex A is another worm, reject
	if (a == v1) return 1;
	if (a->data.v.type != NO) return 1;

	d = pair_ver(c);          //find vertex D, the pair vertex of C
	e2 = d->data.v.ein;       //find the eline coming into vertex D
	b = e2->data.l.tail;      //find vetrex B.
	e2_move = d->data.v.eout; //find the eline coming out from vertex D
	v2 = e2_move->data.l.head;//find vertex V2,  which is the head vertex of the eline(e2_move to removed) coming from C

	if (b == d) return 1;
	if (v2 == d) return 1;
	if (a == b) return 1;

	if (a == v2 || a == c || a == d) return 1;

	//Check
	//1.check whether vertex V2 is the other worm
	if (v1->data.v.type == S) { //if v1 is S , v2 has to be T
		if (v2 != t->t_worm) return 1;
	}
	else {                      //if v1 is T, v2 has to be S 
		if (v2 != t->s_worm) return 1;
	}

	//2.if one of the to-be-removed lines happens to be dummy line, reject update
	if (p1->data.l.measure == MEASURE || eline1_move->data.l.measure == MEASURE || e2_move->data.l.measure == MEASURE) return 1;

	//check spins
	int s1, s1_m, s2, s2_m;
	s1 = eline1->data.l.phys.s;
	s1_m = eline1_move->data.l.phys.s;
	s2 = e2->data.l.phys.s;
	s2_m = e2_move->data.l.phys.s;

	if (s1 != s1_m || s2 != s2_m) return 1;

	//caluculate probability denominator, before the update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2)*green_bold(&(t->prf_g),eline1_move)*green_bold(&(t->prf_g),e2_move)*interaction(&(t->prf_w),p1);

	//read out momentum on eline1 and p1
	/*struct momentum q, k1;
	q = p1->data.l.phys.p;
	k1 = eline1->data.l.phys.p;*/

	//printf("A:%d B:%d v1:%d v2:%d\n", a->data.id, b->data.id, v1->data.id, v2->data.id);
	//printf("C:%d D:%d \n", c->data.id, d->data.id);

	//caluculate probability numerator, after the update
	complex <double> ratio_nu = green(a->data.v.phys.t - v1->data.v.phys.t)
		*green(v2->data.v.phys.t - b->data.v.phys.t);

	//calculate the factor for probability, according to the screen type of new pline
	double factor = 1;
	//if screen of new pline is J, factor keeps 1, if is ~W, 1/T
	if (p1->data.l.phys.screen == W) factor = TEM;

	//printf("factor %f\n", factor);

	//probability
	double ratio = abs(ratio_nu / ratio_de)*factor * TEM / 2 * u9 / u10;
	//printf("remove %f\n", ratio);

	//printf("eline1_move:%d e2_move:%d eline1:%d e2:%d p1:%d\n", eline1_move->data.id, e2_move->data.id, eline1->data.id, e2->data.id, p1->data.id);
	if (ratio >= 1) ratio = 1;
	else{//if ratio smaller than 1, probability is ratio, need compare it to random number
		if (rand_01() > ratio) return 1; //do nothing
	}

	//delete 2 elines, 2 vertices and 1 pline
	delete_lin(t, eline1_move);   //delete eline1_move
	delete_lin(t, e2_move);   //delete e2_move
	delete_ver(t, c);         //delete vertex C
	delete_ver(t, d);         //delete vertex D
	delete_lin(t, p1);        //delete p1

	//connect all the lines
	connect_lin(eline1, a, v1);
	connect_lin(e2, v2, b);

	//printf("ein of V1:%d eout of V1:%d pline:%d\n", v1->data.v.ein->data.id, v1->data.v.eout->data.id, v1->data.v.pend->data.id);
	//change momentum excess
	t->delta = momentum_excess_ver(t->s_worm);


	return OK;
}
//===========================================================
int Dress(struct diag *t){
	/*The dress update starts from random selection of vertex A and identification of vertex B and C, linked
	to it by propagator lines, if B=A, the update is rejected(we donot care one of the vertices is worm).The
	proposal is to add new vertices d (between B and A)and e(between C and A), and link them with W interaction
	line with momentum q. The new time variables are generated from the probability density.*/
	if (diag_order(t) >= Max_order) return 1; //when diagram has up-limit order, reject the updates which will increase order
    if (check_irreducibility(t) == NOTOK) return 1;
	
	struct node *a, *b, *c, *eline1, *e2;
	a = rand_ver(&t->vertex); //select a vertex at random

	//two elines connecting vertex A
	eline1 = a->data.v.ein;  //ein
	e2 = a->data.v.eout; //eout

	b = eline1->data.l.tail; //Vertex B, tail ver of eline1
	c = e2->data.l.head; //Vertex C, head ver of e2

	if (b == a) return 1;//if B=A(C=A), vertex A forms a closed loop, reject update

	//diagram order before the update
	int n = diag_order(t);

	//caluculate probability denominator, before the update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2);

	//propose update
	struct node *d, *e, *eline1_new, *e2_new, *p_new;
	int spin_1, spin_2;

	//read out the spin components
	spin_1 = eline1->data.l.phys.s; //spin on eline1
	spin_2 = e2->data.l.phys.s; //spin on e2

	//create two new elines and two new vertices
	eline1_new = create_lin(t, ELECTRON, REGULAR, spin_1);
	e2_new = create_lin(t, ELECTRON, REGULAR, spin_2);
	d = create_ver(t);
	e = create_ver(t);
	p_new = create_lin(t, PHONON, REGULAR, W);

	//generate the new time on vertex D randomly
	d->data.v.phys.t = BETA* rand_01();

	//generate the new time on new vertex E, and it's different from that on D
	e->data.v.phys.t = BETA* rand_01();

	//connection, 
	connect_lin(eline1, d, b);
	connect_lin(eline1_new, a, d);
	connect_lin(e2, c, e);
	connect_lin(e2_new, e, a);
	connect_lin(p_new, e, d);

	//assign real space position to 2 new vertices
	d->data.v.phys.r = e->data.v.phys.r = a->data.v.phys.r;

	//caluculate probability numerator, after the update
	complex <double> ratio_nu = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2)*green_bold(&(t->prf_g),eline1_new)*green_bold(&(t->prf_g),e2_new)*interaction(&(t->prf_w),p_new);

	/*//count the number of vertex can be undress
	//int n_undress = vertex_undress(t);*/

	//calculate the factor for probability, according to the screen type of new pline
	double factor = 1.0 *n / (n + 1) / TEM / TEM; 


	//probability
	double ratio = abs(ratio_nu / ratio_de)*factor*u12 / u11;
	//printf("Probability %f\n", ratio);

	//if ratio equal or larger than 1, probability is 1, don't need to compare to random number
	if (ratio >= 1) ratio = 1;
	else{//if ratio smaller than 1, probability is ratio, need compare it to random number
		if (rand_01() > ratio) { //if probability is too small, remove all modification
			delete_lin(t, eline1_new);
			delete_lin(t, e2_new);

			delete_ver(t, d);
			delete_ver(t, e);

			delete_lin(t, p_new);

			connect_lin(eline1, a, b);
			connect_lin(e2, c, a);

			return 1;
		}
	}

	//change momentum 
	p_new->data.l.phys.p = assign_ran_momentum();
	eline1_new->data.l.phys.p = momentum_subtract(eline1->data.l.phys.p, p_new->data.l.phys.p);
	e2_new->data.l.phys.p = momentum_subtract(e2->data.l.phys.p, p_new->data.l.phys.p);

	return OK;
}
//===========================================================
int Undress(struct diag *t){
	/*In "undress", we select vertex A at random, identify vertices D, B ,C and E, using links along the
	propagator lines, and verify that the topology of lines and their parameters are consistent with the
	dressed vertex configuration. If the D-E line is not of the diagonal W type or one of the propagators
	D-A and E-C is a dummy line, the update is rejected.*/
	if (check_irreducibility(t) == NOTOK) return 1;

	struct node *a, *b, *c, *eline1, *e2;
	struct node *d, *e, *eline1_move, *e2_move, *p_move;

	a = rand_ver(&t->vertex); //select a vertex at random

	eline1_move = a->data.v.ein;  //find two lines to be removed
	e2_move = a->data.v.eout;
	d = eline1_move->data.l.tail; //find two vertices D and E
	e = e2_move->data.l.head;


	if (d == a || e == a) return 1; //if selected vertex is on a bubble, reject
	if (d == e) return 1; //if d,e,a are on a bubble, reject

	if (d->data.v.pend != e->data.v.pend) return 1; //if vertex D and E are not directly connected by a pline, reject

	//read out the pline, which connects D and E
	p_move = d->data.v.pend;

	if (p_move->data.l.phys.screen != W) return 1;  //if the pline is not a W type, reject

	//find vertices B and C
	eline1 = d->data.v.ein;
	e2 = e->data.v.eout;
	b = eline1->data.l.tail;
	c = e2->data.l.head;

	//if (b == c) return 1;
	if (b == a || c == a) return 1;//if B or C is same as A, reject

	if (b == d) return 1;
	if (c == d) return 1;

	//printf("ID of vertex D:%d and vertex E:%d\n", d->data.id, e->data.id);
	//printf("ID of vertex B:%d and vertex C:%d vertex A:%d\n", b->data.id, c->data.id, a->data.id);


	//2.if one of the to-be-removed lines happens to be dummy line, reject update
	if (p_move->data.l.measure == MEASURE || eline1_move->data.l.measure == MEASURE
		|| e2_move->data.l.measure == MEASURE) return 1;

	//check spins, if not diagonal, reject
	int s1, s1_m, s2, s2_m;
	s1 = eline1->data.l.phys.s;
	s1_m = eline1_move->data.l.phys.s;
	s2 = e2->data.l.phys.s;
	s2_m = e2_move->data.l.phys.s;

	if (s1 != s1_m || s2 != s2_m) return 1;

	//diagram order
	int n = diag_order(t);

	/*//count number of vertex can be undressed
	//int n_undress = vertex_undress(t);//*/

	//caluculate probability denominator, before the update
	complex <double> ratio_de = green_bold(&(t->prf_g),eline1)*green_bold(&(t->prf_g),e2)*green_bold(&(t->prf_g),eline1_move)*green_bold(&(t->prf_g),e2_move)*interaction(&(t->prf_w),p_move);


	//caluculate probability numerator, after the update
	complex <double> ratio_nu = green(a->data.v.phys.t - b->data.v.phys.t)
		*green(c->data.v.phys.t - a->data.v.phys.t);

	//calculate the factor for probability, according to the screen type of new pline
	double factor = 1.0 * n / (n - 1) *TEM *TEM;

	//probability
	double ratio = abs(ratio_nu / ratio_de)*factor*u11 / u12;
	//printf(" with Probability %f\n", ratio);

	//if ratio equal or larger than 1, probability is 1, don't need to compare to random number
	if (ratio >= 1) ratio = 1;
	else{//if ratio smaller than 1, probability is ratio, need compare it to random number
		if (rand_01() > ratio) return 1;
	}

	//do update
	//delete 2 elines, 2 vertices and 1 pline
	delete_lin(t, eline1_move);   //delete eline1_move
	delete_lin(t, e2_move);   //delete e2_move
	delete_ver(t, d);         //delete vertex D
	delete_ver(t, e);         //delete vertex E
	delete_lin(t, p_move);    //delete p_move


	//connect lines
	connect_lin(eline1, a, b);
	connect_lin(e2, c, a);

	return OK;
}
//===========================================================
int Recolor(struct diag *t){
	/*A easy and efficient way to CHANGE the SPIN index of propagator lines is to select a random vertex and
	use it to construct a closed lopp by following the propagator lines attached to it. If all propagators in
	the loop have the same spinindex "alpha" it can be changed to "-alpha" with acceptance ratio 1 in the
	absence of external magnetic field.*/
	if (check_irreducibility(t) == NOTOK) return 1;

	struct node *v1, *efirst;
	v1 = rand_ver(&t->vertex); //select a vertex V1 at random
	//printf("v1 id :%d\n", v1->data.id);

	struct node *eline1, *v2;

	efirst = eline1 = v1->data.v.eout;    //eline coming from selected vertex V1
	//if (eline1->data.l.measure == MEASURE) return 1; //if eline is a dummy line, no closed loop

	int spin_0, spin_1;
	spin_0 = eline1->data.l.phys.s;
	//printf("eline1 id :%d\n", eline1->data.id);

	//find the loop originating from vertex V1, and check whether the spin components in the loop are the same
	//and check if the the there is a dummy line in the loop
	do{
		v2 = eline1->data.l.head;
		//printf("v2 id :%d\n", v2->data.id);
		eline1 = v2->data.v.eout;
		spin_1 = eline1->data.l.phys.s;
		if (spin_1 != spin_0) return 1;  //check spin and dummy line
	} while (v2 != v1);

	//change the spin in the loop to  the opposite of spin_0
	if (spin_0 == UP) efirst->data.l.phys.s = DOWN;
	else efirst->data.l.phys.s = UP;


	eline1 = efirst;
	//change the other spins
	do{
		v2 = eline1->data.l.head;
		eline1 = v2->data.v.eout;
		if (spin_0 == DOWN) eline1->data.l.phys.s = UP;
		else eline1->data.l.phys.s = DOWN;
	} while (v2 != v1);

	//in the absence of external potential, accepatance ratio is 1
	double ratio = 1;

	return OK;
}
//===========================================================
int Move_T(struct diag *t){
	/*This seld-complementary update is designed to sample time variables of the diagram without changing its order
	and topology. The proposal is to select one of the vertices at random: let it be vertex A and update its imaginary
	time variable from t to t', using probability density distribution. The interaction line attached to A can not be
	of the J type, whether it is the physical or unphysical F line does not matter.*/
	
	if (check_irreducibility(t) == NOTOK) return 1;

	struct node *a, *p;
	a = rand_ver(&t->vertex); //select a vertex V1(A) at random

	p = a->data.v.pend;       //find the interaction line connecting to selected vertex
	if (p->data.l.phys.screen == J) return 1;  //if it is J type, reject update

	//identifier of whether the vertex A is on a bubble, 0:general and  1:bubble
	int general_or_bubble = 0;

	//check whether the vertex A is on a bubble

	//read out the eline coming in and out from vertex A
	struct node *ein, *eout;
	ein = a->data.v.ein;
	eout = a->data.v.eout;

	if (ein == eout) general_or_bubble = 1;//if vertex A on a bubble, identifier is 1

	//read out the original time on Vertex A
	double time = a->data.v.phys.t;

	//calculate probability
	double ratio = 1.0;
	double factor = 1.0* TEM / TEM;

	struct node *d;
	d = pair_ver(a);//the pair of A, denoted as vertex D
	
	
	//calculate the main part of ratio, depending on whether A is on a bubble or not
	if (general_or_bubble == 1 && d->data.v.type != NO)  {//calculate probability when A on a bubble
		
		//before update
		complex <double> ratio_de = f_worm(&(t->prf_f),p);
		
		//change the time component of selected vertex
		a->data.v.phys.t = rand_time();
		complex <double> ratio_nu = f_worm(&(t->prf_f),p);

		ratio = factor*abs(ratio_nu/ratio_de);
		//printf("move-Tratio: %g\n",ratio);

		if (ratio >= 1) ratio = 1;
		else{
			
			if (rand_01() > ratio) {
				//undo the update
				a->data.v.phys.t = time;
				return 1; //probability is too small, 
			}
		}
		
	}
	else{//when A is not on a bubble


		//caluculate probability denominator1(eline part), before the update
		complex <double> ratio_deline1 = green_bold(&(t->prf_g),ein)*green_bold(&(t->prf_g),eout);

		//caluculate probability denominator2(pline part), before the update
		complex <double> ratio_de2 = interaction(&(t->prf_w),p);

		//change the time component on the selected vertex
		//if a and d happen to be nearest site, make sure that new time isn't same as that on Vertex D
		if (abs(a->data.v.phys.r - d->data.v.phys.r) != 0){
			do{
				a->data.v.phys.t = rand_time();
			} while (fabs(a->data.v.phys.t - d->data.v.phys.t) < 1E-4);
		}

		else a->data.v.phys.t = rand_time();

		//caluculate probability numerator1(eline part), before the update
		complex <double> ratio_nu1 = green_bold(&(t->prf_g),ein)*green_bold(&(t->prf_g),eout);

		//caluculate probability numerator2(pline part), before the update
		complex <double> ratio_nu2 = interaction(&(t->prf_w),p);

		//probability
		ratio = abs(ratio_nu1 / ratio_deline1)*abs(ratio_nu2 / ratio_de2)*factor;

		if (ratio >= 1) ratio = 1;
		else{
			if (rand_01() > ratio) {//if ratio is too small, remove the update

				a->data.v.phys.t = time;//set the time component back to original value
				return 1;
			}
		}

	}
	return OK;
}
//===========================================================
