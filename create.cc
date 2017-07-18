/*Functions to create/delete lines in diagrams*/

#define MAXID 100000000 //Limit of ID number



double rand_time();

int idiag(struct diag *t);

int create_tadpole_ladder_nonphys(struct diag *t);
int create_tadpole_ladder(struct diag *t);
int create_tadpole_ladder_phy(struct diag *t);
int create_bubble_phy(struct diag *t);
int create_tadpole_ladder_ir(struct diag *t);
int create_tadpole_ladder_phy_s1(struct diag *t);
int startpoint(struct diag *t);
int firstdiag(struct diag *t);

struct node *create_lin(struct diag *t, int type, int measure, int spin);
struct node *create_ver(struct diag *t);
struct node *create_p_rand(struct diag *t);
struct node *create_e_rand(struct diag *t);

int delete_lin(struct diag *t, struct node *e);
int delete_ver(struct diag *t, struct node *v);

int connect_lin(struct node *line, struct node *head_ver, struct node *tail_ver);

void prof_xerror(const char* buff);
struct profile_multi prof_create_multi(int dim, int nbin[],double min[],double max[]);
complex <double> assign_complx (double re, double im);


//====================================================
float rand_01() {//return random number between 0 to 1
	return rand() / (float)RAND_MAX;
}
//==========================================================
double rand_time(){
	/*This function returns a double between 0 and BETA. Uniform distribution of time interval*/
	double rtn = BETA*rand_01();
	//rtn = BETA*(2 * rand_01() - 1);
	return rtn;
}
//==========================================================

//==========================================================
int idiag(struct diag *t){
	/*This function initializes an empty diagram */

	//Initializes the number of elements
	t->n_el = 0;
	t->n_pl = 0;
	t->n_ver = 0;
	t->n_ml = 0;
	t->n_wo = 0;

	//initialize ids of lines and vertices
	t->id_lin = 0;
	t->id_ver = 0;

	//Initializes the identifers of worms and measuring line, these special line and vertices
	t->meas_line = NULL;

	t->s_worm = NULL;
	t->t_worm = NULL;

	initialize_list(&t->line);   //Initializes the lines, no node in the list.
	initialize_list(&t->vertex); //Initializes the vertices

	//Initializes the excess momentum on I worm
	t->delta = momentum_init(0.0, 0.0, 0.0);

	t->parity = NOTDEFINE;
	
	//create two profiles to store the information of bold propagator and the screening interaction
	int nft[1] = {N_t};
	double min[1] = {Tmin};
	double max[1] = {BETA};
	
	t->prf_g = prof_create_multi(1, nft, min, max);//profile which saves the data of green's function G(t)
	
	int nft_w[2] = {N,N_t};
	double min_w[2] = {XL,Tmin};
	double max_w[2] = {XU,BETA};
	
	t->prf_w = prof_create_multi(2, nft_w, min_w, max_w);//profile which saves the profile of ~W(r,t)
	t->prf_f = prof_create_multi(2, nft_w, min_w, max_w);//profile which saves the profile of F(r,t)

	return 1;
}

//Create a tadpole with momentum conserved
//==========================================================
int create_tadpole_ladder_nonphys(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *e1, *e2, *e3, *e4, *p1, *p2, *v1, *v2, *v3, *v4;

	e1 = create_lin(t, ELECTRON, MEASURE, UP);    //create a eline e1, make it a measuring line
	e2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line
	e3 = create_lin(t, ELECTRON, REGULAR, DOWN);  //create a eline e3, make it a regular line
	e4 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e4, make it a regular line
	p1 = create_lin(t, PHONON, REGULAR, W);       //create a pline p1, make it a regular line
	p2 = create_lin(t, PHONON, REGULAR, W);       //create a pline p2, make it a regular line

	v1 = create_ver(t);                           //create 4 vertices
	v2 = create_ver(t);
	v3 = create_ver(t);
	v4 = create_ver(t);

	//connect all lines and vertices
	connect_lin(e1, v1, v3);
	connect_lin(e2, v4, v2);
	connect_lin(e3, v3, v1);
	connect_lin(e4, v2, v4);
	connect_lin(p1, v1, v2);
	connect_lin(p2, v4, v3);

	//send momentum to every line
	e1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	e2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	p1->data.l.phys.p = assign_ran_momentum();  //send random momenta to p1

	//from momenta conservation, assign momentum to the other 3 lines
	p2->data.l.phys.p = p1->data.l.phys.p;
	e3->data.l.phys.p = assign_ran_momentum();
	e4->data.l.phys.p = momentum_add(p1->data.l.phys.p, e2->data.l.phys.p);


	//assign random time to two vertices
	v1->data.v.phys.t = 0.2;//rand_time();
	v2->data.v.phys.t = 0.6;
	v3->data.v.phys.t = 0.1;//rand_time();
	v4->data.v.phys.t = 0.4;

	//assign space position to vertices
	v1->data.v.phys.r = v3->data.v.phys.r = 0;
	v2->data.v.phys.r = v4->data.v.phys.r = 1;

	v1->data.v.type = S;
	v3->data.v.type = T;

	t->n_wo = 2;
	t->s_worm = v1;
	t->t_worm = v3;

	t->delta = momentum_excess_ver(v1);


	return 0;
};
//==========================================================
int create_tadpole_ladder_phy_s(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *e1, *e2, *e3, *e4, *p1, *p2, *v1, *v2, *v3, *v4;

	//all propagators have same spin index
	e1 = create_lin(t, ELECTRON, MEASURE, UP);    //create a eline e1, make it a measuring line
	e2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line
	e3 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e3, make it a regular line
	e4 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e3, make it a regular line
	p1 = create_lin(t, PHONON, REGULAR, J);       //create a pline p1, make it a regualr line
	p2 = create_lin(t, PHONON, REGULAR, J);       //create a pline p2, make it a regular line
	v1 = create_ver(t);                           //create 4 vertices
	v2 = create_ver(t);
	v3 = create_ver(t);
	v4 = create_ver(t);

	//connect all lines and vertices
	connect_lin(e1, v1, v3);
	connect_lin(e2, v4, v2);
	connect_lin(e3, v3, v1);
	connect_lin(e4, v2, v4);
	connect_lin(p1, v1, v2);
	connect_lin(p2, v4, v3);

	//send momentum to every line
	e1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	e2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	p1->data.l.phys.p = assign_ran_momentum();  //send random momenta to p1

	//from momenta conservation, assign momentum to the other 3 lines
	p2->data.l.phys.p = p1->data.l.phys.p;
	e3->data.l.phys.p = momentum_add(e1->data.l.phys.p, p1->data.l.phys.p);
	e4->data.l.phys.p = momentum_add(e2->data.l.phys.p, p1->data.l.phys.p);

	//assign random time to two vertices
	v1->data.v.phys.t = 0.2;//rand_time();
	v2->data.v.phys.t = 0.2;
	v3->data.v.phys.t = 0.3;//rand_time();
	v4->data.v.phys.t = 0.3;

	//assign space position to vertices
	v1->data.v.phys.r = v3->data.v.phys.r = 0;
	v2->data.v.phys.r = v4->data.v.phys.r = 1;

	return 0;
};
//==========================================================
int create_tadpole_ladder_phy_t(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *e1, *e2, *e3, *e4, *p1, *p2, *v1, *v2, *v3, *v4;

	//all propagators have same spin index
	e1 = create_lin(t, ELECTRON, MEASURE, UP);    //create a eline e1, make it a measuring line
	e2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line
	e3 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e3, make it a regular line
	e4 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e3, make it a regular line

	//one pline has W screen type
	p1 = create_lin(t, PHONON, REGULAR, J);       //create a pline p1, make it a regualr line
	p2 = create_lin(t, PHONON, REGULAR, J);       //create a pline p2, make it a regular line
	v1 = create_ver(t);                           //create 4 vertices
	v2 = create_ver(t);
	v3 = create_ver(t);
	v4 = create_ver(t);

	//connect all lines and vertices
	connect_lin(e1, v1, v3);
	connect_lin(e2, v4, v2);
	connect_lin(e3, v3, v1);
	connect_lin(e4, v2, v4);
	connect_lin(p1, v1, v2);
	connect_lin(p2, v4, v3);

	//send momentum to every line
	e1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	e2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	p1->data.l.phys.p = assign_ran_momentum();  //send random momenta to p1

	//from momenta conservation, assign momentum to the other 3 lines
	p2->data.l.phys.p = p1->data.l.phys.p;
	e3->data.l.phys.p = momentum_add(e1->data.l.phys.p, p1->data.l.phys.p);
	e4->data.l.phys.p = momentum_add(e2->data.l.phys.p, p1->data.l.phys.p);

	//assign random time to two vertices
	v1->data.v.phys.t = rand_time();
	v2->data.v.phys.t = rand_time();
	v3->data.v.phys.t = v1->data.v.phys.t;
	v4->data.v.phys.t = v2->data.v.phys.t;

	//assign space position to vertices
	v1->data.v.phys.r = v3->data.v.phys.r = 0;
	v2->data.v.phys.r = v4->data.v.phys.r = 1;

	return 0;
};
//==========================================================
int create_tadpole_ladder_ir(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *e1, *e2, *e3, *e4, *p1, *p2, *v1, *v2, *v3, *v4;

	e1 = create_lin(t, ELECTRON, MEASURE, UP);    //create a eline e1, make it a measuring line
	e2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line
	e3 = create_lin(t, ELECTRON, REGULAR, UP);  //create a eline e3, make it a regular line
	e4 = create_lin(t, ELECTRON, REGULAR, UP);  //create a eline e3, make it a regular line
	p1 = create_lin(t, PHONON, REGULAR, J);       //create a pline p1, make it a regualr line
	p2 = create_lin(t, PHONON, REGULAR, J);       //create a pline p2, make it a regular line
	v1 = create_ver(t);                           //create 4 vertices
	v2 = create_ver(t);
	v3 = create_ver(t);
	v4 = create_ver(t);

	//connect all lines and vertices
	connect_lin(e1, v1, v3);
	connect_lin(e2, v4, v2);
	connect_lin(e3, v3, v1);
	connect_lin(e4, v2, v4);
	connect_lin(p1, v1, v2);
	connect_lin(p2, v3, v4);

	//send momentum to every line
	e1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	e2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	p1->data.l.phys.p = assign_ran_momentum();  //send random momenta to p1

	//from momenta conservation, assign momentum to the other 3 lines
	p2->data.l.phys.p = momentum_factor(-1.0, p1->data.l.phys.p);
	e3->data.l.phys.p = momentum_add(e1->data.l.phys.p, p1->data.l.phys.p);
	e4->data.l.phys.p = momentum_add(e2->data.l.phys.p, p1->data.l.phys.p);

	//assign random time to two vertices
	v1->data.v.phys.t = rand_time();
	v2->data.v.phys.t = v1->data.v.phys.t;
	v3->data.v.phys.t = rand_time();
	v4->data.v.phys.t = v3->data.v.phys.t;

	//assign space position to vertices
	v1->data.v.phys.r = v3->data.v.phys.r = 0;
	v2->data.v.phys.r = v4->data.v.phys.r = 1;

	return 0;
};
//==========================================================
int startpoint(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *e1, *e2, *p1, *v1, *v2;

	e1 = create_lin(t, ELECTRON, MEASURE, UP);    //create a eline e1, make it a measuring line
	e2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line
	p1 = create_lin(t, PHONON, REGULAR, J);       //create a pline p1, make it a regualr line

	v1 = create_ver(t);                           //create 2 vertices
	v2 = create_ver(t);

	connect_lin(e1, v1, v1);
	connect_lin(e2, v2, v2);

	connect_lin(p1, v1, v2);

	//send momentum to every line
	e1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	e2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	p1->data.l.phys.p = momentum_init(0, 0, 0);

	//assign random time to two vertices
	v1->data.v.phys.t = rand_time();
	v2->data.v.phys.t = v1->data.v.phys.t;
	//v2->data.v.phys.t = rand_time();
	//printf("t1:%f t2:%f\n",v1->data.v.phys.t,v2->data.v.phys.t);

	//assign space position to vertices
	v1->data.v.phys.r = 0;
	v2->data.v.phys.r = 1;

	t->n_wo = 0;
	t->s_worm = NULL;
	t->t_worm = NULL;

	t->delta.x = 0.0;
	t->delta.y = 0.0;
	t->delta.z = 0.0;

	return 0;
};
//==========================================================
int firstdiag(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *eline1, *eline2, *p1, *v1, *v2;

	eline1 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e1, make it a measuring line
	eline2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line
	p1 = create_lin(t, PHONON, MEASURE, W);       //create a pline p1, make it a regualr line

	v1 = create_ver(t);                           //create 2 vertices
	v2 = create_ver(t);

	connect_lin(eline1, v2, v1);
	connect_lin(eline2, v1, v2);

	connect_lin(p1, v1, v2);

	//send momentum to every line
	eline1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	eline2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	
	p1->data.l.phys.p = momentum_subtract(eline1->data.l.phys.p, eline2->data.l.phys.p) ; //momentum_init(0, 0, 0);

	//assign random time to two vertices
	v1->data.v.phys.t = rand_time();
	//v2->data.v.phys.t = v1->data.v.phys.t;
	v2->data.v.phys.t = rand_time();
	//printf("t1:%f t2:%f\n",v1->data.v.phys.t,v2->data.v.phys.t);

	//assign space position to vertices
	v1->data.v.phys.r = 0;
	v2->data.v.phys.r = 0;
	
	t->n_wo = 0;
	t->s_worm = NULL;
	t->t_worm = NULL;

	t->delta.x = 0.0;
	t->delta.y = 0.0;
	t->delta.z = 0.0;

	/*t->n_wo = 2;
	t->s_worm = v1;
	t->t_worm = v2;
	
	v1->data.v.type = S;
	v2->data.v.type = T;
	

	//t->delta.x = 0.0;
	//t->delta.y = 0.0;
	//t->delta.z = 0.0;
	t->delta = p1->data.l.phys.p;*/

	return 0;
};
//==========================================================
int create_bubble_phy(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *e1, *e2, *p1, *v1, *v2;

	e1 = create_lin(t, ELECTRON, MEASURE, UP);    //create a eline e1, make it a measuring line
	e2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line

	p1 = create_lin(t, PHONON, REGULAR, W);       //create a pline p1, make it a regualr line

	v1 = create_ver(t);                           //create 2 vertices
	v2 = create_ver(t);


	//connect all lines and vertices
	connect_lin(e1, v1, v1);
	connect_lin(e2, v2, v2);

	connect_lin(p1, v1, v2);


	//send momentum to every line
	e1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	e2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	p1->data.l.phys.p = assign_ran_momentum();  //send random momenta to p1

	//assign random time to two vertices
	v1->data.v.phys.t = rand_time();
	v2->data.v.phys.t = rand_time();

	//assign space position to vertices
	v1->data.v.phys.r = v2->data.v.phys.r = 0;
	
	t->n_wo = 0;
	t->s_worm = NULL;
	t->t_worm = NULL;

	t->delta.x = 0.0;
	t->delta.y = 0.0;
	t->delta.z = 0.0;

	return 0;
};
//==========================================================
int create_tadpole_ladder_phy_s1(struct diag *t){
	/*This function takes an empty diag and make it into a double tadpole diagram*/
	struct node *e1, *e2, *e3, *e4, *p1, *p2, *v1, *v2, *v3, *v4;

	//all propagators have same spin index
	e1 = create_lin(t, ELECTRON, MEASURE, UP);    //create a eline e1, make it a measuring line
	e2 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e2, make it a regular line
	e3 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e3, make it a regular line
	e4 = create_lin(t, ELECTRON, REGULAR, UP);    //create a eline e3, make it a regular line
	p1 = create_lin(t, PHONON, REGULAR, W);       //create a pline p1, make it a regualr line
	p2 = create_lin(t, PHONON, REGULAR, W);       //create a pline p2, make it a regular line
	v1 = create_ver(t);                           //create 4 vertices
	v2 = create_ver(t);
	v3 = create_ver(t);
	v4 = create_ver(t);

	//connect all lines and vertices
	connect_lin(e1, v1, v3);
	connect_lin(e2, v4, v2);
	connect_lin(e3, v3, v1);
	connect_lin(e4, v2, v4);
	connect_lin(p1, v1, v2);
	connect_lin(p2, v4, v3);

	//send momentum to every line
	e1->data.l.phys.p = assign_ran_momentum();  //send random momenta to e1
	e2->data.l.phys.p = assign_ran_momentum();  //send random momenta to e2
	p1->data.l.phys.p = assign_ran_momentum();  //send random momenta to p1

	//from momenta conservation, assign momentum to the other 3 lines
	p2->data.l.phys.p = p1->data.l.phys.p;
	e3->data.l.phys.p = momentum_add(e1->data.l.phys.p, p1->data.l.phys.p);
	e4->data.l.phys.p = momentum_add(e2->data.l.phys.p, p1->data.l.phys.p);

	//assign random time to two vertices
	v1->data.v.phys.t = 0.2;//rand_time();
	v2->data.v.phys.t = 0.2;
	v3->data.v.phys.t = 0.3;//rand_time();
	v4->data.v.phys.t = 0.3;

	//assign space position to vertices
	v1->data.v.phys.r = v3->data.v.phys.r = 0;
	v2->data.v.phys.r = v4->data.v.phys.r = 1;

	return 0;
};
//==========================================================

//Create lines and vertices
//==========================================================
struct node *create_lin(struct diag *t, int type, int measure, int spin_or_screen){
	/*
	This function creates an line in the diagram provided as an argument. measure indicates if it's a measuring line, up_down is spin.
	It returns the new node index, where the line was created.
	*/
	struct node *rtn;  //declare a node structure 
	//rtn = (struct node *)malloc(sizeof(struct node)); //call malloc to give a address to the pionter
	//this is wrong, because in the makenode, we already assign a meomery to rtn

	makenode(&t->line); //create a new node on line list, which the current pointer is pointing to right now, after creating this new node
	rtn = t->line.current;//node pointer= current pointer in the list

	//assign id to line
	rtn->data.id = t->id_lin + 1;
	t->id_lin++;

	if (t->id_lin > MAXID) { //reset the ID numbers

		t->id_lin = 0;
		t->id_ver = 0;
		//for lines list
		reset_current(&t->line);
		do{
			t->line.current->data.id = t->id_lin++;
			t->id_lin = t->line.current->data.id;
		} while (next_one(&t->line) == OK);


		//for vertex list
		reset_current(&t->vertex);
		do{
			t->vertex.current->data.id = t->id_ver++;
			t->id_ver = t->vertex.current->data.id;
		} while (next_one(&t->vertex) == OK);

	};

	//assign print to line
	rtn->data.print = 0;

	//assign line type to line
	if (type == ELECTRON){ //if the new line is electron
		rtn->data.l.type = ELECTRON;
		t->n_el++;
		//spin index of electron
		if (spin_or_screen == UP) rtn->data.l.phys.s = UP;
		else rtn->data.l.phys.s = DOWN;
		//screen part bad, an electron does not have screen
		//rtn->data.l.phys.screen = BAD;/*NULL*/
	}
	else{ //then this line is a interaction line
		rtn->data.l.type = PHONON;
		t->n_pl++;

		//spin part bad, an interaction does not have spin
		rtn->data.l.phys.s = NONE;

		//interaction type 
		switch (spin_or_screen){
		case J:  rtn->data.l.phys.screen = J;
			break;
		case W:  rtn->data.l.phys.screen = W;
			break;
		/*case F:  rtn->data.l.phys.screen = F;
			break;*/
		}
	}

	//assign the measure or not to line , both interaction line and propagator
	if (measure == MEASURE) {
		rtn->data.l.measure = MEASURE;
		t->n_ml++;           //increments the number of measuring line
		t->meas_line = rtn;  //pointer of measuring line
	}
	else rtn->data.l.measure = REGULAR;

	return rtn;
}
//==========================================================
struct node *create_ver(struct diag *t){
	/*
	This function creates an line in the diagram provided as an argument. measure indicates if it's a measuring line, up_down is spin.
	It returns the new node pointer, where the line was created.
	*/
	struct node *rtn;  //declare a node structure 
	//rtn = (struct node *)malloc(sizeof(struct node)); //call malloc to give a address to the pionter


	makenode(&t->vertex);   //create a new node on line list, which the current pointer is pointing to it
	rtn = t->vertex.current;//node pointer= current pointer in the list
	t->vertex.current->data.v.type = NO; //create a normal vertex

	//assign id to vertex
	rtn->data.id = t->id_ver + 1;
	t->id_ver++;

	t->n_ver++;

	//assign print to vertex
	rtn->data.print = 0;

	return rtn;
}
//==========================================================
struct node *create_p_rand(struct diag *t){
	/*This function create a regular interaction line, whether it is "J or W "type, is decided at random*/
	struct node *rtn;  //declare a node structure 
	int ran = 0;       //declare a integer
	ran = rand() % 2;  //generate a number between 0 and 1 at random

	if (ran == 0) rtn = create_lin(t, PHONON, REGULAR, J);  //J type interaction
	else rtn = create_lin(t, PHONON, REGULAR, W);           //W type interaction

	return rtn;
}
//==========================================================
struct node *create_e_rand(struct diag *t){
	/*This function create a regular propagator, the spin type "Up or Down", is decided at random*/
	struct node *rtn;  //declare a node structure 
	int ran = 0;       //declare a integer
	ran = rand() % 2;  //generate a number between 0 and 1 at random

	if (0 == ran) rtn = create_lin(t, ELECTRON, REGULAR, UP);  //spin Up
	else  rtn = create_lin(t, ELECTRON, REGULAR, DOWN);        //spin Down

	return rtn;
}
//==========================================================
//Delete lines and vertices
//==========================================================
int delete_lin(struct diag *t, struct node *e){
	/*This function delete a line node from the line list in the diagram structure, according to its node id.*/

	t->line.current = e;
	if (t->line.current->data.l.type == ELECTRON) t->n_el--; //decrease number of eline
	else t->n_pl--;                                          //decrease number of pline
	delete_current(&t->line); //delete the current node


	return OK;

}
//==========================================================
int delete_ver(struct diag *t, struct node *v){
	/*This function delete a vertex node from the vertex list in the diagram structure,
	according to its  vertex node pointer.*/
	//int id = v->data.id;         //read out the node id
	//delete_node(&t->vertex, id); //delete a vertex
	t->vertex.current = v;
	delete_current(&t->vertex); //delete the current nod
	t->n_ver--;

	return OK;
}
//==========================================================

/*Function to connect lines*/

int connect_lin(struct node *line, struct node *head_ver, struct node *tail_ver) {
	/*This function connect the lines in line list to intended vetices, to make sure the connection is consistant,
	we connect lines from both sides*/

	//Connects the line to the vertives
	line->data.l.head = head_ver;
	line->data.l.tail = tail_ver;

	//Connects the vertices to line
	if (line->data.l.type == ELECTRON){
		//Connects the electron end of vertices to the eline
		head_ver->data.v.ein = line;
		tail_ver->data.v.eout = line;
	}
	else{
		//Connects the pend of vertices to the pline
		head_ver->data.v.pend = line;
		tail_ver->data.v.pend = line;
		head_ver->data.v.p_in_out = IN; //this vertex has this pline coming into it
		tail_ver->data.v.p_in_out = OUT;//this vertex has this pline coming out off it
	}

	return 0;
}
//==========================================================

/*Create a multiple_dimension profile*/
//==========================================================
void prof_xerror(const char* buff) {
  /*Prints an error message and quits*/
  cerr << buff <<"\n";
  cerr << "Profile function error report\n";
  cerr << "Terminating now\n";
  exit(0);
}
//==============================================================
//==========================================================
struct profile_multi prof_create_multi(int dim, int nbin[],double min[],double max[]) {
  /*This function creates and initializes a profile structure 
    of dimension dim with n[] bins in each direction and parameter in 
    the range from min[] to max[] in each of the dim directions. */
  struct profile_multi rtn;
  rtn.dim=dim; //Sets the dimension
  rtn.ntot=1; //Total number of bins in the grid
  
  for(int i=0;i<dim;i++) {//loop through all dimensions
    if(nbin[i]<0) {
      prof_xerror("In prof_create, the number of bins must be >0.");
      
    }
    if(max[i]<=min[i]) {
      prof_xerror("In prof_reate, we need to have max>min.");
      
    }
    rtn.ntot=rtn.ntot*nbin[i]; //Build the product of the numbers of bins
  }
  
  //Declare the arrays for the grid parameters
  rtn.max = new double[dim]; //arraysize =  dimensions of the profile
  rtn.min = new double[dim];
  rtn.nbin = new int[dim];
  
  //Initialize these arrays
  for(int i=0;i < dim;i++) {
    rtn.max[i] = max[i];
    rtn.min[i] = min[i];
    rtn.nbin[i] = nbin[i];
  }

  rtn.outrange=0;
  rtn.nentries=0;  //Initialize the total number of entries
  
  //Declare the array for the data
  rtn.c = new complex <double>[rtn.ntot];
  rtn.c2 = new complex <double>[rtn.ntot];
 //rtn.c=new struct cmplx[rtn.ntot]; //Mean field
 //rtn.c2=new struct cmplx[rtn.ntot];//Mean squared field
  rtn.n=new int[rtn.ntot];          //Entries in that bin
  
  //Initialize the data arrays
  
  for(int i=0;i<rtn.ntot;i++) {
    rtn.c[i] = assign_complx(0.0,0.0);
    rtn.c2[i] = assign_complx(0.0,0.0);
    rtn.n[i]=0;
  }
  
  return rtn;
}
//==============================================================
complex <double> assign_complx (double re, double im){
	/*This function give the real and imagniary parts of a complex number*/ 
	
	complex <double> rtn(re,im);
	
	return rtn;
}
//===============================================================