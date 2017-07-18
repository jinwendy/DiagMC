//Functions to do the test, 



int test_lin(struct diag *t);
int test_ver(struct diag *t);
int test_diag(struct diag *t);

double test_spin_conservation(struct node *v);
double convert_spin(int s);
void test_detached(struct diag *t);

//==========================================================
void errorquit(const char* buff) {
	/*Prints an error message and quits*/

	fprintf(stderr, "ERROR:\n");
	fprintf(stderr, "%s\n", buff);
	fprintf(stderr, "QUITTING NOW\n");
	exit(0);
}
//==========================================================
//Test functions, check diagram
//==========================================================
int test_lin(struct diag *t){
	//this function tests whether every lines connected to correct vertex
	int counter = 0;          //counter how many lines we run this test
	int ml_counter = 0;       //measuring line counter
	reset_current(&t->line);  //set the current node to be the first node on the line list
	//check all node on line list
	for (int i = 0; i < t->line.n_nodes; i++){
		struct node *vend1; vend1 = t->line.current->data.l.head;
		struct node *vend2; vend2 = t->line.current->data.l.tail;

		if (NULL == vend1 || NULL == vend2){
			printf("line head/tail is NULL\n");
			printf("problem_line:%d \n", t->line.current->data.id);
			return 1;
		}

		//for electron line
		if (t->line.current->data.l.type == ELECTRON){//this is a propagator

			//check connection
			if ((vend1->data.v.ein == vend2->data.v.eout) && vend2->data.v.eout == t->line.current);
			else{
				printf("%d %d %d\n", vend1->data.v.ein->data.id, vend2->data.v.eout->data.id, t->line.current->data.id);
				printf("problem_eline:%d \n", t->line.current->data.id);
				errorquit("Error with eline connection");//if wrong, error
			}
			
			//electron connect different position
			if (vend1->data.v.phys.r== vend2->data.v.phys.r);//two ends vertices should have same position
			else{
				printf("%d %d %d\n", vend1->data.v.ein->data.id, vend2->data.v.eout->data.id, t->line.current->data.id);
				printf("problem_eline:%d \n", t->line.current->data.id);
				errorquit("Error with eline's end vertex having different position");//if wrong, error
			}
			//check electron property
			/*if (t->line.current->data.l.phys.screen == J || t->line.current->data.l.phys.screen == W){
				printf("problem_eline:%d \n", t->line.current->data.id);
				errorquit("eline has a screen component");//if wrong, error 
			}*/
		}

		//for interaction line
		else{//this is a interaction line
			//1.check connection
			if ((vend1->data.v.pend == vend2->data.v.pend) && vend2->data.v.pend == t->line.current);
			else{
				printf("problem_pline:%d \n", t->line.current->data.id);
				errorquit("Error with pline connection");//error message
			}

			//2. check the type of interaction line is consistent

			//this pline is has worm on its end
			/*if (t->line.current->data.l.phys.p_f == F){
				if (vend1->data.v.type == NO && vend2->data.v.type == NO){
					printf("problem_pline:%d \n", t->line.current->data.id);
					errorquit("Error with type F, no worm at its end");//error message
				}
			}*/
			
			//this pline is J type
			if (t->line.current->data.l.phys.screen == J){
				if (fabs(time_diff(t->line.current)) > 1E-4) {
					printf("time differenct:%f\n", fabs(time_diff(t->line.current)));
					printf("problem_pline:%d \n", t->line.current->data.id);
					errorquit("Error with type J, not a instant interaction");//error message
				}
				
				if (position_diff(t->line.current) == 0) {
					printf("position differenct:%d\n", position_diff(t->line.current));
					printf("problem_pline:%d \n", t->line.current->data.id);
					errorquit("Error with type J, not a interaction between different site");//error message
				}
				/*if (vend1->data.v.type != NO || vend2->data.v.type != NO){
					printf("problem_pline:%d \n", t->line.current->data.id);
					errorquit("Error with type J, has worm on its end");//error message
				}*/
			}

			//this pline is W type
			else {
				if (fabs(time_diff(t->line.current)) < 1E-4 && abs(vend1->data.v.phys.r - vend1->data.v.phys.r) == 1) {
					printf("problem_pline:%d \n", t->line.current->data.id);
					errorquit("Error with type W, is a instant interaction");//error message
				}
				/*if (vend1->data.v.type != NO || vend2->data.v.type != NO){
					printf("problem_pline:%d \n", t->line.current->data.id);
					errorquit("Error with type W, has worm on its end");//error message
				}*/
			}

		}

		//check special lines, namely the measuring line
		if (t->line.current->data.l.measure == MEASURE){
			ml_counter++;
			if (t->meas_line != t->line.current) errorquit("Inconsistant pointer of measuring line");
			
			//2.check if there is a detached loop on one end of measuring pline
			/*if (t->line.current->data.l.type == PHONON){
				if (vend1->data.v.ein == vend1->data.v.eout || vend2->data.v.ein == vend2->data.v.eout){
					printf("problem_pline:%d \n", t->line.current->data.id);
					printf("diagram order : %d\n",diag_order(t));
				    errorquit("Detached fermion loop");//error message
				}
			}*/
			
		}
		//check complete, set the current pointer to be the next node on the list
		next_one(&t->line);
		counter++; //increase counter by one
	}

	if (ml_counter != t->n_ml){//number of measuring line we count,after the loop
		errorquit("Inconsistant number of measuring line");
	}
	return counter; //returns line-element/node number
}
//----------------------------------------------------------
int test_ver(struct diag *t){//this function tests if vertex connect to correct eline and pline
	int counter = 0;
	int wo_counter = 0;
	reset_current(&t->vertex);  //set the current node to be the first node on the vertex list
	//check all node on vertex list
	for (int i = 0; i < t->vertex.n_nodes; i++){
		struct node *e1; e1 = t->vertex.current->data.v.ein;
		struct node *e2; e2 = t->vertex.current->data.v.eout;
		struct node *p; p = t->vertex.current->data.v.pend;

		if (NULL == e1 || NULL == e2 || NULL == p){
			printf("vertex ein/eout/pend is NULL\n");
			printf("problem_vertex:%d \n", t->vertex.current->data.id);
			return 1;
		}

		//check connection
		if ((e1->data.l.head == e2->data.l.tail) && (p->data.l.head == t->vertex.current || p->data.l.tail == t->vertex.current) && (e2->data.l.tail == t->vertex.current));
		else {
			printf("problem_vertex:%d \n", t->vertex.current->data.id);
			errorquit("Error with vertices: connection");//if wrong, error
		}

		//check position of vertex

		//check whether momentum is conserved, only under physical sector
		if (t->vertex.current->data.v.type == NO){//Regular vertex, not a worm, check momentum
			struct momentum total = momentum_init(0, 0, 0);
			total = momentum_excess_ver(t->vertex.current);
			if (
				(fabs(total.x) < 1E-3) && (fabs(total.y) < 1E-3) && (fabs(total.z) < 1E-3)){
				//printf("vertex id conserved %d\n", j);
			}
			else{
				printf("%g %g %g\n", total.x, total.y, total.z);
				printf("vertex id non-conserved: %d\n", t->vertex.current->data.id);
				errorquit("non-conserved momentum with active vertex");
			}
		}

		//check special vertices, namely the worms, number and ID
		if (t->vertex.current->data.v.type == S){
			wo_counter++;
			if (t->s_worm != t->vertex.current) errorquit("Inconsistant pointer of S-worm");
			if (test_spin_conservation(t->vertex.current) != 1.0)
				errorquit("Inconsistant spin excess of S-worm");
		}
		if (t->vertex.current->data.v.type == T){
			wo_counter++;
			if (t->t_worm != t->vertex.current) errorquit("Inconsistant pointer of T-worm");
			if (test_spin_conservation(t->vertex.current) != -1.0)
				errorquit("Inconsistant spin excess of T-worm");
		}

		//printf("vertex position: %d\n", t->vertex.current->data.v.phys.r);
		//test position component of a vertex, in two-spin system, it can be 0 or 1
		if (t->vertex.current->data.v.phys.r != 0){
			if (t->vertex.current->data.v.phys.r != 1){
				errorquit("Invalid vertex position");
			}

		}
		//set current pointer to be the next node on the list
		next_one(&t->vertex);
		counter++;
	}
	//printf("wo_counter in test_worm:%d\n", wo_counter);
	if (wo_counter != t->n_wo){
		errorquit("Inconsistant number of worms");
	}
	return counter;//return active vertex number
}
//----------------------------------------------------------
int test_diag(struct diag *t){
	//this function test the structure of the diagram
	//test_lin(t);
	//test_ver(t);
	if (t->n_ml != 1) errorquit("number of measuring line is not 1"); //n_ml should always be one for any meaningful diagram
	if (test_lin(t) != t->n_el + t->n_pl) errorquit("Inconsistent number of lines");//counter not equal to number of lines
	if (test_ver(t) != t->n_ver) errorquit("Inconsistent number of vertices");

	//printf("Test passed!\n");
	return 0;
}
//==========================================================
double test_spin_conservation(struct node *v){
	/*Check the spin conservation of a vertex.*/
	struct node *v_pair, *e1_in, *e1_out, *e2_in, *e2_out;
	double rtn = 0;

	v_pair = pair_ver(v);

	e1_in = v->data.v.ein;
	e1_out = v->data.v.eout;
	e2_in = v_pair->data.v.ein;
	e2_out = v_pair->data.v.eout;

	int s1_in, s1_out, s2_in, s2_out; //spin components on four elines
	s1_in = e1_in->data.l.phys.s;
	s1_out = e1_out->data.l.phys.s;
	s2_in = e2_in->data.l.phys.s;
	s2_out = e2_out->data.l.phys.s;

	double s_total = (convert_spin(s1_in) + convert_spin(s2_in))
		- (convert_spin(s1_out) + convert_spin(s2_out));

	return s_total;
}
//==========================================================
double convert_spin(int s){
	/*This function convert_spin to numeric value.*/
	double rtn = 0;
	if (s == UP) rtn = 0.5;
	else rtn = -0.5;
	return rtn;
}
//==========================================================
void test_detached(struct diag *t){
	//if (t->n_wo != 0) return;
	reset_current(&t->line);  //set the current node to be the first node on the line list
	for (int i = 0; i < t->line.n_nodes; i++){
		struct node *vend1; vend1 = t->line.current->data.l.head;
		struct node *vend2; vend2 = t->line.current->data.l.tail;
		if (t->line.current->data.l.measure == MEASURE && t->line.current->data.l.type == PHONON){
			if (vend1->data.v.ein == vend1->data.v.eout || vend2->data.v.ein == vend2->data.v.eout){
				printf("problem_pline:%d \n", t->line.current->data.id);
				printf("diagram order : %d\n",diag_order(t));
				diag_print(t);
			    errorquit("Detached fermion loop");//error message
			}
		}
		next_one(&t->line);
	}
	
	
	return;
}
//==========================================================