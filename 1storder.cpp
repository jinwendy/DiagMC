/*This file munipulate 1st order diagrams*/

#define OK      2     //indicates that updates been applied
#define NOTOK   0
#define G 1           //indicates self-energy sector
#define P 2           //indicates polarization sector

int fake_hartree(struct diag *t);
int hartree_or_g(struct diag *t);
int hartree(struct diag *t);
int ggdiag(struct diag *t);



//====================================================
int fake_hartree(struct diag *t){
	/*This function tests whether a diagram is Hartree type with W-type interaction line.*/
	if (sector(t) != G) return 1;
	if (diag_order(t) != 1) return 1;
	if (t->meas_line->data.l.head != t->meas_line->data.l.tail) return 1; //if no bubble, exit

	reset_current(&t->line);
	for (int i = 0; i < t->line.n_nodes; i++){
		if ((t->line.current->data.l.type == PHONON) & (t->line.current->data.l.measure == REGULAR)){
			if (t->line.current->data.l.phys.screen == W) {
				//printf("hartree\n");
				return OK;
			}
		}
		next_one(&t->line);
	}
	return NOTOK;
}
//====================================================
int hartree_or_g(struct diag *t){
	/*This function tests whether a 1st order diagram is Hartree type with J-type interaction line /or with ~W line.*/
	if (diag_order(t) != 1) return 1;

	if (sector(t) == G){//if measuring line is eline
		reset_current(&t->line);
		do{
			next_one(&t->line);
		} while (t->line.current->data.l.type != PHONON);
		if (t->line.current->data.l.phys.screen == J) return OK;
	}
	else{//measuring line is pline
	  if (count(t) == 1){//only one close loop, to exclude two bubbles
		  if (t->meas_line->data.l.phys.screen == W){
		  	struct node *vhead, *vtail; 
		  	vhead = t->meas_line->data.l.head; //read out two vertices
		  	vtail = t->meas_line->data.l.tail;
		    if(vhead->data.v.phys.r == vtail->data.v.phys.r && vhead->data.v.ein->data.l.phys.s == vtail->data.v.ein->data.l.phys.s) 
				return OK;//only if two vertices have the same space coordinate and two propagators have the same spin
		  }
		}
	}
       

	return NOTOK;
}
//====================================================
int hartree(struct diag *t){
	/*This function tests whether a 1st order diagram is Hartree type with J-type interaction line /or with ~W line.*/
	if (diag_order(t) != 1) return 1;
	if (sector(t) != G) return 1;
	if (count(t) != 1) {
		printf("number of loop is not 1\n");
		return 1; //number of fermi loop
	}
		
	
	reset_current(&t->line);
	do{
		next_one(&t->line);
	} while (t->line.current->data.l.type != PHONON);
	if (t->line.current->data.l.phys.screen == J) return OK;

	else return NOTOK;
}
//====================================================
int ggdiag(struct diag *t){
	/*This function tests whether a 1st order diagram is Hartree type with J-type interaction line /or with ~W line.*/
	if (diag_order(t) != 1) return 1;
	if (sector(t) != P) return 1;
	if (count(t) != 1) {
		printf("number of loop is not 1\n");
		return 1; //number of fermi loop
	}
		
	reset_current(&t->line);
	do{
		next_one(&t->line);
	} while (t->line.current->data.l.type != PHONON);
	
	if (t->line.current->data.l.phys.screen == W) {
	  	struct node *vhead, *vtail; 
	  	vhead = t->meas_line->data.l.head; //read out two vertices
	  	vtail = t->meas_line->data.l.tail;
	    if(vhead->data.v.phys.r == vtail->data.v.phys.r && vhead->data.v.ein->data.l.phys.s == vtail->data.v.ein->data.l.phys.s) 
			return OK;
		else {
			printf("not ok\n");
			return NOTOK;
		}
			
	}
	else return NOTOK;
}
//====================================================