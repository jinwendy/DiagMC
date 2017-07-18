#define PI 3.14159265358979  //Pi

/*Functions declarations*/
double B_zone(double x);
struct momentum back_B(struct momentum p);
struct momentum momentum_add(struct momentum p1, struct momentum p2);
struct momentum momentum_subtract(struct momentum p1, struct momentum p2);
double momentum_dotproduct(struct momentum p1, struct momentum p2);
struct momentum assign_momentum(struct momentum p1);
struct momentum assign_ran_momentum();
struct momentum momentum_init(double x, double y, double z);
struct momentum momentum_excess_ver(struct node *ver);
void print_momentum(struct momentum p);
struct momentum momentum_factor(double f, struct momentum p);
struct momentum momenta_mul(int f1, struct momentum p1, int f2, struct momentum p2);


//Verctor operations: mementum
//==========================================================
double B_zone(double x){
	/*This function set one component of the momentum back to 1st Brillouin zone */
	double ratio, fractpart, intpart;
	ratio = x / (2 * PI);

	/*double modf(double x, double *integer), returns the fraction component(part after the decimal),
	and sets integer to the integer component.*/
	fractpart = modf(ratio, &intpart);

	//printf("Integral part = %lf\n", intpart);
	//printf("Fraction Part = %lf \n", fractpart);

	x = x - intpart * 2 * PI;
	//if x+Tn, x and n have different sign, the x got above will lay outside of 1st Brillouin zone,we need set x to Bzone
	if (fabs(x) > PI){
		if (x > 0) x = x - 2 * PI;
		else x = x + 2 * PI;
	}

	return x;
}
//-------------------------------------------------------------------------------
struct momentum back_B(struct momentum p){
	/*This function set the momentum back to 1st Brillouin zone */
	struct momentum rtn;
	rtn.x = B_zone(p.x);
	rtn.y = B_zone(p.y);
	rtn.z = B_zone(p.z);

	return rtn;
}
//--------------------------------------------------------------------------------
struct momentum momentum_add(struct momentum p1, struct momentum p2){
	/*Add two momentum a and b*/
	struct momentum rtn;
	rtn.x = p1.x + p2.x;
	rtn.y = p1.y + p2.y;
	rtn.z = p1.z + p2.z;
	rtn = back_B(rtn);//set product momentum back to 1st Brillioun zone
	return rtn;
}
//-------------------------------------------------------------------------------
struct momentum momentum_subtract(struct momentum p1, struct momentum p2){
	/*Subtract momentum p2 from p1*/
	struct momentum rtn;
	rtn.x = p1.x - p2.x;
	rtn.y = p1.y - p2.y;
	rtn.z = p1.z - p2.z;
	rtn = back_B(rtn);//set product momentum back to 1st Brillioun zone

	return rtn;
}
//-------------------------------------------------------------------------------
double momentum_dotproduct(struct momentum p1, struct momentum p2){
	/*Dot product of two momentum p1 and p2, returns a double number*/
	double rtn = 0.0;
	rtn = rtn + p1.x*p2.x;
	rtn = rtn + p1.y*p2.y;
	rtn = rtn + p1.z*p2.z;
	return rtn;
}
//-------------------------------------------------------------------------------
struct momentum momentum_init(double x, double y, double z) {
	/*Initialize a vector momentum from its cartesian coordinates*/
	struct momentum rtn;
	rtn.x = x;
	rtn.y = y;
	rtn.z = z;
	return rtn;
}
//-------------------------------------------------------------------------------
struct momentum assign_momentum(struct momentum p1){
	/*Return a momentum structure which is equivalent to p1*/
	struct momentum rtn;
	rtn.x = p1.x;
	rtn.y = p1.y;
	rtn.z = p1.z;
	return rtn;
}
//-------------------------------------------------------------------------------
struct momentum  assign_ran_momentum(){
	/*assign a random momentum (-pi to pi) structure p1*/
	struct momentum rtn;
	int dim = DIM - 1;//dim in space = Dim - time
	rtn.x = 0.0;
	rtn.y = 0.0;
	rtn.z = 0.0;

	/*if (dim == 1){
		rtn.x = PI*(2 * rand_01() - 1);
	}
	else if (dim == 2){
		rtn.x = PI*(2 * rand_01() - 1);
		rtn.y = PI*(2 * rand_01() - 1);
	}
	else if (dim == 3){
		rtn.x = PI*(2 * rand_01() - 1);
		rtn.y = PI*(2 * rand_01() - 1);
		rtn.z = PI*(2 * rand_01() - 1);
	}
	else errorquit("Space dimension exceedes 3");*/
	
	rtn.x = PI*(2 * rand_01() - 1);
	rtn.y = PI*(2 * rand_01() - 1);
	rtn.z = PI*(2 * rand_01() - 1);
	return rtn;
}
//-------------------------------------------------------------------------------
struct momentum momentum_excess_ver(struct node *ver){
	/*this function calculate the momentum excess on a vertex,momenta in - momenta out, and returns this momentum excess*/
	struct momentum rtn; //return a momentum structure

	//read out three pointers of three lines on this vertex
	struct node *e1; e1 = ver->data.v.ein;
	struct node *e2; e2 = ver->data.v.eout;
	struct node *p; p = ver->data.v.pend;

	//read out the momentum on each lines
	struct momentum p1 = e1->data.l.phys.p;
	struct momentum p2 = e2->data.l.phys.p;
	struct momentum p3 = p->data.l.phys.p;

	//read out the line direction of interaction line
	int p_in_or = ver->data.v.p_in_out;

	//calculate momentum excess on this vertex
	rtn.x = p1.x - p2.x + p_in_or*p3.x;
	rtn.y = p1.y - p2.y + p_in_or*p3.y;
	rtn.z = p1.z - p2.z + p_in_or*p3.z;

	rtn = back_B(rtn);//set product momentum back to 1st Brillioun zone
	return rtn;
}
//-------------------------------------------------------------------------------
void print_momentum(struct momentum p){
	printf("x: %g,y: %g,z: %g\n", p.x, p.y, p.z);
}
//-------------------------------------------------------------------------------
struct momentum momentum_factor(double f, struct momentum p){
	/*multiply momentum by a scalar.*/
	struct momentum rtn;
	rtn.x = f*p.x;
	rtn.y = f*p.y;
	rtn.z = f*p.z;

	rtn = back_B(rtn);//set product momentum back to 1st Brillioun zone

	return rtn;
}
//-------------------------------------------------------------------------------
struct momentum momenta_mul(int f1, struct momentum p1, int f2, struct momentum p2){
	struct momentum rtn;
	rtn.x = f1*p1.x + f2*p2.x;
	rtn.y = f1*p1.y + f2*p2.y;
	rtn.z = f1*p1.z + f2*p2.z;
	rtn = back_B(rtn);//set product momentum back to 1st Brillioun zone
	return rtn;
}
//===============================================================================