
//----------------------------------------------------------------------------------------------------
/*
This project, do simulation by use of linked list. save information of lines/ vertices to each node of the two list.
*/
//--------------------------------------------------------------------------------------------------

#define PROBA 1.0     //Probability
#define OK      2     //indicates that updates been applied
#define NOTOK   0
#define NOTDEFINE -1  //use for initialize ID
#define TEM   1.0     //temperature
#define BETA  1.0    //inverse temperature
#define N_t   32   
#define N    2        //number of total lattice sites
#define Z    1.0      //number of neighbouring sites

#define DIM  2        //Dimension of the space : time + position
#define XL -0.5       //the lower limit of x-coordinate, for two sites only
#define XU  1.5       //the upper limit of x-coordinate, for two sites only
#define Tmin 0.0      //The lower limit fo t-domain

#define J0   1.0      //instant  coupling
#define W0   0.01      //delayed  coupling strength 
#define F0   0.05      //coupling in unphysical sector

#define SUM  0.25     //The factor 1/4 to satisfied the sum rule

#define Max_order 4   //the limit order of diagram
#define Max_meas 3   //the limit order of the diagram we measure


enum line_type { ELECTRON, PHONON };          //Distinguishes electron line and phonon line.
enum ver_worm { NO, S, T };                   //Determines  whether this vertex is regular or not.

/*In this project, there is only one measuring line, it can be a propogator or a interaction line.*/
enum line_measure { MEASURE, REGULAR };       //Determines whether this line is regular or not.
enum pline { J, W };                          //Determines the interaciton type, nearest neighbouring or not.
enum spin { NONE, UP = -1, DOWN = 1 };        //Spin type of propogators, up or down for eline. +/-1.
enum pend { EN, IN = 1, OUT = -1 };           //Type of pend of a vertex, whether this pline comes in or out.

//enum p_f {NORMAL,F};//Determine the interaction line whether has worm on one of its end

//enum space { RT, QW };          //Distinguishes the cooridinate space for a profile

struct momentum {//momentum component in 3d
	double x;      //momentum in x direction
	double y;      //momentum in y direction
	double z;      //momentum in z direction
};
struct line_prop {//spin + momentum + interaction type, properties of a line
	enum spin s;       //spin type, when it's a eline
	enum pline screen; //interaction type, when it's a pline
	struct momentum p; //momentum on line
};

struct ver_prop {//time + position in real space, properties of a vertex
	double t;    //time
	int r;       //position of vertex, in 1d, mutiple of lattice spacing
};

struct line{//structure of lines, store the information of a line,
	enum line_type type;        //line is a electron or phonon
	enum line_measure measure;  //whether it is a measuring line or not

	struct node *head; //pointer to to the vertex_node, which the line's head connects to 
	struct node *tail; //pointer to to the vertex_node, which the line's tail connects to

	struct line_prop phys;//line properties
	//int print;         //used for print out information of lines
};

struct vertex {//structure of vertices, includes all information of a vertex
	enum ver_worm type;   //vertex type: regular(NO), S worm or T worm
	struct node *ein;     //pointer to the node_line, which line coming into this vertex
	struct node *eout;    //pointer to the node_line, which line coming out this vertex
	struct node *pend;    //pointer to the node_line, which interaction lines connects to this vertex
	enum pend p_in_out;   //pend, whether this pline comes into the vertex or out of that 
	struct ver_prop phys; //vertices property
};

struct node_data { //Actual node data structure,it stores the information of a line or a vertex
	int id;        //Identification number
	int print;     //used for print out information of nodes
	union {        /*This allows you to use one type of data or the other. The actual memory space used is the largest of the different members of the union*/
		struct line l;
		struct vertex v;
	};
};

struct node {
	struct node_data data;  //information in this node
	struct node *prev;      // pointer to the previous node, which is a line/vertex
	struct node *next;      // pointer to the next node,line/vertex
};

struct list { //structure implementing the linked list of lines/vertices
	struct node *first;    //pointer to the first node on the list
	struct node *last;     //pointer to the last node on the list
	struct node *current;  //pointer to the current node 
	int n_nodes;           //number of nodes in this list
};

struct profile_multi {
  int dim;             //Dimension of the profile
  double *min;         //Lower boundary of the parameter
  double *max;         //Upper boundary of the parameter
  int *nbin;           //Number of bin subdividing the parameter
  int ntot;            //Total number of bins
  int outrange;        //Number of entries out of range
  //struct cmplx *c;   //Array of accumulated values
  //struct cmplx *c2;  //Array of accumulated squared values
  complex <double> *c;
  complex <double> *c2;
  int *n;              //Array of the number of entries per bin 
  int nentries;        //Total number of entries
};

struct diag {//diagram structure
	int n_el;   //number of electron lines
	int n_pl;   //number of phonon lines
	int n_ml;   //number of measuring line

	int n_ver;  //number of vertices
	int n_wo;   //number of worms

	int id_lin; //numbering the most recent assigned id of line node
	int id_ver; //numbering the most recent assigned id of vertex node

	struct node *meas_line; //pointer to the node of the measuring eline, only one measuring line

	struct node *s_worm;    //pointer to the node of S worm
	struct node *t_worm;    //pointer to the node of T worm

	struct momentum delta;  //momentum excess on S worm
	//two lists of lines and vertices
	struct list line;
	struct list vertex;

	int parity;//parity of number of fermion loops, even (0) or odd (1)
	
	struct profile_multi prf_g;//the profile of bold propagator, need to be updated after each iteration
	struct profile_multi prf_w;//the profile of ~W, need to be updated after each iteration
	struct profile_multi prf_f;//the profile of F(in the nonphysical sector), need to be updated after each iteration
	
};

//=============================================================================
//=============================================================================

//Declare functions

float rand_01();

void diag_print(struct diag *t);
void errorquit(const char* buff);

struct node *rand_eline(struct list *l);
struct node *rand_pline(struct list *l);
struct node *rand_pline_meas(struct list *l);
struct node *rand_ver(struct list *l);
struct node *rand_meas_line(struct list *l,int s);
int rand_worm_tobe(struct node *ver);
struct node *pair_ver(struct node *ver);
int tildeW_present(struct list *l);
int apply_worm_rule(struct node *ver);
int test_worm_rule(struct node *ver);
struct node *select_worm(struct diag *t);

int vertex_undress(struct diag *t);

int change_pline_type(struct node *pline);
int position_diff(struct node *pline);
int tell_end_worm (struct node *pline);

int diag_order(struct diag *t);
int count(struct diag *t);
void assign_parity(struct diag *t);
void change_parity(struct diag *t);

int sector(struct diag *t);

double time_diff(struct node *line);

int check_irreducibility(struct diag *t);
int check_irreducibility_old(struct diag *t);
int simulation(struct diag *t);

int simulation_iteration(struct diag *t);
int simulation_t(struct diag *t);


struct profile first_iteraion(struct diag *t);
struct profile next_iteraion(struct diag *t, struct profile h);
struct profile last_iteraion(struct diag *t, struct profile h);

int pro_updates();

int simulation_8(struct diag *t);
int simu_w_nonequalp(struct diag *t);

int print_information(struct diag *t);


int balance_create_delete(struct diag *t);
int balance_createH_deleteH(struct diag *t);
int balance_insert_remove(struct diag *t);
int balance_dress_undress(struct diag *t);

int c_d_m_m(struct diag *t);





#define DECOMP_OK    1
#define DECOMP_NOTOK 0
#define DECOMP_PI 3.14159265358979


complex <double> inverse_FT(struct ft_comp *ft, double t);
void print_cmplx(complex <double> c);
complex <double> inFT_diag(struct diag *t, int n);




#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr 
#define PROFTWOPI 6.28318530717959

void times_profile(struct profile *h, double d);

double gg_value(struct profile *h);

//double gaussdev(void);
//double expdev(void);

int simulation_test(struct diag *t);
int print_momenta(struct diag *t);
int simulation_irr(struct diag *t);



