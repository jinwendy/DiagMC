/*Here are functions to manipulates the linked list*/
//#define OK      2
//#define NOTOK   0
//#define NOTDEFINE -1  //use for initialize ID

int initialize_list(struct list *l);
int makenode(struct list *l);
//int makenode(struct list *l, char *line_ver);
int reset_current(struct list *l);
int next_one(struct list *l);
int delete_current(struct list *l);
int clearlist(struct list *l);
int insert_after_node(struct list *l, int n);

int find_node(struct list *l, int id);
int nth_node(struct list *l, int nth);
int delete_node(struct list *l, int id);

/*Functions mulnipulate the list*/
//=============================================================================
int initialize_list(struct list *l){
	/*Initializae a linked list so it can start being used. A call to this function should be found after the declaration of any linked list and
	probably only there. */
	l->n_nodes = 0;   //There are no nodes yet.
	l->first = NULL;  //There is no first node.
	l->last = NULL;   //There is no last node.
	l->current = NULL;//There is no current node. 
	return OK;
}
//-----------------------------------------------------------------------------
int initialize_node(struct node *n){

	n->data.id = NOTDEFINE;
	return OK;
}
//-----------------------------------------------------------------------------
int makenode(struct list *l){
	/*Creates a node. The new node is the new tail of the linked list. The current node of the list is the freshly created node*/

	//Allocate memory for the new node
	struct node *newnode;
	if ((newnode = (struct node *)malloc(sizeof(struct node))) == NULL) {
		fprintf(stderr, "\n\n Memory unavailable for more nodes\n"); //standard error stream, "stderr", which is used for error messages and diagnostics issued by the program.
		//declared in the header file stdio.h.
		fprintf(stderr, "Exiting from Code :");
		exit(0);
	}

	if (l->current == NULL) {   //if there is no single node in the list, we create the very first one,
		l->first = newnode;     //The new node is the head
		l->last = newnode;     //The new node is the tail
		newnode->prev = NULL;  //There is no node before the current node
		newnode->next = NULL;  //There is no node after the current node
	}
	else {               //If there are already nodes in the list, we create a new node AT THE END OF THE LIST. 
		newnode->prev = l->last; //The node preceding the newnode will be the present tail node
		int id_last = 0;
		newnode->next = NULL;    //The new node will not have a next one
		newnode->prev->next = newnode; /*The node preceding the newnode will be
									   followed by the new node. */
		l->last = newnode;       //The list tail will be the new node
	}

	l->current = newnode;  //The new node is the current node
	l->n_nodes++;          //increase # of nodes by 1

	return OK;
}
//-----------------------------------------------------------------------------
int insert_after_node(struct list *l, int n) {
	/* This function inserts a node in linked list *l after the n-th node of
	the list. If the list is empty, it prints an error message onto the
	standard error and exit. If n < 1 or n > one plus the number of nodes
	in the list, it prints an error message onto the standard error and
	exit. When the node is created, the current node is made to be the
	new node.*/

	return OK;
}
//-----------------------------------------------------------------------------
int reset_current(struct list *l){
	/*Sets the head node of the list as the current node*/
	l->current = l->first;
	return OK;
}
//-----------------------------------------------------------------------------
int next_one(struct list *l){
	/*Change the current node to being the next one in the list*/
	if (l->current == NULL) return NOTOK;  //There are no nodes in the list
	if (l->current->next == NULL) return NOTOK;  //The current node is the tail
	else l->current = l->current->next; /*Set the node following the current
										node as the new current node. */
	return OK;
}
//-----------------------------------------------------------------------------
int delete_current(struct list *l){
	/*Delete the current node*/
	struct node *tmp;
	tmp = l->current; //Stores the current node pointer to a temporary variable

	//Make sure there is a node to delete
	if (tmp == NULL) return NOTOK;

	//If current node is the head and the tail. The list will be empty.
	if (tmp->prev == NULL && tmp->next == NULL) {
		l->first = NULL;    //No head
		l->last = NULL;     //No tail
		l->current = NULL;  //No current
		free(tmp);          //Release the memory used for the node being deleted
		l->n_nodes--;       //Decrement the number of nodes in the list. 
		return OK;
	}

	//If current node is the list head
	if (tmp->prev == NULL) {
		l->first = tmp->next; /*The new head will be the node after the one
							  being deleted*/
		l->current = tmp->next;/*The current node will be the one after the
							   one being deleted*/
		l->first->prev = NULL;  /*The new list head is not preceded by a node*/
		free(tmp);        //Release the memory used for the node being deleted
		l->n_nodes--;     //Decrement the number of nodes in the list. 
		return OK;
	}

	//If current node is the list tail
	if (tmp->next == NULL) {
		l->last = tmp->prev; /*The new tail is the node preceding the node
							 being deleted*/
		l->current = tmp->prev;/*The new current node is the one preceding the
							   node being deleted*/
		l->last->next = NULL;  //The new tail is not followed by a node
		free(tmp);        //Release the memory used for the node being deleted
		l->n_nodes--;     //Decrement the number of nodes in the list. 
		return OK;
	}

	//If we get here, the curent node has one before and one after
	tmp->next->prev = tmp->prev; /*The node preceding the node following the
								 node being deleted will be the node preceding
								 the node being deleted*/
	tmp->prev->next = tmp->next; /*The node following the node that precede the
								 node being deleted will be the node that
								 follows the node being deleted. */
	l->current = tmp->next; /*The current node is set to be the node following
							the node being deleted. */
	free(tmp);        //Release the memory used for the node being deleted
	l->n_nodes--;     //Decrement the number of nodes in the list. 
	//printf("delete\n");
	return OK;
}
//-----------------------------------------------------------------------------
int clearlist(struct list *l){
	/*Deletes all the nodes in the linked list*/
	if (l->n_nodes == 0) return OK; //Nothing to delete
	reset_current(l); //Move to the head of the lit
	do{
		delete_current(l);  //Delete nodes 
	} while (l->n_nodes>0); //as long as there are nodes to delete. 
	return OK;
}
//-----------------------------------------------------------------------------
int find_node(struct list *l, int id){
	/*Finding node according to its id.
	This function finds the node with ID "id" on the list, namely set the current pointer to this node.*/

	reset_current(l);  //set current pointer to the first node in the list
	if (l->current->data.id != id){
		do {               //set the current pointer to the node with ID "id"
			next_one(l);
		} while (l->current->data.id != id);
	}

	return OK;
}
//-----------------------------------------------------------------------------
int nth_node(struct list *l, int nth){
	/*Find node according to its order on the list.
	This function could help us to find the nth node on the list, namely set the current pointer to nth node.*/
	reset_current(l);
	for (int i = 1; i < nth; i++){
		next_one(l);
	}
	return OK;
}
//----------------------------------------------------------------------------
int delete_node(struct list *l, int id){
	/*This function delete the idth node in the list.(the list could be a line/vertex list)*/
	find_node(l, id);   // set the current pointer to the id_th node of the list
	delete_current(l);  //delete the current node

	return OK;
}
//=============================================================================
//=============================================================================