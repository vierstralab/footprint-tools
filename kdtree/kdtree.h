#ifndef __KDTREE_H__
#define __KDTREE_H__


typedef struct kdnode
{
	double *pos;
	int dir;

	struct kdnode *left, *right;
} kdnode_t;

typedef struct kdtree
{
	int dim;
	struct kdnode *root;
} kdtree_t;

typedef struct res_node
{
	struct kdnode *item;
	double dist;
} res_node_t;


typedef struct rheap
{
	int size;
	int capacity;
	struct res_node *heaparr;
} rheap_t;


//KDTREE
#define SQ(x)			((x) * (x))

typedef double (*distfunc_t)(const double*, const double*, int);

double dist_max_norm(const double *x, const double *y, int dim);

kdtree_t* kd_create(int dim);
void kd_free_recursive(kdnode_t *node);
void kd_free(kdtree_t *tree);

/*
void kd_swap_node(kdnode_t *a, kdnode_t *b, int dim);
kdnode_t* kd_find_median_quickselect(kdnode_t *start, kdnode_t *end, int dir, int dim);
kdnode_t* kd_build_tree(kdnode_t *node, int n, int dir, int dim);
*/

static int kd_insert_recursive(kdnode_t **nptr, const double *pos, int dir, int dim);
int kd_insert(kdtree_t *tree, const double *pos);

static int kd_nearest_range_recursive(kdnode_t *node, const double *pos, double range, distfunc_t distfunc, int dim, rheap_t *res);
rheap_t* kd_nearest_range(kdtree_t *tree, const double *pos, double range);

static int kd_nearest_n_recursive(kdnode_t *node, const double *pos, double range, distfunc_t distfunc, int dim, int n, rheap_t* res);;
rheap_t* kd_nearest_n(kdtree_t *tree, const double *pos, double range, int k);


//MAX HEAP
#define LCHILD(x) (2*x) + 1
#define RCHILD(x) (2*x) + 2
#define PARENT(x) (x-1)/2
#define HEAP_INIT_SIZE 500

rheap_t* rheap_init(void);
void rheap_free(rheap_t *heap);

void max_heapify(res_node_t* heaparr, int index, int size);
void rheap_push(rheap_t *heap, kdnode_t *node, double dist);
void rheap_remove(rheap_t *heap);
res_node_t* rheap_max_element(rheap_t *heap);




#endif
