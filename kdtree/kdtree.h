#ifndef __KDTREE_H__
#define __KDTREE_H__


struct kdnode;
struct kdtree;

typedef struct kdnode kdnode_t;
typedef struct kdtree kdtree_t;

struct res_node;
typedef struct res_node res_node_t;

struct rheap;
typedef struct rheap rheap_t;

//KDTREE
#define SQ(x)			((x) * (x))

kdtree_t* kd_create(int dim);

/*
void kd_swap_node(kdnode_t *a, kdnode_t *b, int dim);
kdnode_t* kd_find_median_quickselect(kdnode_t *start, kdnode_t *end, int dir, int dim);
kdnode_t* kd_build_tree(kdnode_t *node, int n, int dir, int dim);
*/

static int kd_insert_recursive(kdnode_t **nptr, const double *pos, int dir, int dim);
int kd_insert(kdtree_t *tree, const double *pos);

static int kd_nearest_range_recursive(kdnode_t *node, const double *pos, double range, rheap_t *res, int dim);
rheap_t* kd_nearest_range(kdtree_t *tree, const double *pos, double range);

static int kd_nearest_n_recursive(kdnode_t *node, const double *pos, double range, int n, rheap_t* res, int dim);
rheap_t* kd_nearest_n(kdtree_t *tree, const double *pos, double range, int k);


//MAX HEAP
#define LCHILD(x) (2*x) + 1
#define RCHILD(x) (2*x) + 2
#define PARENT(x) (x-1)/2
#define HEAP_INIT_SIZE 500

rheap_t* rheap_init(void);

void max_heapify(res_node_t* heaparr, int index, int size);
void rheap_push(rheap_t *heap, kdnode_t *node, double dist_sq);
void rheap_remove(rheap_t *heap);
res_node_t* rheap_max_element(rheap_t *heap);




#endif
