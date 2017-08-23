#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "kdtree.h"


struct kdnode
{
	double *pos;
	int dir;

	struct kdnode *left, *right;
};

struct kdtree
{
	int dim;
	struct kdnode *root;
	struct kdhyperrect *rect;
};

struct res_node
{
	struct kdnode *item;
	double dist_sq;
};


struct rheap
{
	int size;
	int capacity;
	struct res_node *heaparr;
};



kdtree_t* kd_create(int dim)
{
	kdtree_t *tree;

	if(!(tree = malloc(sizeof *tree))) {
		return 0;
	}

	tree->dim = dim;
	tree->root = 0;
	tree->rect = 0;

	return tree;
}

static int kd_insert_recursive(kdnode_t **nptr, const double *pos, int dir, int dim)
{
	int new_dir;
	kdnode_t *node;

	if (!*nptr)
	{
		if (!(node = malloc(sizeof *node)))
		{
			return -1;
		}
		if (!(node->pos = malloc(dim * sizeof *node->pos)))
		{
			free(node);
			return -1;
		}
		memcpy(node->pos, pos, dim * sizeof *node->pos);
		node->dir = dir;
		node->left = 0;
		node->right = 0;
		*nptr = node;
		return 0;
	}

	node = *nptr;
	new_dir = (node->dir + 1) % dim;
	if (pos[node->dir] < node->pos[node->dir])
	{
		return kd_insert_recursive(&(*nptr)->left, pos, new_dir, dim);
	} else {
		return kd_insert_recursive(&(*nptr)->right, pos, new_dir, dim);
	}
}

int kd_insert(kdtree_t *tree, const double *pos)
{
	if (kd_insert_recursive(&tree->root, pos, 0, tree->dim))
	{
		return -1;
	}

	return 0;
}

/*
void kd_swap_node(kdnode_t *a, kdnode_t *b, int dim)
{
	double *tmp = (double*)malloc(dim * sizeof *a->pos);
	memcpy(tmp, a->pos, dim * sizeof *a->pos);
	memcpy(a->pos, b->pos, dim * sizeof *b->pos);
	memcpy(b->pos, tmp, dim * sizeof *tmp);
}

kdnode_t* kd_find_median_quickselect(kdnode_t *start, kdnode_t *end, int dir, int dim)
{
	if (end <= start) return 0;
	if (end == start + 1) return start;

	kdnode_t *p, *store, *m = start + (end - start) / 2;
	double pivot;

	while(1)
	{
		pivot = m->pos[dir];

		kd_swap_node(m, end - 1, dim);
		for(store = p = start; p < end; p++)
		{
			if (p->pos[dir] < pivot)
			{
				if (p != store)
				{
					kd_swap_node(p, store, dim);
				} 
				store++;
			}
		}
		kd_swap_node(store, end - 1, dim);

		if (store->pos[dir] == m->pos[dir])
		{
			return m;
		}

		if (store > m)
		{
			end = store;
		} else {
			start = store;
		}
	}

}

kdnode_t* kd_build_tree(kdnode_t *node, int n, int dir, int dim)
{
	kdnode_t *m; //median node
	
	if (!n) return 0;

	if ((m = kd_find_median_quickselect(node, node + n, dir, dim)))
	{
		dir = (dir + 1) % dim;
		m->left = kd_build_tree(node, m - node, dir, dim);
		m->right = kd_build_tree(m+1, node + n - (m+1), dir, dim);
	}
	return m;
}
*/

static int kd_nearest_range_recursive(kdnode_t *node, const double *pos, double range, rheap_t *res, int dim)
{
	double dist_sq, dx;
	int i, ret, ret_summed = 0;

	if (!node) return 0;

	dist_sq = 0;
	for (i=0; i<dim; i++)
	{
		dist_sq += SQ(node->pos[i] - pos[i]);
	}

	if (dist_sq <= SQ(range))
	{	
		rheap_push(res, node, dist_sq);
		ret_summed = 1;
	}

	dx = pos[node->dir] - node->pos[node->dir];

	ret = kd_nearest_range_recursive(dx <= 0.0 ? node->left : node->right, pos, range, res, dim);

	if (ret >= 0 && fabs(dx) < range)
	{
		ret_summed += ret;
		ret = kd_nearest_range_recursive(dx <= 0.0 ? node->right : node->left, pos, range, res, dim);
	}

	if (ret == -1) return -1;

	ret_summed += ret;

	return ret_summed;
}

rheap_t* kd_nearest_range(kdtree_t *tree, const double *pos, double range)
{
	rheap_t *res = rheap_init();

	kd_nearest_range_recursive(tree->root, pos, range, res, tree->dim);

	return res;
}

static int kd_nearest_n_recursive(kdnode_t *node, const double *pos, double range, int k, rheap_t* res, int dim)
{
	double dist_sq, dx;
	int i, ret, ret_summed = 0;

	if(!node) return 0;

	dist_sq = 0;
	for(i=0; i<dim; i++)
	{
		dist_sq += SQ(node->pos[i] - pos[i]);
	}

	if(dist_sq <= SQ(range))
	{
		if(res->size >= k)
		{
			res_node_t *farthest = rheap_max_element(res);
			if(farthest->dist_sq > dist_sq)
			{
				rheap_remove(res);
				rheap_push(res, node, dist_sq);
				ret_summed = 1;

				range = sqrt(dist_sq);
			}
		} else {
			rheap_push(res, node, dist_sq);
			ret_summed = 1;
		}

	}

	dx = pos[node->dir] - node->pos[node->dir];

	ret = kd_nearest_n_recursive(dx <= 0.0 ? node->left : node->right, pos, range, k, res, dim);

	if (ret >= 0 && fabs(dx) < range)
	{
		ret_summed += ret;
		ret = kd_nearest_n_recursive(dx <= 0.0 ? node->right : node->left, pos, range, k, res, dim);
	}

	if (ret == -1) return -1;

	ret_summed += ret;

	return ret_summed;
}


rheap_t* kd_nearest_n(kdtree_t *tree, const double *pos, double range, int k)
{
	rheap_t *res = rheap_init();

	kd_nearest_n_recursive(tree->root, pos, SQ(range), k, res, tree->dim);

	return res;
}


rheap_t* rheap_init()
{
	rheap_t* heap = (rheap_t*)malloc(sizeof(rheap_t));
	heap->size = 0;
	heap->capacity = HEAP_INIT_SIZE;

	heap->heaparr = (res_node_t*)malloc(sizeof(res_node_t) * heap->capacity);
	return heap;
}

void max_heapify(res_node_t* heaparr, int index, int size)
{
	int left, right, largest;
	res_node_t temp;

	left = LCHILD(index);
	right = RCHILD(index);
	largest = index;

	if (left <= size && heaparr[left].dist_sq > heaparr[largest].dist_sq)
	{
		largest = left;
	}
	if (right <= size && heaparr[right].dist_sq > heaparr[largest].dist_sq)
	{
		largest = right;
	}

	if (largest != index)
	{
		temp = heaparr[index];
		heaparr[index] = heaparr[largest];
		heaparr[largest] = temp;
		max_heapify(heaparr, largest, size);
	}
}

void rheap_push(rheap_t *heap, kdnode_t *node, double dist_sq)
{
	int index, parent;

	if (heap->size == heap->capacity)
	{
		heap->capacity += 1;
		heap->heaparr = realloc(heap->heaparr, sizeof(res_node_t) * heap->capacity);
	}

	index = heap->size++;

	for(;index;index = parent)
	{
		parent = PARENT(index);
		if (heap->heaparr[parent].dist_sq >= dist_sq) break;
		heap->heaparr[index] = heap->heaparr[parent];
	}

	heap->heaparr[index].item = node;
	heap->heaparr[index].dist_sq = dist_sq;
}

void rheap_remove(rheap_t *heap)
{
	res_node_t temp = heap->heaparr[--heap->size];

	if ((heap->size <= (heap->capacity+2)) && (heap->capacity > HEAP_INIT_SIZE))
	{
		heap->capacity--;
		heap->heaparr = realloc(heap->heaparr, sizeof(res_node_t) * heap->capacity);
		//TODO: check mem allocation here
		//assert(heap->heaparr);
	}

	
	heap->heaparr[0] = temp;

	max_heapify(heap->heaparr, 0, heap->size);	
}

res_node_t* rheap_max_element(rheap_t *heap)
{
	return &heap->heaparr[0];
}

static double rd( void ) {
	return (double)rand()/RAND_MAX * 20.0 - 10.0;
}


int main(int argc, char **argv)
{
	int i;

	srand( time(0) );

	double buf[2];

	kdtree_t *tree = kd_create(2);

	for(i = 0; i < 297; i++) {
		buf[0] = rd();
		buf[1] = rd();

		kd_insert(tree, buf);
	}

	buf[0] = 4.5;
	buf[1] = 4.3;

	rheap_t *res = kd_nearest_n(tree, buf, 50, 3);

	res_node_t *resnode;
	
	fprintf(stderr, "results size: %d\n", res->size);

	for(i=0; i<3; i++) {
		resnode = rheap_max_element(res);
		fprintf(stderr, "pt: %f, %f\n", resnode->item->pos[0], resnode->item->pos[1]);
		rheap_remove(res);
	}

	return 0;
}
