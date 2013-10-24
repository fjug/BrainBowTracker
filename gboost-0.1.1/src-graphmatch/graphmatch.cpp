
/* graphmatch.cpp - MEX interface to VFlib2
 *
 * Original Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 1st August 2006
 *
 * Modifying Author: Florian Jug <florian.jug@gmail.com>
 *
 * Usage:
 *
 * [count, matches] = graphmatch (subg1, g2, return_at_most, directed, eq1, eq2);
 *
 * Finds all occurences of subgraph subg1 in the larger graph g2.
 * Optionally stores all found matches.
 *
 ***
 * Input
 *
 * subg1, g2: structures defining graphs, of this layout:
 *    g.nodelabels: (n,1) discrete integer labels [L_1 ; L_2 ; ... ; L_n];
 *    g.edges: (m,2) edges, [from to] at each line:
 *       [e_1_{from} e_1_{to} edgelabel_1 ; ... ; e_m_{from} e_m_{to} edgelabel_m]
 *       The node indices go from 1 to n.
 *
 *    subg1 needs to be smaller or equally sized as g2, nodelabels and edges
 *    need to be uint32's.
 *
 * (optional) return_at_most: The maximum number of mappings returned. 
 *    Default: 0 (meaning: return ALL possible mappings)
 *
 * (optional) directed: If zero, undirected graphs are matched, if non zero,
 *    the graphs are assumed to be directed graphs (default: 1).
 *
 * (optional) eq1: maximal nominal difference within the first byte of given 
 *    edge or vertex labels (makes only sense if FPs are encoded by those labels).
 *
 * (optional) eq2: maximal nominal difference within the 2nd to 4th byte of given 
 *    edge or vertex labels (makes only sense if FPs are encoded by those labels).
 ***
 * Output
 *
 * count (double): The number of times the subg1 graph appears in g2.
 *
 * (optional) matches: (count,K) uint32 matrix, where count is the number of
 *    matches, one match per row; and K is the number of nodes in subg1.  In
 *    each row the node indices of g2 are stored such that the i'th column
 *    corresponds to the i'th node in subg1.  If there is a 5 as the 3rd
 *    element, then node number 5 from g2 is matched with node number 3 from
 *    subg1.
 */

// #define DEBUG 1
unsigned char EQ_VOLUME_ABS = 0;
unsigned char EQ_COLOR_ABS  = 0;

#include <assert.h>

/* VFlib2 include files
 */
#include <argraph.h>
#include <argedit.h>
#include <ull_sub_state.h>
#include <vf2_sub_state.h>
#include <vf2_mono_state.h>
#include <match.h>
#include <vector>
#include <cmath>

/* Matlab MEX interface include files
 */
#include <mex.h>

/* Maximum number of labels at each edge.
 */
#define	ELABEL_MAX	256
#define	ELABEL_END	0xffffffff

/*** Local functions and classes
 */
typedef struct mv_chain {
	struct mv_chain* next;
	mxArray* row;
} mv_chain;

typedef struct {
	unsigned int count;	/* Count of the number of matches */
	int do_collection;	/* Collect all positive match permutations */
	int max_results;	/* How many mappings to retrieve at max */
	mv_chain* last;		/* Pointer to match permutations */
} mv_usr;

// If directed is 0: graph is undirected, directed is 1: directed graph
// (default)
static ARGraph<unsigned int, unsigned int> *
convert_mat_to_graph (const mxArray* par, int directed);

static bool
matched_visitor (int n, node_id ni1[], node_id ni2[], void *usr);

class UIntComparator
	: public AttrComparator
{
	public:
		virtual bool compatible (void* pa, void* pb)
		{
			unsigned int a = *((unsigned int *) (&pa));
			unsigned int b = *((unsigned int *) (&pb));
			return (a == b);
		}
};

class UIntPtrComparator
	: public AttrComparator
{
	public:
		virtual bool compatible (void* pa, void* pb)
		{
			unsigned int* a = (unsigned int*) pa;
			unsigned int* b = (unsigned int*) pb;

			for (unsigned int n1 = 0 ; a[n1] != ELABEL_END ; ++n1)
				for (unsigned int n2 = 0 ; b[n2] != ELABEL_END ; ++n2)
					if (a[n1] == b[n2])
						return (true);

			return (false);
		}
};

class UIntFingerprintComparator
	: public AttrComparator
{
	public:
		virtual bool compatible (void* pa, void* pb)
		{
			unsigned int a = *((unsigned int *) (&pa));
			unsigned int b = *((unsigned int *) (&pb));

            if (sizeof(a) != 4) 
            {
                mexPrintf ("WARNING: on your system unsigned int is NOT 4 bytes long!\n");
            }
            
            // '0' matches always to anything else!
            if (a==0 || b==0) {
#ifdef	DEBUG
                printf("V -- UBIQUITOUS FIT\n");
#endif
                return true;
            }
            
            // read out individual bytes
            // bytes encode: <Voxels><avg red><avg green><avg blue>
            unsigned int mask = 0xff;
            unsigned char a_b = (a & mask); 
            unsigned char b_b = (b & mask); 
            
            unsigned char a_g = (a>>8) & mask; 
            unsigned char b_g = (b>>8) & mask; 
            
            unsigned char a_r = (a>>16) & mask; 
            unsigned char b_r = (b>>16) & mask; 
            
            unsigned char a_v = (a>>24) & mask;  
            unsigned char b_v = (b>>24) & mask; 
            
			bool ret = ( (std::abs(a_v-b_v)<=EQ_VOLUME_ABS) &&
                       (std::abs(a_r-b_r)<=EQ_COLOR_ABS)  &&
                       (std::abs(a_g-b_g)<=EQ_COLOR_ABS)  &&
                       (std::abs(a_b-b_b)<=EQ_COLOR_ABS)  );
            if (1 || ret)
//                 printf("V -- %d:\t%d --> (%d,%d,%d,%d)\t\t%d --> (%d,%d,%d,%d)\n",ret,a, a_v, a_r, a_g, a_b,b, b_v, b_r, b_g, b_b);
// #ifdef	DEBUG
//             printf("VERTEX: std::abs(a_v-b_v)=%d\t\t%d\n",std::abs(a_v-b_v),std::abs(a_v-b_v)<=EQ_VOLUME_ABS);
//             printf("VERTEX: std::abs(a_r-b_r)=%d\t\t%d\n",std::abs(a_r-b_r),std::abs(a_r-b_r)<=EQ_COLOR_ABS);
//             printf("VERTEX: std::abs(a_g-b_g)=%d\t\t%d\n",std::abs(a_g-b_g),std::abs(a_g-b_g)<=EQ_COLOR_ABS);
//             printf("VERTEX: std::abs(a_b-b_b)=%d\t\t%d\n",std::abs(a_b-b_b),std::abs(a_b-b_b)<=EQ_COLOR_ABS);
//             printf("VERTEX: ===> %d\n",ret);
// #endif
            return ret;
		}
};

class UIntPtrFingerprintComparator
	: public AttrComparator
{
	public:
		virtual bool compatible (void* pa, void* pb)
		{
			unsigned int* arr = (unsigned int*) pa;
			unsigned int* brr = (unsigned int*) pb;

			for (unsigned int n1 = 0 ; arr[n1] != ELABEL_END ; ++n1) {
				for (unsigned int n2 = 0 ; brr[n2] != ELABEL_END ; ++n2) {
                    // TODO: yes, this should of course use the UIntFPcomperator...
                    
                    unsigned int a = arr[n1];
                    unsigned int b = brr[n2];

                    // '0' matches always to anything else!
                    if (a==0 || b==0) {
#ifdef	DEBUG
                        printf("E -- UBIQUITOUS FIT\n");
#endif
                        return true;
                    }

                    // read out individual bytes
                    // bytes encode: <Voxels><avg red><avg green><avg blue>
                    unsigned int mask = 0xff;
                    unsigned char a_b = (a & mask); 
                    unsigned char b_b = (b & mask); 

                    unsigned char a_g = (a>>8) & mask; 
                    unsigned char b_g = (b>>8) & mask; 

                    unsigned char a_r = (a>>16) & mask; 
                    unsigned char b_r = (b>>16) & mask; 

                    unsigned char a_v = (a>>24) & mask;  
                    unsigned char b_v = (b>>24) & mask; 

                    bool ret = ( (std::abs(a_v-b_v)<=EQ_VOLUME_ABS) &&
                               (std::abs(a_r-b_r)<=EQ_COLOR_ABS)  &&
                               (std::abs(a_g-b_g)<=EQ_COLOR_ABS)  &&
                               (std::abs(a_b-b_b)<=EQ_COLOR_ABS)  );
                    if (1 || ret) 
//                         printf("E -- %d:\t%d --> (%d,%d,%d,%d)\t\t%d --> (%d,%d,%d,%d)\n",ret,a, a_v, a_r, a_g, a_b,b, b_v, b_r, b_g, b_b);
// #ifdef	DEBUG
//                     printf("EDGE: std::abs(a_v-b_v)=%d\t\t%d\n",std::abs(a_v-b_v),std::abs(a_v-b_v)<=EQ_VOLUME_ABS);
//                     printf("EDGE: std::abs(a_r-b_r)=%d\t\t%d\n",std::abs(a_r-b_r),std::abs(a_r-b_r)<=EQ_COLOR_ABS);
//                     printf("EDGE: std::abs(a_g-b_g)=%d\t\t%d\n",std::abs(a_g-b_g),std::abs(a_g-b_g)<=EQ_COLOR_ABS);
//                     printf("EDGE: std::abs(a_b-b_b)=%d\t\t%d\n",std::abs(a_b-b_b),std::abs(a_b-b_b)<=EQ_COLOR_ABS);
//                     printf("EDGE: ===> %d\n",ret);
// #endif       
                    if ( ret ) {
                        return (true);
                    }
                }
            }
			return (false);
		}
};


static std::vector<unsigned int*> elist_ptr;


/* [count, associations] = graphmatch (subg1, g2);
 */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	elist_ptr.clear ();

	if (nlhs > 2) {
		mexPrintf ("Too many left hand side parameters.\n");
		return;
	}

	if (nrhs < 2 || nrhs > 6) {
		mexPrintf ("Wrong number of right hand side parameters.\n");
		return;
	}

	if (mxIsStruct (prhs[0]) == 0 || mxIsStruct (prhs[1]) == 0) {
		mexPrintf ("subg1 and g2 must be graph structures.\n");
		return;
	}
    
	if (nrhs >= 4) {
		EQ_COLOR_ABS = (unsigned char) mxGetScalar (prhs[5]);
    }

    if (nrhs >= 4) {
		EQ_VOLUME_ABS = (unsigned char) mxGetScalar (prhs[4]);
    }

	unsigned int directed = 1;
	if (nrhs >= 4) {
		directed = (unsigned int) mxGetScalar (prhs[3]);
    }

   	unsigned int return_at_most = 0;
	if (nrhs >= 3) {
		return_at_most = (unsigned int) mxGetScalar (prhs[2]);
    }

	/* Convert Matlab structures to VFlib graphs.
	 */
	ARGraph<unsigned int, unsigned int>* subg1 = convert_mat_to_graph (prhs[0],directed);
	ARGraph<unsigned int, unsigned int>* g2 = convert_mat_to_graph (prhs[1],directed);
	if (subg1 == NULL || g2 == NULL) {
		mexPrintf ("Error parsing subg1 or g2 structure.\n");
		return;
	}

	/* Do the graph matching.
	 */
    bool matchUsingFPs = true;
    if ( !matchUsingFPs ) {
        subg1->SetNodeComparator (new UIntComparator ());
        subg1->SetEdgeComparator (new UIntPtrComparator ());
        g2->SetNodeComparator (new UIntComparator ());
        g2->SetEdgeComparator (new UIntPtrComparator ());
    } else {
        subg1->SetNodeComparator (new UIntFingerprintComparator ());
        subg1->SetEdgeComparator (new UIntPtrFingerprintComparator ());
        g2->SetNodeComparator (new UIntFingerprintComparator ());
        g2->SetEdgeComparator (new UIntPtrFingerprintComparator ());
    }

	/* Do the matching.
	 */
	//UllSubState s0(subg1, g2);	// graph-subgraph isomorphism
	//VF2SubState s0(subg1, g2);	// graph-subgraph isomorphism
	VF2MonoState s0(subg1, g2);	// monomorphism (YES, you want that!)
	mv_usr mvs = { 0, 
                   (nlhs >= 2) ? 1 : 0,
		           return_at_most, 
                   NULL 
                 };
    
//  int n;
//  node_id ni1[12000], ni2[12000];
//	if (match (&s0, &n, ni1, ni2) )
//      matched_visitor(n, ni1, ni2, 0);
    match (&s0, matched_visitor, &mvs);

	delete subg1;
	delete g2;

	plhs[0] = mxCreateDoubleScalar (mvs.count);

	if (nlhs >= 2 && mvs.count == 0) {
		plhs[1] = mxCreateNumericMatrix (0, 0, mxUINT32_CLASS, 0);
	} else if (nlhs >= 2) {
		unsigned int cols = mxGetN (mvs.last->row);
		plhs[1] = mxCreateNumericMatrix (mvs.count, cols, mxUINT32_CLASS, 0);
		uint32_T* ass = (uint32_T*) mxGetPr (plhs[1]);

		int drow = mvs.count - 1;
		for (mv_chain* cur = mvs.last ; cur != NULL ; ) {
			uint32_T* row = (uint32_T*) mxGetPr (cur->row);

			/* Copy from each single row into the large matrix, reverse row
			 * order.
			 */
			for (int i = 0 ; i < cols ; ++i)
				ass[i*mvs.count+drow] = row[i];

			mv_chain* tbd = cur;
			cur = cur->next;

			mxDestroyArray (tbd->row);
			free (tbd);

			drow -= 1;
		}
	}

	for (unsigned int n = 0 ; n < elist_ptr.size () ; ++n)
		free (elist_ptr[n]);
	elist_ptr.clear ();
}


/* Produce an attributed relational graph from a Matlab structure.
 */
static ARGraph<unsigned int, unsigned int> *convert_mat_to_graph (const mxArray* par, int directed) {
#ifdef	DEBUG
    FILE * fp;
    fp = fopen ("scheiss.log", "a+");
    fprintf (fp, "muh1\n");fflush(fp);
#endif
    
	mxArray* nodelabels = mxGetField (par, 0, "nodelabels");
	mxArray* edges = mxGetField (par, 0, "edges");

	if (nodelabels == NULL || edges == NULL) {
		mexPrintf ("Graph structures are missing the \"nodelabels\" "
			"or \"edges\" fields.\n");

		return (NULL);
	}

	if (mxIsUint32 (nodelabels) == 0 || mxIsUint32 (edges) == 0) {
		mexPrintf ("\"nodelabels\" or \"edges\" not of type Uint32.\n");

		return (NULL);
	}

	/* Create graph
	 */
//     fprintf (fp, "muh2\n");fflush(fp);
	ARGEdit ed;
	uint32_T* nmat = (uint32_T*) mxGetPr (nodelabels);
	unsigned int nodes = mxGetM (nodelabels);
	for (unsigned int i = 0 ; i < nodes ; ++i) {
		ed.InsertNode ((void *) nmat[i]);
#ifdef	DEBUG
		mexPrintf ("inserting node %d with label %d\n", i+1, nmat[i]);
        fprintf (fp,"inserting node %d with label %d\n",i+1, nmat[i]);fflush(fp);
#endif
	}

	/* emat is column organised.
	 */
	uint32_T* emat = (uint32_T*) mxGetPr (edges);
	unsigned int edgecount = mxGetM (edges);

	for (unsigned int i = 0 ; i < edgecount ; ++i) {
		unsigned int from = emat[i+0*edgecount];
		unsigned int to = emat[i+1*edgecount];
		unsigned int elabel = emat[i+2*edgecount];

		/* Collect all labels
		 */
		unsigned int elabel_set[ELABEL_MAX];
		unsigned int n;
		for (n = 0 ; (i + n) < edgecount ; ++n) {
			if (from == emat[(i+n)+0*edgecount] &&
				to == emat[(i+n)+1*edgecount])
			{
				elabel_set[n] = emat[(i+n)+2*edgecount];
			} else
				break;
		}
		i += n-1;

		elabel_set[n] = ELABEL_END;
		n += 1;
		unsigned int* elabel_set_p = (unsigned int*) calloc (n, sizeof (unsigned int));
		for (unsigned int k = 0 ; k < n ; ++k)
			elabel_set_p[k] = elabel_set[k];

		/* Append to the global junkyard list.
		 */
		elist_ptr.push_back (elabel_set_p);

#ifdef	DEBUG
		mexPrintf ("inserting edge from %d to %d, label %d (in total %d) \n",
			from, to, elabel, n-1);
        fprintf (fp,"inserting edge from %d to %d, label %d (in total %d) \n",
			from, to, elabel, n-1);fflush(fp);
#endif
		ed.InsertEdge (from-1, to-1, (void *) elabel_set_p);
		if (directed == 0)
			ed.InsertEdge (to-1, from-1, (void *) elabel_set_p);
	}

	ARGraph<unsigned int, unsigned int>* gr =
		new ARGraph<unsigned int, unsigned int> (&ed);

	return (gr);
}


/* Callback function: called whenever the subgraph is successfully matched.
 * Return value false means "continue searching", true means "stop".
 */
static bool matched_visitor (int n, node_id ni1[], node_id ni2[], void *usr) {
    mv_usr* mv = (mv_usr*) usr;
    
#ifdef DEBUG
    mexPrintf ("MATCHED_VISITOR CALLED\n");    
	for (int i = 0 ; i < n ; ++i)
		mexPrintf ("%4d ", ni1[i]+1);
	mexPrintf ("\n");
	for (int i = 0 ; i < n ; ++i)
		mexPrintf ("%4d ", ni2[i]+1);
	mexPrintf ("\n");
#endif

	mv->count += 1;
    // in case we do not want to collect the mappings...
	if (mv->do_collection == 0) {    
		return (mv->max_results==0 || mv->count < mv->max_results ? false : true);
    }
    
	// in case we do want to cellect the mappings...
	mv_chain* c = (mv_chain*) calloc (1, sizeof (mv_chain));

	// XXX: Check that ni1 is always in the right order.
	for (int i = 0 ; i < (n-1) ; ++i) {
		assert (ni1[i] <= ni1[i+1]);
	}

	c->row = mxCreateNumericMatrix(1, n, mxUINT32_CLASS, 0);
	uint32_T* rowp = (uint32_T*) mxGetPr (c->row);
	for (int i = 0 ; i < n ; ++i)
		rowp[i] = ni2[i]+1;	// Fix index +1, conforming to Matlab convention

	// Link.
	c->next = mv->last;
	mv->last = c;

#ifdef DEBUG
	printf ("matched\n");
#endif
	return (mv->max_results==0 || mv->count < mv->max_results ? false : true);
}


