/* QBPO.h */
/*
	Version 1.32

    Copyright 2006-2008 Vladimir Kolmogorov (vnk@ist.ac.at).
    This software can be used for research purposes only.
	This software or its derivatives must not be publicly distributed
	without a prior consent from the author (Vladimir Kolmogorov).

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	http://pub.ist.ac.at/~vnk/software.html

	/////////////////////////////////////////////////////////////////////////////

	Software for minimizing energy functions of the form
	E(x_1, ..., x_n) = \sum_i Ei(x_i) + \sum_{ij} Eij(x_i,x_j)
	where x_i are binary labels (0 or 1).

	Terms Eij can be submodular or supermodular, so in general the task is NP-hard.
	The software produces a *partial* labeling: each node is labeled as either 0,1 or
	``unknown'' (represented by negative numbers). 
	This labeling is guaranteed to be a part of at least one optimal solution.

	The following techniques are implemented:

	1. Basic roof duality algorithm ("QPBO"):

		P. L. Hammer, P. Hansen, and B. Simeone. 
		Roof duality, complementation and persistency in quadratic 0-1 optimization. 
		Mathematical Programming, 28:121–155, 1984.

		E. Boros, P. L. Hammer, and X. Sun.
		Network flows and minimization of quadratic pseudo-Boolean functions. 
		Technical Report RRR 17-1991, RUTCOR Research Report, May 1991.

	2. "Probe" technique:

		E. Boros, P. L. Hammer, and G. Tavares
		Preprocessing of Unconstrained Quadratic Binary Optimization
		Technical Report RRR 10-2006, RUTCOR Research Report, April 2006.

	with implementational details described in

		C. Rother, V. Kolmogorov, V. Lempitsky, and M. Szummer
		Optimizing binary MRFs via extended roof duality
		CVPR 2007.

	3. QPBOI ("Improve") technique:

		C. Rother, V. Kolmogorov, V. Lempitsky, and M. Szummer
		Optimizing binary MRFs via extended roof duality
		CVPR 2007.

	The maxflow algorithm used is from

		Y. Boykov, V. Kolmogorov
		An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision
		PAMI, 26(9):1124-1137, September 2004.

	Functions Improve() and Probe() reuse search trees as described in

		"Efficiently Solving Dynamic Markov Random Fields Using Graph Cuts."
		Pushmeet Kohli and Philip H.S. Torr
		International Conference on Computer Vision (ICCV), 2005

	*************************************************************************************

	Example usage: minimize energy E(x,y) = 2*x + 3*(y+1) + (x+1)*(y+2), where x,y \in {0,1}.

	#include <stdio.h>
	#include "QPBO.h"

	int main()
	{
		typedef int REAL;
		QPBO<REAL>* q;

		q = new QPBO<REAL>(2, 1); // max number of nodes & edges
		q->AddNode(2); // add two nodes

		q->AddUnaryTerm(0, 0, 2); // add term 2*x
		q->AddUnaryTerm(1, 3, 6); // add term 3*(y+1)
		q->AddPairwiseTerm(0, 1, 2, 3, 4, 6); // add term (x+1)*(y+2)

		q->Solve();
		q->ComputeWeakPersistencies();

		int x = q->GetLabel(0);
		int y = q->GetLabel(1);
		printf("Solution: x=%d, y=%d\n", x, y);

		return 0;
	}

	*************************************************************************************
*/

#ifndef __QPBO_H__
#define __QPBO_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "block.h"

// NOTE: in UNIX, use -DNDEBUG flag to suppress assertions!!!
#include <assert.h>
#define user_assert assert // used for checking user input
#define code_assert assert // used for checking algorithm's correctness

// #define user_assert(ignore)((void) 0) 
// #define code_assert(ignore)((void) 0) 



// REAL: can be int, float, double.
// Current instantiations are in instances.inc
// NOTE: WITH FLOATING POINT NUMBERS ERRORS CAN ACCUMULATE.
// IT IS STRONGLY ADVISABLE TO USE INTEGERS!!! (IT IS ALSO *MUCH* FASTER).
template <typename REAL> class QPBO
{
public:
	typedef int NodeId;
	typedef int EdgeId;

	// Constructor. 
	// The first argument gives an estimate of the maximum number of nodes that can be added
	// to the graph, and the second argument is an estimate of the maximum number of edges.
	// The last (optional) argument is the pointer to the function which will be called 
	// if an error occurs; an error message is passed to this function. 
	// If this argument is omitted, exit(1) will be called.
	//
	// IMPORTANT: 
	// 1. It is possible to add more nodes to the graph than node_num_max 
	// (and node_num_max can be zero). However, if the count is exceeded, then 
	// the internal memory is reallocated (increased by 50%) which is expensive. 
	// Also, temporarily the amount of allocated memory would be more than twice than needed.
	// Similarly for edges.
	// 
	// 2. If Probe() is used with option=1 or option=2, then it is advisable to specify
	// a larger value of edge_num_max (e.g. twice the number of edges in the original energy).
	QPBO(int node_num_max, int edge_num_max, void (*err_function)(const char *) = NULL);
	// Copy constructor
	QPBO(QPBO<REAL>& q);

	// Destructor
	~QPBO();

	// Save current reparameterisation of the energy to a text file. (Note: possibly twice the energy is saved).
	// Returns true if success, false otherwise.
	bool Save(char* filename);
	// Load energy from a text file. Current terms of the energy (if any) are destroyed.
	// Type identifier in the file (int/float/double) should match the type QPBO::REAL.
	// Returns true if success, false otherwise.
	bool Load(char* filename);

	// Removes all nodes and edges. 
	// After that functions AddNode(), AddUnaryTerm(), AddPairwiseTerm() must be called again. 
	//
	// Advantage compared to deleting QPBO and allocating it again:
	// no calls to delete/new (which could be quite slow).
	void Reset();

	int GetMaxEdgeNum(); // returns the number of edges for which the memory is allocated. 
	void SetMaxEdgeNum(int num); // If num > edge_num_max then memory for edges is reallocated. Important for Probe() with option=1,2.

	///////////////////////////////////////////////////////////////

	// Adds node(s) to the graph. By default, one node is added (num=1); then first call returns 0, second call returns 1, and so on. 
	// If num>1, then several nodes are added, and NodeId of the first one is returned.
	// IMPORTANT: see note about the constructor 
	NodeId AddNode(int num = 1);

	// Adds unary term Ei(x_i) to the energy function with cost values Ei(0)=E0, Ei(1)=E1.
	// Can be called multiple times for each node.
	void AddUnaryTerm(NodeId i, REAL E0, REAL E1);

	// Adds pairwise term Eij(x_i, x_j) with cost values E00, E01, E10, E11.
	// IMPORTANT: see note about the constructor 
	EdgeId AddPairwiseTerm(NodeId i, NodeId j, REAL E00, REAL E01, REAL E10, REAL E11);

	// This function modifies an already existing pairwise term.
	void AddPairwiseTerm(EdgeId e, NodeId i, NodeId j, REAL E00, REAL E01, REAL E10, REAL E11);

	// If AddPairwiseTerm(i,j,...) has been called twice for some pairs of nodes,
	// then MergeParallelEdges() must be called before calling Solve()/Probe()/Improve().
	void MergeParallelEdges();

	///////////////////////////////////////////////////////////////

	// Returns 0 or 1, if the node is labeled, and a negative number otherwise.
	// Can be called after Solve()/ComputeWeakPersistencies()/Probe()/Improve().
	int GetLabel(NodeId i);

	// Sets label for node i. 
	// Can be called before Stitch()/Probe()/Improve().
	void SetLabel(NodeId i, int label);

	///////////////////////////////////////////////////////////////
	// Read node & edge information.
	// Note: NodeId's are consecutive integers 0,1,...,GetNodeNum()-1.
	// However, EdgeId's are not necessarily consecutive.
	// The list of EdgeId's can be obtained as follows:
	//   QPBO<int>* q;
	//   QPBO<int>::EdgeId e;
	//   ...
	//   for (e=q->GetNextEdgeId(-1); e>=0; e=q->GetNextEdgeId(e))
	//   {
	//       ...
	//   }
	int GetNodeNum();
	EdgeId GetNextEdgeId(EdgeId e);

	// Read current reparameterization. Cost values are multiplied by 2 in the returned result.
	void GetTwiceUnaryTerm(NodeId i, REAL& E0, REAL& E1);
	void GetTwicePairwiseTerm(EdgeId e, /*output*/ NodeId& i, NodeId& j, REAL& E00, REAL& E01, REAL& E10, REAL& E11);

	///////////////////////////////////////////////////////////////

	// Return energy/lower bound.
	// NOTE: in the current implementation Probe() may add constants to the energy
	// during transormations, so after Probe() the energy/lower bound would be shifted by some offset.

	// option == 0: returns 2 times the energy of internally stored solution which would be
	//              returned by GetLabel(). Negative values (unknown) are treated as 0. 
	// option == 1: returns 2 times the energy of solution set by the user (via SetLabel()).
	REAL ComputeTwiceEnergy(int option = 0);
	// labeling must be an array of size nodeNum. Values other than 1 are treated as 0.
	REAL ComputeTwiceEnergy(int* labeling);
	// returns the lower bound defined by current reparameterizaion.
	REAL ComputeTwiceLowerBound();





	///////////////////////////////////////////////////////////////
	//                   Basic QPBO algorithm                    //
	///////////////////////////////////////////////////////////////

	// Runs QPBO. After calling Solve(), use GetLabel(i) to get label of node i.
	// Solve() produces a STRONGLY PERSISTENT LABELING. It means, in particular,
	// that if GetLabel(i)>=0 (i.e. node i is labeled) then x_i == GetLabel(i) for ALL global minima x.
	void Solve();

	// Can only be called immediately after Solve()/Probe() (and before any modifications are made to the energy).
	// Computes WEAKLY PERSISTENT LABELING. Use GetLabel() to read the result.
	// NOTE: if the energy is submodular, then ComputeWeakPersistences() will label all nodes (in general, this is not necessarily true for Solve()).
	void ComputeWeakPersistencies();

	// GetRegion()/Stitch():
	// ComputeWeakPersistencies() also splits pixels into regions (``strongly connected components'') U^0, U^1, ..., U^k as described in
	//
	//         A. Billionnet and B. Jaumard. 
	//         A decomposition method for minimizing quadratic pseudoboolean functions. 
	//         Operation Research Letters, 8:161–163, 1989.	
	//
	//     For a review see also 
	//
	//         V. Kolmogorov, C. Rother
	//         Minimizing non-submodular functions with graph cuts - a review
	//         Technical report MSR-TR-2006-100, July 2006. To appear in PAMI.
	//
	//     Nodes in U^0 are labeled, nodes in U^1, ..., U^k are unlabeled.
	//     (To find out to what region node i belongs, call GetRegion(i)).
	//     The user can use these regions as follows:
	//      -- For each r=1..k, compute somehow minimum x^r of the energy corresponding to region U^r.
	//         This energy can be obtained by calling GetPairwiseTerm() for edges inside the region.
	//         (There are no unary terms). Note that computing the global minimum is NP-hard;
	//	       it is up to the user to decide how to solve this problem.
	//      -- Set the labeling by calling SetLabel().
	//      -- Call Stitch(). It will compute a complete global minimum (in linear time).
	//      -- Call GetLabel() for nodes in U^1, ..., U^k to read new solution.
	//      Note that if the user can provides approximate rather than global minima x^r, then the stitching
	//      can still be done but the result is not guaranteed to be a *global* minimum.
	//
	// GetRegion()/Stitch() can be called only immediately after ComputeWeakPersistencies().
	// NOTE: Stitch() changes the stored energy!
	void Stitch();
	int GetRegion(NodeId i); // Returns a nonegative number which identifies the region. 0 corresponds to U^0.
	                         // The numbers are not necessarily consecutive (i.e. some number may be missed).
	                         // The maximum possible number is 2*nodeNum-5.

	//////////////////////////////////////////////////////////
	//                   QPBO extensions                    //
	//////////////////////////////////////////////////////////

	// Tries to improve the labeling provided by the user (via SetLabel()).
	// The new labeling is guaranteed to have the same or smaller energy than the input labeling.
	//
	// The procedure is as follows:
	//   1. Run QBPO
	//   2. Go through nodes in the order order_array[0], ..., order_array[N-1].
	//      If a node is unlabeled, fix it to the label provided by the user and run QBPO again.
	//   3. For remaining unlabeled nodes run set their labels to values provided by the user.
	//      (If order_array[] contains all nodes, then there should be no unlabeled nodes in step 3).
	//
	// New labeling can be obtained via GetLabel(). (The procedure also calls SetLabel() with
	// new labels, so Improve() can be called again). Returns true if success 
	// (i.e. the labeling has changed and, thus, the energy has decreased), and false otherwise.
	//
	// If array fixed_pixels of size nodeNum is provided, then it is set as follows:
	// fixed_nodes[i] = 1 if node i was fixed during Improve(), and false otherwise.
	// order_array and fixed_pixels can point to the same array.
	bool Improve(int N, int* order_array, int* fixed_nodes = NULL);

	// Calls the function above with random permutation of nodes.
	// The user should initialize the seed before the first call (using srand()).
	// NOTE: IF THE CURRENT ITERATION IS UNSUCCESSFUL, THE NEXT
	// ITERATION MAY STILL BE SUCCESFULL SINCE A DIFFERENT PERMUTATION WILL BE USED.
	// A typical number of iterations could be e.g. 10-100.
	bool Improve();

	struct ProbeOptions
	{
		ProbeOptions()
			: directed_constraints(2),
			  weak_persistencies(0),
			  C(100000),
			  order_array(NULL),
			  order_seed(0),
			  dilation(3),
			  callback_fn(NULL)
		{
		}

		int directed_constraints; // 0: directed constraints are added only for existing edges
		                          // 1: all possible directed constraints are added, if there is sufficient space for edges (as specified by edge_num_max; see SetEdgeNumMax() function)
		                          // 2: all possible directed constraints are added. If necessary, new memory for edges is allocated.
		int weak_persistencies; // 0: use only strong persistency
		                        // 1: use weak persistency in the main loop (but not for probing operations)

		REAL C; // Large constant used inside Probe() for enforcing directed constraints. 
		        // Note: small value may increase the number of iterations, large value may cause overflow.

		int* order_array; // if array of size nodeNum() is provided, then nodes are tested in the order order_array[0], order_array[1], ...
		unsigned int order_seed; // used only if order_array == NULL:
		                         // 0: default order (0,1,...,nodeNum()-1) is used.
		                         // otherwise: random permutation with random seed 'order_seed' is used.
		int dilation; // determines order of processing nodes (see Rother et al. CVPR'07):
		              // d<0:  one iteration tests all unlabeled nodes (i.e. fixes them to 0 and 1). 
		              // d>=0: nodes within distance d from successful nodes are tested in the next iteration.

		bool (*callback_fn)(int unlabeled_num); // if callback_fn!=NULL, then after every testing a node Probe calls callback_fn();
		                                        // unlabeled_num is the current number of remaining nodes in the energy.
		                                        // If callback_fn returns true then Probe() terminates.
	};

	// Fixes some nodes to 0 or 1, contracts other nodes. These transformations
	// do not change global minima. The internally stored energy is modified accordingly.
	// (In particular, the new energy may have a different number of nodes and edges).
	//
	// Nodes of the old energy are associated with nodes of the new energy 
	// (possibly with inversion: 0<-->1). This association is returned
	// in array mapping as follows:
	//   If old node i corresponds to new node j with inversion x (x=0,1) then 
	//     mapping[i] = 2*j + x.
	//
	// If y is a global minimum of the new energy, then solution x defined by
	//     x[i] = (y[mapping[i]/2] + mapping[i]) % 2
	// is a global minimum of the original energy.
	//
	// Node 0 of the new energy is guaranteed to have optimal label 0 (y[0]=0),
	// therefore if mapping[i] < 2 then this is the optimal label for node i.
	//
	// Before calling Probe() you can call SetLabel() to set an input labeling x0.
	// During the procedure this labeling is transformed. The new labeling y0 can
	// be read via GetLabel() after Probe() (P+I method - see Rother et al, CVPR'07).
	void Probe(int* mapping, ProbeOptions& option);

	// If Probe() is called two times, then mappings mapping0 and mapping1 produced by the first and
	// second run can be combined using MergeMappings. Array mapping0 is updated accordingly.
	static void MergeMappings(int nodeNum0, int* mapping0, int* mapping1);


	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////










/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
	
private:
	// internal variables and functions

	struct Arc;

	struct Node
	{
		Arc		*first;		// first outcoming Arc

		Node	*next;		// pointer to the next active Node
							// (or to itself if it is the last Node in the list)

		unsigned int is_sink : 1;	// flag showing whether the node is in the source or in the sink tree (if parent!=NULL)
		unsigned int is_marked : 1;	// set by mark_node()
		unsigned int is_in_changed_list : 1; // set by maxflow if the node is added to changed_list
		unsigned int is_removed : 1; // 1 means that the node is removed (for node[0][...])

		int	         label : 2;
		int	         label_after_fix0 : 2;
		int	         label_after_fix1 : 2;

		unsigned int list_flag : 2; // used in Probe() and Improve()

		unsigned int user_label : 1; // set by calling SetLabel()

		union
		{
			struct
			{
				// used inside maxflow algorithm
				int		TS;			// timestamp showing when DIST was computed
				int		DIST;		// distance to the terminal
				Arc		*parent;	// Node's parent
			};
			struct
			{
				int		region;
				Node	*dfs_parent;
				Arc		*dfs_current;
			};
		};

		REAL		tr_cap;		// if tr_cap > 0 then tr_cap is residual capacity of the Arc SOURCE->Node
								// otherwise         -tr_cap is residual capacity of the Arc Node->SINK 
	};

	struct Arc
	{
		Node		*head;		// Node the Arc points to
		Arc			*next;		// next Arc with the same originating Node
		Arc			*sister;	// reverse Arc

		REAL		r_cap;		// residual capacity
	};

	struct nodeptr
	{
		Node    	*ptr;
		nodeptr		*next;
	};
	static const int NODEPTR_BLOCK_SIZE = 128;

	Node	*nodes[2], *node_last[2], *node_max[2]; // node_last[k] = nodes[k]+node_num
			                                        // node_max[k] = nodes[k]+node_num_max
			                                        // nodes[1] = node_max[0]
	Arc		*arcs[2], *arc_max[2]; // arc_max[k] = arcs[k]+2*edge_num_max
			                       // arcs[1] = arc_max[0]

	Arc*	first_free; // list of empty spaces for edges. 
	void InitFreeList();

	int		node_num;
	int		node_shift; // = node_num_max*sizeof(Node)
	int		arc_shift; // = 2*edge_num_max*sizeof(Arc)

	DBlock<nodeptr>		*nodeptr_block;

	void	(*error_function)(const char *);	// this function is called if a error occurs,
										// with a corresponding error message
										// (or exit(1) is called if it's NULL)

	REAL	zero_energy; // energy of solution (0,...,0)

	// reusing trees & list of changed pixels
	int					maxflow_iteration; // counter
	bool				keep_changed_list;
	Block<Node*>		*changed_list;

	/////////////////////////////////////////////////////////////////////////

	void get_type_information(const char*& type_name, const char*& type_format);

	void reallocate_nodes(int node_num_max_new);
	void reallocate_arcs(int arc_num_max_new);

	int	stage; // 0: maxflow is solved only for nodes in [nodes[0],node_last[0]).
		       //    Arcs corresponding to supermodular edges are present in arcs[0] and arcs[1],
		       //    but nodes do not point to them.
		       // 1: maxflow is solved for the entire graph.
	bool all_edges_submodular;
	void TransformToSecondStage(bool copy_trees);

	static void ComputeWeights(REAL A, REAL B, REAL C, REAL D, REAL& ci, REAL& cj, REAL& cij, REAL& cji);
	bool IsNode0(Node* i) { return (i<nodes[1]); }
	Node* GetMate0(Node* i) { code_assert(i< nodes[1]); return (Node*)((char*)i + node_shift); }
	Node* GetMate1(Node* i) { code_assert(i>=nodes[1]); return (Node*)((char*)i - node_shift); }
	Node* GetMate(Node* i) { return IsNode0(i) ? GetMate0(i) : GetMate1(i); }
	bool IsArc0(Arc* a) { return (a<arcs[1]); }
	Arc* GetMate0(Arc* a) { code_assert(a< arcs[1]); return (Arc*)((char*)a + arc_shift); }
	Arc* GetMate1(Arc* a) { code_assert(a>=arcs[1]); return (Arc*)((char*)a - arc_shift); }
	Arc* GetMate(Arc* a) { return IsArc0(a) ? GetMate0(a) : GetMate1(a); }

	ProbeOptions probe_options;
	bool user_terminated;
	bool Probe(int* mapping); // Probe(int*,ProbeOptions&) iteratively calls Probe(int*)

	void TestRelaxedSymmetry(); // debug function

	REAL DetermineSaturation(Node* i);
	void AddUnaryTerm(Node* i, REAL E0, REAL E1);
	void FixNode(Node* i, int x); // fix i to label x. there must hold IsNode0(i).
	void ContractNodes(Node* i, Node* j, int swap); // there must hold IsNode0(i) && IsNode0(j) && (swap==0 || swap==1)
	                                                // enforces constraint i->label = (j->label + swap) mod 2
	                                                // i is kept, all arcs from j are deleted.
	int MergeParallelEdges(Arc* a1, Arc* a2); // there must hold (a1->sister->head == a2->sister->head) && IsNode0(a1->sister->head) &&
	                                          //                 (a1->head == a2->head || a1->head = GetMate(a2->head))
	                                          // returns 0 if a1 is removed, 1 otherwise
	bool AddDirectedConstraint0(Arc* a, int xi, int xj); // functions return true if the energy was changed.
	bool AddDirectedConstraint1(Arc* a, int xi, int xj); // ...0 checks whether submodurality needs to be swapped, ...1 preserves submodularity.
	void AddDirectedConstraint(Node* i, Node* j, int xi, int xj); // adds new edge. first_free must not be NULL.
	void AllocateNewEnergy(int* mapping);


	static void ComputeRandomPermutation(int N, int* permutation);

	struct FixNodeInfo { Node* i; REAL INFTY; };
	Block<FixNodeInfo>* fix_node_info_list;

	/////////////////////////////////////////////////////////////////////////

	Node				*queue_first[2], *queue_last[2];	// list of active nodes
	nodeptr				*orphan_first, *orphan_last;		// list of pointers to orphans
	int					TIME;								// monotonically increasing global counter

	/////////////////////////////////////////////////////////////////////////

	// functions for processing active list
	void set_active(Node *i);
	Node *next_active();

	// functions for processing orphans list
	void set_orphan_front(Node* i); // add to the beginning of the list
	void set_orphan_rear(Node* i);  // add to the end of the list

	void mark_node(Node* i);
	void add_to_changed_list(Node* i);

	void maxflow(bool reuse_trees = false, bool keep_changed_list = false);
	void maxflow_init();             // called if reuse_trees == false
	void maxflow_reuse_trees_init(); // called if reuse_trees == true
	void augment(Arc *middle_arc);
	void process_source_orphan(Node *i);
	void process_sink_orphan(Node *i);

	int what_segment(Node* i, int default_segm = 0);

	void test_consistency(Node* current_node=NULL); // debug function
};











///////////////////////////////////////
// Implementation - inline functions //
///////////////////////////////////////



template <typename REAL> 
	inline typename QPBO<REAL>::NodeId QPBO<REAL>::AddNode(int num)
{
	user_assert(num >= 0);

	if (node_last[0] + num > node_max[0]) 
	{
		int node_num_max = node_shift / sizeof(Node);
		node_num_max += node_num_max / 2;
		if (node_num_max < (int)(node_last[0] + num - nodes[0]) + 1) node_num_max = (int)(node_last[0] + num - nodes[0]) + 1;
		reallocate_nodes(node_num_max);
	}

	memset(node_last[0], 0, num*sizeof(Node));
	NodeId i = node_num;
	node_num += num;
	node_last[0] += num;

	if (stage)
	{
		memset(node_last[1], 0, num*sizeof(Node));
		node_last[1] += num;
	}

	return i;
}

template <typename REAL> 
	inline void QPBO<REAL>::AddUnaryTerm(NodeId i, REAL E0, REAL E1)
{
	user_assert(i >= 0 && i < node_num);

	nodes[0][i].tr_cap += E1 - E0;
	if (stage) nodes[1][i].tr_cap -= E1 - E0;

	zero_energy += E0;
}

template <typename REAL> 
	inline void QPBO<REAL>::AddUnaryTerm(Node* i, REAL E0, REAL E1)
{
	code_assert(i >= nodes[0] && i<node_last[0]);

	i->tr_cap += E1 - E0;
	if (stage) GetMate0(i)->tr_cap -= E1 - E0;

	zero_energy += E0;
}

template <typename REAL> 
	inline int QPBO<REAL>::what_segment(Node* i, int default_segm)
{
	if (i->parent)
	{
		return (i->is_sink) ? 1 : 0;
	}
	else
	{
		return default_segm;
	}
}

template <typename REAL> 
	inline void QPBO<REAL>::mark_node(Node* i)
{
	if (!i->next)
	{
		/* it's not in the list yet */
		if (queue_last[1]) queue_last[1] -> next = i;
		else               queue_first[1]        = i;
		queue_last[1] = i;
		i -> next = i;
	}
	i->is_marked = 1;
}

template <typename REAL> 
	inline int QPBO<REAL>::GetLabel(NodeId i)
{
	user_assert(i >= 0 && i < node_num);

	return nodes[0][i].label;
}

template <typename REAL> 
	inline int QPBO<REAL>::GetRegion(NodeId i)
{
	user_assert(i >= 0 && i < node_num);
	user_assert(stage == 1);

	return nodes[0][i].region;
}

template <typename REAL> 
	inline void QPBO<REAL>::SetLabel(NodeId i, int label)
{
	user_assert(i >= 0 && i < node_num);

	nodes[0][i].user_label = label;
}

template <typename REAL> 
	inline void QPBO<REAL>::GetTwiceUnaryTerm(NodeId i, REAL& E0, REAL& E1)
{
	user_assert(i >= 0 && i < node_num);

	E0 = 0; 
	if (stage == 0) E1 = 2*nodes[0][i].tr_cap;
	else            E1 = nodes[0][i].tr_cap - nodes[1][i].tr_cap;
}

template <typename REAL> 
	inline void QPBO<REAL>::GetTwicePairwiseTerm(EdgeId e, NodeId& _i, NodeId& _j, REAL& E00, REAL& E01, REAL& E10, REAL& E11)
{
	user_assert(e >= 0 && arcs[0][2*e].sister);

	Arc* a;
	Arc* a_mate;
	if (IsNode0(arcs[0][2*e+1].head))
	{
		a = &arcs[0][2*e];
		a_mate = &arcs[1][2*e];
	}
	else
	{
		a = &arcs[1][2*e+1];
		a_mate = &arcs[0][2*e+1];
	}
	Node* i = a->sister->head;
	Node* j = a->head;
	_i = (int)(i - nodes[0]);

	if (IsNode0(j))
	{
		E00 = E11 = 0;
		if (stage == 0) { E01 = 2*a->r_cap; E10 = 2*a->sister->r_cap; }
		else            { E01 = a->r_cap + a_mate->r_cap; E10 = a->sister->r_cap + a_mate->sister->r_cap; }
		_j = (int)(j - nodes[0]);
	}
	else
	{
		E01 = E10 = 0;
		if (stage == 0) { E00 = 2*a->r_cap; E11 = 2*a->sister->r_cap; }
		else            { E00 = a->r_cap + a_mate->r_cap; E11 = a->sister->r_cap + a_mate->sister->r_cap; }
		_j = (int)(j - nodes[1]);
	}
}

template <typename REAL> 
	inline int QPBO<REAL>::GetNodeNum() 
{
	return (int)(node_last[0] - nodes[0]); 
}

template <typename REAL> 
	inline typename QPBO<REAL>::EdgeId QPBO<REAL>::GetNextEdgeId(EdgeId e) 
{
	Arc* a;
	for (a=&arcs[0][2*(++e)]; a<arc_max[0]; a+=2, e++)
	{ 
		if (a->sister) return e;
	}
	return -1;
}

template <typename REAL> 
	inline int QPBO<REAL>::GetMaxEdgeNum() 
{
	return (int)(arc_max[0]-arcs[0])/2;
}


template <typename REAL> 
	inline void QPBO<REAL>::ComputeWeights(
	REAL A, REAL B, REAL C, REAL D, // input - E00=A, E01=B, E10=C, E11=D
	REAL& ci, REAL& cj, REAL& cij, REAL& cji // output - edge weights
	)
{
	/* 
	E = A A  +  0   B-A
		D D     C-D 0
	Add edges for the first term
	*/
	ci = D - A;
	B -= A; C -= D;

	/* now need to represent
	0 B
	C 0
	*/

	if (B < 0)
	{
		/* Write it as
		B B  +  -B 0  +  0   0
		0 0     -B 0     B+C 0
		*/
		ci += -B; /* first term */
		cj = B; /* second term */
		cji = B+C; /* third term */
		cij = 0;
	}
	else if (C < 0)
	{
		/* Write it as
		-C -C  +  C 0  +  0 B+C
			0  0     C 0     0 0
		*/
		ci += C; /* first term */
		cj = -C; /* second term */
		cij = B+C; /* third term */
		cji = 0;
	}
	else /* B >= 0, C >= 0 */
	{
		cj = 0;
		cij = B;
		cji = C;
	}
}

/*
	special constants for node->parent
*/
#define QPBO_MAXFLOW_TERMINAL ( (Arc *) 1 )		/* to terminal */
#define QPBO_MAXFLOW_ORPHAN   ( (Arc *) 2 )		/* orphan */


#include "instances.inc"

/* QPBO.cpp */

template <typename REAL> 
	QPBO<REAL>::QPBO(int node_num_max, int edge_num_max, void (*err_function)(const char *))
	: node_num(0),
	  nodeptr_block(NULL),
	  changed_list(NULL),
	  fix_node_info_list(NULL),
	  stage(0),
	  all_edges_submodular(true),
	  error_function(err_function),
	  zero_energy(0)
{
	node_num_max += 4;
	if (node_num_max < 16) node_num_max = 16;
	if (edge_num_max < 16) edge_num_max = 16;

	nodes[0] = (Node*) malloc(2*node_num_max*sizeof(Node));
	arcs[0] = (Arc*) malloc(4*edge_num_max*sizeof(Arc));
	if (!nodes[0] || !arcs[0]) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	node_last[0] = nodes[0];
	node_max[0] = nodes[1] = node_last[1] = nodes[0] + node_num_max;
	node_max[1] = nodes[1] + node_num_max;
	node_shift = node_num_max*sizeof(Node);

	arc_max[0] = arcs[1] = arcs[0] + 2*edge_num_max;
	arc_max[1] = arcs[1] + 2*edge_num_max;
	arc_shift = 2*edge_num_max*sizeof(Arc);

	maxflow_iteration = 0;

	memset(arcs[0], 0, 2*arc_shift);
	InitFreeList();
}

template <typename REAL> 
	void QPBO<REAL>::InitFreeList()
{
	Arc* a;
	Arc* a_last_free;

	first_free = a_last_free = NULL;
	for (a=arcs[0]; a<arc_max[0]; a+=2)
	if (!a->sister)
	{
		if (a_last_free) a_last_free->next = a;
		else        first_free = a;
		a_last_free = a;
	}
	if (a_last_free) a_last_free->next = NULL;
}

template <typename REAL> 
	QPBO<REAL>::QPBO(QPBO<REAL>& q)
	: node_num(q.node_num),
	  nodeptr_block(NULL),
	  changed_list(NULL),
	  fix_node_info_list(NULL),
	  stage(q.stage),
	  all_edges_submodular(q.all_edges_submodular),
	  error_function(q.error_function),
	  zero_energy(q.zero_energy)
{
	int node_num_max = q.node_shift/sizeof(Node);
	int arc_num_max = (int)(q.arc_max[0] - q.arcs[0]);
	Node* i;
	Arc* a;

	nodes[0] = (Node*) malloc(2*node_num_max*sizeof(Node));
	arcs[0] = (Arc*) malloc(2*arc_num_max*sizeof(Arc));
	if (!nodes[0] || !arcs[0]) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	node_last[0] = nodes[0] + node_num;
	node_max[0] = nodes[1] = nodes[0] + node_num_max;
	node_last[1] = nodes[1] + node_num;
	node_max[1] = nodes[1] + node_num_max;
	node_shift = node_num_max*sizeof(Node);

	arc_max[0] = arcs[1] = arcs[0] + arc_num_max;
	arc_max[1] = arcs[1] + arc_num_max;
	arc_shift = arc_num_max*sizeof(Arc);

	maxflow_iteration = 0;

	memcpy(nodes[0], q.nodes[0], 2*node_num_max*sizeof(Node));
	memcpy(arcs[0], q.arcs[0], 2*arc_num_max*sizeof(Arc));

	for (i=nodes[0]; i<node_last[stage]; i++)
	{
		if (i==node_last[0]) i = nodes[1];
		if (i->first) i->first = (Arc*) ((char*)i->first + (((char*) arcs[0]) - ((char*) q.arcs[0])));
	}

	for (a=arcs[0]; a<arc_max[stage]; a++)
	{
		if (a == arc_max[0]) a = arcs[1];
		if (a->sister)
		{
			a->head              = (Node*) ((char*)a->head   + (((char*) nodes[0]) - ((char*) q.nodes[0])));
			if (a->next) a->next = (Arc*)  ((char*)a->next   + (((char*) arcs[0])  - ((char*) q.arcs[0])));
			a->sister            = (Arc*)  ((char*)a->sister + (((char*) arcs[0])  - ((char*) q.arcs[0])));
		}
	}

	InitFreeList();
}

template <typename REAL> 
	QPBO<REAL>::~QPBO()
{
	if (nodeptr_block) 
	{ 
		delete nodeptr_block; 
		nodeptr_block = NULL; 
	}
	if (changed_list)
	{
		delete changed_list;
		changed_list = NULL;
	}
	if (fix_node_info_list)
	{
		delete fix_node_info_list;
		fix_node_info_list = NULL;
	}
	free(nodes[0]);
	free(arcs[0]);
}

template <typename REAL> 
	void QPBO<REAL>::Reset()
{
	node_last[0] = nodes[0];
	node_last[1] = nodes[1];
	node_num = 0;

	if (nodeptr_block) 
	{ 
		delete nodeptr_block; 
		nodeptr_block = NULL; 
	}
	if (changed_list)
	{
		delete changed_list;
		changed_list = NULL;
	}
	if (fix_node_info_list)
	{
		delete fix_node_info_list;
		fix_node_info_list = NULL;
	}

	maxflow_iteration = 0;
	zero_energy = 0;

	stage = 0;
	all_edges_submodular = true;

	memset(arcs[0], 0, 2*arc_shift);
	InitFreeList();
}

template <typename REAL> 
	void QPBO<REAL>::reallocate_nodes(int node_num_max_new)
{
	code_assert(node_num_max_new > node_shift/((int)sizeof(Node)));
	Node* nodes_old[2] = { nodes[0], nodes[1] };

	int node_num_max = node_num_max_new;
	nodes[0] = (Node*) realloc(nodes_old[0], 2*node_num_max*sizeof(Node));
	if (!nodes[0]) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	node_shift = node_num_max*sizeof(Node);
	node_last[0] = nodes[0] + node_num;
	node_max[0] = nodes[1] = nodes[0] + node_num_max;
	node_last[1] = nodes[1] + node_num;
	node_max[1] = nodes[1] + node_num_max;
	if (stage)
	{
		memmove(nodes[1], (char*)nodes[0] + ((char*)nodes_old[1] - (char*)nodes_old[0]), node_num*sizeof(Node));
	}

	Arc* a;
	for (a=arcs[0]; a<arc_max[stage]; a++)
	{
		if (a->sister)
		{
			int k = (a->head < nodes_old[1]) ? 0 : 1;
			a->head = (Node*) ((char*)a->head + (((char*) nodes[k]) - ((char*) nodes_old[k])));
		}
	}
}

template <typename REAL> 
	void QPBO<REAL>::reallocate_arcs(int arc_num_max_new)
{
	int arc_num_max_old = (int)(arc_max[0] - arcs[0]);
	int arc_num_max = arc_num_max_new; if (arc_num_max & 1) arc_num_max ++;
	code_assert(arc_num_max > arc_num_max_old);
	Arc* arcs_old[2] = { arcs[0], arcs[1] };

	arcs[0] = (Arc*) realloc(arcs_old[0], 2*arc_num_max*sizeof(Arc));
	if (!arcs[0]) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	arc_shift = arc_num_max*sizeof(Arc);
	arc_max[0] = arcs[1] = arcs[0] + arc_num_max;
	arc_max[1] = arcs[1] + arc_num_max;

	if (stage)
	{
		memmove(arcs[1], arcs[0]+arc_num_max_old, arc_num_max_old*sizeof(Arc));
		memset(arcs[0]+arc_num_max_old, 0, (arc_num_max-arc_num_max_old)*sizeof(Arc));
		memset(arcs[1]+arc_num_max_old, 0, (arc_num_max-arc_num_max_old)*sizeof(Arc));
	}
	else
	{
		memset(arcs[0]+arc_num_max_old, 0, (2*arc_num_max-arc_num_max_old)*sizeof(Arc));
	}

	Node* i;
	Arc* a;
	for (i=nodes[0]; i<node_last[stage]; i++)
	{
		if (i==node_last[0]) i = nodes[1];

		if (i->first) 
		{
			int k = (i->first < arcs_old[1]) ? 0 : 1;
			i->first = (Arc*) ((char*)i->first + (((char*) arcs[k]) - ((char*) arcs_old[k])));
		}
	}
	for (a=arcs[0]; a<arc_max[stage]; a++)
	{
		if (a->sister)
		{
			if (a->next) 
			{
				int k = (a->next < arcs_old[1]) ? 0 : 1;
				a->next = (Arc*) ((char*)a->next + (((char*) arcs[k]) - ((char*) arcs_old[k])));
			}
			int k = (a->sister < arcs_old[1]) ? 0 : 1;
			a->sister = (Arc*) ((char*)a->sister + (((char*) arcs[k]) - ((char*) arcs_old[k])));
		}
	}

	InitFreeList();
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

template <typename REAL> 
	bool QPBO<REAL>::Save(char* filename)
{
	int e;
	int edge_num = 0;
	for (e=GetNextEdgeId(-1); e>=0; e=GetNextEdgeId(e)) edge_num ++;

	FILE* fp;
	REAL E0, E1, E00, E01, E10, E11;
	int i, j;
	const char* type_name;
	const char* type_format;
	char FORMAT_LINE[64];
	int factor = (stage == 0) ? 2 : 1;

	get_type_information(type_name, type_format);

	fp = fopen(filename, "w");
	if (!fp) return false;

	fprintf(fp, "nodes=%d\n", GetNodeNum());
	fprintf(fp, "edges=%d\n", edge_num);
	fprintf(fp, "labels=2\n");
	fprintf(fp, "type=%s\n", type_name);
	fprintf(fp, "\n");

	sprintf(FORMAT_LINE, "n %%d %%%s %%%s\n", type_format, type_format);
	for (i=0; i<GetNodeNum(); i++)
	{
		GetTwiceUnaryTerm(i, E0, E1);
		REAL delta = (E0 < E1) ? E0 : E1;
		fprintf(fp, FORMAT_LINE, i, (E0-delta)/factor, (E1-delta)/factor);
	}
	sprintf(FORMAT_LINE, "e %%d %%d %%%s %%%s %%%s %%%s\n", type_format, type_format, type_format, type_format);
	for (e=GetNextEdgeId(-1); e>=0; e=GetNextEdgeId(e))
	{
		GetTwicePairwiseTerm(e, i, j, E00, E01, E10, E11);
		fprintf(fp, FORMAT_LINE, i, j, E00/factor, E01/factor, E10/factor, E11/factor);
	}
	fclose(fp);
	return true;
}

template <typename REAL> 
	bool QPBO<REAL>::Load(char* filename)
{
	Reset();

	FILE* fp;
	REAL E0, E1, E00, E01, E10, E11;
	int i, j;
	const char* type_name;
	const char* type_format;
	char LINE[256], FORMAT_LINE_NODE[64], FORMAT_LINE_EDGE[64];
	int NODE_NUM, EDGE_NUM, K;

	get_type_information(type_name, type_format);

	fp = fopen(filename, "r");
	if (!fp) { printf("Cannot open %s\n", filename); return false; }

	if (fscanf(fp, "nodes=%d\n", &NODE_NUM) != 1) { printf("%s: wrong format\n", filename); fclose(fp); return false; }
	if (fscanf(fp, "edges=%d\n", &EDGE_NUM) != 1) { printf("%s: wrong format\n", filename); fclose(fp); return false; }
	if (fscanf(fp, "labels=%d\n", &K) != 1) { printf("%s: wrong format\n", filename); fclose(fp); return false; }
	if (K != 2) { printf("%s: wrong number of labels\n", filename); fclose(fp); return false; }
	if (fscanf(fp, "type=%10s\n", LINE) != 1) { printf("%s: wrong format\n", filename); fclose(fp); return false; }
	if (strcmp(LINE, type_name)) { printf("%s: type REAL mismatch\n", filename); fclose(fp); return false; }

	AddNode(NODE_NUM+4);
	node_num -= 4;
	node_last[0] -= 4;
	node_last[1] -= 4;

	sprintf(FORMAT_LINE_NODE, "n %%d %%%s %%%s\n", type_format, type_format);
	sprintf(FORMAT_LINE_EDGE, "e %%d %%d %%%s %%%s %%%s %%%s\n", type_format, type_format, type_format, type_format);
	while (fgets(LINE, sizeof(LINE), fp))
	{
		if (sscanf(LINE, FORMAT_LINE_EDGE, &i, &j, &E00, &E01, &E10, &E11) == 6)
		{
			if (i<0 || i>=NODE_NUM || j<0 || j>=NODE_NUM || i==j) { printf("%s: wrong format\n", filename); fclose(fp); return false; }
			AddPairwiseTerm(i, j, E00, E01, E10, E11);
		}
		else if (sscanf(LINE, FORMAT_LINE_NODE, &i, &E0, &E1) == 3)
		{
			if (i<0 || i>=NODE_NUM) { printf("%s: wrong format\n", filename); fclose(fp); return false; }
			AddUnaryTerm(i, E0, E1);
		}
	}

	fclose(fp);
	return true;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#define SET_SISTERS(a, a_rev)    (a)->sister = (a_rev); (a_rev)->sister = (a);
#define SET_FROM(a, i)           (a)->next = (i)->first; (i)->first = (a);
#define REMOVE_FROM(a, i)        if ((i)->first==(a)) (i)->first=(a)->next;\
								 else { Arc* a_TMP_REMOVE_FROM; for (a_TMP_REMOVE_FROM=i->first; ; a_TMP_REMOVE_FROM=a_TMP_REMOVE_FROM->next)\
								                 if (a_TMP_REMOVE_FROM->next==(a)) { a_TMP_REMOVE_FROM->next=(a)->next; break; } }
#define SET_TO(a, j)             (a)->head = (j);


template <typename REAL> 
	typename QPBO<REAL>::EdgeId QPBO<REAL>::AddPairwiseTerm(NodeId _i, NodeId _j, REAL E00, REAL E01, REAL E10, REAL E11)
{
	user_assert(_i >= 0 && _i < node_num);
	user_assert(_j >= 0 && _j < node_num);
	user_assert(_i != _j);

	REAL ci, cj, cij, cji;

	if (!first_free) 
	{
		reallocate_arcs(2*(GetMaxEdgeNum() + GetMaxEdgeNum()/2));
	}

	EdgeId e = (int)(first_free - arcs[IsArc0(first_free) ? 0 : 1])/2;
	first_free = first_free->next;

	if (stage == 0)
	{
		Arc *a, *a_rev;
		a     = &arcs[0][2*e];
		a_rev = &arcs[0][2*e+1];

		Node* i = nodes[0] + _i;
		Node* j = nodes[0] + _j;

		if (E01 + E10 >= E00 + E11)
		{
			ComputeWeights(E00, E01, E10, E11, ci, cj, cij, cji);

			SET_TO(a, j);
			SET_FROM(a,     i);
			SET_FROM(a_rev, j);

			j->tr_cap += cj;
		}
		else
		{
			all_edges_submodular = false;
			ComputeWeights(E01, E00, E11, E10, ci, cj, cij, cji);

			SET_TO(a, GetMate0(j));
			a->next = NULL;
			a_rev->next = NULL;

			j->tr_cap -= cj;
		}

		SET_SISTERS(a, a_rev);
		SET_TO(a_rev, i);

		i->tr_cap += ci;
		a->r_cap = cij;
		a_rev->r_cap = cji;
	}
	else
	{
		Arc *a[2], *a_rev[2];
		a[0]     = &arcs[0][2*e];
		a_rev[0] = &arcs[0][2*e+1];
		a[1]     = &arcs[1][2*e];
		a_rev[1] = &arcs[1][2*e+1];

		Node* i[2] = { nodes[0] + _i, nodes[1] + _i };
		Node* j[2];

		if (E01 + E10 >= E00 + E11)
		{
			j[0] = nodes[0] + _j; j[1] = nodes[1] + _j;
			ComputeWeights(E00, E01, E10, E11, ci, cj, cij, cji);
		}
		else
		{
			j[1] = nodes[0] + _j; j[0] = nodes[1] + _j;
			ComputeWeights(E01, E00, E11, E10, ci, cj, cij, cji);
		}

		SET_SISTERS(a[0], a_rev[0]);
		SET_SISTERS(a[1], a_rev[1]);

		SET_TO(a[0],     j[0]);
		SET_TO(a_rev[0], i[0]);
		SET_TO(a[1],     i[1]);
		SET_TO(a_rev[1], j[1]);

		SET_FROM(a[0],     i[0]);
		SET_FROM(a_rev[0], j[0]);
		SET_FROM(a[1],     j[1]);
		SET_FROM(a_rev[1], i[1]);

		i[0]->tr_cap += ci; i[1]->tr_cap -= ci;
		j[0]->tr_cap += cj; j[1]->tr_cap -= cj;
		a[0]->r_cap = a[1]->r_cap = cij;
		a_rev[0]->r_cap = a_rev[1]->r_cap = cji;
	}

	zero_energy += E00;

	return e;
}

template <typename REAL> 
	void QPBO<REAL>::AddPairwiseTerm(EdgeId e, NodeId _i, NodeId _j, REAL E00, REAL E01, REAL E10, REAL E11)
{
	user_assert(e >= 0 && arcs[0][2*e].sister);
	user_assert(arcs[0][2*e].head==&nodes[0][_i] || arcs[0][2*e].head==&nodes[1][_i] || arcs[0][2*e].head==&nodes[0][_j] || arcs[0][2*e].head==&nodes[1][_j]);
	user_assert(arcs[0][2*e+1].head==&nodes[0][_i] || arcs[0][2*e+1].head==&nodes[1][_i] || arcs[0][2*e+1].head==&nodes[0][_j] || arcs[0][2*e+1].head==&nodes[1][_j]);
	user_assert(_i != _j);

	REAL delta, ci, cj, cij, cji;

	if (stage == 0)
	{
		Arc* a = &arcs[0][2*e];
		Arc* a_rev = &arcs[0][2*e+1];
		code_assert(a->sister==a_rev && a->sister==a_rev);

		Node* i = a_rev->head;
		Node* j = a->head;
		code_assert(IsNode0(i));
		if (i != &nodes[0][_i]) { delta = E01; E01 = E10; E10 = delta; }
		if (IsNode0(j))
		{
			ComputeWeights(E00, E01, E10, E11, ci, cj, cij, cji);
			
			i->tr_cap += ci;
			j->tr_cap += cj;
			a->r_cap += cij;
			a_rev->r_cap += cji;

			if (a->r_cap < 0)
			{
				delta = a->r_cap;
				a->r_cap = 0;
				a_rev->r_cap += delta;
				i->tr_cap -= delta;
				j->tr_cap += delta;
			}
			if (a_rev->r_cap < 0)
			{
				delta = a_rev->r_cap;
				a_rev->r_cap = 0;
				a->r_cap += delta;
				j->tr_cap -= delta;
				i->tr_cap += delta;
			}

			if (a->r_cap < 0)
			{
				all_edges_submodular = false;
				REMOVE_FROM(a, i);
				REMOVE_FROM(a_rev, j);
				SET_TO(a, GetMate0(j));

				delta = a->r_cap;
				i->tr_cap -= delta;
				a->r_cap = -delta;
			}
		}
		else
		{
			j = GetMate1(j);
			ComputeWeights(E01, E00, E11, E10, ci, cj, cij, cji);
			
			i->tr_cap += ci;
			j->tr_cap -= cj;
			a->r_cap += cij;
			a_rev->r_cap += cji;

			if (a->r_cap < 0)
			{
				delta = a->r_cap;
				a->r_cap = 0;
				a_rev->r_cap += delta;
				i->tr_cap -= delta;
				j->tr_cap -= delta;
			}
			if (a_rev->r_cap < 0)
			{
				delta = a_rev->r_cap;
				a_rev->r_cap = 0;
				a->r_cap += delta;
				j->tr_cap += delta;
				i->tr_cap += delta;
			}

			if (a->r_cap < 0)
			{
				SET_FROM(a, i);
				SET_FROM(a_rev, j);
				SET_TO(a, j);

				delta = a->r_cap;
				i->tr_cap -= delta;
				a->r_cap = -delta;
			}
		}
	}
	else
	{
		Arc* a[2] = { &arcs[0][2*e], &arcs[1][2*e] };
		Arc* a_rev[2] = { &arcs[0][2*e+1], &arcs[1][2*e+1] };
		code_assert(a[0]->sister==a_rev[0] && a[1]->sister==a_rev[1] && a[0]==a_rev[0]->sister && a[1]==a_rev[1]->sister);

		Node* i[2] = { a_rev[0]->head, a[1]->head };
		Node* j[2] = { a[0]->head, a_rev[1]->head };
		int k = IsNode0(i[0]) ? 0 : 1;
		if (i[k] != &nodes[0][_i]) { delta = E01; E01 = E10; E10 = delta; }
		if (IsNode0(j[k]))
		{ 
			ComputeWeights(E00, E01, E10, E11, ci, cj, cij, cji);
		}
		else
		{ 
			ComputeWeights(E01, E00, E11, E10, ci, cj, cij, cji);
		};

		// make sure that a[0]->r_cap == a[1]->r_cap and a_rev[0]->r_cap == a_rev[1]->r_cap by pushing flow
		delta = a[1]->r_cap - a[0]->r_cap;
		//a[1]->r_cap -= delta;   // don't do the subtraction - later we'll set explicitly a[1]->r_cap = a[0]->r_cap
		//a[1]->sister->r_cap += delta;
		a_rev[1]->head->tr_cap -= delta;
		a[1]->head->tr_cap     += delta;

		i[0]->tr_cap += ci; i[1]->tr_cap -= ci;
		j[0]->tr_cap += cj; j[1]->tr_cap -= cj;
		a[0]->r_cap += cij;
		a_rev[0]->r_cap += cji;

		if (a[0]->r_cap < 0)
		{
			delta = a[0]->r_cap;
			a[0]->r_cap = 0;
			a_rev[0]->r_cap += delta;
			i[0]->tr_cap -= delta; i[1]->tr_cap += delta;
			j[0]->tr_cap += delta; j[1]->tr_cap -= delta;
		}
		if (a_rev[0]->r_cap < 0)
		{
			delta = a_rev[0]->r_cap;
			a_rev[0]->r_cap = 0;
			a[0]->r_cap += delta;
			j[0]->tr_cap -= delta; j[1]->tr_cap += delta;
			i[0]->tr_cap += delta; i[1]->tr_cap -= delta;
		}

		if (a[0]->r_cap < 0)
		{
			// need to swap submodular <-> supermodular
			SET_TO(a[0], j[1]);
			SET_TO(a_rev[1], j[0]);
			REMOVE_FROM(a_rev[0], j[0]);
			SET_FROM(a_rev[0], j[1]);
			REMOVE_FROM(a[1], j[1]);
			SET_FROM(a[1], j[0]);

			delta = a[0]->r_cap;
			i[0]->tr_cap -= delta; i[1]->tr_cap += delta;
			a[0]->r_cap = -delta;
		}

		a[1]->r_cap = a[0]->r_cap;
		a_rev[1]->r_cap = a_rev[0]->r_cap;
	}

	zero_energy += E00;
}

template <typename REAL> 
	void QPBO<REAL>::TransformToSecondStage(bool copy_trees)
{
	// add non-submodular edges
	Node* i[2];
	Node* j[2];
	Arc* a[2];

	memset(nodes[1], 0, node_num*sizeof(Node));
	node_last[1] = nodes[1] + node_num;

	if (!copy_trees)
	{
		for (i[0]=nodes[0], i[1]=nodes[1]; i[0]<node_last[0]; i[0]++, i[1]++)
		{
			i[1]->first = NULL;
			i[1]->tr_cap = -i[0]->tr_cap;
		}

		for (a[0]=arcs[0], a[1]=arcs[1]; a[0]<arc_max[0]; a[0]+=2, a[1]+=2)
		{
			if (!a[0]->sister) continue;

			code_assert(IsNode0(a[0]->sister->head));
			SET_SISTERS(a[1], a[1]+1);
			if (IsNode0(a[0]->head))
			{
				i[1] = GetMate0(a[0]->sister->head);
				j[1] = GetMate0(a[0]->head);

				SET_FROM(a[1],         j[1]);
				SET_FROM(a[1]->sister, i[1]);
				SET_TO(a[1],         i[1]);
				SET_TO(a[1]->sister, j[1]);
			}
			else
			{
				i[0] = a[0]->sister->head;
				i[1] = GetMate0(i[0]);
				j[1] = a[0]->head;
				j[0] = GetMate1(j[1]);

				SET_FROM(a[0],         i[0]);
				SET_FROM(a[0]->sister, j[1]);
				SET_FROM(a[1],         j[0]);
				SET_FROM(a[1]->sister, i[1]);
				SET_TO(a[1],         i[1]);
				SET_TO(a[1]->sister, j[0]);
			}
			a[1]->r_cap = a[0]->r_cap;
			a[1]->sister->r_cap = a[0]->sister->r_cap;
		}
	}
	else
	{
		for (i[0]=nodes[0], i[1]=nodes[1]; i[0]<node_last[0]; i[0]++, i[1]++)
		{
			i[1]->first = NULL;
			i[1]->tr_cap = -i[0]->tr_cap;
			i[1]->is_sink = i[0]->is_sink ^ 1;
			i[1]->DIST = i[0]->DIST;
			i[1]->TS = i[0]->TS;

			if (i[0]->parent == NULL || i[0]->parent == QPBO_MAXFLOW_TERMINAL) i[1]->parent = i[0]->parent;
			else i[1]->parent = GetMate0(i[0]->parent->sister);
		}

		for (a[0]=arcs[0], a[1]=arcs[1]; a[0]<arc_max[0]; a[0]+=2, a[1]+=2)
		{
			if (!a[0]->sister) continue;

			code_assert(IsNode0(a[0]->sister->head));
			SET_SISTERS(a[1], a[1]+1);
			if (IsNode0(a[0]->head))
			{
				i[1] = GetMate0(a[0]->sister->head);
				j[1] = GetMate0(a[0]->head);

				SET_FROM(a[1],         j[1]);
				SET_FROM(a[1]->sister, i[1]);
				SET_TO(a[1],         i[1]);
				SET_TO(a[1]->sister, j[1]);
			}
			else
			{
				i[0] = a[0]->sister->head;
				i[1] = GetMate0(i[0]);
				j[1] = a[0]->head;
				j[0] = GetMate1(j[1]);

				SET_FROM(a[0],         i[0]);
				SET_FROM(a[0]->sister, j[1]);
				SET_FROM(a[1],         j[0]);
				SET_FROM(a[1]->sister, i[1]);
				SET_TO(a[1],         i[1]);
				SET_TO(a[1]->sister, j[0]);

				mark_node(i[0]);
				mark_node(i[1]);
				mark_node(j[0]);
				mark_node(j[1]);
			}
			a[1]->r_cap = a[0]->r_cap;
			a[1]->sister->r_cap = a[0]->sister->r_cap;
		}
	}

	stage = 1;
}

template <typename REAL>
	void QPBO<REAL>::MergeParallelEdges()
{
	if (stage == 0) TransformToSecondStage(false);
	Node* i;
	Node* j;
	Arc* a;
	Arc* a_next;

	for (i=nodes[0]; i<node_last[0]; i++)
	{
		for (a=i->first; a; a=a->next)
		{
			j = a->head;
			if (!IsNode0(j)) j = GetMate1(j);
			j->parent = a;
		}
		for (a=i->first; a; a=a_next)
		{
			a_next = a->next;
			j = a->head;
			if (!IsNode0(j)) j = GetMate1(j);
			if (j->parent == a) continue;
			if (MergeParallelEdges(j->parent, a)==0) 
			{
				j->parent = a;
				a_next = a->next;
			}
		}
	}
}

template <typename REAL>
	void QPBO<REAL>::Solve()
{
	Node* i;

	maxflow();

	if (stage == 0)
	{
		if (all_edges_submodular)
		{
			for (i=nodes[0]; i<node_last[0]; i++)
			{
				i->label = what_segment(i);
			}
			return;
		}

		TransformToSecondStage(true);
		maxflow(true);
	}

	for (i=nodes[0]; i<node_last[0]; i++)
	{
		i->label = what_segment(i);
		if (i->label == what_segment(GetMate0(i))) i->label = -1;
	}
}

template <typename REAL>
	REAL QPBO<REAL>::ComputeTwiceEnergy(int option)
{
	REAL E = 2*zero_energy, E1[2], E2[2][2];
	int i, j, e;
	int xi, xj;

	for (i=0; i<GetNodeNum(); i++)
	{
		GetTwiceUnaryTerm(i, E1[0], E1[1]);
		if (option == 0) xi = (nodes[0][i].label < 0) ? 0 : nodes[0][i].label;
		else             xi = nodes[0][i].user_label;
		code_assert(xi==0 || xi==1);
		E += E1[xi] - E1[0];
	}
	for (e=GetNextEdgeId(-1); e>=0; e=GetNextEdgeId(e))
	{
		GetTwicePairwiseTerm(e, i, j, E2[0][0], E2[0][1], E2[1][0], E2[1][1]);
		if (option == 0)
		{
			xi = (nodes[0][i].label < 0) ? 0 : nodes[0][i].label;
			xj = (nodes[0][j].label < 0) ? 0 : nodes[0][j].label;
		}
		else
		{
			xi = nodes[0][i].user_label;
			xj = nodes[0][j].user_label;
		}
		E += E2[xi][xj] - E2[0][0];
	}
	return E;
}

template <typename REAL>
	REAL QPBO<REAL>::ComputeTwiceEnergy(int* solution)
{
	REAL E = 2*zero_energy, E1[2], E2[2][2];
	int i, j, e;

	for (i=0; i<GetNodeNum(); i++)
	{
		GetTwiceUnaryTerm(i, E1[0], E1[1]);
		if (solution[i] == 1) E += E1[1];
	}
	for (e=GetNextEdgeId(-1); e>=0; e=GetNextEdgeId(e))
	{
		GetTwicePairwiseTerm(e, i, j, E2[0][0], E2[0][1], E2[1][0], E2[1][1]);
		E += E2[(solution[i] == 1) ? 1 : 0][(solution[j] == 1) ? 1 : 0] - E2[0][0];
	}
	return E;
}

template <typename REAL>
	REAL QPBO<REAL>::ComputeTwiceLowerBound()
{
	REAL lowerBound = 2*zero_energy, E0, E1, E00, E01, E10, E11;
	int i, j, e;

	for (i=0; i<GetNodeNum(); i++)
	{
		GetTwiceUnaryTerm(i, E0, E1);
		if (E0 > E1) lowerBound += E1 - E0;
	}
	for (e=GetNextEdgeId(-1); e>=0; e=GetNextEdgeId(e))
	{
		GetTwicePairwiseTerm(e, i, j, E00, E01, E10, E11);
		lowerBound -= E00;
	}

	return lowerBound;
}

template <typename REAL>
	void QPBO<REAL>::TestRelaxedSymmetry()
{
	Node* i;
	Arc* a;
	REAL c1, c2;

	if (stage == 0) return;

	for (i=nodes[0]; i<node_last[0]; i++)
	{
		if (i->is_removed) continue;
		c1 = i->tr_cap;
		for (a=i->first; a; a=a->next) c1 += a->sister->r_cap;
		c2 = -GetMate0(i)->tr_cap;
		for (a=GetMate0(i)->first; a; a=a->next) c2 += a->r_cap;
		if (c1 != c2)
		{
			code_assert(0);
			exit(1);
		}
	}
}

/* QPBO_extra.cpp */

template <typename REAL>
	void QPBO<REAL>::ComputeRandomPermutation(int N, int* permutation)
{
	int i, j, k;
	for (i=0; i<N; i++)
	{
		permutation[i] = i;
	}
	for (i=0; i<N-1; i++)
	{
		j = i + (int)((rand()/(1.0+(double)RAND_MAX))*(N-i));
		if (j>N-1) j = N-1;
		k = permutation[j]; permutation[j] = permutation[i]; permutation[i] = k;
	}
}

template <typename REAL>
	void QPBO<REAL>::MergeMappings(int nodeNum0, int* mapping0, int* mapping1)
{
	int i;
	for (i=0; i<nodeNum0; i++)
	{
		int j = mapping0[i] / 2;
		int k = mapping1[j] / 2;
		mapping0[i] = 2*k + ((mapping0[i] + mapping1[j]) % 2);
	}
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#define SET_SISTERS(a, a_rev)    (a)->sister = (a_rev); (a_rev)->sister = (a);
#define SET_FROM(a, i)           (a)->next = (i)->first; (i)->first = (a);
#define REMOVE_FROM(a, i)        if ((i)->first==(a)) (i)->first=(a)->next;\
								 else { Arc* a_TMP_REMOVE_FROM; for (a_TMP_REMOVE_FROM=i->first; ; a_TMP_REMOVE_FROM=a_TMP_REMOVE_FROM->next)\
								                 if (a_TMP_REMOVE_FROM->next==(a)) { a_TMP_REMOVE_FROM->next=(a)->next; break; } }
#define SET_TO(a, j)             (a)->head = (j);


template <typename REAL>
	inline void QPBO<REAL>::FixNode(Node* i, int x)
{
	Node* _i[2] = { i, GetMate0(i) };
	Arc* a;
	Arc* a_next;

	for (a=_i[x]->first; a; a=a->next)
	{
		mark_node(a->head);
		a->head->tr_cap += a->r_cap;
		REMOVE_FROM(a->sister, a->head);
		a->sister->sister = NULL;
		a->sister = NULL;
	}
	for (a=_i[1-x]->first; a; a=a_next)
	{
		mark_node(a->head);
		a->head->tr_cap -= a->sister->r_cap;
		REMOVE_FROM(a->sister, a->head);
		a->sister->sister = NULL;
		a->sister = NULL;

		a_next = a->next;
		a->next = first_free;
		first_free = a;
	}
	_i[0]->first = _i[1]->first = NULL;
}

template <typename REAL>
	inline void QPBO<REAL>::ContractNodes(Node* i, Node* j, int swap)
{
	code_assert(IsNode0(i) && IsNode0(j) && swap>=0 && swap<=1);

	Node* _i[2] = { i, GetMate0(i) };
	Node* _j[2];
	Arc* a;
	Arc* a_selfloop = NULL;
	int x;

	if (swap == 0) { _j[0] = j; _j[1] = GetMate0(j); }
	else           { _j[1] = j; _j[0] = GetMate0(j); }

	_i[0]->tr_cap += _j[0]->tr_cap;
	_i[1]->tr_cap += _j[1]->tr_cap;
	for (x=0; x<2; x++)
	{
		Arc* a_next;
		for (a=_j[x]->first; a; a=a_next)
		{
			mark_node(a->head);
			a_next = a->next;
			if (a->head == _i[x])
			{
				REMOVE_FROM(a->sister, _i[x]);
				a->sister->sister = NULL;
				a->sister = NULL;
				a_selfloop = a;
			}
			else if (a->head == _i[1-x])
			{
				REMOVE_FROM(a->sister, _i[1-x]);
				_i[x]->tr_cap   -= a->r_cap;
				_i[1-x]->tr_cap += a->r_cap;
				a->sister->sister = NULL;
				a->sister = NULL;
			}
			else
			{
				SET_FROM(a, _i[x]);
				SET_TO(a->sister, _i[x]);
			}
		}
	}
	_j[0]->first = _j[1]->first = NULL;

	if (a_selfloop)
	{
		a_selfloop->next = first_free;
		first_free = a_selfloop;
	}
}

template <typename REAL>
	int QPBO<REAL>::MergeParallelEdges(Arc* a1, Arc* a2)
{
	code_assert(a1->sister->head == a2->sister->head && IsNode0(a1->sister->head));
	code_assert(a1->head == a2->head || a1->head == GetMate(a2->head));

	REAL delta;
	int x;
	Node* _i[2];
	Node* _j[2];
	Arc* _a1[2] = { a1, GetMate(a1) };
	Arc* _a2[2] = { a2, GetMate(a2) };
	_i[0] = a1->sister->head; _i[1] = GetMate0(_i[0]);

	if (a1->head == a2->head)
	{
		a1->r_cap += a2->r_cap;
		a1->sister->r_cap += a2->sister->r_cap;
		_a1[1]->r_cap += _a2[1]->r_cap;
		_a1[1]->sister->r_cap += _a2[1]->sister->r_cap;
		x = 1;

		_j[0] = a1->head;
		_j[1] = GetMate(_j[0]);
	}
	else
	{
		code_assert(a1->head == GetMate(a2->head));
		// make sure that _a1[0]->r_cap == _a1[1]->r_cap and _a1[0]->sister->r_cap == _a1[1]->sister->r_cap by pushing flow
		delta = _a1[1]->r_cap - _a1[0]->r_cap;
		//_a1[1]->r_cap -= delta;   // don't do the subtraction - later we'll set explicitly _a1[1]->r_cap = _a1[0]->r_cap
		//_a1[1]->sister->r_cap += delta;
		_a1[1]->sister->head->tr_cap -= delta;
		_a1[1]->head->tr_cap         += delta;
		// same for a2
		delta = _a2[1]->r_cap - _a2[0]->r_cap;
		//_a2[1]->r_cap -= delta;   // don't do the subtraction - later we'll set explicitly _a2[1]->r_cap = _a2[0]->r_cap
		//_a2[1]->sister->r_cap += delta;
		_a2[1]->sister->head->tr_cap -= delta;
		_a2[1]->head->tr_cap         += delta;

		if (a1->r_cap + a1->sister->r_cap >= a2->r_cap + a2->sister->r_cap) x = 1;
		else // swap a1 <-> a2
		{		
			Arc* tmp;
			tmp = a1; a1 = a2; a2 = tmp; 
			_a1[0] = a1; 
			_a2[0] = a2;
			tmp = _a1[1]; _a1[1] = _a2[1]; _a2[1] = tmp;
			x = 0;
		}

		_j[0] = a1->head;
		_j[1] = a2->head;

		REAL ci, cj, cij, cji;
		ci = a2->sister->r_cap - a2->r_cap;
		cj = 0;
		cij = -a2->r_cap;
		cji = -a2->sister->r_cap;

		_i[0]->tr_cap += ci; _i[1]->tr_cap -= ci;
		_j[0]->tr_cap += cj; _j[1]->tr_cap -= cj;
		_a1[0]->r_cap += cij;
		_a1[0]->sister->r_cap += cji;

		if (_a1[0]->r_cap < 0)
		{
			delta = _a1[0]->r_cap;
			_a1[0]->r_cap = 0;
			_a1[0]->sister->r_cap += delta;
			_i[0]->tr_cap -= delta; _i[1]->tr_cap += delta;
			_j[0]->tr_cap += delta; _j[1]->tr_cap -= delta;
		}
		if (_a1[0]->sister->r_cap < 0)
		{
			delta = _a1[0]->sister->r_cap;
			_a1[0]->sister->r_cap = 0;
			_a1[0]->r_cap += delta;
			_j[0]->tr_cap -= delta; _j[1]->tr_cap += delta;
			_i[0]->tr_cap += delta; _i[1]->tr_cap -= delta;
		}

		code_assert (_a1[0]->r_cap >= 0);

		_a1[1]->r_cap = _a1[0]->r_cap;
		_a1[1]->sister->r_cap = _a1[0]->sister->r_cap;
	}

	REMOVE_FROM(a2, _i[0]);
	REMOVE_FROM(a2->sister, a2->head);
	REMOVE_FROM(_a2[1], _a2[1]->sister->head);
	REMOVE_FROM(_a2[1]->sister, _i[1]);
	a2->sister->sister = NULL;
	a2->sister = NULL;
	_a2[1]->sister->sister = NULL;
	_a2[1]->sister = NULL;

	_a2[1]->next = first_free;
	first_free = _a2[1];

	return x;
}

template <typename REAL>
	inline REAL QPBO<REAL>::DetermineSaturation(Node* i)
{
	Arc* a;
	REAL c1 = -i->tr_cap;
	REAL c2 = i->tr_cap;
	for (a=i->first; a; a=a->next)
	{
		c1 += a->r_cap;
		c2 += a->sister->r_cap;
	}

	return (c1 > c2) ? c1 : c2;
}

template <typename REAL>
	inline void QPBO<REAL>::AddDirectedConstraint(Node* i, Node* j, int xi, int xj)
{
	code_assert(first_free && IsNode0(i) && IsNode0(j) && i!=j);

	int e = ((int)(first_free - arcs[IsArc0(first_free) ? 0 : 1])) & (~1);
	first_free = first_free->next;

	Arc* _a[2] = { &arcs[0][e], &arcs[1][e] };
	Node* _i[2] = { i, GetMate0(i) };
	Node* _j[2];
	if (xi==xj) { _j[0] = j; _j[1] = GetMate0(j); }
	else        { _j[1] = j; _j[0] = GetMate0(j); }

	SET_SISTERS(_a[0], _a[0]+1);
	SET_SISTERS(_a[1], _a[1]+1);

	SET_FROM(_a[0], _i[0]);
	SET_TO(_a[0], _j[0]);
	SET_FROM(_a[0]->sister, _j[0]);
	SET_TO(_a[0]->sister, _i[0]);
	SET_FROM(_a[1], _j[1]);
	SET_TO(_a[1], _i[1]);
	SET_FROM(_a[1]->sister, _i[1]);
	SET_TO(_a[1]->sister, _j[1]);

	if (xi==0) { _a[0]->r_cap = probe_options.C; _a[0]->sister->r_cap = 0; }
	else       { _a[0]->r_cap = 0; _a[0]->sister->r_cap = probe_options.C; }

	_a[1]->r_cap = _a[0]->r_cap;
	_a[1]->sister->r_cap = _a[0]->sister->r_cap;
}
/*
template <typename REAL>
	inline bool QPBO<REAL>::AddDirectedConstraint(Arc* a, int xi, int xj)
{
	Node* i = a->sister->head;
	Node* j = a->head;
	Node* _i[2] = { i, GetMate0(i) };
	Node* _j[2];
	Arc* _a[2] = { a, GetMate(a) };
	int x; // 0 if current edge is submodular, 1 otherwise
	REAL delta;

	_j[0] = j;
	if (IsNode0(j)) { _j[1] = GetMate0(j); x = 0; }
	else            { _j[1] = GetMate1(j); x = 1; }

	if ((xi + xj + x)%2 == 0)
	{
		// easy case - graph structure doesn't need to be changed
		if (xi == 0)
		{
			if (a->r_cap + _a[1]->r_cap >= 2*probe_options.C) return false;
			mark_node(_j[0]);
			mark_node(_j[1]);
			a->r_cap += 2*probe_options.C; 
			_a[1]->r_cap += 2*probe_options.C;
			return true;
		}
		else
		{
			if (a->sister->r_cap + _a[1]->sister->r_cap >= 2*probe_options.C) return false;
			mark_node(_j[0]);
			mark_node(_j[1]);
			a->sister->r_cap += 2*probe_options.C; 
			_a[1]->sister->r_cap += 2*probe_options.C;
			return true;
		}
	}

	mark_node(_j[0]);
	mark_node(_j[1]);

	// make sure that _a[0]->r_cap == _a[1]->r_cap and _a[0]->sister->r_cap == _a[1]->sister->r_cap by pushing flow
	delta = _a[1]->r_cap - _a[0]->r_cap;
	//_a[1]->r_cap -= delta;   // don't do the subtraction - later we'll set explicitly _a[1]->r_cap = _a[0]->r_cap
	//_a[1]->sister->r_cap += delta;
	_a[1]->sister->head->tr_cap -= delta;
	_a[1]->head->tr_cap         += delta;

	SET_TO(_a[0], _j[1]);
	SET_TO(_a[1]->sister, _j[0]);
	REMOVE_FROM(_a[0]->sister, _j[0]);
	SET_FROM(_a[0]->sister, _j[1]);
	REMOVE_FROM(_a[1], _j[1]);
	SET_FROM(_a[1], _j[0]);

	i->tr_cap += a->sister->r_cap - a->r_cap; _i[1]->tr_cap -= a->sister->r_cap - a->r_cap;
	a->r_cap = -a->r_cap;

	if (xi == 0) a->r_cap += 2*probe_options.C; 
	else         a->sister->r_cap += 2*probe_options.C;

	if (a->r_cap < 0)
	{
		delta = a->r_cap;
		a->r_cap = 0;
		a->sister->r_cap += delta;
		i->tr_cap -= delta; _i[1]->tr_cap += delta;
		_j[1]->tr_cap += delta; j->tr_cap -= delta;
	}
	if (a->sister->r_cap < 0)
	{
		delta = a->sister->r_cap;
		a->sister->r_cap = 0;
		a->r_cap += delta;
		_j[1]->tr_cap -= delta; j->tr_cap += delta;
		i->tr_cap += delta; _i[1]->tr_cap -= delta;
	}

	_a[1]->r_cap = a->r_cap;
	_a[1]->sister->r_cap = a->sister->r_cap;

	return true;
}
*/
template <typename REAL>
	inline bool QPBO<REAL>::AddDirectedConstraint0(Arc* a, int xi, int xj)
{
	Node* i = a->sister->head;
	Node* j = a->head;
	Node* _i[2] = { i, GetMate0(i) };
	Node* _j[2];
	Arc* _a[2] = { a, GetMate(a) };
	int x; // 0 if current edge is submodular, 1 otherwise
	REAL delta;

	_j[0] = j;
	if (IsNode0(j)) { _j[1] = GetMate0(j); x = 0; }
	else            { _j[1] = GetMate1(j); x = 1; }

	if ((xi + xj + x)%2 == 0)
	{
		// easy case - graph structure doesn't need to be changed
		if (a->r_cap + a->sister->r_cap + _a[1]->r_cap + _a[1]->sister->r_cap >= 2*probe_options.C) return false;
		mark_node(_j[0]);
		mark_node(_j[1]);
		if (xi == 0)
		{
			a->r_cap += probe_options.C; 
			_a[1]->r_cap += probe_options.C;
		}
		else
		{
			a->sister->r_cap += probe_options.C; 
			_a[1]->sister->r_cap += probe_options.C;
		}
		return true;
	}

	mark_node(_j[0]);
	mark_node(_j[1]);

	// make sure that _a[0]->r_cap == _a[1]->r_cap and _a[0]->sister->r_cap == _a[1]->sister->r_cap by pushing flow
	delta = _a[1]->r_cap - _a[0]->r_cap;
	//_a[1]->r_cap -= delta;   // don't do the subtraction - later we'll set explicitly _a[1]->r_cap = _a[0]->r_cap
	//_a[1]->sister->r_cap += delta;
	_a[1]->sister->head->tr_cap -= delta;
	_a[1]->head->tr_cap         += delta;

	SET_TO(_a[0], _j[1]);
	SET_TO(_a[1]->sister, _j[0]);
	REMOVE_FROM(_a[0]->sister, _j[0]);
	SET_FROM(_a[0]->sister, _j[1]);
	REMOVE_FROM(_a[1], _j[1]);
	SET_FROM(_a[1], _j[0]);

	i->tr_cap += a->sister->r_cap - a->r_cap; _i[1]->tr_cap -= a->sister->r_cap - a->r_cap;
	a->r_cap = -a->r_cap;

	if (xi == 0) a->r_cap         += probe_options.C + a->sister->r_cap - a->r_cap;
	else         a->sister->r_cap += probe_options.C + a->sister->r_cap - a->r_cap;

	if (a->r_cap < 0)
	{
		delta = a->r_cap;
		a->r_cap = 0;
		a->sister->r_cap += delta;
		i->tr_cap -= delta; _i[1]->tr_cap += delta;
		_j[1]->tr_cap += delta; j->tr_cap -= delta;
	}
	if (a->sister->r_cap < 0)
	{
		delta = a->sister->r_cap;
		a->sister->r_cap = 0;
		a->r_cap += delta;
		_j[1]->tr_cap -= delta; j->tr_cap += delta;
		i->tr_cap += delta; _i[1]->tr_cap -= delta;
	}

	_a[1]->r_cap = a->r_cap;
	_a[1]->sister->r_cap = a->sister->r_cap;

	return true;
}

template <typename REAL>
	inline bool QPBO<REAL>::AddDirectedConstraint1(Arc* a, int xi, int xj)
{
	Node* j = a->head;
	Node* _j[2];
	Arc* _a[2] = { a, GetMate(a) };
	int x; // 0 if current edge is submodular, 1 otherwise

	_j[0] = j;
	if (IsNode0(j)) { _j[1] = GetMate0(j); x = 0; }
	else            { _j[1] = GetMate1(j); x = 1; }

	code_assert((xi + xj + x)%2 == 0);

	if (xi == 0)
	{
		if (a->r_cap > 0 && _a[1]->r_cap > 0) return false;
		mark_node(_j[0]);
		mark_node(_j[1]);
		a->r_cap += probe_options.C; 
		_a[1]->r_cap += probe_options.C;
		return true;
	}
	else
	{
		if (a->sister->r_cap > 0 && _a[1]->sister->r_cap > 0) return false;
		mark_node(_j[0]);
		mark_node(_j[1]);
		a->sister->r_cap += probe_options.C; 
		_a[1]->sister->r_cap += probe_options.C;
		return true;
	}
}

template <typename REAL>
	void QPBO<REAL>::AllocateNewEnergy(int* mapping)
{
	int i_index, j_index;
	int nodeNumOld = GetNodeNum();
	int nodeNumNew = 1;
	int edgeNumOld = GetMaxEdgeNum();
	int e;
	Node* i;

	////////////////////////////////////////////////////////////////
	for (i_index=0, i=nodes[0]; i_index<nodeNumOld; i_index++, i++)
	{
		if (mapping[i_index] < 0)
		{
			mapping[i_index] = 2*nodeNumNew + ((i->user_label) % 2);
			nodeNumNew ++;
		}
		else if (mapping[i_index]>=2) mapping[i_index] = -mapping[i_index];
	}

	////////////////////////////////////////////////////////////////
	code_assert(nodes[0] + nodeNumNew <= node_max[0]);
	// Reset:
	node_last[0] = nodes[0];
	node_last[1] = nodes[1];
	node_num = 0;

	if (nodeptr_block) 
	{ 
		delete nodeptr_block; 
		nodeptr_block = NULL; 
	}
	if (changed_list)
	{
		delete changed_list;
		changed_list = NULL;
	}
	if (fix_node_info_list)
	{
		delete fix_node_info_list;
		fix_node_info_list = NULL;
	}

	maxflow_iteration = 0;
	zero_energy = 0;

	stage = 0;
	all_edges_submodular = true;
	////////////////////////////////////////////////////////////////


	AddNode(nodeNumNew);
	AddUnaryTerm(0, (REAL)0, (REAL)1);

	i = nodes[0];
	i->user_label = i->label = 0;
	for (i_index=0; i_index<nodeNumOld; i_index++)
	{
		if (mapping[i_index] >= 2) 
		{
			i = nodes[0] + (mapping[i_index]/2);
			i->user_label = i->label = mapping[i_index] & 1;
			mapping[i_index] &= ~1;
		}
	}
	////////////////////////////////////////////////////////////////
	for (i_index=0; i_index<nodeNumOld; i_index++)
	if (mapping[i_index] < 0)
	{
		int y[2];
		int x = 0;
		j_index = i_index;
		do
		{
			x = (x - mapping[j_index]) % 2;
			j_index = (-mapping[j_index])/2 - 1;
		} while (mapping[j_index] < 0);
		y[x] = mapping[j_index];
		y[1-x] = mapping[j_index] ^ 1;
		
		x = 0;
		j_index = i_index;
		do
		{
			int m_old = mapping[j_index];
			mapping[j_index] = y[x];
			x = (x - m_old) % 2;
			j_index = (-m_old)/2 - 1;
		} while (mapping[j_index] < 0);
	}

	////////////////////////////////////////////////////////////////
	int edgeNumNew = 0;
	for (e=0; e<edgeNumOld; e++)
	{
		if ( arcs[0][2*e].sister )
		{
			Arc* a;
			Arc* a_mate;
			if (IsNode0(arcs[0][2*e].sister->head))
			{
				a = &arcs[0][2*e];
				a_mate = &arcs[1][2*e];
			}
			else
			{
				a = &arcs[1][2*e+1];
				a_mate = &arcs[0][2*e+1];
			}
			i_index = mapping[(int)(a->sister->head - nodes[0])] / 2;
			code_assert(i_index > 0 && i_index < nodeNumNew);

			first_free = &arcs[0][2*edgeNumNew++];
			if (IsNode0(a->head))
			{
				j_index = mapping[(int)(a->head - nodes[0])] / 2;
				code_assert(j_index > 0 && j_index < nodeNumNew);
				AddPairwiseTerm(i_index, j_index, 
					0, a->r_cap+a_mate->r_cap, a->sister->r_cap+a_mate->sister->r_cap, 0);
			}
			else
			{
				j_index = mapping[(int)(a->head - nodes[1])] / 2;
				code_assert(j_index > 0 && j_index < nodeNumNew);
				AddPairwiseTerm(i_index, j_index, 
					a->r_cap+a_mate->r_cap, 0, 0, a->sister->r_cap+a_mate->sister->r_cap);
			}
		}
	}

	first_free = &arcs[0][2*edgeNumNew];
	memset(first_free, 0, (int)((char*)arc_max[0] - (char*)first_free));
	InitFreeList();
}

const int LIST_NUM = 4;
struct List // contains LIST_NUM lists containing integers 0,1,...,num-1. In the beginning, each integer is in list 1.
{
	List(int _num, int* order)
		: num(_num)
	{
		int i;
		prev = new int[num+LIST_NUM]; prev += LIST_NUM;
		next = new int[num+LIST_NUM]; next += LIST_NUM;
		if (order)
		{
			for (i=0; i<num; i++)
			{
				prev[order[i]] = (i==0)     ? -1 : order[i-1];
				next[order[i]] = (i==num-1) ? -1 : order[i+1];
			}
			prev[-1] = order[num-1];
			next[-1] = order[0];
		}
		else
		{
			for (i=0; i<num; i++)
			{
				prev[i] = i-1;
				next[i] = i+1;
			}
			prev[-1] = num-1;
			next[-1] = 0;
			next[num-1] = -1;
		}
		for (i=2; i<=LIST_NUM; i++)
		{
			prev[-i] = -i;
			next[-i] = -i;
		}
	}

	~List()
	{
		delete [] (prev - LIST_NUM);
		delete [] (next - LIST_NUM);
	}
	// i must be in the list
	void Remove(int i)
	{
		next[prev[i]] = next[i];
		prev[next[i]] = prev[i];
	}
	void Move(int i, int r_to) // moves node i to list r_to
	{
		code_assert (r_to>=1 && r_to<=LIST_NUM);
		Remove(i);
		next[i] = -r_to;
		prev[i] = prev[-r_to];
		next[prev[-r_to]] = i;
		prev[-r_to] = i;
	}
	void MoveList(int r_from, int r_to) // append list r_from to list r_to. List r_from becomes empty.
	{
		code_assert (r_from>=1 && r_from<=LIST_NUM && r_to>=1 && r_to<=LIST_NUM && r_from!=r_to);
		if (next[-r_from] < 0) return; // list r_from is empty
		prev[next[-r_from]] = prev[-r_to];
		next[prev[-r_to]] = next[-r_from];
		prev[-r_to] = prev[-r_from];
		next[prev[-r_from]] = -r_to;
		prev[-r_from] = next[-r_from] = -r_from;
	}
	// i must be in the list
	// (or -r, in which case the first element of list r is returned). 
	// Returns -1 if no more elements.
	int GetNext(int i) { return next[i]; }

private:
	int num;
	int* next;
	int* prev;
};

template <typename REAL>
	void QPBO<REAL>::SetMaxEdgeNum(int num)
{
	if (num > GetMaxEdgeNum()) reallocate_arcs(2*num);
}

template <typename REAL>
	bool QPBO<REAL>::Probe(int* mapping)
{
	int i_index, i_index_next, j_index;
	Node* i;
	Node* j;
	Node* _i[2];
	Arc* a;
	Arc* a_next;
	int x;
	Node** ptr;
	bool is_enough_memory = true;
	int unlabeled_num;

	if (probe_options.order_array) probe_options.order_seed = 0;
	if (probe_options.order_seed != 0)
	{
		srand(probe_options.order_seed);
		probe_options.order_array = new int[GetNodeNum()];
		ComputeRandomPermutation(GetNodeNum(), probe_options.order_array);
	}
	List list(GetNodeNum(), probe_options.order_array);
	if (probe_options.order_seed != 0)
	{
		delete [] probe_options.order_array;
	}

	if (node_last[0] == node_max[0]) reallocate_nodes((int)(node_last[0] - nodes[0]) + 1);

	all_edges_submodular = false;
	Solve();

	int MASK_CURRENT = 1;
	int MASK_NEXT = 2;

	unlabeled_num = 0;
	for (i=nodes[0], i_index=0; i<node_last[0]; i++, i_index++)
	{
		if (i->label >= 0)
		{
			FixNode(i, i->label);
			mapping[i_index] = i->label;
			i->is_removed = 1;
			list.Remove(i_index);
		}
		else
		{
			i->is_removed = 0;
			i->list_flag = MASK_CURRENT;
			mapping[i_index] = -1;
			unlabeled_num ++;
		}
		i->label_after_fix0 = -1;
		i->label_after_fix1 = -1;
	}

	maxflow();

	// INVARIANTS: 
	//    node i_index is removed <=> mapping[i_index] >= 0 <=> nodes[0][i_index].is_removed == 1
	//    edge e is removed <=> Arc::sister does not point to the correct arc for at least one out of the 4 arcs


	// Four lists are used: 1,2,3,4. In the beginning of each iteration
	// lists 3 and 4 are empty. current_list is set to the lowest non-empty list (1 or 2).
	// After the iteration current_list becomes empty, and its former nodes are moved
	// either to list 3 (if the node or its neighbor has changed) or to list 4.
	//
	// Invariants during the iteration: 
	//   - i->list_flag == MASK_CURRENT              => i is in current_list
	//   - i->list_flag == MASK_CURRENT & MASK_NEXT  => i is in current_list, after processing should be moved to list 3
	//   - i->list_flag ==                MASK_NEXT  => i is in list 3
	//   - i->list_flag == 0                         => i is not in current_list and not in list 3
	//
	// After the iteration, list 3 is eroded probe_dilation-1 times (it acquired nodes from lists 2 and 4).
	// It is then moved to list 1 (which is now empty) and list 4 is added to list 2.

	while ( 1 )
	{
		bool success = true;

		// Try fixing nodes to 0 and 1
		int current_list = (list.GetNext(-1) >= 0) ? 1 : 2;
		for (i_index=list.GetNext(-current_list) ; i_index>=0; i_index=i_index_next)
		{
			i_index_next = list.GetNext(i_index);
			i = nodes[0]+i_index;
			code_assert (!i->is_removed && i->label<0);

			_i[0] = i;
			_i[1] = GetMate0(i);
			bool is_changed = false;

			REAL INFTY0 = DetermineSaturation(i);
			REAL INFTY1 = DetermineSaturation(_i[1]);
			REAL INFTY = ((INFTY0 > INFTY1) ? INFTY0 : INFTY1) + 1;


			// fix to 0, run maxflow
			mark_node(i);
			mark_node(_i[1]);
			AddUnaryTerm(i, 0, INFTY);
			maxflow(true, true);
			if (what_segment(i)!=0 || what_segment(_i[1])!=1)
			{
				printf("Error in Probe()! Perhaps, overflow due to large capacities?\n");
			}
			for (ptr=changed_list->ScanFirst(); ptr; ptr=changed_list->ScanNext())
			{
				j = *ptr;
				code_assert(!j->is_removed && j->label<0);
				j->label_after_fix0 = what_segment(j);
				if (j->label_after_fix0 == what_segment(GetMate0(j))) j->label_after_fix0 = -1;
				else if (i->user_label == 0) j->user_label = j->label_after_fix0;
			}

			// fix to 1, run maxflow
			mark_node(i);
			mark_node(_i[1]);
			AddUnaryTerm(i, 0, -2*INFTY);
			maxflow(true, true);
			if (what_segment(i)!=1 || what_segment(_i[1])!=0)
			{
				printf("Error in Probe()! Perhaps, overflow due to large capacities?\n");
			}
			mark_node(i);
			mark_node(_i[1]);
			AddUnaryTerm(i, 0, INFTY);

			// go through changed pixels
			bool need_to_merge = false;

			for (ptr=changed_list->ScanFirst(); ptr; ptr=changed_list->ScanNext())
			{
				j = *ptr;
				code_assert(!j->is_removed && j->label<0);
				j->label_after_fix1 = what_segment(j);
				if (j->label_after_fix1 == what_segment(GetMate0(j))) j->label_after_fix1 = -1;
				else if (i->user_label == 1) j->user_label = j->label_after_fix1;

				if (i == j || j->label_after_fix0 < 0 || j->label_after_fix1 < 0) continue;

				j_index = (int)(j - nodes[0]);

				is_changed = true;
				if (j->label_after_fix0 == j->label_after_fix1)
				{
					// fix j
					FixNode(j, j->label_after_fix0);
					mapping[j_index] = j->label_after_fix0;
				}
				else
				{
					// contract i and j
					ContractNodes(i, j, j->label_after_fix0);
					mapping[j_index] = 2*i_index + 2 + j->label_after_fix0;
					need_to_merge = true;
				}
				j->is_removed = 1;
				if (i_index_next == j_index) i_index_next = list.GetNext(j_index);
				list.Remove(j_index);
				unlabeled_num --;
			}

			if (need_to_merge) 
			{
				// merge parallel edges incident to i
				for (a=i->first; a; a=a->next) // mark neighbor nodes
				{
					j_index = (int)(a->head - nodes[IsNode0(a->head) ? 0 : 1]);
					mapping[j_index] = (int)(a - arcs[0]);
				}
				for (a=i->first; a; a=a_next)
				{
					a_next = a->next;
					j_index = (int)(a->head - nodes[IsNode0(a->head) ? 0 : 1]);
					Arc* a2 = &arcs[0][mapping[j_index]];
					if (a2 == a) continue;
					mark_node(a->head);
					mark_node(GetMate(a->head));
					if (MergeParallelEdges(a2, a)==0) 
					{
						mapping[j_index] = (int)(a - arcs[0]);
						a_next = a->next;
					}
				}
				for (a=i->first; a; a=a->next)
				{
					j_index = (int)(a->head - nodes[IsNode0(a->head) ? 0 : 1]);
					mapping[j_index] = -1;
				}
			}

			// add directed links for neighbors of i
			for (a=i->first; a; a=a->next)
			{
				j = a->head;
				int label[2];
				if (IsNode0(j))
				{
					label[0] = j->label_after_fix0;
					label[1] = j->label_after_fix1;
				}
				else
				{
					label[0] = GetMate1(j)->label_after_fix0;
					label[1] = GetMate1(j)->label_after_fix1;
				}
				for (x=0; x<2; x++)
				{
					if (label[x]>=0 && label[1-x]<0)
					{
						if (AddDirectedConstraint0(a, x, label[x])) is_changed = true;
					}
				}
			}
			maxflow(true, true);
			mark_node(i);
			mark_node(_i[1]);
			for (a=i->first; a; a=a->next)
			{
				j = a->head;
				int label[2];
				if (IsNode0(j))
				{
					label[0] = j->label_after_fix0; j->label_after_fix0 = -1;
					label[1] = j->label_after_fix1; j->label_after_fix1 = -1;
				}
				else
				{
					label[0] = GetMate1(j)->label_after_fix0; GetMate1(j)->label_after_fix0 = -1;
					label[1] = GetMate1(j)->label_after_fix1; GetMate1(j)->label_after_fix1 = -1;
				}
				for (x=0; x<2; x++)
				{
					if (label[x]>=0 && label[1-x]<0)
					{
						if (AddDirectedConstraint1(a, x, label[x])) is_changed = true;
					}
				}
			}

			// add directed constraints for nodes which are not neighbors
			if (probe_options.directed_constraints!=0)
			{
				for (ptr=changed_list->ScanFirst(); ptr; ptr=changed_list->ScanNext())
				{
					j = *ptr;
					int x, y;

					if (j->is_removed) continue;
					if (i == j) continue;
					if      (j->label_after_fix0 >= 0 && j->label_after_fix1 < 0) { x = 0; y = j->label_after_fix0; }
					else if (j->label_after_fix1 >= 0 && j->label_after_fix0 < 0) { x = 1; y = j->label_after_fix1; }
					else continue;

					if (first_free)
					{
						AddDirectedConstraint(i, j, x, y);
						is_changed = true;
					}
					else
					{
						is_enough_memory = false;
						break;
					}
				}
			}

			if (probe_options.dilation >= 0)
			{
				// if is_changed, add i and its neighbors to list 3, otherwise add i to list 2 (unless it's already there)
				if (is_changed)
				{
					i->list_flag |= MASK_NEXT; 
					if (probe_options.dilation >= 1)
					{
						for (a=i->first; a; a=a->next)
						{
							j = a->head;
							if (!IsNode0(j)) j = GetMate1(j);
							if (!(j->list_flag & MASK_NEXT))
							{
								j->list_flag |= MASK_NEXT;
								if (!(j->list_flag & MASK_CURRENT))
								{
									j_index = (int)(j-nodes[0]);
									code_assert(j_index != i_index_next);
									list.Move(j_index, 3);
								}
							}
						}
					}
				}
				code_assert (i->list_flag & MASK_CURRENT);
				i->list_flag &= ~MASK_CURRENT;
				list.Move(i_index, (i->list_flag & MASK_NEXT) ? 3 : 4);
			}

			if (is_changed) success = false;

			// after fixes and contractions run maxflow, check whether more nodes have become labeled
			maxflow(true, true);

			for (ptr=changed_list->ScanFirst(); ptr; ptr=changed_list->ScanNext())
			{
				j = *ptr;

				j->is_in_changed_list = 0;
				j->label_after_fix0 = -1;
				j->label_after_fix1 = -1;

				if (j->is_removed) continue;

				j->label = what_segment(j);
				if (j->label == what_segment(GetMate0(j))) j->label = -1;

				if (j->label >= 0)
				{
					j_index = (int)(j - nodes[0]);
					FixNode(j, j->label);
					mapping[j_index] = j->label;
					j->is_removed = 1;
					if (i_index_next == j_index) i_index_next = list.GetNext(j_index);
					list.Remove(j_index);
					unlabeled_num --;
				}
			}
			changed_list->Reset();

			if (probe_options.callback_fn)
			{
				if (probe_options.callback_fn(unlabeled_num))
				{
					user_terminated = true;
					AllocateNewEnergy(mapping);
					return is_enough_memory;
				}
			}
		}

		if (probe_options.dilation < 0)
		{
			if (success) break;
			else         continue;
		}

		if (current_list == 2 && success) break;

		code_assert(list.GetNext(-1) == -1);
		list.MoveList(4, 2);
		if (list.GetNext(-3) < 0)
		{
			// list 3 is empty
			for (i_index=list.GetNext(-2); i_index>=0; i_index=list.GetNext(i_index))
			{
				i = nodes[0]+i_index;
				i->list_flag = MASK_CURRENT;
			}
		}
		else
		{
			int MASK_TMP = MASK_CURRENT; MASK_CURRENT = MASK_NEXT; MASK_NEXT = MASK_TMP;
			int r;
			for (r=1; r<probe_options.dilation; r++)
			{
				for (i_index=list.GetNext(-3); i_index>=0; i_index=list.GetNext(i_index))
				{
					i = nodes[0]+i_index;
					for (a=i->first; a; a=a->next)
					{
						j = a->head;
						if (!IsNode0(j)) j = GetMate1(j);
						if (!(j->list_flag & MASK_CURRENT))
						{
							j->list_flag = MASK_CURRENT;
							j_index = (int)(j-nodes[0]);
							code_assert(j_index != i_index_next);
							list.Move(j_index, 4);
						}
					}
				}
				list.MoveList(3, 1);
				list.MoveList(4, 3);
			}
			list.MoveList(3, 1);
		}
	}

	// almost done
	AllocateNewEnergy(mapping);
	Solve();

	return is_enough_memory;
}


template <typename REAL>
	void QPBO<REAL>::Probe(int* mapping, ProbeOptions& options)
{
	int nodeNum0 = GetNodeNum();
	bool is_enough_memory;
	user_terminated = false;

	memcpy(&probe_options, &options, sizeof(ProbeOptions));

	is_enough_memory = Probe(mapping);

	while ( 1 )
	{
		if (user_terminated) break;

		bool success = true;
		if ( probe_options.weak_persistencies )
		{
			int i;
			ComputeWeakPersistencies();
			for (i=1; i<GetNodeNum(); i++)
			{
				int ki = GetLabel(i);
				if (ki >= 0)
				{
					AddUnaryTerm(i, 0, (REAL)(1-2*ki));
					success = false;
				}
			}
		}

		if (probe_options.directed_constraints == 2 && !is_enough_memory)
		{
			SetMaxEdgeNum(GetMaxEdgeNum() + GetMaxEdgeNum()/2 + 1);
		}
		else
		{
			if (success) break;
		}

		int* mapping1 = new int[GetNodeNum()];
		is_enough_memory = Probe(mapping1);
		MergeMappings(nodeNum0, mapping, mapping1);
		delete [] mapping1;
	}
}

template <typename REAL>
	bool QPBO<REAL>::Improve(int N, int* order_array, int* fixed_nodes)
{
	int p, i_index;
	Node* i;
	Node* _i[2];
	FixNodeInfo* ptr;

	if (!fix_node_info_list)
	{
		fix_node_info_list = new Block<FixNodeInfo>(128);
	}

	maxflow();
	if (stage == 0)
	{
		TransformToSecondStage(true);
		maxflow(true);
	}

	for (p=0; p<N; p++)
	{
		i_index = order_array[p];
		code_assert(i_index>=0 && i_index<GetNodeNum());

		i = _i[0] = &nodes[0][i_index];
		_i[1] = &nodes[1][i_index];

		i->label = what_segment(i);
		if (i->label != what_segment(GetMate0(i)))
		{
			continue;
		}

		REAL INFTY0 = DetermineSaturation(i);
		REAL INFTY1 = DetermineSaturation(_i[1]);
		REAL INFTY = ((INFTY0 > INFTY1) ? INFTY0 : INFTY1) + 1;

		if (i->user_label == 1) INFTY = -INFTY;

		ptr = fix_node_info_list->New();
		ptr->i = i;
		ptr->INFTY = INFTY;

		mark_node(i);
		mark_node(GetMate0(i));
		AddUnaryTerm(i, 0, INFTY);
		maxflow(true);
	}

	if (fixed_nodes) memset(fixed_nodes, 0, GetNodeNum()*sizeof(int));
	for (ptr=fix_node_info_list->ScanFirst(); ptr; ptr=fix_node_info_list->ScanNext())
	{
		AddUnaryTerm(ptr->i, 0, -ptr->INFTY);
		if (fixed_nodes) fixed_nodes[(int)(ptr->i - nodes[0])] = 1;
	}
	fix_node_info_list->Reset();

	bool success = false;
	for (i=nodes[0]; i<node_last[0]; i++)
	{
		i->label = what_segment(i);
		if (i->label == what_segment(GetMate0(i))) i->label = i->user_label;
		else if (i->label != (int)i->user_label) 
		{
			success = true;
			i->user_label = (unsigned int)i->label;
		}
	}

	return success;
}

template <typename REAL>
	bool QPBO<REAL>::Improve()
{
	int* permutation = new int[node_num];
	ComputeRandomPermutation(node_num, permutation);
	bool success = Improve(node_num, permutation);
	delete [] permutation;
	return success;
}


/* QPBO_maxflow.cpp */
#define INFINITE_D ((int)(((unsigned)-1)/2))		/* infinite distance to the terminal */

/***********************************************************************/

/*
	Functions for processing active list.
	i->next points to the next node in the list
	(or to i, if i is the last node in the list).
	If i->next is NULL iff i is not in the list.

	There are two queues. Active nodes are added
	to the end of the second queue and read from
	the front of the first queue. If the first queue
	is empty, it is replaced by the second queue
	(and the second queue becomes empty).
*/


template <typename REAL> 
	inline void QPBO<REAL>::set_active(Node *i)
{
	if (!i->next)
	{
		/* it's not in the list yet */
		if (queue_last[1]) queue_last[1] -> next = i;
		else               queue_first[1]        = i;
		queue_last[1] = i;
		i -> next = i;
	}
}

/*
	Returns the next active node.
	If it is connected to the sink, it stays in the list,
	otherwise it is removed from the list
*/
template <typename REAL> 
	inline typename QPBO<REAL>::Node* QPBO<REAL>::next_active()
{
	Node *i;

	while ( 1 )
	{
		if (!(i=queue_first[0]))
		{
			queue_first[0] = i = queue_first[1];
			queue_last[0]  = queue_last[1];
			queue_first[1] = NULL;
			queue_last[1]  = NULL;
			if (!i) return NULL;
		}

		/* remove it from the active list */
		if (i->next == i) queue_first[0] = queue_last[0] = NULL;
		else              queue_first[0] = i -> next;
		i -> next = NULL;

		/* a node in the list is active iff it has a parent */
		if (i->parent) return i;
	}
}

/***********************************************************************/

template <typename REAL> 
	inline void QPBO<REAL>::set_orphan_front(Node *i)
{
	nodeptr *np;
	i -> parent = QPBO_MAXFLOW_ORPHAN;
	np = nodeptr_block -> New();
	np -> ptr = i;
	np -> next = orphan_first;
	orphan_first = np;
}

template <typename REAL> 
	inline void QPBO<REAL>::set_orphan_rear(Node *i)
{
	nodeptr *np;
	i -> parent = QPBO_MAXFLOW_ORPHAN;
	np = nodeptr_block -> New();
	np -> ptr = i;
	if (orphan_last) orphan_last -> next = np;
	else             orphan_first        = np;
	orphan_last = np;
	np -> next = NULL;
}

/***********************************************************************/

template <typename REAL> 
	inline void QPBO<REAL>::add_to_changed_list(Node *i)
{
	if (keep_changed_list)
	{
		if (!IsNode0(i)) i = GetMate1(i);
		if (!i->is_in_changed_list)
		{
			Node** ptr = changed_list->New();
			*ptr = i;;
			i->is_in_changed_list = true;
		}
	}
}

/***********************************************************************/

template <typename REAL> 
	void QPBO<REAL>::maxflow_init()
{
	Node *i;

	queue_first[0] = queue_last[0] = NULL;
	queue_first[1] = queue_last[1] = NULL;
	orphan_first = NULL;

	TIME = 0;

	for (i=nodes[0]; i<node_last[stage]; i++)
	{
		if (i==node_last[0]) i = nodes[1];

		i -> next = NULL;
		i -> is_marked = 0;
		i -> is_in_changed_list = 0;
		i -> TS = TIME;
		if (i->tr_cap > 0)
		{
			/* i is connected to the source */
			i -> is_sink = 0;
			i -> parent = QPBO_MAXFLOW_TERMINAL;
			set_active(i);
			i -> DIST = 1;
		}
		else if (i->tr_cap < 0)
		{
			/* i is connected to the sink */
			i -> is_sink = 1;
			i -> parent = QPBO_MAXFLOW_TERMINAL;
			set_active(i);
			i -> DIST = 1;
		}
		else
		{
			i -> parent = NULL;
		}
	}
}

template <typename REAL> 
	void QPBO<REAL>::maxflow_reuse_trees_init()
{
	Node* i;
	Node* j;
	Node* queue = queue_first[1];
	Arc* a;
	nodeptr* np;

	queue_first[0] = queue_last[0] = NULL;
	queue_first[1] = queue_last[1] = NULL;
	orphan_first = orphan_last = NULL;

	TIME ++;

	while ((i=queue))
	{
		queue = i->next;
		if (queue == i) queue = NULL;
		if (IsNode0(i))
		{
			if (i->is_removed) continue;
		}
		else
		{
			if (GetMate1(i)->is_removed) continue;
		}
		i->next = NULL;
		i->is_marked = 0;
		set_active(i);

		if (i->tr_cap == 0)
		{
			if (i->parent) set_orphan_rear(i);
			continue;
		}

		if (i->tr_cap > 0)
		{
			if (!i->parent || i->is_sink)
			{
				i->is_sink = 0;
				for (a=i->first; a; a=a->next)
				{
					j = a->head;
					if (!j->is_marked)
					{
						if (j->parent == a->sister) set_orphan_rear(j);
						if (j->parent && j->is_sink && a->r_cap > 0) set_active(j);
					}
				}
				add_to_changed_list(i);
			}
		}
		else
		{
			if (!i->parent || !i->is_sink)
			{
				i->is_sink = 1;
				for (a=i->first; a; a=a->next)
				{
					j = a->head;
					if (!j->is_marked)
					{
						if (j->parent == a->sister) set_orphan_rear(j);
						if (j->parent && !j->is_sink && a->sister->r_cap > 0) set_active(j);
					}
				}
				add_to_changed_list(i);
			}
		}
		i->parent = QPBO_MAXFLOW_TERMINAL;
		i -> TS = TIME;
		i -> DIST = 1;
	}

	code_assert(stage == 1);
	//test_consistency();

	/* adoption */
	while ((np=orphan_first))
	{
		orphan_first = np -> next;
		i = np -> ptr;
		nodeptr_block -> Delete(np);
		if (!orphan_first) orphan_last = NULL;
		if (i->is_sink) process_sink_orphan(i);
		else            process_source_orphan(i);
	}
	/* adoption end */

	//test_consistency();
}

template <typename REAL> 
	void QPBO<REAL>::augment(Arc *middle_arc)
{
	Node *i;
	Arc *a;
	REAL bottleneck;


	/* 1. Finding bottleneck capacity */
	/* 1a - the source tree */
	bottleneck = middle_arc -> r_cap;
	for (i=middle_arc->sister->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == QPBO_MAXFLOW_TERMINAL) break;
		if (bottleneck > a->sister->r_cap) bottleneck = a -> sister -> r_cap;
	}
	if (bottleneck > i->tr_cap) bottleneck = i -> tr_cap;
	/* 1b - the sink tree */
	for (i=middle_arc->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == QPBO_MAXFLOW_TERMINAL) break;
		if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
	}
	if (bottleneck > - i->tr_cap) bottleneck = - i -> tr_cap;


	/* 2. Augmenting */
	/* 2a - the source tree */
	middle_arc -> sister -> r_cap += bottleneck;
	middle_arc -> r_cap -= bottleneck;
	for (i=middle_arc->sister->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == QPBO_MAXFLOW_TERMINAL) break;
		a -> r_cap += bottleneck;
		a -> sister -> r_cap -= bottleneck;
		if (!a->sister->r_cap)
		{
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	i -> tr_cap -= bottleneck;
	if (!i->tr_cap)
	{
		set_orphan_front(i); // add i to the beginning of the adoption list
	}
	/* 2b - the sink tree */
	for (i=middle_arc->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == QPBO_MAXFLOW_TERMINAL) break;
		a -> sister -> r_cap += bottleneck;
		a -> r_cap -= bottleneck;
		if (!a->r_cap)
		{
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	i -> tr_cap += bottleneck;
	if (!i->tr_cap)
	{
		set_orphan_front(i); // add i to the beginning of the adoption list
	}
}

/***********************************************************************/

template <typename REAL> 
	void QPBO<REAL>::process_source_orphan(Node *i)
{
	Node *j;
	Arc *a0, *a0_min = NULL, *a;
	int d, d_min = INFINITE_D;

	/* trying to find a new parent */
	for (a0=i->first; a0; a0=a0->next)
	if (a0->sister->r_cap)
	{
		j = a0 -> head;
		if (!j->is_sink && (a=j->parent))
		{
			/* checking the origin of j */
			d = 0;
			while ( 1 )
			{
				if (j->TS == TIME)
				{
					d += j -> DIST;
					break;
				}
				a = j -> parent;
				d ++;
				if (a==QPBO_MAXFLOW_TERMINAL)
				{
					j -> TS = TIME;
					j -> DIST = 1;
					break;
				}
				if (a==QPBO_MAXFLOW_ORPHAN) { d = INFINITE_D; break; }
				j = a -> head;
			}
			if (d<INFINITE_D) /* j originates from the source - done */
			{
				if (d<d_min)
				{
					a0_min = a0;
					d_min = d;
				}
				/* set marks along the path */
				for (j=a0->head; j->TS!=TIME; j=j->parent->head)
				{
					j -> TS = TIME;
					j -> DIST = d --;
				}
			}
		}
	}

	if (i->parent = a0_min)
	{
		i -> TS = TIME;
		i -> DIST = d_min + 1;
	}
	else
	{
		/* no parent is found */
		add_to_changed_list(i);

		/* process neighbors */
		for (a0=i->first; a0; a0=a0->next)
		{
			j = a0 -> head;
			if (!j->is_sink && (a=j->parent))
			{
				if (a0->sister->r_cap) set_active(j);
				if (a!=QPBO_MAXFLOW_TERMINAL && a!=QPBO_MAXFLOW_ORPHAN && a->head==i)
				{
					set_orphan_rear(j); // add j to the end of the adoption list
				}
			}
		}
	}
}

template <typename REAL> 
	void QPBO<REAL>::process_sink_orphan(Node *i)
{
	Node *j;
	Arc *a0, *a0_min = NULL, *a;
	int d, d_min = INFINITE_D;

	/* trying to find a new parent */
	for (a0=i->first; a0; a0=a0->next)
	if (a0->r_cap)
	{
		j = a0 -> head;
		if ((a=j->parent) && j->is_sink)
		{
			/* checking the origin of j */
			d = 0;
			while ( 1 )
			{
				if (j->TS == TIME)
				{
					d += j -> DIST;
					break;
				}
				a = j -> parent;
				d ++;
				if (a==QPBO_MAXFLOW_TERMINAL)
				{
					j -> TS = TIME;
					j -> DIST = 1;
					break;
				}
				if (a==QPBO_MAXFLOW_ORPHAN) { d = INFINITE_D; break; }
				j = a -> head;
			}
			if (d<INFINITE_D) /* j originates from the sink - done */
			{
				if (d<d_min)
				{
					a0_min = a0;
					d_min = d;
				}
				/* set marks along the path */
				for (j=a0->head; j->TS!=TIME; j=j->parent->head)
				{
					j -> TS = TIME;
					j -> DIST = d --;
				}
			}
		}
	}

	if (i->parent = a0_min)
	{
		i -> TS = TIME;
		i -> DIST = d_min + 1;
	}
	else
	{
		/* no parent is found */
		add_to_changed_list(i);

		/* process neighbors */
		for (a0=i->first; a0; a0=a0->next)
		{
			j = a0 -> head;
			if ((a=j->parent) && j->is_sink)
			{
				if (a0->r_cap) set_active(j);
				if (a!=QPBO_MAXFLOW_TERMINAL && a!=QPBO_MAXFLOW_ORPHAN && a->head==i)
				{
					set_orphan_rear(j); // add j to the end of the adoption list
				}
			}
		}
	}
}

/***********************************************************************/

template <typename REAL> 
	void QPBO<REAL>::maxflow(bool reuse_trees, bool _keep_changed_list)
{
	Node *i, *j, *current_node = NULL;
	Arc *a;
	nodeptr *np, *np_next;

	if (!nodeptr_block)
	{
		nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);
	}

	if (maxflow_iteration == 0)
	{
		reuse_trees = false;
		_keep_changed_list = false;
	}

	keep_changed_list = _keep_changed_list;
	if (keep_changed_list)
	{
		if (!changed_list) changed_list = new Block<Node*>(NODEPTR_BLOCK_SIZE, error_function);
	}

	if (reuse_trees) maxflow_reuse_trees_init();
	else             maxflow_init();

	// main loop
	while ( 1 )
	{
		// test_consistency(current_node);

		if ((i=current_node))
		{
			i -> next = NULL; /* remove active flag */
			if (!i->parent) i = NULL;
		}
		if (!i)
		{
			if (!(i = next_active())) break;
		}

		/* growth */
		if (!i->is_sink)
		{
			/* grow source tree */
			for (a=i->first; a; a=a->next)
			if (a->r_cap)
			{
				j = a -> head;
				if (!j->parent)
				{
					j -> is_sink = 0;
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
					set_active(j);
					add_to_changed_list(j);
				}
				else if (j->is_sink) break;
				else if (j->TS <= i->TS &&
				         j->DIST > i->DIST)
				{
					/* heuristic - trying to make the distance from j to the source shorter */
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
				}
			}
		}
		else
		{
			/* grow sink tree */
			for (a=i->first; a; a=a->next)
			if (a->sister->r_cap)
			{
				j = a -> head;
				if (!j->parent)
				{
					j -> is_sink = 1;
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
					set_active(j);
					add_to_changed_list(j);
				}
				else if (!j->is_sink) { a = a -> sister; break; }
				else if (j->TS <= i->TS &&
				         j->DIST > i->DIST)
				{
					/* heuristic - trying to make the distance from j to the sink shorter */
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
				}
			}
		}

		TIME ++;

		if (a)
		{
			i -> next = i; /* set active flag */
			current_node = i;

			/* augmentation */
			augment(a);
			/* augmentation end */

			/* adoption */
			while ((np=orphan_first))
			{
				np_next = np -> next;
				np -> next = NULL;

				while ((np=orphan_first))
				{
					orphan_first = np -> next;
					i = np -> ptr;
					nodeptr_block -> Delete(np);
					if (!orphan_first) orphan_last = NULL;
					if (i->is_sink) process_sink_orphan(i);
					else            process_source_orphan(i);
				}

				orphan_first = np_next;
			}
			/* adoption end */
		}
		else current_node = NULL;
	}
	// test_consistency();

	if (!reuse_trees || (maxflow_iteration % 64) == 0)
	{
		delete nodeptr_block; 
		nodeptr_block = NULL; 
	}

	maxflow_iteration ++;
}

/***********************************************************************/


template <typename REAL> 
	void QPBO<REAL>::test_consistency(Node* current_node)
{
	Node *i;
	Arc *a;
	int r;
	int num1 = 0, num2 = 0;

	// test whether all nodes i with i->next!=NULL are indeed in the queue
	for (i=nodes[0]; i<node_last[stage]; i++)
	{
		if (i==node_last[0]) i = nodes[1];
		if ((IsNode0(i) && i->is_removed) || (!IsNode0(i) && GetMate1(i)->is_removed))
		{
			code_assert(i->first == NULL);
			continue;
		}

		if (i->next || i==current_node) num1 ++;
	}
	for (r=0; r<3; r++)
	{
		i = (r == 2) ? current_node : queue_first[r];
		if (i)
		for ( ; ; i=i->next)
		{
			code_assert((IsNode0(i) && !i->is_removed) || (!IsNode0(i) && !GetMate1(i)->is_removed));
			num2 ++;
			if (i->next == i)
			{
				if (r<2) code_assert(i == queue_last[r]);
				else     code_assert(i == current_node);
				break;
			}
		}
	}
	code_assert(num1 == num2);

	for (i=nodes[0]; i<node_last[stage]; i++)
	{
		if (i==node_last[0]) i = nodes[1];
		if ((IsNode0(i) && i->is_removed) || (!IsNode0(i) && GetMate1(i)->is_removed)) continue;

		// test whether all edges in seach trees are non-saturated
		if (i->parent == NULL) {}
		else if (i->parent == QPBO_MAXFLOW_ORPHAN) {}
		else if (i->parent == QPBO_MAXFLOW_TERMINAL)
		{
			if (!i->is_sink) code_assert(i->tr_cap > 0);
			else             code_assert(i->tr_cap < 0);
		}
		else
		{
			if (!i->is_sink) code_assert (i->parent->sister->r_cap > 0);
			else             code_assert (i->parent->r_cap > 0);
		}
		// test whether passive nodes in search trees have neighbors in
		// a different tree through non-saturated edges
		if (i->parent && !i->next)
		{
			if (!i->is_sink)
			{
				code_assert(i->tr_cap >= 0);
				for (a=i->first; a; a=a->next)
				{
					if (a->r_cap > 0) code_assert(a->head->parent && !a->head->is_sink);
				}
			}
			else
			{
				code_assert(i->tr_cap <= 0);
				for (a=i->first; a; a=a->next)
				{
					if (a->sister->r_cap > 0) code_assert(a->head->parent && a->head->is_sink);
				}
			}
		}
		// test marking invariants
		if (i->parent && i->parent!=QPBO_MAXFLOW_ORPHAN && i->parent!=QPBO_MAXFLOW_TERMINAL)
		{
			code_assert(i->TS <= i->parent->head->TS);
			if (i->TS == i->parent->head->TS) code_assert(i->DIST > i->parent->head->DIST);
		}
	}
}

/* QPBO_postprocessing.cpp */
template <typename REAL>
	void QPBO<REAL>::ComputeWeakPersistencies()
{
	if (stage == 0) return;

	Node* i;
	Node* j;
	Node* stack = NULL;
	int component;

	for (i=nodes[0]; i<node_last[0]; i++)
	{
		code_assert(i->label>=-1 && i->label<=1);

		Node* i1 = GetMate0(i);

		if (i->label >= 0)
		{
			i->dfs_parent = i;
			i1->dfs_parent = i1;
			i->region = i1->region = 0;
		}
		else
		{
			i->dfs_parent = i1->dfs_parent = NULL;
			i->region = i1->region = -1;
		}
	}

	// first DFS
	for (i=nodes[0]; i<node_last[1]; i++)
	{
		if (i == node_last[0]) i = nodes[1];
		if (i->dfs_parent) continue;

		// DFS starting from i
		i->dfs_parent = i;
		i->dfs_current = i->first;
		while ( 1 )
		{
			if (!i->dfs_current)
			{
				i->next = stack;
				stack = i;

				if (i->dfs_parent == i) break;
				i = i->dfs_parent;
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j = i->dfs_current->head;
			if (!(i->dfs_current->r_cap>0) || j->dfs_parent)
			{
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j->dfs_parent = i;
			i = j;
			i->dfs_current = i->first;
		}
	}

	// second DFS
	component = 0;
	while ( stack )
	{
		i = stack;
		stack = i->next;
		if (i->region > 0) continue;

		i->region = ++ component;
		i->dfs_parent = i;
		i->dfs_current = i->first;
		while ( 1 )
		{
			if (!i->dfs_current)
			{
				if (i->dfs_parent == i) break;
				i = i->dfs_parent;
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j = i->dfs_current->head;
			if (!(i->dfs_current->sister->r_cap>0) || j->region>=0)
			{
				i->dfs_current = i->dfs_current->next;
				continue;
			}

			j->dfs_parent = i;
			i = j;
			i->dfs_current = i->first;
			i->region = component;
		}
	}

	// assigning labels
	for (i=nodes[0]; i<node_last[0]; i++)
	{
		if (i->label < 0)
		{
			code_assert(i->region > 0);
			if      (i->region > GetMate0(i)->region) { i->label = 0; i->region = 0; }
			else if (i->region < GetMate0(i)->region) { i->label = 1; i->region = 0; }
		}
		else code_assert(i->region == 0);
	}
}

template <typename REAL>
	void QPBO<REAL>::Stitch()
{
	if (stage == 0) return;

	Node* i;
	Node* i_mate;
	Node* j;
	Arc* a;
	Arc* a_mate;

	for (a=arcs[0], a_mate=arcs[1]; a<arc_max[0]; a++, a_mate++)
	if (a->sister)
	{
		a->r_cap = a_mate->r_cap = a->r_cap + a_mate->r_cap;

		i = a->sister->head;
		j = a->head;

		if (i->region==0 || i->region != j->region) continue;
		if (IsNode0(i))
		{
			if (i->user_label != 0) continue;
		}
		else
		{
			if (GetMate1(i)->user_label != 1) continue;
		}
		if (IsNode0(j))
		{
			if (j->user_label != 1) continue;
		}
		else
		{
			if (GetMate1(j)->user_label != 0) continue;
		}

		a->r_cap = a_mate->r_cap = 0;
	}

	for (i=nodes[0], i_mate=nodes[1]; i<node_last[0]; i++, i_mate++)
	{
		i->tr_cap = i->tr_cap - i_mate->tr_cap;
		i_mate->tr_cap = -i->tr_cap;
	}

	ComputeWeakPersistencies();
}


#endif
