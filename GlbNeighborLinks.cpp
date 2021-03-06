
#include "GlbSpatialKdTree.h"

using namespace GlbGlobe;

//两三角形是否共面
static bool towTrianglecoplanar(const GlbGlobe::Vec3&A0,const GlbGlobe::Vec3&A1,const GlbGlobe::Vec3&A2,
								const GlbGlobe::Vec3&B0,const GlbGlobe::Vec3&B1,const GlbGlobe::Vec3&B2)
{
	GlbGlobe::Vec3 U = A1 - A0; //A
	GlbGlobe::Vec3 V = A2 - A0;

	GlbGlobe::Vec3 S = B1 - B0; //P
	GlbGlobe::Vec3 T = B2 - B0;

	GlbGlobe::Vec3 N1 = U ^ V;
	GlbGlobe::Vec3 N2 = S ^ T;

	if((N1 ^ N2).length2() == 0) return false;

	if((B0 - A0) * N1 == 0) return true;

	return false;
}

template <typename T>
static void optimizeRopes(T * ropes[],GlbGlobe::BoundingBox&box)
{
	// Loop through ropes of all faces of node bounding box.
	for ( int i = 0; i < 6; ++i )
	{
		T *rope_node = ropes[i];

		if ( rope_node == NULL )
		{
			continue;
		}

		// Process until leaf node is reached.
		// The optimization - We try to push the ropes down into the tree as far as possible
		// instead of just having the ropes point to the roots of neighboring subtrees.
		while ( !rope_node->is_leaf() )
		{
			GlbGlobe::Axes splitAxis = rope_node->splitEdge->axis;
			double splitValue = rope_node->splitEdge->splitPlanePosition;

			//FLeft = 0x00, FRight, FFront, FBack, FBottom, FTop
			if ( i == FLeft || i == FRight )
			{
				// Case I-A.
				// Handle parallel split plane case.
				if ( splitAxis == X_axis )
				{
					rope_node = ( i == FLeft ) ? rope_node->right : rope_node->left;
				}

				// Case I-B.

				else if ( splitAxis == Y_axis )
				{
					if ( splitValue < ( box._min.y() - KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->right;
					}
					else if ( splitValue > ( box._max.y() + KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->left;
					}
					else
					{
						break;
					}
				}

				// Case I-C.

				// Split plane is Z_AXIS.
				else
				{
					if (splitValue < ( box._min.z() - KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->right;
					}
					else if ( splitValue > ( box._max.z() + KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->left;
					}
					else
					{
						break;
					}
				}
			}

			else if ( i == FFront || i == FBack )
			{
				// Handle parallel split plane case.
				if ( splitAxis == Z_axis )
				{
					rope_node = ( i == FBack ) ? rope_node->right : rope_node->left;
				}

				// Case II-B.

				else if ( splitAxis == X_axis )
				{
					if ( splitValue < ( box._min.x() - KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->right;
					}
					else if ( splitAxis > ( box._max.x() + KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->left;
					}
					else
					{
						break;
					}
				}

				// Split plane is Y_AXIS.
				else
				{
					if ( splitValue < ( box._min.y() - KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->right;
					}
					else if ( splitValue > ( box._max.y() + KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->left;
					}
					else
					{
						break;
					}
				}
			}

			// TOP and BOTTOM.
			else
			{
				// Case III-A.
				// Handle parallel split plane case.
				if ( splitAxis == Y_axis )
				{
					rope_node = ( i == FBottom ) ? rope_node->right : rope_node->left;
				}

				// Case III-B.

				else if ( splitAxis == Z_axis )
				{
					if ( splitValue < ( box._min.z() - KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->right;
					}
					else if ( splitValue > ( box._max.z() + KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->left;
					}
					else
					{
						break;
					}
				}

				// Case III-C.

				// Split plane is X_AXIS.
				else
				{
					if ( splitValue < ( box._min.x() - KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->right;
					}
					else if ( splitValue > ( box._max.x() + KD_TREE_EPSILON ) )
					{
						rope_node = rope_node->left;
					}
					else
					{
						break;
					}
				}
			}
		}
	}
}

/*
	可以参考pdf ,singleRay 为true，创建的neigh link,指向可能为内部节点也可能是叶节点
*/
template <typename T>
static void buildRopeStructure( T *curr_node, T *rs[], bool singleRay)
{
#ifdef KDTREE_NEIGHBORLINKS

	if ( curr_node->is_leaf() )
	{
		for ( unsigned int i = 0; i < 6; ++i )
		{
			curr_node->ropes[i] = rs[i];
		}
	}
	else
	{
		
		if ( singleRay )
		{
			optimizeRopes( rs, curr_node->box );
		}

		GlbGlobe::Faces SL, SR; //FLeft = 0x00, FRight, FFront, FBack, FBottom, FTop
		if ( curr_node->splitEdge->axis == X_axis )
		{
			SL = FLeft;
			SR = FRight;
		}
		else if ( curr_node->splitEdge->axis == Y_axis )
		{
			SL = FBottom;
			SR = FTop;
		}
		// Split plane is Z_AXIS.
		else
		{
			SL = FBack;
			SR = FFront;
		}

		T* RS_left[6];
		T* RS_right[6];
		for ( unsigned int i = 0; i < 6; ++i )
		{
			RS_left[i] = rs[i];
			RS_right[i] = rs[i];
		}

		// Recurse.
		RS_left[SR] = curr_node->right;
		buildRopeStructure( curr_node->left, RS_left, singleRay );

		// Recurse.
		RS_right[SL] = curr_node->left;
		buildRopeStructure( curr_node->right, RS_right, singleRay );
	}

#endif
}


KDTNode* GLbKdTree::FindSingleNeighborLink(KDTNode *node, Faces face, KDTNode *rootNode)
{
	KDTNode *currNode; /* currently accessed node */

	/*if (node->box[face] is coplanar with a face of sceneBox)
		return NULL;*/

	KDTNode *stack[TERMINATION_CRITERIA_D_MAX]; 

	/* search starts from the root node */
	unsigned int stackIndex = 0;
	stack[stackIndex++] = rootNode;

	while(stackIndex != 0)
	{
		currNode = stack[--stackIndex];

		//不是叶节点
		if(!currNode->is_leaf())
		{
			if(currNode->splitEdge->axis == node->splitEdge->axis)
			{
				/*splitting plane 平行于node的一个面*/
				if(node->box.GetExtent(face) < currNode->splitEdge->splitPlanePosition - KD_TREE_EPSILON)
				{
					stack[stackIndex++] = currNode->left;
				}
				else
				{
					if(node->box.GetExtent(face) > currNode->splitEdge->splitPlanePosition - KD_TREE_EPSILON)
					{
						stack[stackIndex++] = currNode->right;
					}
					else
					{
						//splitting plane 在Face里
						if(/*face is min face of node*/)
						{
							stack[stackIndex++] = currNode->left;
						}
						else
						{
							stack[stackIndex++] = currNode->right;
						}
					}
				}
			}
			else
			{
				/*is some other axis so test if it splits the face or not */
				if(node->box._min[currNode->splitEdge->axis] >= currNode->splitEdge->splitPlanePosition)
				{
					stack[stackIndex++] = currNode->right;
				}
				else
				{
					if(node->box._max[currNode->splitEdge->axis] <= currNode->splitEdge->splitPlanePosition)
					{
						stack[stackIndex++] = currNode->left;
					}
				}
			}
		}
	}

	return currNode;
}


void GLbKdTree::BuildNeighborLinksTree(KDTNode *node, Faces face)
{ 
#ifdef KDTREE_NEIGHBORLINKS

	if (node->ropes[face] == NULL) return;
	if (node->ropes[face]->is_leaf()) return;

	/* the neighbor-link points to an interior kd-tree node, */
	/* it will be replaced a neighbor-links tree, it has sense */
	node->ropes[face] = CreateNeighborLinksTree(node, node->ropes[face], face);

#endif
}
/* BuildNeighborLinksTree */

KDTNode* GLbKdTree::CreateNeighborLinksTree(KDTNode *node, KDTNode *subtree, Faces face)
{
	KDTNode *currNode = subtree;

	while(!currNode->is_leaf())
	{
		if(currNode->splitEdge->axis /*is perpendicular to face*/)
		{
			if (/*node->box[face] is to the left of currNode->splitPlane*/1)
			{
				currNode = currNode->left;
			}
			else
			{
				if ( /*node->box[face] is to the right of currNode->splitPlane*/1)
				{
					currNode = currNode->right;
				}
				else
				{
					/* limits are equal make decision regarding to the oposite limit */
					if (/*face is a min face of currNode*/1)
						currNode = currNode->left; /* it was a min limit – go left */
					else
						currNode = currNode->right; /* it was a max limit – go right */
				}
			}
		}
		else
		{
			/* it is some other axis so test if it splits the face or not */
			if (/*node->box is to the right of currNode->splitPlane*/1)
			{
				currNode = currNode->right; /* greater */
			}
			else
			{
				if (/* node->box is to the left of currNode->splitPlane*/1)
				{
					currNode = currNode->left; /* smaller */
				}
				else
				{
					/* this node intersects the face – it must be added to the created neighbor-links tree */
					KDTNode *result = new KDTNode;
					result->splitEdge = currNode->splitEdge;
					result->splitEdge->axis = currNode->splitEdge->axis;
					/* recusively call itself to get whole neighbor-links tree */
					result->left = CreateNeighborLinksTree(node, currNode->left, face);
					result->right = CreateNeighborLinksTree(node, currNode->right, face);
					return result;
				}
			}
		}
	}

	return currNode;
}

void GLbKdTree::BuildRopeStructure()
{
	if(!sahUse)
	{
		KDTNodeM* ropes[6] = { NULL };
		buildRopeStructure(treeRootM, ropes, true);
	}
}
