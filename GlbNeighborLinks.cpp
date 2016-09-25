
#include "GlbSpatialKdTree.h"

using namespace GlbGlobe;

extern const double KD_TREE_EPSILON;

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

template <typename T>
static void buildRopeStructure( T *curr_node, T *rs[], bool is_single_ray_case)
{
    if ( curr_node->is_leaf() )
        {
        for ( unsigned int i = 0; i < 6; ++i )
            {
            curr_node->ropes[i] = rs[i];
            }
        }
    else
        {
            // Only optimize ropes on single ray case.
            // It is not optimal to optimize on packet traversal case.
        if ( is_single_ray_case )
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
        buildRopeStructure( curr_node->left, RS_left, is_single_ray_case );

            // Recurse.
        RS_right[SL] = curr_node->left;
        buildRopeStructure( curr_node->right, RS_right, is_single_ray_case );
        }
}


KDTNode* GLbKdTree::FindSingleNeighborLink(KDTNode *node, Faces face, KDTNode *rootNode)
{
	KDTNode *currNode; /* currently accessed node */

	/*if (node->box[face] is coplanar with a face of sceneBox)
	return ["No neighbor link exists"];
	*/
	/* stack required for traversal to store far child nodes */
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
				 //if(node)
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
