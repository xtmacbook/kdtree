
#include "GlbSpatialKdTree.h"

using namespace GlbGlobe;
///Users/glp/Documents/3Dirty/OpenSceneGraph-3.2.3/examples/osgpdf/kdtree/kdtree


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
	KDTNode* ropes[6] = { NULL };
}
