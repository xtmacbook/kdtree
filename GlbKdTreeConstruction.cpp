#include "GlbSpatialKdTree.h"

#include <algorithm>

using namespace GlbGlobe;

const bool USE_TIGHT_FITTING_BOUNDING_BOXES = false;
const unsigned int INTINFINITY = std::numeric_limits<unsigned int >::max();



bool GLbKdTree::ConstructKdTree(unsigned int num_tris,const BoundingBox&bounds)
{
	if(sahUse)
	{
		treeRoot = ConstructTreeSAHSplit(num_tris,bounds);
	}
	else
	{
		treeRootM = ConstructTreeMedianSpaceSplit(num_tris,bounds);
	}
	return true;
}

KDTNodeM* GLbKdTree::ConstructTreeMedianSpaceSplit(unsigned int num_tris,unsigned int *tri_indices,
	const BoundingBox& bounds, unsigned int curr_depth,KDTNodeM*parent)
{
	//创建新的节点
	KDTNodeM * currentNode = new KDTNodeM();

#ifndef KDTREE_NEIGHBORLINKS 
	 currentNode->parent = parent;
#endif

	currentNode->parent = parent;

	currentNode->num_tris = num_tris;

	if ( USE_TIGHT_FITTING_BOUNDING_BOXES ) 
	{
		currentNode->box = computeTightFittingBoundingBox( num_tris ,tri_indices);
	}
	else 
	{
		currentNode->box = bounds;
	}

	if(num_tris <= TERMINATION_CRITERIA_N_MAX || curr_depth >= TERMINATION_CRITERIA_D_MAX)
	{

		currentNode->tri_indices = tri_indices;

		return currentNode;
	}

	const Vec3& boundMax = bounds._max;
	const Vec3& boundMin = bounds._min;

	//设置splitting plane 轴
	double xlength = boundMax.x() - boundMin.x();
	double ylength = boundMax.y() - boundMin.y();
	double zlength = boundMax.z() - boundMin.z();

	//splitting plane轴向
	Axes bestAxes = ( xlength > ylength && xlength > zlength ) ? X_axis :
		( ylength > zlength ? Y_axis : Z_axis);
	currentNode->splitEdge->axis = bestAxes;

	//对最长轴进行中点划分
	double median_val = 0.0;
	BoundingBox left_bbox  = bounds;
	BoundingBox right_bbox = bounds;

	if ( bestAxes == X_axis ) 
	{
		median_val = boundMin.x() + ( ( boundMax.x() - boundMin.x() ) / 2.0f );
		left_bbox._max.x() = median_val;
		right_bbox._min.x() = median_val;
	}
	else if ( bestAxes == Y_axis ) 
	{
		median_val = boundMin.y() + ( ( boundMax.y() - boundMin.y() ) / 2.0f );
		left_bbox._max.y() = median_val;
		right_bbox._min.y() = median_val;
	}
	else 
	{
		median_val = boundMin.z() + ( ( boundMax.z() - boundMin.z() ) / 2.0f );
		left_bbox._max.z() = median_val;
		right_bbox._min.z() = median_val;
	}
	//设置splitting plane位置
	currentNode->splitEdge->splitPlanePosition = median_val;


	unsigned int *temp_left_tri_indices = new unsigned int[num_tris];
	unsigned int *temp_right_tri_indices = new unsigned int[num_tris];

	unsigned int left_tri_count = 0, right_tri_count = 0;
	double min_tri_val, max_tri_val;

	for ( unsigned int i = 0; i < num_tris; ++i ) 
	{
		unsigned triIndice = tri_indices[i];

		Triangle&tri = meshTriangles[triIndice];

		Vec3& v0 = tri.vertex[0];
		Vec3& v1 = tri.vertex[1];
		Vec3& v2 = tri.vertex[2];

		// 获取splitting plane最大最小值

		if(bestAxes == X_axis)
		{
			min_tri_val = ( v0.x() < v1.x() && v0.x() < v2.x() ) ? v0.x() : ( v1.x() < v2.x() ? v1.x() : v2.x() );
			max_tri_val = ( v0.x() > v1.x() && v0.x() > v2.x() ) ? v0.x() : ( v1.x() > v2.x() ? v1.x() : v2.x() );
		}
		else if(bestAxes == Y_axis)
		{
			min_tri_val = ( v0.y() < v1.y() && v0.y() < v2.y() ) ? v0.y() : ( v1.y() < v2.y() ? v1.y() : v2.y() );
			max_tri_val = ( v0.y() > v1.y() && v0.y() > v2.y() ) ? v0.y() : ( v1.y() > v2.y() ? v1.y() : v2.y() );
		}
		else
		{
			min_tri_val = ( v0.z() < v1.z() && v0.z() < v2.z() ) ? v0.z() : ( v1.z() < v2.z() ? v1.z() : v2.z() );
			max_tri_val = ( v0.z() > v1.z() && v0.z() > v2.z() ) ? v0.z() : ( v1.z() > v2.z() ? v1.z() : v2.z() );
		}

		// 重置
		if ( min_tri_val < median_val )
		{
			temp_left_tri_indices[i] = tri_indices[i];
			++left_tri_count;
		}
		else 
		{
			temp_left_tri_indices[i] = INTINFINITY;
		}

		// Update temp_right_tri_indices.
		if ( max_tri_val >= median_val ) 
		{
			temp_right_tri_indices[i] = tri_indices[i];
			++right_tri_count;
		}
		else 
		{
			temp_right_tri_indices[i] = INTINFINITY;
		}
	}

	//为左右子节点分配索引容器
	unsigned int *left_tri_indices = new unsigned int[left_tri_count];
	unsigned int *right_tri_indices = new unsigned int[right_tri_count];

	// Populate lists of triangle indices.
	unsigned int left_index = 0, right_index = 0;
	for (unsigned int i = 0; i < num_tris; ++i )
	{
		if ( temp_left_tri_indices[i] != INTINFINITY )
		{
			left_tri_indices[left_index] = temp_left_tri_indices[i];
			++left_index;
		}
		if ( temp_right_tri_indices[i] != INTINFINITY )
		{
			right_tri_indices[right_index] = temp_right_tri_indices[i];
			++right_index;
		}
	}

	// 删除临时索引
	delete[] temp_left_tri_indices;
	delete[] temp_right_tri_indices;

	delete[] tri_indices;
	tri_indices = NULL;
	// 递归子节点
	currentNode->left = ConstructTreeMedianSpaceSplit( left_tri_count, left_tri_indices, left_bbox, curr_depth + 1,currentNode );
	currentNode->right = ConstructTreeMedianSpaceSplit( right_tri_count, right_tri_indices, right_bbox, curr_depth + 1 ,currentNode);


	return currentNode;
}

KDTNodeM* GLbKdTree::ConstructTreeMedianSpaceSplit(unsigned int num_tris,
	const BoundingBox& bounds)
{
	unsigned int * tris_indics = new unsigned int[num_tris];

	for(unsigned int i = 0;i < num_tris;i++)
	{
		tris_indics[i] = i;
	}

	treeRootM =  ConstructTreeMedianSpaceSplit(num_tris,tris_indics,bounds,0,NULL);

	return treeRootM;
}


KDTNode* GLbKdTree::ConstructTreeSAHSplit(unsigned int num_tris,const BoundingBox& bounds)
{

	// boxEdgeList列表
	//分别保存XYZ三个坐标轴信息
	vv_BoxEdge boxEdgeList(3);

	// 构建模型三角形面片索引
	std::vector<int> triangleIndices;
	for (unsigned int i = 0; i < num_tris; i++) 
	{
		triangleIndices.push_back(i);
	}

	//保存Start end两个点
	unsigned int n = triangleIndices.size();
	for (unsigned int i=0;i<3;i++) 
	{
		boxEdgeList[i].resize(2*n);
	}

	// 初始化edges列表
	for (unsigned int i = 0; i < 3; i++) 
	{
		for (unsigned int j = 0; j < n; j++) 
		{
			boxEdgeList[i][j*2]   = BoxEdge(meshTriangles[triangleIndices[j]].bound._min[i], triangleIndices[j], START, (Axes )i);
			boxEdgeList[i][j*2+1] = BoxEdge(meshTriangles[triangleIndices[j]].bound._max[i], triangleIndices[j], END, (Axes )i);
		}
	}

	// boxEdgeList排序
	sort(boxEdgeList[0].begin(), boxEdgeList[0].end());
	sort(boxEdgeList[1].begin(), boxEdgeList[1].end());
	sort(boxEdgeList[2].begin(), boxEdgeList[2].end());

	//构建树结构
	treeRoot = buildTree_boxEdges(bounds,boxEdgeList,TERMINATION_CRITERIA_D_MAX,NULL);

	return treeRoot;
}

KDTNode * GLbKdTree::buildTree_boxEdges(const BoundingBox& nodeExtent, vv_BoxEdge& boxEdgeList,
	int maxDepth,KDTNode*parent)
{
	//创建节点
	KDTNode * newNode = new KDTNode();

#ifndef KDTREE_NEIGHBORLINKS
	newNode->parent = parent;
#endif

	unsigned int triangles = boxEdgeList[0].size()/2;

	if (maxDepth == 0 || triangles == 0) 
	{
		newNode->left = NULL;
		newNode->right = NULL;
		newNode->box = nodeExtent;
		newNode->triangleIndices = new std::vector<int>();

		for (v_BoxEdge::const_iterator I=boxEdgeList[0].begin(),
			E=boxEdgeList[0].end(); I!=E; I++) 
		{
			if ((*I).edgeType == START) 
			{
				newNode->triangleIndices->push_back((*I).triangleIndex);
			}
		}
		return newNode;
	} 
	else 
	{
		double SAH_best = triangles * sah.m_Ci;
		
		//最佳分割线
		BoxEdge *bestEdge = NULL;

		//SAH 获取最佳的分割线
		for (unsigned int i=0;i<3;i++) 
		{
			unsigned int nA = 0, nB = triangles;
			for (unsigned int j=0;j<boxEdgeList[i].size(); j++) 
			{

				BoxEdge edge = boxEdgeList[i][j];
				if (edge.edgeType == END)
				{
					nB--;
				}

				double SAH_now = sah(nodeExtent, i, nA, nB, edge.splitPlanePosition);
				//遍历得出最小值
				if (SAH_now < SAH_best) 
				{
					SAH_best = SAH_now;
					bestEdge = &boxEdgeList[i][j];
				}

				if (edge.edgeType == START)
				{
					nA++;
				}
			}
		}

		// 此时分割效率不大
		if (!bestEdge)
		{
			newNode->right = NULL;
			newNode->left = NULL;
			newNode->box = nodeExtent;
			newNode->triangleIndices = new std::vector<int>();

			for (v_BoxEdge::const_iterator I=boxEdgeList[0].begin(),
				E=boxEdgeList[0].end(); I!=E; I++)
			{
				if ((*I).edgeType == START) 
				{
					newNode->triangleIndices->push_back((*I).triangleIndex);
				}
			}
			return newNode;
		}

		std::vector<char> membership(triangleSize, 0);
		vv_BoxEdge left(3), right(3); //左右分支boxEdge

		unsigned int left_s = 0, right_s = 0;
		v_BoxEdge::const_iterator I = boxEdgeList[bestEdge->axis].begin(),
			E = boxEdgeList[bestEdge->axis].end();

		//遍历 分开左右子节点
		for (; (&(*I)) != bestEdge; I++) 
		{
			BoxEdge edge = *I;
			if (edge.edgeType == START)
			{
				membership[edge.triangleIndex] += 1;
				left_s++;  //左子树++
			}
		}

		for (++I; I != E; I++)
		{
			BoxEdge edge = *I;
			if (edge.edgeType == END)
			{
				membership[edge.triangleIndex] += 2;
				right_s++; //右子树++
			}
		}

		/////如果分割线将三角形分割了，则三角形放在左右子树下
		for (unsigned int i=0;i<3;i++) 
		{
			for(v_BoxEdge::const_iterator I=boxEdgeList[i].begin(),
				E=boxEdgeList[i].end(); I!=E; I++)
			{
				BoxEdge edge = *I;
				if (membership[edge.triangleIndex] & 1) 
				{
					left[i].push_back(edge);
				}
				if (membership[edge.triangleIndex] & 2) 
				{
					right[i].push_back(edge);
				}
			}
		}

		//左右子树的AABB
		BoundingBox leftNodeExtent(nodeExtent), rightNodeExtent(nodeExtent);
		leftNodeExtent._max[bestEdge->axis] = bestEdge->splitPlanePosition;
		rightNodeExtent._min[bestEdge->axis] = bestEdge->splitPlanePosition;


		newNode->box = nodeExtent;
		//newNode->splitEdge = bestEdge;
		newNode->splitEdge->axis = bestEdge->axis;
		newNode->splitEdge->edgeType = bestEdge->edgeType;
		newNode->splitEdge->splitPlanePosition = bestEdge->splitPlanePosition;
		newNode->splitEdge->triangleIndex = bestEdge->triangleIndex;

		newNode->left = buildTree_boxEdges(leftNodeExtent, left, maxDepth-1,newNode);
		newNode->right = buildTree_boxEdges(rightNodeExtent, right, maxDepth-1,newNode);

		return newNode;
	}

	return newNode; //
}