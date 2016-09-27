
#include "GlbSAH.h"

#include <limits>

using namespace GlbGlobe;

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

double GlbSAH::SA(const osg::BoundingBox&V)
{

	return  2 * (V._max - V._min).x() * (V._max - V._min).y()
		+ 2 * (V._max - V._min).x() * (V._max - V._min).z()
		+ 2 * (V._max - V._min).y() * (V._max - V._min).z();
}

double GlbSAH::P_Vsub(const osg::BoundingBox&Vsub,const osg::BoundingBox&V)
{
	return SA(Vsub)/SA(V);
}

GlbSAH::GlbSAH(double Ct, double Ci, double emptyBonus)
	: m_Ct(Ct), m_Ci(Ci), m_emptyBonus(emptyBonus) {}

double GlbSAH::calculateSAH(const osg::BoundingBox &nodeExtent, char axis,
	unsigned int  nA, unsigned int  nB, double position) const
{
	if (position < nodeExtent._min[axis] || position > nodeExtent._max[axis])
	{
		// outside the bounding box
		return std::numeric_limits<double>::max();
	} 
	else
	{
		double delta[3] = 
		{ nodeExtent._max[0] - nodeExtent._min[0],
		nodeExtent._max[1] - nodeExtent._min[1],
		nodeExtent._max[2] - nodeExtent._min[2] };

		double oneOverTotalSurfaceArea = 0.5f /
			(delta[0]*delta[1] + delta[1]*delta[2] + delta[2]*delta[0]);

		int otherAxis[2] = { (axis+1)%3 , (axis+2)%3 };

		double crossSectionArea = delta[otherAxis[0]]*delta[otherAxis[1]];

		if (nA == 0 && position > nodeExtent._min[axis]) 
		{
			double pB = 2*(crossSectionArea + (nodeExtent._max[axis] -
				position) *
				(delta[otherAxis[0]] +
				delta[otherAxis[1]]))
				* oneOverTotalSurfaceArea;

			return this->m_Ct + this->m_Ci*pB*nB*(1.0f-this->m_emptyBonus);

		}
		else if (nB == 0 && position < nodeExtent._max[axis]) 
		{
			double pA = 2*(crossSectionArea +
				(position - nodeExtent._min[axis])*
				(delta[otherAxis[0]] + delta[otherAxis[1]]))*
				oneOverTotalSurfaceArea;

			return this->m_Ct + this->m_Ci*pA*(nA+nB)*(1.0f-this->m_emptyBonus);

		} 
		else
		{
			//左子树 进行遍历的概率
			double pA = 2*(crossSectionArea +
				(position - nodeExtent._min[axis])*
				(delta[otherAxis[0]] + delta[otherAxis[1]]))*
				oneOverTotalSurfaceArea;

			//右子树 进行遍历的概率
			double pB = 2*(crossSectionArea +
				(nodeExtent._max[axis] - position)*
				(delta[otherAxis[0]] + delta[otherAxis[1]]))*
				oneOverTotalSurfaceArea;

			/*
			//总共消耗时间(当前节点的遍历时间 和 左右子节点的消耗时间<该节点的遍历、子节点的消耗时间>)
			//遍历时间 、求交时间
			//f(b) = LSA(b) * NL(b) + RSA(b) * (N - NL) 
			b:代表splitting plane的位置（0,1)
			LSA:左子树的节点面积,
			RSA:右子树的节点面积,
			NL:左子树的节点个数
			*/
			return this->m_Ct + this->m_Ci*(pA*nA + pB*nB);

		}
	}
}

double GlbSAH::operator()(const osg::BoundingBox &nodeExtent, char axis,
	unsigned int  nA, unsigned int  nB, double position) const
{
	//计算cost time
	return calculateSAH(nodeExtent, axis, nA, nB, position);
}
