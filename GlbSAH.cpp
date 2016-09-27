
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
			//������ ���б����ĸ���
			double pA = 2*(crossSectionArea +
				(position - nodeExtent._min[axis])*
				(delta[otherAxis[0]] + delta[otherAxis[1]]))*
				oneOverTotalSurfaceArea;

			//������ ���б����ĸ���
			double pB = 2*(crossSectionArea +
				(nodeExtent._max[axis] - position)*
				(delta[otherAxis[0]] + delta[otherAxis[1]]))*
				oneOverTotalSurfaceArea;

			/*
			//�ܹ�����ʱ��(��ǰ�ڵ�ı���ʱ�� �� �����ӽڵ������ʱ��<�ýڵ�ı������ӽڵ������ʱ��>)
			//����ʱ�� ����ʱ��
			//f(b) = LSA(b) * NL(b) + RSA(b) * (N - NL) 
			b:����splitting plane��λ�ã�0,1)
			LSA:�������Ľڵ����,
			RSA:�������Ľڵ����,
			NL:�������Ľڵ����
			*/
			return this->m_Ct + this->m_Ci*(pA*nA + pB*nB);

		}
	}
}

double GlbSAH::operator()(const osg::BoundingBox &nodeExtent, char axis,
	unsigned int  nA, unsigned int  nB, double position) const
{
	//����cost time
	return calculateSAH(nodeExtent, axis, nA, nB, position);
}
