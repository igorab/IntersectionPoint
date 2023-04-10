#pragma once
#include "IntersectionObjects.cpp"

class SegmentSegmentIntersection
{
private:
	Segment segmentAB;
	Segment segmentCD;

	Point A, B, 
		  C, D;

	Vector P;

public:
	
	SegmentSegmentIntersection(Segment _segmentAB, Segment _segmentCD) : segmentAB(_segmentAB), segmentCD(_segmentCD) 
	{ 
		A = segmentAB.getA(); B = segmentAB.getB(); C = segmentCD.getA(); D = segmentCD.getB(); 
	};
	SegmentSegmentIntersection(Point _pointA, Point _pointB, Point _pointC, Point _pointD) :  segmentAB(Segment(_pointA, _pointB)), segmentCD(Segment(_pointC, _pointD)), A(_pointA), B(_pointB), C(_pointB), D(_pointD) {};


	static num Area2(Point _A, Point _B, Point _C)
	{
		num area;

		area = (_B.X - _A.X) * (_C.Y - _A.Y) - (_C.X - _A.Y) * (_B.Y - _A.Y);

		return area;
	}

	static bool Collinear(Point _A, Point _B, Point _C)
	{
		return Area2(_A, _B, _C) >= 0;
	}


	char ParallelInt()
	{
		if (!Collinear(A, B, C))
			return '0';

		return '0';
	}

	char Calc()
	{
		num s, t;
		num num, denom;
		
		char code = '?';

		denom = A.X * (D.Y - C.Y) + B.X * (C.Y - D.Y) + D.X * (B.Y - A.Y) + C.X * (A.Y - B.Y);

		if (denom == 0)
		{
			return ParallelInt();
		}

		num  = A.X * (D.Y - C.Y) + C.X * (A.Y - D.Y) + D.X * (C.Y - A.Y);
			
		if (num == 0.0 || num == denom)
			code = 'v';

		s = num / denom;

		num = - (A.X * (C.Y - B.Y) + B.X * (A.Y - C.Y) + C.X * (B.Y - A.Y));
		if (num == 0.0 || num == denom)
			code = 'v';

		t = num / denom;

		if (0. < s && s < 1. && 0. < t && t < 1.)
			code = '1';
		else if (0. > s || s > 1. || 0. > t || t > 1.)
			code = '0';

		Vector AB = Vector(B) - Vector(A);
		P = A + s * AB;

		return code;
	}


};