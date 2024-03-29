#include "IntersectionObjects.cpp"

enum class IntersectionType {	
	// Check plane: 	
	p, //  Segment lies wholly within the plane  
	q, //  Q from-point is on the plane 
	r, //  Second R endpoint is on the plane 
	p0, // Segment parallel - lies strictly to one side or the other of the plane. 
	p1, // The segment intersects the plane

	//  Check  segment / triangle:
	v,  // segment includes a vertex of triangle 
	e,  // segment includes a point of an edge of triangle 
	f,  // segment includes a point in the face of triangle. 
	f0, // segment lays outside the triangle 

	//  Check point / triangle:
	V, // p coincides with a Vertex of T. 
	E, // p is in the relative interior of an Edge of T.
	F, // p is in the relative interior of a Face of T.
	T0, // point lays outside

	// Check segment/segment
	s_e, // segments overlap
	s_v,  // trough vertex
	s_1,  // intersects properly
	s_0	 // do not intersect
};

const char* txt_result[] = {"The segment lies wholly within the plane", 
					 "First Q point is on the plane",   
					 "Second R endpoint is on the plane",
				     "Segment parallel to triangle",
					 "Segment intersects the plane",
					 "Segment includes a vertex of triangle",
					 "Segment includes a point in the relative interior of an edge of triangle",
					 "Segment includes a point in the relative interior of a face of triangle",
					 "Segment lays outside the triangle",
					 "Point coincides with a Vertex of triangle",
					 "Point is in the relative interior of an Edge of triangle",
					 "Point is in the relative interior of a Face of triangle",
					 "Point lays outside triangle",
					 "Segment collinearly overlap",
					 "Segment endpoint is on the other segment",
				     "Segments intersect properly",
					 "Segments do not intersect"
					 };


ostream& operator<< (std::ostream& out, const Point& point)
{
	out << "Point(" << point.X << ", " << point.Y << ", " << point.Z << ')';
	return out;
};

istream& operator>> (istream& in, Point& point)
{
	in >> point.X >> point.Y >> point.Z;
	return in;
};


class SegmentSegmentIntersection
{
private:
	Segment segmentAB;
	Segment segmentCD;

	Point A, B,
		C, D;

	Vector P;
public:

	Point getP() { return P.getPoint(); }

	SegmentSegmentIntersection(Segment _segmentAB, Segment _segmentCD) : segmentAB(_segmentAB), segmentCD(_segmentCD)
	{
		A = segmentAB.getA(); B = segmentAB.getB(); 
		C = segmentCD.getA(); D = segmentCD.getB();
	};
	SegmentSegmentIntersection(Point _pointA, Point _pointB, Point _pointC, Point _pointD) : segmentAB(Segment(_pointA, _pointB)), segmentCD(Segment(_pointC, _pointD)), A(_pointA), B(_pointB), C(_pointC), D(_pointD) {};

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

	static bool Between(Point _A, Point _B, Point _C)
	{
		bool ret = false;

		if (_A.X != _B.X)
		{
			ret = (_A.X <= _C.X && _C.X <= _B.X || _A.X >= _C.X && _C.X >= _B.X);
		}
		else
		{
			ret = (_A.Y <= _C.Y && _C.Y <= _B.Y || _A.Y >= _C.Y && _C.Y >= _B.Y);
		}

		return ret;
	}

	IntersectionType ParallelInt()
	{
		IntersectionType ret = IntersectionType::s_e;

		if (!Collinear(A, B, C))
			return IntersectionType::s_0;

		if (Between(A, B, C))
		{
			P = Point(C);
			return ret;
		}
		else if (Between(A, B, D))
		{
			P = Point(D);
			return ret;
		}
		else if (Between(C, D, A))
		{
			P = Point(A);
			return ret;
		}
		else if (Between(C, D, B))
		{
			P = Point(B);
			return ret;
		}

		return ret;
	}

	IntersectionType Calc()
	{
		num s, t; //coeff in line parametric equation
		
		IntersectionType intype = IntersectionType::s_0;
	
		num denom = A.X * (D.Y - C.Y) + B.X * (C.Y - D.Y) + D.X * (B.Y - A.Y) + C.X * (A.Y - B.Y);

		if (abs(denom) < EPSILON) // denom == 0
		{
			return ParallelInt();
		}

		num num_s = A.X * (D.Y - C.Y) + C.X * (A.Y - D.Y) + D.X * (C.Y - A.Y);

		if (abs(num_s) < EPSILON || abs(num_s - denom) < EPSILON )  //num_s == 0.0 || num_s == denom
			intype = IntersectionType::s_v;

		s = num_s / denom;

		num num_t = -(A.X * (C.Y - B.Y) + B.X * (A.Y - C.Y) + C.X * (B.Y - A.Y));

		if (abs(num_t) < EPSILON || abs(num_t - denom) < EPSILON) // num_t == 0.0 || num_t == denom
			intype = IntersectionType::s_v;

		t = num_t / denom;

		if (0. < s && s < 1. && 0. < t && t < 1.)
			intype = IntersectionType::s_1;
		else if (0. > s || s > 1. || 0. > t || t > 1.)
			intype = IntersectionType::s_0;

		Vector AB = Vector(B) - Vector(A);
		P = A + s * AB;

		return intype;
	}
};

class SegmentTriangleIntersection
{
private:
	Segment* segment;
	Triangle* triangle;

	Point iPoint;

public:

	SegmentTriangleIntersection(Segment* _segment, Triangle* _triangle)
	{
		segment = _segment;
		triangle = _triangle;
		iPoint = Point(0, 0, 0);
	}

	SegmentTriangleIntersection(Point _from, Point _to, Point _A, Point _B, Point _C)
	{
		segment = new Segment(_from, _to);
		triangle = new Triangle(_A, _B, _C);
		iPoint = Point(0, 0, 0);
	}

	~SegmentTriangleIntersection()
	{
		delete segment;
		delete triangle;
	}

	Point* getIntersectionPoint()
	{ 
		return &iPoint; 
	}

	
	static int AreaSign(Point _A, Point _B, Point _C)
	{
		num area2;

		area2 = (_B.X - _A.X) * (_C.Y - _A.Y) - (_C.X - _A.X) * (_B.Y - _A.Y);

		if (area2 >= 0.5) 
			return 1;
		else if (area2 <= -0.5) 
			return -1;
		
		return 0;
	}

	static int VolumeSign(Vector v_A, Vector v_B, Vector v_C, Vector v_P )
	{
		num vol;						
		Vector v_dA = v_A - v_P;
		Vector v_dB = v_B - v_P;
		Vector v_dC = v_C - v_P;

		vol = (v_dA * (v_dB ^ v_dC));

		if (vol >= 0.5)
			return 1;
		if (vol <= -0.5)
			return -1;
		return 0;
	}

	// ���������� �� ����� ��������� ����������
	int PlaneCoefficients(double* D)
	{
		int i;
		double t;
		double biggest = 0.0;
		int max_XYZ = 0;

		Vector v_N = triangle->norm();

		Vector v_A(triangle->getA());

		*D = (v_A * v_N);

		num N[] = {v_N.X, v_N.Y, v_N.Z};

		for (i = 0; i < dim3D; i++) 
		{
			t = fabs(N[i]);

			if (t > biggest) 
			{
				biggest = t;
				max_XYZ = i;
			}
		}
		return max_XYZ;
	}
	

	IntersectionType IntersectionTriangle2D(Point _projPointXY, Triangle _triangleXY)
	{
		int area0_sign, area1_sign, area2_sign; 
		Point A_2D = _triangleXY.getA();
		Point B_2D = _triangleXY.getB();
		Point C_2D = _triangleXY.getC();

		area0_sign = AreaSign(_projPointXY, A_2D, B_2D);
		area1_sign = AreaSign(_projPointXY, B_2D, C_2D);
		area2_sign = AreaSign(_projPointXY, C_2D, A_2D);

		if (area0_sign == 0 && area1_sign == 0 && area2_sign == 0)
			exit(EXIT_FAILURE);

		if (area0_sign == 0 && area1_sign == 0 || area0_sign == 0 && area2_sign == 0 || area1_sign == 0 && area2_sign == 0)
			return IntersectionType::V; //coincides with a Vertex of T. 

		if ((area0_sign == 0 && area1_sign > 0 && area2_sign > 0) || (area1_sign == 0 && area0_sign > 0 && area2_sign > 0) || (area2_sign == 0 && area0_sign > 0 && area1_sign > 0))
			return IntersectionType::E; // in the relative interior of an Edge of T.

		if ((area0_sign == 0 && area1_sign < 0 && area2_sign > 0) || (area1_sign == 0 && area0_sign < 0 && area2_sign > 0) || (area2_sign == 0 && area0_sign < 0 && area1_sign > 0))
			return IntersectionType::E; //in the relative interior of an Edge of T.

		if ((area0_sign > 0 && area1_sign > 0 && area2_sign > 0) || (area0_sign < 0 && area1_sign < 0 && area2_sign < 0))
			return IntersectionType::F; //in the relative interior of a Face of T.
				
		return IntersectionType::T0;
	}

	IntersectionType IntersectionTriangle3D(Point point, int m_XYZ)
	{
		Point projPointXY;
			
		if (m_XYZ == 0) // X
		{
			projPointXY = Point(point.Y, point.Z, 0);
		}
		else if (m_XYZ == 1) // Y
		{
			projPointXY = Point(point.X, point.Z, 0);
		}
		else if (m_XYZ == 2) // Z
		{
			projPointXY = Point(point.X, point.Y, 0);
		}

		Triangle triangleProjected = triangle->Triangle_Proj2D(m_XYZ);

		IntersectionType itypePoint =  IntersectionTriangle2D(projPointXY, triangleProjected);

		return itypePoint;
	}
		
	IntersectionType SegmentPlaneIntersection(int* _m_XYZ)
	{
		Vector v_q  = segment->Q;
		Vector v_r  = segment->R;
		Vector v_qr = segment->dir;
		Vector v_N  = triangle->norm();

		Vector v_P;
		double D = 0;		
		double num, denom, t;
		int i;

		*_m_XYZ = PlaneCoefficients(&D);

		num = D - (v_q * v_N);
		denom = (v_qr * v_N);

		if (abs(denom) >= EPSILON) //denom != 0
		{
			t = num / denom;
		}
		else
		{
			if (abs(num) < EPSILON) //num == 0
				return IntersectionType::p;
			else
				return IntersectionType::p0;
		}
		 
		v_P = v_q + t * v_qr;

		iPoint = v_P.getPoint(); // ����������� � ����������
		
		if (t > EPSILON && t < 1.0-EPSILON) // t>0 && t<1
			return IntersectionType::p1;

		if (abs(num) <= EPSILON) // num == 0
			return IntersectionType::q;

		if (abs(num - denom) <= EPSILON ) // num == denom
			return IntersectionType::r;
				
		return IntersectionType::p0;
	}

	IntersectionType SegmentTriangleCross()
	{
		int vol0, vol1, vol2;

		Vector v_A = triangle->getA(), 
			   v_B = triangle->getB(), 
			   v_C = triangle->getC();


		vol0 = VolumeSign(segment->Q, v_A, v_B, segment->R);
		vol1 = VolumeSign(segment->Q, v_B, v_C, segment->R);
		vol2 = VolumeSign(segment->Q, v_C, v_A, segment->R);

		if (vol0 == 0 && vol1 == 0 && vol2 == 0)
			exit(EXIT_FAILURE);

		//Same sign: segment intersects interior of triangle. 
		if ((vol0 > 0 && vol1 > 0 && vol2 > 0) || (vol0 < 0 && vol1 < 0 && vol2 < 0))
			return IntersectionType::f;

		//Opposite sign: no intersection between segment and triangle. 
		if ((vol0 > 0 || vol1 > 0 || vol2 > 0) && (vol0 < 0 || vol1 < 0 || vol2 < 0))
			return IntersectionType::f0;

		//Two zeros: segment intersects vertex. 
		if ((vol0 == 0 && vol1 == 0) || (vol0 == 0 && vol2 == 0) || (vol1 == 0 && vol2 == 0))
			return IntersectionType::v;

		//One zero : segment intersects edge
		if (vol0 == 0 || vol1 == 0 || vol2 == 0)
			return IntersectionType::e;

		return IntersectionType::f0;
	}

	// lies entirely in the plane 
	// 
	// ***possible many points of intersection: 2 edges of triangle or lies on edge. I show the first founded point
	//
	IntersectionType InPlane()
	{
		Point A = triangle->getA();
		Point B = triangle->getB();
		Point C = triangle->getC();

		IntersectionType itype;

		Segment side1(A, B);
		SegmentSegmentIntersection ssI_1(*segment, side1);
		itype = ssI_1.Calc();
		if (itype != IntersectionType::s_0)
		{
			iPoint = ssI_1.getP();
			return IntersectionType::p;
		}
		
		Segment side2(B, C);
		SegmentSegmentIntersection ssI_2(*segment, side2);
		itype = ssI_2.Calc();
		if (itype != IntersectionType::s_0)
		{
			iPoint = ssI_2.getP();
			return IntersectionType::p;
		}

		Segment side3(C, A);
		SegmentSegmentIntersection ssI_3(*segment, side3);
		itype = ssI_3.Calc();
		if (itype != IntersectionType::s_0)
		{
			iPoint = ssI_3.getP();
			return IntersectionType::p;
		}

		if (itype == IntersectionType::s_0)
		{
			return IntersectionType::f0;
		}

		return IntersectionType::p;
	}

	IntersectionType IntersectionCalculate()
	{		
		int max_XYZ = 0;
		Vector v_P;
		 
		IntersectionType itypePlane = SegmentPlaneIntersection(&max_XYZ);
		IntersectionType itypeTriangle = IntersectionType::T0;

		if (itypePlane == IntersectionType::q)
		{
			Point Q = segment->getA();
			itypeTriangle = IntersectionTriangle3D(Q, max_XYZ);
		}
		else if (itypePlane == IntersectionType::r)
		{
			Point R = segment->getB();
			itypeTriangle = IntersectionTriangle3D(R, max_XYZ);
		}
		else if (itypePlane == IntersectionType::p)
		{
			itypeTriangle = InPlane();
		}
		else if (itypePlane == IntersectionType::p0)
		{
			itypeTriangle = itypePlane;
		}
		else if (itypePlane == IntersectionType::p1)
		{
			itypeTriangle = SegmentTriangleCross();
		}

		return itypeTriangle;					
	}
	
	// no use in this task, only to compare results
	//Moeller�Trumbore intersection algorithm 
	static bool RayIntersectsTriangle(Segment* segment,
		Triangle* inTriangle,
		Vector& outIntersectionPoint)
	{
		Vector rayOrigin = segment->getA();
		Vector rayVector = segment->getB();


		Vector vertex0 = inTriangle->getA();
		Vector vertex1 = inTriangle->getB();
		Vector vertex2 = inTriangle->getC();

		Vector edge1, edge2, h, s, q;

		num a, f, u, v;
		edge1 = vertex1 - vertex0;
		edge2 = vertex2 - vertex0;

		h = (rayVector ^ edge2);

		a = (edge1 * h);

		if (a > -EPSILON && a < EPSILON)
			return false;    // This ray is parallel to this triangle.

		f = 1.0 / a;
		s = rayOrigin - vertex0;

		u = f * (s * h);

		if (u < 0.0 || u > 1.0)
			return false;

		q = (s ^ edge1);

		v = f * (rayVector * q);

		if (v < 0.0 || u + v > 1.0f)
			return false;

		// At this stage we can compute t to find out where the intersection point is on the line.
		num t = f * (edge2 * q);

		if (t > EPSILON) // ray intersection
		{
			outIntersectionPoint = rayOrigin + (t * rayVector);
			return true;
		}
		else
		{
			// This means that there is a line intersection but not a ray intersection.
			return false;
		}
	}

	static bool is_ray_cross_triangle(Segment* _segment, Triangle* _triangle)
	{
		Vector intersectionPoint;
		
		RayIntersectsTriangle(_segment, _triangle, intersectionPoint);

		std::cout << "Intersection: " << intersectionPoint.getPoint() << std::endl;

		return true;
	}
};

#include "IntersectionTest.cpp"

int main()
{		
	const Point pA(vA), pB(vB), pC(vC);	
	const Point pntFrom(fromA), pntTo(toB);

	Segment* segment = new Segment(pntFrom, pntTo);
	Triangle* triangle = new Triangle(pA, pB, pC);
	SegmentTriangleIntersection* segmentTriangleIntersection = new SegmentTriangleIntersection(segment, triangle);

	if (segment->len() == 0 || triangle->area() == 0)
	{		
		cout << "incorrect size";
	}
	else
	{			
		IntersectionType inResult = segmentTriangleIntersection->IntersectionCalculate();
							
		if (inResult == IntersectionType::F || inResult == IntersectionType::E || inResult == IntersectionType::V || 
			inResult == IntersectionType::f || inResult == IntersectionType::e || inResult == IntersectionType::v ||
			inResult == IntersectionType::p)
		{
			cout << "Result: triangle intersection." << endl;
			cout << "Intersection " << *segmentTriangleIntersection->getIntersectionPoint() << endl ;
		}
		else
		{
			cout << "Result: no intersection " << endl;
		}

		cout << endl << "Description: " << txt_result[(int)inResult] << endl;				
	}	
	delete segmentTriangleIntersection;

	return 0;
}

