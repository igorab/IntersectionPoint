#pragma once
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
	T0 // point lays outside
	
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
					 "Point lays outside triangle"
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

	// определ€ем на какую плоскость проецируем
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

		for (i = 0; i < DIM; i++) 
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

		iPoint = v_P.getPoint(); // пересечение с плоскостью
		
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

	//lies entirely in the plane 
	IntersectionType InPlane()
	{
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
	//MoellerЦTrumbore intersection algorithm 
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

//Test 1, IntersectionType::F 
// конец отрезка попал в плоскоть треугольника
#define vA 2, 0, 1
#define vB 0, 1, 0
#define vC 0, 0, 3 

#define fromA 0, 0, 0
#define toB 0.6, 0.6, 0.6 

//Test 2, IntersectionType::E
// конец отрезка попал в грань
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 1, 1, 1 

#define fromA 0, 0, 0
#define toB 0.5, 0.5, 0

//Test 3, IntersectionType::V
// конец отрезка попал в вершину
#define vA 3, 0, 0
#define vB 0, 1, 0
#define vC 1, 1, 1 

#define fromA 0, 0, 0
#define toB 1., 1., 1.

//Test 4, IntersectionType::F 
// начало отрезка находитс€ в площади треугольника
#define vA 2, 0, 1
#define vB 0, 1, 0
#define vC 0, 0, 3 

#define fromA 0.6, 0.6, 0.6
#define toB 1., 1., 1.


//Test 4, IntersectionType::V 
// начало отрезка находитс€ в вершине треугольника
#define vA 2, 0, 1
#define vB 0, 1, 0
#define vC 0, 0, 3 

#define fromA 2, 0., 1
#define toB 5, 5, 6

//Test 5, IntersectionType::V 
// отрезок пересекает вершину треугольника
#define vA 2, 1, 1
#define vB 1, 4, -1
#define vC 0.5, 0.5, 7 

#define fromA 0, 0., 0
#define toB 10, 10, 10

//Test 6, IntersectionType::f 
// отрезок пересекает плоскость треугольника
#define vA 2, 1, 1
#define vB 1, 4, -1
#define vC 0.5, 0.5, 7 

#define fromA 10, 10., 10
#define toB 0, 0, 0

//Test 7, IntersectionType::v 
// отрезок пересекает вершину треугольника
#define vA 2, 1, 1
#define vB 1, 4, -1
#define vC 0.5, 0.5, 7 

#define fromA 0, -1, -1
#define toB 4, 3, 3

//Test 8, IntersectionType::e 
// отрезок пересекает грань треугольника
#define vA 0, 0, 0
#define vB 3, 5, 1
#define vC 0, 0, 7 

#define fromA -2, -2, 0
#define toB 2, 2, 4

//Test 9, 
// нет пересечени€
#define vA 0, 0, 0
#define vB 3, 5, 1
#define vC 0, 0, 7 

#define fromA 4, 6, 3
#define toB 9, 12, 24

//Test 10, 
// параллельно плоскости
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 0, 0, 7 

#define fromA 1, 0, 0
#define toB 1, 4, 5


//Test 11, 
// лежит в плоскости
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 0, 0, 3 

#define fromA 0, 0, 0
#define toB 0, 4, 5



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
			inResult == IntersectionType::f || inResult == IntersectionType::e || inResult == IntersectionType::v)
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

