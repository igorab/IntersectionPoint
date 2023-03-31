// �������� ��� �� �++��� ����������� ����� ����������� ������� �
// ������������ � 3D ������������.������� ����� ������������ ������,
// ����������� ����� ������������ ���� ���� �����.

#include <iostream>

using namespace std;
typedef float num;
const int DIM = 3;
const float EPSILON = 0.0000001;

typedef int tPointi[DIM];
typedef double tPointd[DIM];


enum class IntersectionType{
	p, //The segment lies wholly within the plane.  
	q, //The (first) q endpoint is on the plane (but not 'p'). 
	r, //The (second) r endpoint is on the plane (but not 'p'). 
	p0, //The segment lies strictly to one side or the other of the plane. 
	p1  //The segment intersects the plane, and none of {p, q, r} hold. 
};


struct IShape
{
	virtual num len() = 0;
	virtual num area() = 0;
};

// ����� � ������������
struct Point : IShape
{
	num X = 0;
	num Y = 0;
	num Z = 0;
	
	explicit Point(num _X = 0, num _Y = 0, num _Z = 0) : X(_X), Y(_Y), Z(_Z) {}

	friend ostream& operator<< (ostream& out, const Point& point);

	friend istream& operator>> (istream& in, Point& point);

	num* getTPoint()
	{
		static num Tp[DIM] = { X, Y, Z };
		return Tp;
	}

	num len() { return 0; }
	num area() { return 0; }
};

// ������
struct Vector
{
	num X, Y, Z;

	Vector() : X(0), Y(0), Z(0) {}
	Vector(num _X, num _Y, num _Z) : X(_X), Y(_Y), Z(_Z) {}
	Vector(Point _point) {
		X = _point.X;
		Y = _point.Y;
		Z = _point.Z;
	}
	Vector(Point _from, Point _to) {
		X = _to.X - _from.X;
		Y = _to.Y - _from.Y;
		Z = _to.Z - _from.Z;
	}

	Point getPoint() { return Point(X, Y, Z); }

	// ������ �������
	num len() const {
		num sumV = X * X + Y * Y + Z * Z;
		return sqrt(sumV);
	}

	// ������������ ������� �� �����
	Vector mult3(const num C)
	{
		Vector vc(C * X, C * Y, C * Z);
		return vc;
	}

	//�������� ��������
	void sum3(const Vector _V, int sign)
	{
		X += sign * _V.X;
		Y += sign * _V.Y;
		Z += sign * _V.Z;
	}

	num dot3(const Vector _b)
	{
		return X * _b.X + Y * _b.Y + Z * _b.Z;
	}

	//��������� ������������
	Vector cross3(const Vector _V)
	{
		Vector vNorm(Y * _V.Z - Z * _V.Y,
			Z * _V.X - X * _V.Z,
			X * _V.Y - Y * _V.X);
		return vNorm;
	}

	num det3(Vector a, Vector _V)
	{
		num det;
		det = a.X * (Y * _V.Z - Z * _V.Y) +
			a.Y * (Z * _V.X - X * _V.Z) +
			a.Z * (X * _V.Y - Y * _V.X);
		return det;
	}

	friend Vector operator+ (Vector const _A, Vector const _B) { return Vector(_A.X + _B.X, _A.Y + _B.Y, _A.Z + _B.Z); }
	friend Vector operator- (Vector const _A, Vector const _B) { return Vector(_A.X - _B.X, _A.Y - _B.Y, _A.Z - _B.Z); }
	// ��������� ������������
	friend num operator*(Vector const _a, Vector const _b) { return _a.X * _b.X + _a.Y * _b.Y + _a.Z * _b.Z; } 
	// ��������� �� �����
	friend Vector operator*(num const _c, Vector const _V) { 
		Vector V(_V);
		Vector c_V = V.mult3(_c);
		return c_V;
	} 
	//��������� ������������
	friend Vector operator^(Vector const _a, Vector const _b) {
		Vector A(_a);
		Vector normV = A.cross3(_b);		
		return normV;
	}
};

// �������
class Segment :IShape
{
private:
	Point L_A;
	Point L_B;
public:

	Point getA() const { return L_A; };
	Point getB() const { return L_B; };

	Vector Q;
	Vector R;

	// ������������ ������
	Vector dir;

	explicit Segment(Point _A, Point _B) : L_A(_A), L_B(_B) {
		dir = Vector(L_A, L_B);
		Q = Vector(L_A);
		R = Vector(L_B);
	}

	num len() { return dir.len(); }
	num area() { return 0; }
};


//    �����������
class Triangle : IShape
{
private:
	Point P_A;
	Point P_B;
	Point P_C;

public:

	Vector V_AB;
	Vector V_AC;
	Vector V_BC;
	
	explicit Triangle( Point _A, Point _B, Point _C) : P_A(_A), P_B(_B), P_C(_C)
	{
		V_AB = Vector(P_A, P_B);
		V_AC = Vector(P_A, P_C);		
		V_BC = Vector(P_B, P_C);
	}
	
	Point getA() const { return P_A; };
	Point getB() const { return P_B; };
	Point getC() const { return P_C; };

	// ������� ��������� ��������� ����� 3 �����
	// 
	//  ������ �������
	Vector norm()
	{
		Vector V_N = V_AB ^ V_AC;
		return V_N;
	}

	friend ostream& operator<<(ostream& out, Triangle& _T) {
		out << _T.P_A << " " << _T.P_B << " " << _T.P_C << endl;
		return out;
	}

	Triangle Triangle_Proj(int _O_XYZ)
	{
		Point _A, _B, _C;
		Triangle tr(_A, _B, _C);
		
		if (_O_XYZ == 0)
		{
			_A = Point(0, P_A.Y, P_A.Z);
			_B = Point(0, P_B.Y, P_B.Z);
			_C = Point(0, P_C.Y, P_C.Z);
			tr = Triangle(_A, _B, _C);
		}
		else if (_O_XYZ == 1)
		{
			_A = Point(P_A.X, 0, P_A.Z);
			_B = Point(P_B.X, 0, P_B.Z);
			_C = Point(P_C.X, 0, P_C.Z);
			tr = Triangle(_A, _B, _C);
		}
		else if (_O_XYZ == 2)
		{
			_A = Point(P_A.X, P_A.Y, 0);
			_B = Point(P_B.X, P_B.Y, 0);
			_C = Point(P_C.X, P_C.Y, 0);
			tr = Triangle(_A, _B, _C);
		}

		return tr;
	}

	num len() { return V_AB.len() + V_AC.len() + V_BC.len(); }
	num area() { return (V_AB ^ V_AC).len() / 2 ; }
};

ostream& operator<< (std::ostream& out, const Point& point)
{
	out << "Point(" << point.X << ", " << point.Y << ", " << point.Z << ')';
	return out;
}

istream& operator>> (istream& in, Point& point)
{
	in >> point.X >> point.Y >> point.Z ;
	return in;
}

bool RayIntersectsTriangle(Segment *segment,
						   Triangle* inTriangle,
						   Vector& outIntersectionPoint)
{	
	Vector rayOrigin = segment->getA() ;
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


	static num AreaSign(Point A, Point B, Point C)
	{
		num area2;

		area2 = (B.X - A.X) * (C.Y - A.Y) - (C.X - A.X) * (B.Y - A.Y);

		if (area2 > 0.5) 
			return 1;
		else if (area2 < 0.5) 
			return -1;
		else 
			return 0;
	}

	int PlaneCoefficients(double* D)
	{
		int i;
		double t;
		double biggest = 0.0;
		int m = 0;

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
				m = i;
			}
		}
		return m;
	}
	
	char IntersectionTriangle2D(Point pp, Triangle _triangle)
	{
		int area0, area1, area2;
		
		Point Tp[] = { _triangle.getA(), _triangle.getB(), _triangle.getC()};

		area0 = AreaSign(pp, Tp[0], Tp[1]);
		area1 = AreaSign(pp, Tp[1], Tp[2]);
		area2 = AreaSign(pp, Tp[2], Tp[0]);

		if ((area0 == 0 && area1 > 0 && area2 > 0) || (area1 == 0 && area0 > 0 && area2 > 0) || (area2 == 0 && area0 > 0 && area1 > 0))
			return 'E';

		if ((area0 == 0 && area1 < 0 && area2 > 0) || (area1 == 0 && area0 < 0 && area2 > 0) || (area2 == 0 && area0 < 0 && area1 > 0))
			return 'E';

		if ((area0 > 0 && area1 > 0 && area2 > 0) || (area0 < 0 && area1 < 0 && area2 < 0))
			return 'F';

		if (area0 == 0 && area1 == 0 && area2 == 0)
			exit(EXIT_FAILURE);

		if (area0 == 0 && area1 == 0 || area0 == 0 && area2 == 0 || area1 == 0 && area2 == 0)
			return 'V';

		return '0';
	}

	char IntersectionTriangle3D(Point point, int m_XYZ)
	{
		Point projectedPoint;
		int i, j, k;
		num* p;
		num pp[DIM];
		p = point.getTPoint();

		Point tr[] = { triangle->getA(), triangle->getB(), triangle->getC() };

		Point projectedTriangle[DIM];

		j = 0;
		for (i = 0; i < DIM; i++)
		{
			if (i != m_XYZ)
			{
				//pp[j] = p[i];
				
				projectedTriangle[i] = tr[i];

				j++;
			}
		}

		Triangle triangleProjected = triangle->Triangle_Proj(m_XYZ);

		return IntersectionTriangle2D(projectedPoint, triangleProjected);
	}
		
	char SegmentPlaneIntersection(int* _m_XYZ)
	{
		Vector v_q  = segment->Q;
		Vector v_r  = segment->R;
		Vector v_rq = segment->dir;
		Vector v_N  = triangle->norm();

		Vector v_P;
		double D = 0;		
		double num, denom, t;
		int i;

		*_m_XYZ = PlaneCoefficients(&D);

		num = D - (v_q * v_N);
		denom = (v_rq * v_N);

		if (denom != 0)
		{
			t = num / denom;
		}
		else
		{
			return (num == 0) ? 'p' : '0';
		}
		 
		v_P = v_q + t * v_rq;

		iPoint = v_P.getPoint(); // ����������� � ����������

		if (t > 0 && t < 1)
			return '1';

		if (num == 0)
			return 'q';

		if (num == denom)
			return 'r';

		return '0';
	}

	char SegmentTriangleCross()
	{
		int volO, voll, vol2;

		return '0';
	}

	//lies entirely in the plane 
	char InPlane()
	{
		return '0';
	}

	char IntersectionCalculate()
	{
		int code;
		int m_XYZ;
		Vector v_P;

		code = SegmentPlaneIntersection(&m_XYZ);

		if (code == 'q')
		{
			return IntersectionTriangle3D(segment->getA(), m_XYZ);
		}
		else if (code == 'r')
		{
			return IntersectionTriangle3D(segment->getB(), m_XYZ);
		}
		else if (code == 'p')
		{
			return InPlane();
		}
		else 
		{
			return SegmentTriangleCross();
		}			
	}
		
	
	// ����� ���������� ����������� ?
	static bool is_ray_cross_triangle(Segment* _segment, Triangle* _triangle)
	{
		Vector intersectionPoint;
		
		RayIntersectsTriangle(_segment, _triangle, intersectionPoint);

		std::cout << "Intersection: " << intersectionPoint.getPoint() << std::endl;

		return true;
	}

};


// triangle
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 0, 0, 1 
// segment
#define fromA 0, 0, 0
#define toB 1, 1, -1 

int main()
{	
	system("chcp 1251");

	const Point pA(vA), pB(vB), pC(vC);	
	const Point pntFrom(fromA), pntTo(toB);

	Segment* segment = new Segment(pntFrom, pntTo);
	Triangle* triangle = new Triangle(pA, pB, pC);

	if (segment->len() == 0 || triangle->area() == 0)
	{
		cout << "incorrect size";
	}
	else
	{	
		SegmentTriangleIntersection* segmentTriangleIntersection = new SegmentTriangleIntersection(segment, triangle);

		segmentTriangleIntersection->IntersectionCalculate();
		
		SegmentTriangleIntersection::is_ray_cross_triangle(segment, triangle);

		std::cout << "Intersection: " << *segmentTriangleIntersection->getIntersectionPoint() << std::endl;

		delete segmentTriangleIntersection;


		//
	}
	//delete triangle;
	//delete segment;

	return 0;
}

