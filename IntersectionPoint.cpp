// �������� ��� �� �++��� ����������� ����� ����������� ������� �
// ������������ � 3D ������������.������� ����� ������������ ������,
// ����������� ����� ������������ ���� ���� �����.
//
//https://www.youtube.com/watch?v=RIFXebcuryc&list=PLtNPgSbW9TX7acrQa2LeBAMGxO5WRAVsz&index=58


//https://question-it.com/questions/3961314/3d-peresechenie-mezhdu-segmentom-i-treugolnikom

//https://www.geeksforgeeks.org/equation-of-a-line-in-3d/

//https://algocode.ru/page/c-23-geometry/

//https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d


// ��������� � 
// https://ru.wikipedia.org/wiki/OpenGL_Mathematics

#include <iostream>

using namespace std;
typedef int num;

struct IShape
{
	virtual num len() = 0;
};


struct Point : IShape
{
	num X = 0;
	num Y = 0;
	num Z = 0;

	explicit Point(num _X = 0, num _Y = 0, num _Z = 0) : X(_X), Y(_Y), Z(_Z) {}

	friend ostream& operator<< (ostream& out, const Point& point);

	friend istream& operator>> (istream& in, Point& point);

	num len() { return 0; }
};

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

	// ������ �������
	num len() const {
		num sumV = X * X + Y * Y + Z * Z;
		return sqrt(sumV);
	}

	// ������������ ������� �� �����
	Vector mult3(const num C)
	{
		Vector vc(X,Y,Z);

		vc.X = X * C;
		vc.Y = Y * C;
		vc.Z = Z * C;

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
		Vector vNorm(0, 0, 0);

		vNorm.X = Y * _V.Z - Z * _V.Y;
		vNorm.Y = Z * _V.X - X * _V.Z;
		vNorm.Z = X * _V.Y - Y * _V.X;

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
	friend num operator*(Vector const _a, Vector const _b) { return _a.X * _b.X + _a.Y * _b.Y + _a.Z * _b.Z; } // ��������� ������������
	friend num operator*(num const _c, Vector const _V) { return _c * _V.X + _c * _V.Y + _c * _V.Z; } // ��������� �� �����

	//��������� ������������
	friend Vector operator^(Vector _a, Vector _b) {
		Vector normV(0, 0, 0);
		normV.X = _a.Y * _b.Z - _a.Z * _b.Y;
		normV.Y = _a.X * _b.Z - _a.Z * _b.X;
		normV.Z = _a.X * _b.Y - _a.Y * _b.X;
		return normV;
	}

};

class VectorMath
{
public:

	double* cross3(double* _X, double* _Y)
	{
		double* v = (double*)malloc(3 * sizeof(double));

		v[0] = _X[1] * _Y[2] - _X[2] * _Y[1];
		v[1] = _X[2] * _Y[0] - _X[0] * _Y[2];
		v[2] = _X[0] * _Y[1] - _X[1] * _Y[0];

		return v;
	}

	//��������� ������������
	double dotProduct3(double* _X, double* _Y)
	{
		double val;

		val = _X[0] * _Y[0] + _X[1] * _Y[1] + _X[2] * _Y[2];

		return val;
	}

	//������ �������
	double norma3(double* val)
	{
		double sumV = 0;

		for (int i = 0; i < 3; i++)
		{
			sumV += pow(val[i], 2);
		}

		return sqrt(sumV);
	}

	// ������������ ������� �� �����
	double* mult3(double* _X, double C)
	{
		double* V = (double*)malloc(3 * sizeof(double));

		for (int i = 0; i < 3; i++)
		{
			V[i] = _X[i] * C;
		}

		return V;
	}

	//�������� ��������
	double* sum3(double* _X, double* _Y, int sign)
	{
		double* v = (double*)malloc(3 * sizeof(double));

		for (int i = 0; i < 3; i++)
		{
			v[i] = _X[i] + sign * _Y[i];
		}

		return v;
	}

	double determinant(double* a, double* _X, double* _Y)
	{
		double det;

		det = a[0] * (_X[1] * _Y[2] - _X[2] * _Y[1]) +
			a[1] * (_X[2] * _Y[0] - _X[0] * _Y[2]) +
			a[2] * (_X[0] * _Y[1] - _X[1] * _Y[0]);

		return det;
	}
	
	//��������������� ������������ 	
	static num v_cross_product(Vector v1, Vector v2)
	{
		num ret;

		ret = v1.cross3(v2).X;

		return ret;
	}

};

/**
 https://radioprog.ru/post/1240 
 */
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

/**
 * \brief �������
 * Finding Parametric Equations Passing Through Two Points
 *  https://www.youtube.com/watch?v=NXazSzbK6n8
 * 
 */
class Segment :IShape
{
private:
	Point L_A;
	Point L_B;
public:

	Point getA() const { return L_A; };
	Point getB() const { return L_B; };

	// ������������ ������
	Vector dir;

	explicit Segment(Point _A, Point _B) : L_A(_A), L_B(_B) {
		dir = Vector(L_A, L_B);
	}

	num len() { return dir.len(); }
};


struct r {
	double x, y;
	r() {}
	r(int _x, int _y) { x = _x, y = _y; }
};

struct VGeom
{
	double len(r a) { return sqrt(a.x * a.x + a.y * a.y); }
	
	friend r operator+ (r const a, r const b) { return r(a.x + b.x, a.y + b.y); }

	friend r operator- (r a, r b) { return r(a.x - b.x, a.y - b.y); }

	// ��������� ������������
	friend int operator*(r a, r b) { return a.x * b.x + a.y * b.y; }
	//��������� ������������
	friend int operator^(r a, r b) { return a.x * b.y - b.x * a.y; }

	friend istream& operator>>(istream& in, r& p) {
		in >> p.x >> p.y;
		return in;
	}

	friend ostream& operator<<(ostream& out, r& p) {
		out << p.x << " " << p.y << endl;
		return out;
	}	
};


/**
 *  \brief �����������
 */
class Triangle: IShape
{
private:
	Point P_A;
	Point P_B;
	Point P_C;

public:

	Vector V_AB;
	Vector V_AC;
	
	explicit Triangle(const Point _A, Point _B, Point _C) : P_A(_A), P_B(_B), P_C(_C) 
	{
		V_AB = Vector(P_A, P_B);
		V_AC = Vector(P_A, P_C);
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
		out << _T.P_A << " " << _T.P_B <<  " " << _T.P_C <<  endl;
		return out;
	}

	num len() { return 0; }
};

/**
 * \brief ����������� �� ����� ������������?
 * \param _point 
 * \param _triangle 
 * \return 
 */
bool hit_into_triangle(Point _point, Triangle _triangle)
{
	bool ret = false;

	Point vertex1;
	Point vertex2;
	Point vertex3;

	Vector vector1(vertex1, _point) ;
	Vector vector2(vertex2, _point) ;
	Vector vector3(vertex3, _point) ;

	
	num product_1 = VectorMath::v_cross_product(static_cast<Vector>(vertex1), Vector(vertex1, vertex2) ); // use static cast ?


	return ret;
}

/// <summary>
/// ����� ���������
/// </summary>
num SignedVolume(Vector A, Vector B, Vector C, Vector D)
{	
	num signedVol = 1/6 * (((B - A) ^ (C - A)) * (D - A));

	return signedVol;
}


bool RayIntersectsTriangle(Vector rayOrigin,
						   Vector rayVector,
						   Triangle* inTriangle,
						   Vector& outIntersectionPoint)
{
	const float EPSILON = 0.0000001;

	Vector vertex0 = inTriangle->getA();
	Vector vertex1 = inTriangle->getB();
	Vector vertex2 = inTriangle->getC();

	Vector edge1, edge2, h, s, q;

	float a, f, u, v;
	edge1 = vertex1 - vertex0;
	edge2 = vertex2 - vertex0;

	h = rayVector.cross3(edge2);

	a = edge1.dot3(h);

	if (a > -EPSILON && a < EPSILON)
		return false;    // This ray is parallel to this triangle.

	f = 1.0 / a;
	s = rayOrigin - vertex0;

	u = f * s.dot3(h);

	if (u < 0.0 || u > 1.0)
		return false;

	q = s.cross3(edge1);

	v = f * rayVector.dot3(q);

	if (v < 0.0 || u + v > 1.0)
		return false;

	// At this stage we can compute t to find out where the intersection point is on the line.
	float t = f * edge2.dot3(q);

	if (t > EPSILON) // ray intersection
	{
		outIntersectionPoint = rayOrigin + rayVector.mult3(t);
		return true;
	}
	else // This means that there is a line intersection but not a ray intersection.
		return false;
}




// ����� ���������� ����������� ?
bool is_line_cross_triangle(Triangle _triangle, Segment _segment)
{
	//����� p1, p2, p3 ���������� �����������
	Vector p1 = _triangle.getA();
	Vector p2 = _triangle.getB();
	Vector p3 = _triangle.getC();

	Vector q1 = _segment.getA();
	Vector q2 = _segment.getB();



	return true;
}

// �������� ����� �����������
void CrossPoint(Triangle _triangle, Segment _segment)
{			
	Point p1 = _triangle.getA();
	Point p2 = _triangle.getB();
	Point p3 = _triangle.getC();
	
	num t = 1;	
	Point p_t;

	Point q1 = _segment.getA();
	Point q2 = _segment.getB();
	
	Vector d_q = _segment.dir;

	//��������� ������ � ��������������� �����: p (t) = q1 + t * (q2-q1)
	p_t.X = q1.X + t * d_q.X;
	p_t.Y = q1.Y + t * d_q.Y;
	p_t.Z = q1.Z + t * d_q.Z;
	
	//��������� ��������� : 
	// dot(p, N) � dot(p, p1) = 0, ��� N = �����(p2 - p1, p3 - p1)
	Triangle triangle(p1, p2, p3);

	Vector N = triangle.norm();

	//������� p(t) � ��������� ��������� : �����(q1 + t * (q2 - q1), N - p1) = 0

	//�������� t = -dot(q1, N - p1) / dot(q1, q2 - q1)
	Vector v_q1(q1);
	Vector v_p1(p1);
	Vector v_q2_1(q2, q1);
	
	num denom = v_q1 * v_q2_1;

	t = - (v_q1 * (N - v_p1)) / denom;

	// ����� ����������� q1 + t * (q2 - q1)
}


int main()
{		
	const Point pA(1, 0, 0);
	const Point pB(0, 1, 0);
	const Point pC(0, 0, 1);

	Triangle triangle(pA, pB, pC);

	Vector normT = triangle.norm();


	const Point pntFrom(3, 4, 5);
	const Point pntTo(4, 4, 7);
	Segment segment(pntFrom, pntTo);

	CrossPoint(triangle, segment);

	Vector segDir = segment.dir;

	std::cout << triangle << std::endl;


	return 0;
}

