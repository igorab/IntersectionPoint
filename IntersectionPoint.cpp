// Написать код на С++для определения точки пересечения отрезка и
// треугольника в 3D пространстве.Отрезок задан координатами концов,
// треугольник задан координатами всех трех углов.

#include <iostream>

using namespace std;
typedef float num;

const float EPSILON = 0.0000001;

struct IShape
{
	virtual num len() = 0;
};

// точка в пространстве
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

// вектор
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

	// модуль вектора
	num len() const {
		num sumV = X * X + Y * Y + Z * Z;
		return sqrt(sumV);
	}

	// произведение вектора на число
	Vector mult3(const num C)
	{
		Vector vc(C * X, C * Y, C * Z);
		return vc;
	}

	//сложение векторов
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

	//векторное произведение
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
	// скалярное произведение
	friend num operator*(Vector const _a, Vector const _b) { return _a.X * _b.X + _a.Y * _b.Y + _a.Z * _b.Z; } 
	// умножение на число
	friend Vector operator*(num const _c, Vector const _V) { 
		Vector V(_V);
		Vector c_V = V.mult3(_c);
		return c_V;
	} 
	//Векторное произведение
	friend Vector operator^(Vector const _a, Vector const _b) {
		Vector A(_a);
		Vector normV = A.cross3(_b);		
		return normV;
	}
};

// Отрезок
class Segment :IShape
{
private:
	Point L_A;
	Point L_B;
public:

	Point getA() const { return L_A; };
	Point getB() const { return L_B; };

	// направляющий вектор
	Vector dir;

	explicit Segment(Point _A, Point _B) : L_A(_A), L_B(_B) {
		dir = Vector(L_A, L_B);
	}

	num len() { return dir.len(); }
};


//    треугольник
class Triangle : IShape
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

	// описать уравнение плоскости через 3 точки
	// 
	//  вектор нормали
	Vector norm()
	{
		Vector V_N = V_AB ^ V_AC;
		return V_N;
	}

	friend ostream& operator<<(ostream& out, Triangle& _T) {
		out << _T.P_A << " " << _T.P_B << " " << _T.P_C << endl;
		return out;
	}

	num len() { return 0; }
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

	float a, f, u, v;
	edge1 = vertex1 - vertex0;
	edge2 = vertex2 - vertex0;

	h = rayVector.cross3(edge2);

	a = edge1.dot3(h);

	if (a > -EPSILON && a < EPSILON)
		return false;    // This ray is parallel to this triangle.

	f = 1.0 / a;
	s = rayOrigin - vertex0;

	u = f * (s * h);

	if (u < 0.0 || u > 1.0)
		return false;
	
	q = s ^ edge1;

	v = f * rayVector.dot3(q);

	if (v < 0.0 || u + v > 1.0f)
		return false;

	// At this stage we can compute t to find out where the intersection point is on the line.
	num t = f * edge2.dot3(q);

	if (t > EPSILON) // ray intersection
	{		
		outIntersectionPoint = rayOrigin + (rayVector.mult3(t));
		return true;
	}
	else
	{
		// This means that there is a line intersection but not a ray intersection.
		return false;
	}
}


// Линия пересекает треугольник ?
bool is_line_cross_triangle(Segment* _segment, Triangle* _triangle)
{		
	Vector intersectionPoint;

	RayIntersectsTriangle(_segment, _triangle, intersectionPoint);
	
	std::cout << "Intersection: " << intersectionPoint.getPoint() << std::endl;

	return true;
}

// triangle
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 0, 0, 1
// segment
#define fromA 0, 0, 0
#define toB 5, 5, 5

int main()
{		
	const Point pA(vA), pB(vB), pC(vC);	
	//Triangle triangle(pA, pB, pC);
	Triangle* triangle = new Triangle(pA, pB, pC);
	Vector normT = triangle->norm();

	const Point pntFrom(fromA), pntTo(toB);
	
	Segment segment(pntFrom, pntTo);
	Vector segDir = segment.dir;

	is_line_cross_triangle(&segment, triangle);
	
	delete triangle;
	//std::cout << triangle << std::endl;

	return 0;
}

