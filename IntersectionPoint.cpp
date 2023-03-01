// Написать код на С++для определения точки пересечения отрезка и
// треугольника в 3D пространстве.Отрезок задан координатами концов,
// треугольник задан координатами всех трех углов.
//
//https://www.youtube.com/watch?v=RIFXebcuryc&list=PLtNPgSbW9TX7acrQa2LeBAMGxO5WRAVsz&index=58

//https://web-answers.ru/c/peresechenie-linii-i-treugolnika-v-3d.html

//https://question-it.com/questions/3961314/3d-peresechenie-mezhdu-segmentom-i-treugolnikom

//https://www.geeksforgeeks.org/equation-of-a-line-in-3d/

//https://algocode.ru/page/c-23-geometry/


#include <iostream>

using namespace std;
typedef int num;

struct Point
{
	num X = 0;
	num Y = 0;
	num Z = 0;

	explicit Point(num _X = 0, num _Y = 0, num _Z = 0) : X(_X), Y(_Y), Z(_Z) {}

	friend ostream& operator<< (ostream& out, const Point& point);

	friend istream& operator>> (istream& in, Point& point);
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

	// модуль вектора
	num len() const {
		num sumV = X * X + Y * Y + Z * Z;
		return sqrt(sumV);
	}

	// произведение вектора на число
	void mult3(const num C)
	{
		X = X * C;
		Y = Y * C;
		Z = Z * C;
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
	friend num operator*(Vector const _a, Vector const _b) { return _a.X * _b.X + _a.Y * _b.Y + _a.Z * _b.Z; } // скалярное произведение
	friend num operator*(num const _c, Vector const _V) { return _c * _V.X + _c * _V.Y + _c * _V.Z; } // умножение на число

	//Векторное произведение
	friend Vector operator^(Vector _a, Vector _b) {
		Vector normV(0, 0, 0);
		normV.X = _a.Y * _b.Z - _b.Z * _a.Y;
		normV.Y = _a.Z * _b.X - _b.X * _a.Z;
		normV.Z = _a.X * _b.Y - _b.Y * _a.X;
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

	//скалярное произведение
	double dotProduct3(double* _X, double* _Y)
	{
		double val;

		val = _X[0] * _Y[0] + _X[1] * _Y[1] + _X[2] * _Y[2];

		return val;
	}

	//модуль вектора
	double norma3(double* val)
	{
		double sumV = 0;

		for (int i = 0; i < 3; i++)
		{
			sumV += pow(val[i], 2);
		}

		return sqrt(sumV);
	}

	// произведение вектора на число
	double* mult3(double* _X, double C)
	{
		double* V = (double*)malloc(3 * sizeof(double));

		for (int i = 0; i < 3; i++)
		{
			V[i] = _X[i] * C;
		}

		return V;
	}

	//сложение векторов
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
	
	//псевдоскалярное произведение 	
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
 * \brief отрезок
 * Finding Parametric Equations Passing Through Two Points
 *  https://www.youtube.com/watch?v=NXazSzbK6n8
 * 
 */
class Segment
{
private:
	Point L_A;
	Point L_B;
public:

	// направляющий вектор
	Vector dir;

	explicit Segment(Point _A, Point _B) : L_A(_A), L_B(_B) {
		dir = Vector(L_A, L_B);
	}			
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

	// скалярное произведение
	friend int operator*(r a, r b) { return a.x * b.x + a.y * b.y; }
	//Векторное произведение
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
 *  \brief треугольник
 */
class Triangle
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
	
	// описать уравнение плоскости через 3 точки
	// 
	//  вектор нормали
	Vector norm()
	{
		Vector V_N = V_AB ^ V_AC;
		return V_N;
	}

	friend ostream& operator<<(ostream& out, Triangle& _T) {
		out << _T.P_A << " " << _T.P_B <<  " " << _T.P_C <<  endl;
		return out;
	}
};

/**
 * \brief принадлежит ли точка треугольнику?
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

// Линия пересекает треугольник ?
void is_line_cross_triangle()
{
	//пусть p1, p2, p3 обозначают ваш треугольник


}

// Получить точку пересечения
void CrossPoint()
{
	Point p1(0, 0, 0);
	Point p2(1, 0, 0);
	Point p3(1, 1, 1);


	//запишите уравнение прямой в параметрической форме: p (t) = q1 + t * (q2-q1)
	num t = 0;	
	Point p_t;
	Point q1;
	Point q2;
	
	Segment segment(q1, q2);

	p_t.X = q1.X + t * (q2.X - q1.X);
	p_t.Y = q1.Y + t * (q2.Y - q1.Y);
	p_t.Z = q1.Z + t * (q2.Z - q1.Z);
	
	//Запишите уравнение плоскости : точка(p, N) — точка(p, p1) = 0, где N = крест(p2 - p1, p3 - p1)
	Triangle triangle(p1, p2, p3);

	//Введите p(t) в уравнение плоскости : точка(q1 + t * (q2 - q1), N - p1) = 0


	//Выведите t = -dot(q1, N - p1) / dot(q1, q2 - q1)
	Vector v(q1);
	
	Vector v_q1(q1);
	Vector v_q21(q2, q1);

	t = -v.dot3(p1) / (v_q1 * v_q21);

	// Точка пересечения q1 + t * (q2 - q1)
}


int main()
{
	const Point pntFrom(3, 4, 5);
	const Point pntTo(4, 4, 7);

	std::cout << pntFrom << std::endl;

	Point pA(0, 0, 0);
	Point pB(1, 0, 0);
	Point pC(1, 1, 1);

	Triangle triangle(pA, pB, pC);

	Vector normT = triangle.norm();

	Segment segment(pntFrom, pntTo);
	Vector segDir = segment.dir;

	std::cout << triangle << std::endl;


	return 0;
}

