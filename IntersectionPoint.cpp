// Ќаписать код на —++дл€ определени€ точки пересечени€ отрезка и
// треугольника в 3D пространстве.ќтрезок задан координатами концов,
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

static class VectorMath
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

	//скал€рное произведение
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
};


struct Point
{
	num X = 0;
	num Y = 0;
	num Z = 0;

	explicit Point(num _X=0, num _Y=0, num _Z=0) : X(_X), Y(_Y), Z(_Z) {}

	friend std::ostream& operator<< (std::ostream& out, const Point& point);
};

/**
 https://radioprog.ru/post/1240 
 */
std::ostream& operator<< (std::ostream& out, const Point& point)
{
	out << "Point(" << point.X << ", " << point.Y << ", " << point.Z << ')';
	return out;
}

/**
 * \brief отрезок
 */
struct Segment
{
	Point L_A;
	Point L_B;

	explicit Segment(Point _A, Point _B) : L_A(_A), L_B(_B) {}

	/**
	 * \brief Finding Parametric Equations Passing Through Two Points
	 * \return / https://www.youtube.com/watch?v=NXazSzbK6n8
	 */
	num func_AB()
	{
		//Vector direction_V(L_B, L_A);

		//direction_V;
			//num 
	}
};

struct Vector
{
	num X, Y, Z;
	
	Vector(num _X, num _Y, num _Z) : X(_X), Y(_Y), Z(_Z) {}
	
	Vector(Point _point)
	{
		X = _point.X;
		Y = _point.Y;
		Z = _point.Z;
	}

	Vector(Point _from, Point _to)
	{
		X = _to.X - _from.X;
		Y = _to.Y - _from.Y;
		Z = _to.Z - _from.Z;
	}
	

	// модуль вектора
	num len()
	{
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
	
	friend Vector operator+ (Vector const _a, Vector const _b) { return Vector(_a.X + _b.X , _a.Y + _b.Y, _a.Z + _b.Z ); }
	friend Vector operator- (Vector const _a, Vector const _b) { return Vector(_a.X - _b.X, _a.Y - _b.Y, _a.Z - _b.Z); }	
	friend num operator*(Vector const _a, Vector const _b) { return _a.X * _b.X + _a.Y * _b.Y + _a.Z * _b.Z; } // скал€рное произведение

	//¬екторное произведение
	friend Vector operator^(Vector _a, Vector _b) { 
		Vector normV(0,0,0);
		normV.X = _a.Y * _b.Z - _b.Z * _a.Y; 
		normV.Y = _a.Z * _b.X - _b.X * _a.Z; 
		normV.Z = _a.X * _b.Y - _b.Y * _a.X; 
		return normV;
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

	// скал€рное произведение
	friend int operator*(r a, r b) { return a.x * b.x + a.y * b.y; }
	//¬екторное произведение
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
struct Triangle
{
	Point V_A;
	Point V_B;
	Point V_C;

	explicit Triangle(const Point p_A, Point p_B, Point p_C) : V_A(p_A), V_B(p_B), V_C(p_C) {}
	
	// описать уравнение плоскости через 3 точки
};


/**
 * \brief псевдоскал€рное произведение
 * \param v1 вектор 1
 * \param v2 вектор 2
 * \return число
 */
num v_cross_product(Vector v1, Vector v2)
{
	num ret;

	ret = v1.cross3(v2).X;

	return ret;
}

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


	num product_1 = v_cross_product(static_cast<Vector>(vertex1), Vector(vertex1, vertex2) ); // use static cast ?


	return ret;
}

// Ћини€ пересекает треугольник ?
void is_line_cross_triangle()
{
	//пусть p1, p2, p3 обозначают ваш треугольник


}

// ѕолучить точку пересечени€
void CrossPoint()
{
	Point p1(0, 0, 0);
	Point p2(1, 0, 0);
	Point p3(1, 1, 1);


	//запишите уравнение пр€мой в параметрической форме: p (t) = q1 + t * (q2-q1)
	num t = 0;	
	Point p_t;
	Point q1;
	Point q2;
	
	Segment segment(q1, q2);

	p_t.X = q1.X + t * (q2.X - q1.X);
	p_t.Y = q1.Y + t * (q2.Y - q1.Y);
	p_t.Z = q1.Z + t * (q2.Z - q1.Z);
	
	//«апишите уравнение плоскости : точка(p, N) Ч точка(p, p1) = 0, где N = крест(p2 - p1, p3 - p1)
	Triangle triangle(p1, p2, p3);

	//¬ведите p(t) в уравнение плоскости : точка(q1 + t * (q2 - q1), N - p1) = 0


	//¬ыведите t = -dot(q1, N - p1) / dot(q1, q2 - q1)
	Vector v(q1);
	Vector v(p1);
	Vector v_q1(q1);
	Vector v_q21(q2, q1);

	t = -v.dot3(p1) / (v_q1 * v_q21);

	// “очка пересечени€ q1 + t * (q2 - q1)
}


int main()
{
	const Point pnt(3, 4, 5);

	std::cout << pnt << std::endl;

    //std::cout << "Point: X= " << pnt.X << " Y= " << pnt.Y << " Z= " << pnt.Z << std::endl;
	return 0;
}

