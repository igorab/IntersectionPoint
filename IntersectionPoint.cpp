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

//https://www.youtube.com/watch?v=NXazSzbK6n8
	/**
	 * \brief Finding Parametric Equations Passing Through Two Points
	 * \return /
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

	Vector (num _X, num _Y, num _Z): X(_X), Y(_Y), Z(_Z) {}

	// модуль вектора
	double norma3()
	{
		num sumV = 0;
		
		sumV = X*X + Y*Y + Z*Z;
		
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

	//векторное произведение
	Vector cross3(const Vector _V)
	{		
		Vector vNorm(0, 0, 0);

		vNorm.X = Y * _V.Z - Z * _V.Y;
		vNorm.Y = Z * _V.X - X * _V.Z;
		vNorm.Z = X * _V.Y - Y * _V.X;

		return vNorm;
	}

	num determinant(Vector a, Vector _V)
	{
		num det;

		det = a.X * (Y * _V.Z - Z * _V.Y) +
			a.Y * (Z * _V.X - X * _V.Z) +
			a.Z * (X * _V.Y - Y * _V.X);

		return det;
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
	/*
	r operator+ (r a, r b) { return r(a.x + b.x, a.y + b.y); }
	r operator- (r a, r b) { return r(a.x - b.x, a.y - b.y); }
	// скалярное произведение
	int operator*(r a, r b) { return a.x * b.x + a.y * b.y; }
	//Векторное произведение
	int operator^(r a, r b) { return a.x * b.y - b.x * a.y; }
	istream& operator>>(istream& in, r& p) {
		in >> p.x >> p.y;
		return in;
	}

	ostream& operator<<(ostream& out, r& p) {
		out << p.x << " " << p.y << endl;
		return out;
	}
	*/
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
	
};


/**
 * \brief псевдоскалярное произведение
 * \param v1 вектор 1
 * \param v2 вектор 2
 * \return число
 */
num v_cross_product(Vector v1, Vector v2)
{
	num ret;

	ret = v1.X * v2.Y - v2.X * v1.Y;

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


int main()
{
	const Point pnt(3, 4, 5);

	std::cout << pnt << std::endl;

    //std::cout << "Point: X= " << pnt.X << " Y= " << pnt.Y << " Z= " << pnt.Z << std::endl;
	return 0;
}

