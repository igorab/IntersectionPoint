#include <iostream>

//https://www.youtube.com/watch?v=RIFXebcuryc&list=PLtNPgSbW9TX7acrQa2LeBAMGxO5WRAVsz&index=58
//https://question-it.com/questions/3961314/3d-peresechenie-mezhdu-segmentom-i-treugolnikom
//https://www.geeksforgeeks.org/equation-of-a-line-in-3d/
//https://algocode.ru/page/c-23-geometry/
//https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
// проверить с 
// https://ru.wikipedia.org/wiki/OpenGL_Mathematics
/**
 https://radioprog.ru/post/1240
 */
 /**
  * \brief отрезок
  * Finding Parametric Equations Passing Through Two Points
  *  https://www.youtube.com/watch?v=NXazSzbK6n8
  *
  */


  /// <summary>
  /// Объем тетраэдра
  /// </summary>
  /*
num SignedVolume(Vector A, Vector B, Vector C, Vector D)
{
	num signedVol = 1 / 6 * (((B - A) ^ (C - A)) * (D - A));

	return signedVol;
}
*/

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
	/*
	friend istream& operator>>(istream& in, r& p) {
		in >> p.x >> p.y;
		return in;
	}

	friend ostream& operator<<(ostream& out, r& p) {
		out << p.x << " " << p.y << endl;
		return out;
	}
	*/
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

	/*
	//псевдоскалярное произведение 	
	static num v_cross_product(Vector v1, Vector v2)
	{
		num ret;

		ret = v1.cross3(v2).X;

		return ret;
	}
	*/

	/*
	// Получить точку пересечения
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

		//уравнение прямой в параметрической форме: p (t) = q1 + t * (q2-q1)
		p_t.X = q1.X + t * d_q.X;
		p_t.Y = q1.Y + t * d_q.Y;
		p_t.Z = q1.Z + t * d_q.Z;

		//уравнение плоскости : 
		// dot(p, N) — dot(p, p1) = 0, где N = крест(p2 - p1, p3 - p1)
		Triangle triangle(p1, p2, p3);

		Vector N = triangle.norm();

		//Введите p(t) в уравнение плоскости : точка(q1 + t * (q2 - q1), N - p1) = 0

		//Выведите t = -dot(q1, N - p1) / dot(q1, q2 - q1)
		Vector v_q1(q1);
		Vector v_p1(p1);
		Vector v_q2_1(q2, q1);

		num denom = v_q1 * v_q2_1;

		t = -(v_q1 * (N - v_p1)) / denom;

		// Точка пересечения q1 + t * (q2 - q1)
	}
	*/

};