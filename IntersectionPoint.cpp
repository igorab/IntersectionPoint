// Написать код на С++для определения точки пересечения отрезка и
// треугольника в 3D пространстве.Отрезок задан координатами концов,
// треугольник задан координатами всех трех углов.
//
//https://www.youtube.com/watch?v=RIFXebcuryc&list=PLtNPgSbW9TX7acrQa2LeBAMGxO5WRAVsz&index=58

#include <iostream>

typedef int num;

struct Point
{
	num X = 0;
	num Y = 0;
	num Z = 0;

	explicit Point(num _X=0, num _Y=0, num _Z=0) : X(_X), Y(_Y), Z(Z) {}

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


	num product_1 = v_cross_product((Vector)vertex1, Vector(vertex1, vertex2) ); // use static cast ?


	return ret;
}


int main()
{
	const Point pnt(3, 4, 5);

	std::cout << pnt << std::endl;

    //std::cout << "Point: X= " << pnt.X << " Y= " << pnt.Y << " Z= " << pnt.Z << std::endl;
	return 0;
}

