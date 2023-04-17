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
// начало отрезка находится в площади треугольника
#define vA 2, 0, 1
#define vB 0, 1, 0
#define vC 0, 0, 3 

#define fromA 0.6, 0.6, 0.6
#define toB 1., 1., 1.


//Test 4, IntersectionType::V 
// начало отрезка находится в вершине треугольника
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
// нет пересечения
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

//Test 12,
// лежит в плоскости
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 3, 0, 0 

#define fromA 0, 1, 0
#define toB 5, 4, 0

//Test 13,
// лежит в плоскости, параллельно одной из граней
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 3, 0, 0

#define fromA 0, -1, 0
#define toB 0, 4, 0

//Test 14
// лежит в плоскости, не пересекает ни одну из граней
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 3, 0, 0

#define fromA 4, 5, 0
#define toB 6, 7, 0

//Test 15
// лежит в плоскости
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 0, 0, 1

#define fromA 1, 0, 0
#define toB 0, 1, 0

//Test 16
// лежит в плоскости, пересекает параллельно
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 0, 0, 1

#define fromA 0.5, 0, 0.5
#define toB 0, 0.5, 0.5



#pragma endregion