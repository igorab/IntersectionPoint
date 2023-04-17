//Test 1, IntersectionType::F 
// ����� ������� ����� � �������� ������������
#define vA 2, 0, 1
#define vB 0, 1, 0
#define vC 0, 0, 3 

#define fromA 0, 0, 0
#define toB 0.6, 0.6, 0.6 

//Test 2, IntersectionType::E
// ����� ������� ����� � �����
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 1, 1, 1 

#define fromA 0, 0, 0
#define toB 0.5, 0.5, 0

//Test 3, IntersectionType::V
// ����� ������� ����� � �������
#define vA 3, 0, 0
#define vB 0, 1, 0
#define vC 1, 1, 1 

#define fromA 0, 0, 0
#define toB 1., 1., 1.

//Test 4, IntersectionType::F 
// ������ ������� ��������� � ������� ������������
#define vA 2, 0, 1
#define vB 0, 1, 0
#define vC 0, 0, 3 

#define fromA 0.6, 0.6, 0.6
#define toB 1., 1., 1.


//Test 4, IntersectionType::V 
// ������ ������� ��������� � ������� ������������
#define vA 2, 0, 1
#define vB 0, 1, 0
#define vC 0, 0, 3 

#define fromA 2, 0., 1
#define toB 5, 5, 6

//Test 5, IntersectionType::V 
// ������� ���������� ������� ������������
#define vA 2, 1, 1
#define vB 1, 4, -1
#define vC 0.5, 0.5, 7 

#define fromA 0, 0., 0
#define toB 10, 10, 10

//Test 6, IntersectionType::f 
// ������� ���������� ��������� ������������
#define vA 2, 1, 1
#define vB 1, 4, -1
#define vC 0.5, 0.5, 7 

#define fromA 10, 10., 10
#define toB 0, 0, 0

//Test 7, IntersectionType::v 
// ������� ���������� ������� ������������
#define vA 2, 1, 1
#define vB 1, 4, -1
#define vC 0.5, 0.5, 7 

#define fromA 0, -1, -1
#define toB 4, 3, 3

//Test 8, IntersectionType::e 
// ������� ���������� ����� ������������
#define vA 0, 0, 0
#define vB 3, 5, 1
#define vC 0, 0, 7 

#define fromA -2, -2, 0
#define toB 2, 2, 4

//Test 9, 
// ��� �����������
#define vA 0, 0, 0
#define vB 3, 5, 1
#define vC 0, 0, 7 

#define fromA 4, 6, 3
#define toB 9, 12, 24

//Test 10, 
// ����������� ���������
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 0, 0, 7 

#define fromA 1, 0, 0
#define toB 1, 4, 5


//Test 11, 
// ����� � ���������
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 0, 0, 3 

#define fromA 0, 0, 0
#define toB 0, 4, 5

//Test 12,
// ����� � ���������
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 3, 0, 0 

#define fromA 0, 1, 0
#define toB 5, 4, 0

//Test 13,
// ����� � ���������, ����������� ����� �� ������
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 3, 0, 0

#define fromA 0, -1, 0
#define toB 0, 4, 0

//Test 14
// ����� � ���������, �� ���������� �� ���� �� ������
#define vA 0, 0, 0
#define vB 0, 3, 0
#define vC 3, 0, 0

#define fromA 4, 5, 0
#define toB 6, 7, 0

//Test 15
// ����� � ���������
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 0, 0, 1

#define fromA 1, 0, 0
#define toB 0, 1, 0

//Test 16
// ����� � ���������, ���������� �����������
#define vA 1, 0, 0
#define vB 0, 1, 0
#define vC 0, 0, 1

#define fromA 0.5, 0, 0.5
#define toB 0, 0.5, 0.5



#pragma endregion