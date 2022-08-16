//import libraries
#include <cmath>
#include <exception>
//#include <format>
#include <stdexcept>
// import internal header files
#include "Point3D.h"

//using namespace std;
using std::cout;
using std::endl;
using std::cerr;

auto sqr = [](double x){return x * x;};


//constructors
Point3D::Point3D()
	:x(0), y(0), z(0)
{}
Point3D::Point3D(double d)
	: x(d), y(d), z(d)
{}

Point3D::Point3D(double x, double y, double z)
	: x(x), y(y), z(z)
{}

/*methods*/


/*operator overloading*/
Point3D& Point3D::operator +=(const Point3D& rhs)
{
	this->x += rhs.x;
	this->y += rhs.y;
	this->z += rhs.z;

	return (*this);
}

int Point3D::get_clusterID()const
{
	return m_clusterID;
}

void Point3D::set_clusterID(int clusterID)
{
	m_clusterID = clusterID;
}

Point3D& Point3D::operator -=(const Point3D& rhs)
{
	(*this) = (*this) - rhs;

	return (*this);
}

Point3D& Point3D::operator *=(const Point3D& rhs)
{
	(*this) = (*this) * rhs;

	return (*this);
}

Point3D& Point3D::operator /=(const Point3D& rhs)
{
	(*this) = (*this) / rhs;

	return (*this);
}

/*global helper functions*/
// operator << for out stream
std::ostream& operator <<(std::ostream& os, const Point3D& p)
{
	os << p.x << " " << p.y << " " << p.z << " ";
	return os;
}

// since c++20 we can use format header as standard library
// and use it for custom type
//template<>
//struct formatter<Point3D>
//{
//public:
	//constexpr auto parse(auto& contex)
	//{
		//auto iter{ contex.begin() };
		//auto end{ contex.end() };

		//if (iter == end || *iter == '}')
		//{
			//outputType = OutputType::ThreeD;
			//return iter;
		//}

		//switch (*iter)
		//{
			//case 'a': //{:a}
				//outputType = OutputType::OneD;
				//break;
			//case 'b': //{:b}
				//outputType = OutputType::TwoD;
				//break;
			//case 'c':
				//outputType = OutputType::ThreeD;
				//break;
			//default:
				//throw format_error{ "Invalid Point3D format specifier." };
		//}

		//++iter;
		//if (iter != end && *iter != '}')
		//{
			//throw format_error{ "Invalide Point3D format specifier." };
		//}

		//return iter;
	//}

	//auto format(const Point3D& p, auto& contex)
	//{
		//switch (outputType)
		//{
			//case OutputType::OneD : return format_to(contex.out(), "{} ", p.x);
			//case OutputType::TwoD: return format_to(contex.out(), "{} {} ", p.x, p.y);
			//default: return format_to(contex.out(), "{} {} {} ", p.x, p.y, p.z);
		//}
	//}

//private:
	//enum class OutputType{
		//OneD,
		//TwoD,
		//ThreeD
	//};

	//OutputType outputType{ OutputType::ThreeD };

//}; it works for std::format in visual studio 2022, it needs some correction to work with fmt::format

Point3D operator +(const Point3D& lhs, const Point3D& rhs)
{
	Point3D result;
	result.x = lhs.x + rhs.x;
	result.y = lhs.y + rhs.y;
	result.z = lhs.z + rhs.z;

	return result;

}

Point3D operator -(const Point3D& lhs, const Point3D& rhs)
{
	Point3D result;
	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;
	result.z = lhs.z - rhs.z;

	return result;
}

Point3D operator *(const Point3D& lhs, const Point3D& rhs)
{
	Point3D result;
	result.x = lhs.x * rhs.x;
	result.y = lhs.y * rhs.y;
	result.z = lhs.z * rhs.z;

	return result;
}

Point3D operator /(const Point3D& lhs, const Point3D& rhs)
{
	Point3D result;

	if (rhs.x == 0 || rhs.y == 0 || rhs.z == 0)
	{
		throw std::invalid_argument{ "Divided by zero is occured" };
	}

	result.x = lhs.x / rhs.x;
	result.y = lhs.y / rhs.y;
	result.z = lhs.z / rhs.z;

	return result;
}

auto Point3D::operator <=>(const Point3D& rhs) const
{
	return this->length() <=> rhs.length();
}

double Point3D::length() const
{
	return sqrt(sqr(x) + sqr(y) + sqr(z));
}

double Point3D::distance(const Point3D& rhs) const
{
	Point3D p;
	p = (*this) - rhs;
	return  p.length();
}

double Point3D::distanceSqr(const Point3D& rhs) const
{
	Point3D p;
	p = (*this) - rhs;
	return (sqr(p.x) + sqr(p.y) + sqr(p.z));
}

//#define __MAIN__
#ifdef __MAIN__
int main()
{
	Point3D p{ 1, 2, 3 };
	cout << p << endl;

	Point3D p2{ 3, 2, 1 };
	Point3D p3 = p + p2;
	Point3D p4 = p + 2;


	p4 += 3;

	p4 -= 3;

	p4 *= 3;

	p4 /= 3;
	
	cout << p3 << endl;
	cout << p4 << endl;
	cout << fmt::format("length of p4 is {}", p4.length()) << endl;
	cout << fmt::format("the sqr({}) is {}", p4.x, sqr(p4.x)) << endl;
	if(p3 < p4)
	{
		cout << "p3 < p4" << endl;
	}
	else if(p3 > p4)
	{
		cout << "p3 > p4" << endl;
	}
	else
	{
		cout << "p3 == p4" << endl;
	}
	
	cout << fmt::format("the distance between p3 and p4 is {}", p3.distance(p4)) << endl;
	cout << fmt::format("the  sqaure of distance between p3 and p4 is {}", p3.distanceSqr(p4)) << endl;
	return 0;
}

#endif
