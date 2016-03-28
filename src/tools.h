/*
 * tool.h
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <cmath>
#include <cassert>
#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>

class Vector;

class Point {

public:
	Point();
	Point(double x, double y, double z);
	double distance(const Point &p1, const Point &p2);
	void setPoint(double x, double y, double z);
	double getX();
	double getY();
	double getZ();

	friend Vector operator^(const Vector& v1, const Vector& v2);
	friend Vector operator-(const Vector& v);
	friend Vector operator-(const Vector& v1, const Vector& v2);
	friend Vector operator*(const double c, const Vector& v2);
	friend Vector operator*(const Vector& v1, const double c);
	friend Vector operator+(const Vector& v1, const Vector& v2);
	friend double operator*(const Point& p, const Vector& v);
	friend double operator*(const Vector& v1, const Vector& v2);
	friend Vector operator/(const Vector& v1, const double c);
	friend Vector operator-(const Point& p1, const Point& p2);
	friend Vector crossProduct(Vector& v1, Vector& v2);

private:
	double x, y, z;
};

class Vector {

public:
	Vector();
	Vector(double x, double y, double z);
	void setVector(double x, double y, double z);
	double length();
	void normalize();

	double getX();
	double getY();
	double getZ();

	friend Vector operator^(const Vector& v1, const Vector& v2);
	friend Vector operator-(const Vector& v);
	friend Vector operator-(const Vector& v1, const Vector& v2);
	friend Vector operator*(const double c, const Vector& v2);
	friend Vector operator*(const Vector& v1, const double c);
	friend Vector operator+(const Vector& v1, const Vector& v2);
	friend double operator*(const Point& p, const Vector& v);
	friend double operator*(const Vector& v1, const Vector& v2);
	friend Vector operator/(const Vector& v1, const double c);
	friend Vector operator-(const Point& p1, const Point& p2);
	friend Vector crossProduct(Vector& v1, Vector& v2);

private:
	double x, y, z;
};

class Ray {

public:
	Ray();
	Ray(Point origin, Vector dir);
	void normalize();
	void setOrigin(const Point newOrigin);
	void setDir(const Vector newDir);

	Point origin;
	Vector dir;

};

class Render {

	friend class Output;

public:
	Render(int heigh, int width);
	void setPixel(int i, int j, Imf::Rgba *color);

private:
	Imf::Array2D<Imf::Rgba> pixels;
	int width;
	int height;

};

class Output {

public:
	static void renderToFile(const char *fileName, const Render *render);

private:
	static void writeRgba(const char *fileName, const Imf::Rgba *pixels, int width, int height);

};

#endif /* TOOLS_H_ */
