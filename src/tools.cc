/*
 * tools.cc
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#include "tools.h"

/*
 * Point
 */
Point::Point() {
	this->x = 0;
	this->y = 0;
	this->z = 0;
}

Point::Point(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

double Point::distance(const Point &p1, const Point &p2) {
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	double dz = p1.z - p2.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}

void Point::setPoint(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

double Point::getX() {
	return x;
}

double Point::getY() {
	return y;
}

double Point::getZ() {
	return z;
}

/*
 * Vector
 */
Vector::Vector() {
	this->x = 0;
	this->y = 0;
	this->z = 0;
}

Vector::Vector(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

void Vector::setVector(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

double Vector::length() {
	return sqrt(x * x + y * y + z * z);
}

void Vector::normalize() {

	assert (! (this->x == 0.0 && this->y == 0.0 && this->z == 0.0));
	double len = length();

//	if (len == 0. || len == 1.) return;
	x /= len;
	y /= len;
	z /= len;
}

double Vector::getX() {
	return x;
}

double Vector::getY() {
	return y;
}

double Vector::getZ() {
	return z;
}

/*
 * Basic mathematical operations of Vector
 */

Vector operator^ (const Vector &v1, const Vector &v2) {
	Vector res(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
	return res;
}

Vector operator+ (const Vector &v1, const Vector &v2) {
	Vector res(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
	return res;
}

Vector operator- (const Vector &v1, const Vector &v2) {
	Vector res(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
	return res;
}

Vector operator* (const Vector &v, const double c) {
    Vector res(c*v.x, c*v.y, c*v.z);
    return res;
}

Vector operator/ (const Vector &v, const double c) {
    Vector res(v.x/c, v.y/c, v.z/c);
    return res;
}

Vector operator* (const double c, const Vector &v) {
    Vector res(c*v.x, c*v.y, c*v.z);
    return res;
}

double operator* (const Point &p, const Vector &v) {
	return (p.x*v.x + p.y*v.y + p.z*v.z);
}

double operator* (const Vector &v1, const Vector &v2) {
	double res = (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
	return res;
}

Vector operator- (const Vector &v) {
	Vector res(-1.0 * v.x, -1.0 * v.y, -1.0 * v.z);
	return res;
}

Vector operator- (const Point &p1, const Point &p2) {
	Vector res(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
	return res;
}

Vector crossProduct(Vector &v1, Vector &v2) {
	double x = v1.y * v2.z - v1.z * v2.y;
	double y = v1.z * v2.x - v1.x * v2.z;
	double z = v1.x * v2.y - v1.y * v2.x;
	return Vector(x, y, z);
}

/*
 * Ray
 */
Ray::Ray() {}

Ray::Ray(Point origin, Vector dir) {
	this->origin = origin;
	this->dir = dir;
	this->dir.normalize();
}

void Ray::normalize() {
	this->dir.normalize();
}

void Ray::setOrigin(Point newOrigin) {
	this->origin = newOrigin;
}

void Ray::setDir(Vector newDir) {
	this->dir = newDir;
}

//Point Ray::getPoint() {
//	return this->origin;
//}
//
//Vector Ray::getVector() {
//	return this->dir;
//}

/*
 * Render
 */
Render::Render(int height, int width) {
	this->height = height;
	this->width = width;
	pixels.resizeErase(this->height, this->width);
}

void Render::setPixel(int i, int j, Imf::Rgba *color) {
	pixels[i][j] = *color;
}

/*
 * Output
 */
void Output::writeRgba(const char *fileName, const Imf::Rgba *pixels, int width, int height) {
    Imf::RgbaOutputFile file (fileName, width, height, Imf::WRITE_RGBA);
    file.setFrameBuffer (pixels, 1, width);
    file.writePixels (height);
}


void Output::renderToFile(const char *fileName, const Render *render) {
	writeRgba(fileName, &render->pixels[0][0], render->width, render->height);
}
