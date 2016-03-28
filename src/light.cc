/*
 * light.c
 *
 *  Created on: Feb 17, 2016
 *      Author: Rui
 */

#include "light.h"

/*
 * Point light
 */
PointLight::PointLight(Point pos, Imf::Rgba instensity_p) {
	this->pos = pos;
	this->intensity_p = instensity_p;
}

Point PointLight::getPoint() {
	return this->pos;
}

Imf::Rgba PointLight::getIntensity_p() {
	return this->intensity_p;
}

/*
 * Directional light
 */
DirectionalLight::DirectionalLight(Vector dir, Imf::Rgba intensity_d) {
	this->dir = dir;
	this->intensity_d = intensity_d;
}

Vector DirectionalLight::getDir() {
	return this->dir;
}

Imf::Rgba DirectionalLight::getIntensity_d() {
	return this->intensity_d;
}

/*
 * Ambient light
 */
AmbientLight::AmbientLight(Imf::Rgba intensity_a) {
	this->intensity_a = intensity_a;
}

Imf::Rgba AmbientLight::getIntensity_a() {
	return this->intensity_a;
}

/*
 * Area Light
 */
AreaLight::AreaLight(Point center, Vector dir, Vector u_dir, double len, Imf::Rgba intensity_area) {
	this->center = center;
	this->dir = dir;
	this->u_dir = u_dir;
	this->len = len;
	this->intensity_area = intensity_area;
}

Point AreaLight::getCenter() {
	return this->center;
}

Vector AreaLight::getDir() {
	return this->dir;
}

Vector AreaLight::getUDir() {
	return this->u_dir;
}

double AreaLight::getLen() {
	return this->len;
}

Imf::Rgba AreaLight::getIntensity_area() {
	return this->intensity_area;
}
