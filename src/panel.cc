/*
 * panel.cc
 *
 *  Created on: Feb 17, 2016
 *      Author: Rui
 */

#include "panel.h"

Imf::Rgba Panel::black = Imf::Rgba(0,0,0,1);

Panel::Panel(int r_samples, int s_samples) {
	this->black = Imf::Rgba(0,0,0,1);
	this->bvhTree = NULL;
	this->r_samples = r_samples;
	this->s_samples = s_samples;
}

Render* Panel::render(Camera &cam,
					  std::vector< Surface * > &mysurfaces,
					  Lights &lights)
{
	this->bvhTree = new BvhNode(mysurfaces, X_AXIS);
	this->bvhTree->create();

    int width = cam.width;
    int height = cam.height;

    double t = cam.t;
    double b = cam.b;

    double l = cam.l;
    double r = cam.r;

    Render *newRender = new Render(height, width);

    srand((unsigned)time(NULL));
    double rand_n;
    int ray_samples = this->r_samples;

    int totalPixels = width * height;
    std::cout << "EXR size: " << width << " " << height << std::endl;
    int renderedPixels = 0;
    double infinity = std::numeric_limits<double>::max();
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            Imf::Rgba c = black;

            for (int p = 0; p < ray_samples; p++) {
            	for (int q = 0; q < ray_samples; q++) {
            		rand_n = ((double)rand()) / RAND_MAX;
            		double u = l + (r - l) * (i + (p + rand_n) / ray_samples) / width;
            		double v = b + (t - b) * (j + (q + rand_n) / ray_samples) / height;
            		Ray primaryRay = Ray(cam.eye, u * cam.u + v * cam.v + (-1.) * cam.focal_length * cam.w);
            		primaryRay.normalize();
            		Vector color = getColorBvh(primaryRay, 0.0001, infinity, 20);
            		c.r += color.getX();
            		c.g += color.getY();
            		c.b += color.getZ();
            	}
            }

            c.r = c.r / (ray_samples * ray_samples);
            c.g = c.g / (ray_samples * ray_samples);
            c.b = c.b / (ray_samples * ray_samples);

            newRender->setPixel(j, i, &c);
            renderedPixels++;
            double pert = ((double)renderedPixels / (double)totalPixels) * 100.;
            std::cout<<pert<<" % "<<std::endl;
        }
    }
    std::cout << "Processing: 100% " << std::endl;
    deleteTree(bvhTree);
    return newRender;
}

Vector Panel::getColorBvh(Ray ray, double t0, double t1, int recurse_limit)
{
	Vector c = Vector(0, 0, 0);

	if (recurse_limit == 0) {
		return c;
	}

	// rec is for record of normal ray, srec is for record of shadow ray which is useless
	HitRecord rec, srec;
	rec.hits.resize(0);
	rec.hitNums = 0;

	bool isIntersected = this->bvhTree->isIntersected(ray, t0, t1, &rec);

	if (isIntersected) {
//		double infinity = std::numeric_limits<double>::infinity();
		double that = rec.hits[0];
    	Vector n = rec.normal;

    	// get the intersection point
    	Point eye = ray.origin;
    	Point ip = Point(eye.getX() + that * ray.dir.getX(), eye.getY() + that * ray.dir.getY(), eye.getZ() + that * ray.dir.getZ());

    	Vector v = -ray.dir;
    	v.normalize();

    	double lamr = 0.0, lamg = 0.0, lamb = 0.0;
    	double bpr = 0.0, bpg = 0.0, bpb = 0.0;
    	double ambient_r = 0.0, ambient_g = 0.0, ambient_b = 0.0;

    	// Point light
    	for (int y = 0; y < lights.pointLights.size(); y++) {

    		PointLight *pl = lights.pointLights[y];

    		Vector l = pl->getPoint() - ip;
    		// computer the distance between light source and intersection point
    		double t_l_to_ip = l.length();
    		l.normalize();

    		// diffuse shading
    		double lDotNormal = n * l;

    		// shadow ray
    		Ray raytmp = Ray(ip, l);

    		if (!this->bvhTree->isIntersected(raytmp, 0.001, t_l_to_ip, &srec)) {
        		lDotNormal = lDotNormal < 0.0 ? 0.0 : lDotNormal;

        		lamr += rec.diffuse.r * pl->getIntensity_p().r * lDotNormal;
        		lamg += rec.diffuse.g * pl->getIntensity_p().g * lDotNormal;
        		lamb += rec.diffuse.b * pl->getIntensity_p().b * lDotNormal;

        		// Blinn-Phong shading
        		Vector h = v + l;

        		if (h.length() > 0.0) {

        			h.normalize();

        			double exponent_res = n * h;
        			exponent_res = exponent_res < 0.0 ? 0.0 : exponent_res;
        			exponent_res = pow (exponent_res, rec.phong_exp);

        			bpr += rec.specular.r * pl->getIntensity_p().r * exponent_res;
        			bpg += rec.specular.g * pl->getIntensity_p().g * exponent_res;
        			bpb += rec.specular.b * pl->getIntensity_p().b * exponent_res;

        		}
    		}
    	}

    	//Directional lights
    	for (int d = 0; d < lights.directionalLights.size(); d++) {
    		DirectionalLight *dl = lights.directionalLights[d];

    		Vector l = dl->getDir();

    		// computer the distance between light source and intersection point
    		double t_l_to_ip = l.length();

    		l.normalize();

    		Ray raytmp = Ray(ip, l);

    		if (!bvhTree->isIntersected(raytmp, 0.001, t_l_to_ip, &srec)) {

        		// diffuse shading
        		double lDotNormal = n * l;

        		lDotNormal = lDotNormal < 0.0 ? 0.0 : lDotNormal;

        		lamr += rec.diffuse.r * dl->getIntensity_d().r * lDotNormal;
        		lamg += rec.diffuse.g * dl->getIntensity_d().g * lDotNormal;
        		lamb += rec.diffuse.b * dl->getIntensity_d().b * lDotNormal;

        		// Blinn-Phong shading
        		Vector h = v + l;

        		if (h.length() > 0.0) {

        			h.normalize();

        			double exponent_res = n * h;
        			exponent_res = exponent_res < 0.0 ? 0.0 : exponent_res;
        			exponent_res = pow (exponent_res, rec.phong_exp);

        			bpr += rec.specular.r * dl->getIntensity_d().r * exponent_res;
        			bpg += rec.specular.g * dl->getIntensity_d().g * exponent_res;
        			bpb += rec.specular.b * dl->getIntensity_d().b * exponent_res;

        		}
    		}
    	}

        // Ambient light, only does once
    	// and add ambient light to every place except complete missing
    	for (int a = 0; a < lights.ambientLights.size(); a++) {

			ambient_r = rec.diffuse.r * lights.ambientLights[0]->getIntensity_a().r;
			ambient_g = rec.diffuse.g * lights.ambientLights[0]->getIntensity_a().g;
			ambient_b = rec.diffuse.b * lights.ambientLights[0]->getIntensity_a().b;
    	}

    	// Area light
    	int shadow_samples = this->s_samples;

    	for (int ar = 0; ar < lights.areaLights.size(); ar++) {
    		AreaLight *area = lights.areaLights[ar];

    		double ir = area->getIntensity_area().r / (shadow_samples * shadow_samples);
    		double ig = area->getIntensity_area().g / (shadow_samples * shadow_samples);
    		double ib = area->getIntensity_area().b / (shadow_samples * shadow_samples);

    		double len = area->getLen();

    		Point center = area->getCenter();

    		Vector normal = area->getDir();
    		Vector u = area->getUDir();
    		Vector v_here = crossProduct(u, normal);
    		v_here.normalize();

    		for (double sx = 0; sx < shadow_samples; sx++) {
    			for (double sy = 0; sy < shadow_samples; sy++) {
    				double rand_a = ((double) rand()) / RAND_MAX;
    				double rand_b = ((double) rand()) / RAND_MAX;
//    				double a = sx+(interval*rand_a);
//    				double b = sy+(interval*rand_b);

    				double p_x = center.getX() + (rand_a - 0.5) * len * u.getX() +  (rand_b - 0.5) * len * v_here.getX();
    				double p_y = center.getY() + (rand_a - 0.5) * len * u.getY() +  (rand_b - 0.5) * len * v_here.getY();
    				double p_z = center.getZ() + (rand_a - 0.5) * len * u.getZ() +  (rand_b - 0.5) * len * v_here.getZ();
    				Point pos = Point(p_x, p_y, p_z);

        			Vector l = pos - ip;
        			double tmax = l.length();
        			l.normalize();
            		// shadow ray
            		Ray raytmp = Ray(ip, l);

            		double lDotNormal = n * l;
            		if (!this->bvhTree->isIntersected(raytmp, 0.001, tmax, &srec)) {
                		lDotNormal = lDotNormal < 0.0 ? 0.0 : lDotNormal;

                		lamr += rec.diffuse.r * ir * lDotNormal / ((tmax + 1) * (tmax + 1));
                		lamg += rec.diffuse.g * ig * lDotNormal / ((tmax + 1) * (tmax + 1));
                		lamb += rec.diffuse.b * ib * lDotNormal / ((tmax + 1) * (tmax + 1));

                		// Blinn-Phong shading
                		Vector h = v + l;

                		if (h.length() > 0.0) {

                			h.normalize();

                			double exponent_res = n * h;
                			exponent_res = exponent_res < 0.0 ? 0.0 : exponent_res;
                			exponent_res = pow (exponent_res, rec.phong_exp);

                			bpr += rec.specular.r * ir * exponent_res / ((tmax + 1) * (tmax + 1));
                			bpg += rec.specular.g * ig * exponent_res / ((tmax + 1) * (tmax + 1));
                			bpb += rec.specular.b * ib * exponent_res / ((tmax + 1) * (tmax + 1));

                		}
            		}
    			}
    		}
    	}


    	if (rec.ideal.r == 0.0 && rec.ideal.g == 0.0 && rec.ideal.b == 0.0) {
    		return Vector(lamr + bpr + ambient_r, lamg + bpg + ambient_g, lamb + bpb + ambient_b);
    	} else {

    		Vector v = -ray.dir;
    		v.normalize();

    		Vector r = 2 * (n * v) * n - v;
    		r.normalize();

    		Ray reflection_ray = Ray(ip, r);

    		Vector reflection_color = getColorBvh(reflection_ray, 0.001, std::numeric_limits<double>::infinity(), recurse_limit - 1);

    		// I choose 0.9 as the radiance weighted coefficient
    		double reflection_r = rec.ideal.r * reflection_color.getX() * 0.9;
    		double reflection_g = rec.ideal.g * reflection_color.getY() * 0.9;
    		double reflection_b = rec.ideal.b * reflection_color.getZ() * 0.9;

    		return Vector(lamr + bpr + reflection_r, lamg + bpg + reflection_g, lamb + bpb + reflection_b);
    	}


    } else {
    	return c;
    }
}

void Panel::deleteTree(BvhNode *root){
    if(root){
        deleteTree(dynamic_cast<BvhNode * >(root->left));
        deleteTree(dynamic_cast<BvhNode *>(root->right));
        delete root;
    }
}
