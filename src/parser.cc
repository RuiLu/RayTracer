/*
 * parser.cc
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 *      Uni   : rl2784
 *
 *  Reference: from -> Demo Code
 *  		   where -> (1) num_cams
 *  		   			(2)
 */

#include "parser.h"
#include "surface.h"
#include "sphere.h"
#include "triangle.h"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

double Parser::getTokenAsDouble (string inString, int whichToken)
{


    double thisDoubleVal = 0.;  // the return value

    if (whichToken == 0) {
        cerr << "error: the first token on a line is a character!" << endl;
        exit (-1);
    }

    // c++ string class has no super-easy way to tokenize, let's use c's:
    char *cstr = new char [inString.size () + 1];

    strcpy (cstr, inString.c_str());

    char *p = strtok (cstr, " ");
    if (p == 0) {
        cerr << "error: the line has nothing on it!" << endl;
        exit (-1);
    }

    for (int i = 0; i < whichToken; i++) {
        p = strtok (0, " ");
        if (p == 0) {
            cerr << "error: the line is not long enough for your token request!" << endl;
            exit (-1);
        }
    }

    thisDoubleVal = atof (p);

    delete[] cstr;

    return thisDoubleVal;
}

void Parser::parser(const char *filename,
					Camera &cam,
					std::vector< Surface * > &mysurfaces,
					Lights &lights)
{

	ifstream inFile(filename);    // open the file
    string line;

    if (! inFile.is_open ()) {  // if it's not open, error out.
        cerr << "can't open scene file" << endl;
        exit (-1);
    }

    Surface *thisSurface = 0;
    Material *lastMaterialLoaded = 0;
    int num_cams = 0;

    while (! inFile.eof ()) {   // go through every line in the file until finished


        getline (inFile, line); // get the line
        switch (line[0])  {     // we'll decide which command based on the first character

            case 's':
            	{
                // it's a sphere, load in the parameters

					double x, y, z, r;
					x = getTokenAsDouble (line, 1);
					y = getTokenAsDouble (line, 2);
					z = getTokenAsDouble (line, 3);
					r = getTokenAsDouble (line, 4);

					if (lastMaterialLoaded) {
						thisSurface = new Sphere(Point(x, y, z), r);
						thisSurface->setMaterial(*lastMaterialLoaded);
						thisSurface->setType("sphere");
						mysurfaces.push_back(thisSurface);
					}
            	}
                break;

            case 't':   // triangle
            	{
            		double v0_x, v0_y, v0_z;
            		double v1_x, v1_y, v1_z;
            		double v2_x, v2_y, v2_z;
            		v0_x = getTokenAsDouble(line, 1);
            		v0_y = getTokenAsDouble(line, 2);
            		v0_z = getTokenAsDouble(line, 3);
            		v1_x = getTokenAsDouble(line, 4);
            		v1_y = getTokenAsDouble(line, 5);
            		v1_z = getTokenAsDouble(line, 6);
            		v2_x = getTokenAsDouble(line, 7);
            		v2_y = getTokenAsDouble(line, 8);
            		v2_z = getTokenAsDouble(line, 9);

            		if (lastMaterialLoaded) {
            			thisSurface = new Triangle(Point(v0_x, v0_y, v0_z),
            									   Point(v1_x, v1_y, v1_z),
												   Point(v2_x, v2_y, v2_z));
            			thisSurface->setMaterial(*lastMaterialLoaded);
            			thisSurface->setType("triangle");
            			mysurfaces.push_back(thisSurface);
            		}
            	}
                break;

            case 'p':   // plane
            	{

            		double nx, ny, nz, scalar_value;
            		nx = getTokenAsDouble(line, 1);
            		ny = getTokenAsDouble(line, 2);
            		nz = getTokenAsDouble(line, 3);
            		scalar_value = getTokenAsDouble(line, 4);

            		if (lastMaterialLoaded) {
            			thisSurface = new Plane(Vector(nx, ny, nz), scalar_value);
            			thisSurface->setMaterial(*lastMaterialLoaded);
            			thisSurface->setType("plane");
            			mysurfaces.push_back(thisSurface);
            		}

            	}
                break;

            case 'c':    // camera
            	{
					// one trick here: the cameras pixel count (width, height) are integers,
					// so cast them.

					double pos_x = getTokenAsDouble(line, 1);
					double pos_y = getTokenAsDouble(line, 2);
					double pos_z = getTokenAsDouble(line, 3);
					double dir_x = getTokenAsDouble(line, 4);
					double dir_y = getTokenAsDouble(line, 5);
					double dir_z = getTokenAsDouble(line, 6);
					double focal_length = getTokenAsDouble(line, 7);
					double iw = getTokenAsDouble(line, 8);
					double ih = getTokenAsDouble(line, 9);
					int width = getTokenAsDouble(line, 10);
					int height = getTokenAsDouble(line, 11);

					Point eye(pos_x, pos_y, pos_z);
					Vector dir(dir_x, dir_y, dir_z);

					cam.initCamera(eye, dir, iw, ih, focal_length, width, height);

					num_cams++;
            	}
                break;

                //
                // lights:
                //
            case 'l':   // light

                // slightly different from the rest, we need to examine the second param,
                // which is at the third position on the line:
                switch (line[2]) {
                    case 'p':   // point light
                    	{

							double px = getTokenAsDouble(line, 2);
							double py = getTokenAsDouble(line, 3);
							double pz = getTokenAsDouble(line, 4);

							double ir = getTokenAsDouble(line, 5);
							double ig = getTokenAsDouble(line, 6);
							double ib = getTokenAsDouble(line, 7);

							Imf::Rgba intensity_p;
							intensity_p.r = ir;
							intensity_p.g = ig;
							intensity_p.b = ib;
							intensity_p.a = 1.0;

							PointLight *p = new PointLight(Point(px, py, pz), intensity_p);
							lights.pointLights.push_back(p);
                    	}
                        break;
                    case 'd':   // directional light
                    	{
							double vx = getTokenAsDouble(line, 2);
							double vy = getTokenAsDouble(line, 3);
							double vz = getTokenAsDouble(line, 4);

							double ir = getTokenAsDouble(line, 5);
							double ig = getTokenAsDouble(line, 6);
							double ib = getTokenAsDouble(line, 7);

							Imf::Rgba intensity_d;
							intensity_d.r = ir;
							intensity_d.g = ig;
							intensity_d.b = ib;
							intensity_d.a = 1.0;

							DirectionalLight *d = new DirectionalLight(Vector(vx, vy, vz), intensity_d);
							lights.directionalLights.push_back(d);
                    	}
                        break;
                    case 'a':   // ambient light
                    	{
							double ir = getTokenAsDouble(line, 2);
							double ig = getTokenAsDouble(line, 3);
							double ib = getTokenAsDouble(line, 4);

							Imf::Rgba intensity_a;
							intensity_a.r = ir;
							intensity_a.g = ig;
							intensity_a.b = ib;
							intensity_a.a = 1.0;

							AmbientLight *a = new AmbientLight(intensity_a);
							lights.ambientLights.push_back(a);
                    	}
                        break;
                    case 's': {	 // area light
                            double x, y, z;
                            double nx, ny, nz;
                            double ux, uy, uz;
                            double len;
                            double r, g, b;

                            x = getTokenAsDouble(line, 2);
                            y = getTokenAsDouble(line, 3);
                            z = getTokenAsDouble(line, 4);
                            Point center = Point(x, y, z);

                            nx = getTokenAsDouble(line, 5);
                            ny = getTokenAsDouble(line, 6);
                            nz = getTokenAsDouble(line, 7);
                            Vector dir = Vector(nx, ny, nz);

                            ux = getTokenAsDouble(line, 8);
                            uy = getTokenAsDouble(line, 9);
                            uz = getTokenAsDouble(line, 10);
                            Vector u_dir = Vector(ux, uy, uz);

                            len = getTokenAsDouble(line, 11);

                            r = getTokenAsDouble(line, 12);
                            g = getTokenAsDouble(line, 13);
                            b = getTokenAsDouble(line, 14);
                            Imf::Rgba intensity_area;
                            intensity_area.r = r;
                            intensity_area.g = g;
                            intensity_area.b = b;

                            AreaLight *areal = new AreaLight(center, dir, u_dir, len, intensity_area);
                            lights.areaLights.push_back(areal);

                        }
                        break;

                }

                break;

                //
                // materials:
                //
            case 'm': {  // material
                // the trick here: we should keep a pointer to the last material we read in,
                // so we can apply it to any subsequent geometry. Say it's called "lastMaterialLoaded"
                // we might then do something like this:
                //
                //  1. read in the 10 material parameters: dr, dg, db, sr, sg, sb, r, ir, ig, ib
                //  2. call lastMaterialLoaded->setMaterial(dr, dg, db,...);
                //

					Imf::Rgba diffuse, specular, ideal_specular;

					diffuse.r = getTokenAsDouble(line, 1);
					diffuse.g = getTokenAsDouble(line, 2);
					diffuse.b = getTokenAsDouble(line, 3);
					diffuse.a = 1.0;

					specular.r = getTokenAsDouble(line, 4);
					specular.g = getTokenAsDouble(line, 5);
					specular.b = getTokenAsDouble(line, 6);
					specular.a = 1.0;

					double phong_exponent = getTokenAsDouble(line, 7);

					ideal_specular.r = getTokenAsDouble(line, 8);
					ideal_specular.g = getTokenAsDouble(line, 9);
					ideal_specular.b = getTokenAsDouble(line, 10);
					ideal_specular.a = 1.0;

					lastMaterialLoaded = new Material(diffuse, specular, ideal_specular, phong_exponent);
				}
                break;


            case '/':
                // don't do anything, it's a comment
                break;


                //
                // options
                //
            case 'o':   // make your own options if you wish
                break;
                //
                // OBJ files
                //
            case 'w':
            	{
            		string obj_filename = line.substr(2, line.size() - 2);

            		std::vector< int > tris;
            		std::vector< float > verts;
            		read_wavefront_file(obj_filename, tris, verts);

            		double v0_x, v0_y, v0_z;
            		double v1_x, v1_y, v1_z;
            		double v2_x, v2_y, v2_z;

            		for (int i = 0; i < (int)(tris.size() / 3); i++) {
                    	v0_x = verts[3 * tris[3 * i]];
                    	v0_y = verts[3 * tris[3 * i] + 1];
                    	v0_z = verts[3 * tris[3 * i] + 2];
                    	v1_x = verts[3 * tris[3 * i + 1]];
                    	v1_y = verts[3 * tris[3 * i + 1] + 1];
                    	v1_z = verts[3 * tris[3 * i + 1] + 2];
                    	v2_x = verts[3 * tris[3 * i + 2]];
                    	v2_y = verts[3 * tris[3 * i + 2] + 1];
                    	v2_z = verts[3 * tris[3 * i + 2] + 2];

                    	if (lastMaterialLoaded) {
                    		thisSurface = new Triangle(Point(v0_x, v0_y, v0_z),
                    								   Point(v1_x, v1_y, v1_z),
        											   Point(v2_x, v2_y, v2_z));
                    		thisSurface->setMaterial(*lastMaterialLoaded);
                    		thisSurface->setType("triangle");
                    		mysurfaces.push_back(thisSurface);
                    	}
            		}
//            		std::cout<<obj_filename<<" is done..."<<std::endl;
            	}
            	break;
        }

    }
    if (num_cams != 1) {
    	std::cout << "scene file error: exactly ONE camera must be defined." << std::endl;
    }

}


void Parser::read_wavefront_file(string file, std::vector< int > &tris, std::vector< float > &verts)
{

    // clear out the tris and verts vectors:
    tris.clear ();
    verts.clear ();

    ifstream in(file);
    char buffer[1025];
    string cmd;


    for (int line=1; in.good(); line++) {
        in.getline(buffer,1024);
        buffer[in.gcount()]=0;

        cmd="";

        istringstream iss (buffer);

        iss >> cmd;

        if (cmd[0]=='#' || cmd.empty() || cmd[0] == 'g') {
            // ignore comments or blank lines
            continue;
        }
        else if (cmd=="v") {
            // got a vertex:

            // read in the parameters:
            double pa, pb, pc;
            iss >> pa >> pb >> pc;

            verts.push_back (pa);
            verts.push_back (pb);
            verts.push_back (pc);
         }
        else if (cmd=="f") {
            // got a face (triangle)

            // read in the parameters:
            int i, j, k;
            iss >> i >> j >> k;

            // vertex numbers in OBJ files start with 1, but in C++ array
            // indices start with 0, so we're shifting everything down by
            // 1
            tris.push_back (i-1);
            tris.push_back (j-1);
            tris.push_back (k-1);
        }
        else {
            std::cerr << "Parser error: invalid command at line " << line << std::endl;
        }

     }

    in.close();

//    std::cout << "found this many tris, verts: " << tris.size () / 3.0 << "  " << verts.size () / 3.0 << std::endl;

}


