/*
 * parser.h
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstring>
#include <string>
#include "camera.h"


class Parser {
public:
	void parser(const char *filename,
				Camera &cam,
				std::vector< Surface * > &mysurfaces,
				Lights &lights);
	double getTokenAsDouble (std::string inString, int whichToken);
	void read_wavefront_file (
		    string file,
		    std::vector< int > &tris,
		    std::vector< float > &verts);
};


#endif /* PARSER_H_ */
