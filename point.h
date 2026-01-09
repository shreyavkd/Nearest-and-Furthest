#ifndef POINT_H
#define POINT_H

//Defining Point structure with (x,y) coordinates
//and Geom enum that applies the chosen method.

struct Point
{
    double x;
    double y;
};

enum class Geometry
{
    Standard,
    Wraparound
};

#endif