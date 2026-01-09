#ifndef CSV_H
#define CSV_H

#include <string>
#include <vector>
#include "point.h"

std::vector<Point> read_csv(const std::string& filename);

#endif