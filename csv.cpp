#include "csv.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<Point> read_csv(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("Could not open CSV file: " + filename);

    std::vector<Point> pts;
    pts.reserve(1024);

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        std::istringstream ss(line);
        double x, y;
        char comma;

        if (ss >> x >> comma >> y)
            pts.push_back({x, y});
    }

    if (pts.empty())
        throw std::runtime_error("CSV file contained no valid points: " + filename);

    return pts;
}