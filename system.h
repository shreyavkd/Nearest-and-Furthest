#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <string>
#include "point.h"

class System
{
public:
    explicit System(Geometry g);

    //Initialising the grid with points - either randomly
    //or from a CSV.
    void generate_random_points(std::size_t n, unsigned int seed);
    void load_points(const std::vector<Point>& pts);
    std::size_t size() const;

    //4 different algorithms.
    void naive_serial();
    void naive_parallel();
    void fast_serial();
    void fast_parallel();

    //Console and CSV Outputs
    void write_results(const std::string& nearestFile,
                       const std::string& furthestFile) const;

    void print_summary(const std::string& label, double seconds) const;

private:
    Geometry geometry_;
    std::vector<Point> points_;
    std::vector<double> nearest_, furthest_;
    double mean_nearest_, mean_furthest_;

    void allocate_arrays(std::size_t n); //Helper function.

    //Helpers for naive methods.
    void separation(const Point& a, const Point& b, double& dx, double& dy) const;
    double distance(const Point& a, const Point& b) const;
};

#endif