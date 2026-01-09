#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "system.h"
#include "rng.h"

using namespace std;

//System constructor with default values.
System::System(Geometry g): geometry_(g), mean_nearest_(0.0), mean_furthest_(0.0)
{ }

//Uses Mersenne Twister method to generate points uniformly.
void System::generate_random_points(std::size_t n, unsigned int seed)
{
    utils::rand::seedRand(seed);
    points_.clear();
    points_.reserve(n);

    for (std::size_t i = 0; i < n; ++i)
        points_.push_back({ utils::rand::randDouble(0.0,1.0),
                            utils::rand::randDouble(0.0,1.0) });
}

void System::load_points(const std::vector<Point>& pts)
{
    points_ = pts;
}

std::size_t System::size() const
{
    return points_.size();
}

//Creates arrays whose size matches the
//number of points in the system.
void System::allocate_arrays(std::size_t n)
{
    nearest_.assign(n, 0.0);
    furthest_.assign(n, 0.0);
}

//Finds the absolute difference between
//coordinates. The computed differences
//are stored as squares temporarily.
void System::separation(const Point& a, const Point& b, double& dx, double& dy) const
{
    dx = a.x-b.x;
    if (dx<0) dx = -dx;
    dy = a.y-b.y;
    if (dy<0) dy = -dy;

    if (geometry_ == Geometry::Wraparound)
    {
        if (dx>0.5) dx = 1.0 - dx;
        if (dy>0.5) dy = 1.0 - dy;
    }
}

//Computing distance for naive method:
//Returns Standard/Wrap-around distance.
double System::distance(const Point& a, const Point& b) const
{
    double dx, dy;
    separation(a, b, dx, dy);
    return std::sqrt(dx*dx + dy*dy);
}

// NAIVE SERIAL ALGORITHM
//Computes distances for ordered pairs (i,j)
//where i<j making sure each pair is
//visited only once.
void System::naive_serial()
{
    const std::size_t n = points_.size();
    allocate_arrays(n);

    const double INF = std::numeric_limits<double>::infinity();
    double sumN = 0, sumF = 0;

    for (std::size_t i = 0; i < n; ++i)
    {
        double best = INF, worst = 0.0; //Initialising dist_max = 0, and dist_min=infinity as 
                                        //best and worst (to avoid typos).

        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j) continue;
            double d = distance(points_[i], points_[j]);
            if (d < best)  best = d;
            if (d > worst) worst = d;
        }

        nearest_[i] = best;
        furthest_[i] = worst;

        sumN+= best;
        sumF+= worst;
    }

    mean_nearest_  = sumN/n;
    mean_furthest_ = sumF/n;
}

// NAIVE PARALLEL ALGORITHM
//Outer loop is divided among the threads.
//dynamic scheduling makes the program faster
//than static scheduling since it avoids
//idle time.
//reduction is used to avoid race condition.
void System::naive_parallel()
{
    const std::size_t n = points_.size();
    allocate_arrays(n);

    const double INF = std::numeric_limits<double>::infinity();
    double sumN = 0, sumF = 0;

    #pragma omp parallel for reduction(+:sumN,sumF) schedule(dynamic, 1)
    for (long ii = 0; ii < (long)n; ++ii)
    {
        std::size_t i = (std::size_t)ii;

        double best = INF, worst = 0.0;

        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j) continue;
            double d = distance(points_[i], points_[j]);
            if (d < best)  best = d;
            if (d > worst) worst = d;
        }

        nearest_[i] = best;
        furthest_[i] = worst;

        sumN+= best;
        sumF+= worst;
    }

    mean_nearest_  = sumN/n;
    mean_furthest_ = sumF/n;
}

// FAST SERIAL ALGORITHM
//Computes distances for unordered pairs
//and does not implement the square root inside
//the loop, which improves run-time.
void System::fast_serial()
{
    const std::size_t n = points_.size();
    allocate_arrays(n);

    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> nearSq(n, INF), farSq(n, 0.0);

    if (geometry_ == Geometry::Standard)
    {
        for (std::size_t i = 0; i < n-1; ++i)
        {
            double xi = points_[i].x;
            double yi = points_[i].y;

            for (std::size_t j = i+1; j < n; ++j)
            {
                double dx = xi - points_[j].x; if (dx < 0) dx = -dx;
                double dy = yi - points_[j].y; if (dy < 0) dy = -dy;

                double d2 = dx*dx + dy*dy;

                if (d2 < nearSq[i]) nearSq[i] = d2;
                if (d2 < nearSq[j]) nearSq[j] = d2;

                if (d2 > farSq[i]) farSq[i] = d2;
                if (d2 > farSq[j]) farSq[j] = d2;
            }
        }
    }
    else   // Wrap-around Method
    {
        for (std::size_t i = 0; i < n-1; ++i)
        {
            double xi = points_[i].x;
            double yi = points_[i].y;

            for (std::size_t j = i+1; j < n; ++j)
            {
                double dx = xi - points_[j].x; if (dx < 0) dx = -dx;
                double dy = yi - points_[j].y; if (dy < 0) dy = -dy;

                if (dx > 0.5) dx = 1.0 - dx;
                if (dy > 0.5) dy = 1.0 - dy;

                double d2 = dx*dx + dy*dy;

                if (d2 < nearSq[i]) nearSq[i] = d2;
                if (d2 < nearSq[j]) nearSq[j] = d2;

                if (d2 > farSq[i]) farSq[i] = d2;
                if (d2 > farSq[j]) farSq[j] = d2;
            }
        }
    }

    double sumN=0, sumF=0;
    for (std::size_t i=0;i<n;++i)
    {
        nearest_[i] = std::sqrt(nearSq[i]);
        furthest_[i] = std::sqrt(farSq[i]);
        sumN += nearest_[i];
        sumF += furthest_[i];
    }

    mean_nearest_  = sumN / n;
    mean_furthest_ = sumF / n;
}

// FAST PARALLEL ALGORITHM
//Each thread works on a local copy
//of the array and outer loop divides tasks dynamically.
//dynamic scheduling once again improves run-time
//compared to static scheduling by avoiding idle time.
//After parallelisation, the local results are joined
//by taking the minimum and maximum of the results
//for each point i.
void System::fast_parallel()
{
    const std::size_t n = points_.size();
    allocate_arrays(n);

    const double INF = std::numeric_limits<double>::infinity();
    int max_threads = omp_get_max_threads();

    std::vector<std::vector<double>> nearSq(max_threads,
        std::vector<double>(n, INF));
    std::vector<std::vector<double>> farSq(max_threads,
        std::vector<double>(n, 0.0));

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        if (geometry_ == Geometry::Standard)
        {
            #pragma omp for schedule(dynamic, 1)
            for (long ii = 0; ii < (long)n - 1; ++ii)
            {
                std::size_t i = (std::size_t)ii;
                double xi = points_[i].x;
                double yi = points_[i].y;

                for (std::size_t j = i+1; j < n; ++j)
                {
                    double dx = xi - points_[j].x; if (dx < 0) dx = -dx;
                    double dy = yi - points_[j].y; if (dy < 0) dy = -dy;

                    double d2 = dx*dx + dy*dy;

                    if (d2 < nearSq[threadID][i]) nearSq[threadID][i] = d2;
                    if (d2 < nearSq[threadID][j]) nearSq[threadID][j] = d2;

                    if (d2 > farSq[threadID][i]) farSq[threadID][i] = d2;
                    if (d2 > farSq[threadID][j]) farSq[threadID][j] = d2;
                }
            }
        }
        else 
        // wrap-around method
        {
            #pragma omp for schedule(dynamic, 1)
            for (long ii = 0; ii < (long)n - 1; ++ii)
            {
                std::size_t i = (std::size_t)ii;
                double xi = points_[i].x;
                double yi = points_[i].y;

                for (std::size_t j = i+1; j < n; ++j)
                {
                    double dx = xi - points_[j].x; if (dx < 0) dx = -dx;
                    double dy = yi - points_[j].y; if (dy < 0) dy = -dy;

                    if (dx>0.5) dx = 1.0 - dx;
                    if (dy>0.5) dy = 1.0 - dy;

                    double d2 = dx*dx + dy*dy;

                    if (d2 < nearSq[threadID][i]) nearSq[threadID][i] = d2;
                    if (d2 < nearSq[threadID][j]) nearSq[threadID][j] = d2;

                    if (d2 > farSq[threadID][i]) farSq[threadID][i] = d2;
                    if (d2 > farSq[threadID][j]) farSq[threadID][j] = d2;
                }
            }
        }
    }

    double sumN = 0, sumF = 0;

    for (std::size_t i = 0; i < n; ++i)
    {
        double mn = INF, mx = 0.0;

        for (int t=0; t<max_threads; ++t)
        {
            if (nearSq[t][i] < mn) mn = nearSq[t][i];
            if (farSq[t][i]  > mx) mx = farSq[t][i];
        }

        nearest_[i]  = std::sqrt(mn);
        furthest_[i] = std::sqrt(mx);

        sumN+= nearest_[i];
        sumF+= furthest_[i];
    }

    mean_nearest_  = sumN/n;
    mean_furthest_ = sumF/n;
}

//Outputing the CSV files - one for farthest distances
//and one for shortest.
void System::write_results(const std::string& nearestFile,
                           const std::string& furthestFile) const
{
    std::ofstream f1(nearestFile.c_str());
    if (!f1) throw std::runtime_error("Cannot open output: " + nearestFile);
    for (double v : nearest_) f1 << v << "\n";

    std::ofstream f2(furthestFile.c_str());
    if (!f2) throw std::runtime_error("Cannot open output: " + furthestFile);
    for (double v : furthest_) f2 << v << "\n";
}

//Printing summary onto console.
void System::print_summary(const std::string& label, double sec) const
{
    cout<< "Method: " << label << "\n";
    cout<< "Geometry: "
              << (geometry_ == Geometry::Standard ? "Standard" : "Wraparound") << "\n";
    cout<< "Points:   " << points_.size() << "\n";
    cout<< "Nearest mean:  " << mean_nearest_  << "\n";
    cout<< "Farthest mean: " << mean_furthest_ << "\n";
    cout<< "Time: " << sec << " seconds\n";
}