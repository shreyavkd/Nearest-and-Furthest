#include <iostream>
#include <chrono>
#include <string>
#include <omp.h>
#include "csv.h"
#include "system.h"
#include "point.h"

using namespace std;

//Interactive interface for user
//to view the summary and produce csv outputs
//for the selected method, initialisation, and algorithm.

//.................................................................
//IMPORTANT:
//To execute this program, please run the following
//lines on your terminal:
//
//g++ main.cpp system.cpp rng.cpp csv.cpp -O3 -fopenmp -o output
//./output
//
//.................................................................

int main()
{
    cout<< "\nProgram to compute nearest and farthest paths between data points in a 1x1 grid!";

    //I. Distance Measure
    cout<< "\nPlease choose your preferred methods:\n";
    cout<< "\nDistance Measure\nEnter 's' for standard or 'w' for wrap-around.\n";

    char distance_type;
    cin>> distance_type;

    //Setting geometry according to the user's choice.
    Geometry geom = (distance_type == 'w' || distance_type == 'W') ? Geometry::Wraparound : Geometry::Standard;

    System sys(geom);

    //II.  System Generation
    cout<< "\nSystem Generation\nEnter 'r' for random points or 'c' to input from a CSV file.\n";
    char sys_type;
    cin>> sys_type;


    std::string initName;
    bool usingRandom = false;

    //Initalising the system either with random
    //points or points from chosen csv file.
    switch (sys_type)
    {
        case 'r':
        case 'R':
        {
            usingRandom = true;
            initName = "random";

            cout<< "Enter number of points for the system.\n";
            std::size_t n; //size_t is the maximum size n can be - this ensures that we don't restrict the user's choice
            cin>> n;

            cout<< "Enter a seed.\n";
            unsigned int seed;
            cin>> seed;

            sys.generate_random_points(n, seed);
            break;
        }

        case 'c':
        case 'C':
        default:
        {
            usingRandom = false;

            cout<< "Enter 1 or 2 for the preferred CSV file:\n";
            cout<< "1: 100000 locations.csv\n";
            cout << "2: 200000 locations.csv\n";

            int csv_type = 0;
            cin >> csv_type;

            std::string fname;
            if (csv_type == 2)
            {
                fname = "200000 locations.csv";
                initName = "200k";
            }
            else
            {
                fname = "100000 locations.csv";
                initName = "100k";
            }

            sys.load_points(read_csv(fname)); //Reading csv file into the system.

            cout << "Loaded " << sys.size() << " points from " << fname << "\n";
            break;
        }
    }

    //III. Algorithm
    cout<< "\nEnter the respective serial number for your choice of algorithm:\n";
    cout<< "1: Naive serial\n";
    cout<< "2: Naive parallel\n";
    cout<< "3: Fast serial\n";
    cout<< "4: Fast parallel\n";

    int alg_type;
    cin>> alg_type;

    //Variables to store the details to print summary in a neat format.
    std::string algoName;
    std::string methodLabel;

    std::string geomName = (geom == Geometry::Standard) ? "Standard" : "Wraparound";

    //Setting values to the variables for displaying.
    switch (alg_type)
    {
        case 1:
            algoName = "naive_serial";
            methodLabel = "Naive serial";
            break;

        case 2:
            algoName = "naive_parallel";
            methodLabel = "Naive parallel";
            break;

        case 3:
            algoName = "fast_serial";
            methodLabel = "Fast serial";
            break;

        case 4:
            algoName = "fast_parallel";
            methodLabel = "Fast parallel";
            break;

        default:
            std::cerr << "Invalid algorithm selection.\n";
            return EXIT_FAILURE;
    }

    //Naming the output files.
    std::string nearestFile = "nearest_" + geomName + "_" + initName + "_" + algoName + ".csv";
    std::string farthestFile = "farthest_" + geomName + "_" + initName + "_" + algoName + ".csv";

    //Starting the timer to calculate the run-time of the algorithm.
    //double runtime stores the runtime at the end of the algorithm execution.
    double startTime = omp_get_wtime();

    //Running required algorithm.
    switch (alg_type)
    {
        case 1:
            sys.naive_serial();
            break;

        case 2:
            sys.naive_parallel();
            break;

        case 3:
            sys.fast_serial();
            break;

        case 4:
            sys.fast_parallel();
            break;
    }

    double runtime = omp_get_wtime() - startTime;

    sys.write_results(nearestFile, farthestFile);
    sys.print_summary(methodLabel, runtime);

    cout<< "\nOutput file for nearest distances: "<< nearestFile;
    cout<< "\nOutput file for farthest distances: "<< farthestFile;

    return EXIT_SUCCESS;
}