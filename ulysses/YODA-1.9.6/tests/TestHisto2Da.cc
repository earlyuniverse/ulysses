#include "YODA/Histo2D.h"
#include "YODA/Scatter3D.h"
#include "YODA/WriterYODA.h"

#include <cmath>
#include <iostream>
#include <unistd.h>
#include <sys/time.h>

using namespace std;
using namespace YODA;


// A stats printing function.
//
// A very unpolished stats printing function. It prints some stats, sometimes
// looking at the full distribution and sometimes not. A better version is not
// too high on the priority list right now!
void printStats(Histo2D& h, bool full=false){
    cout << "-----------------------------" << endl;
    cout << "LowEdgeX = " << h.xMin() << " HighEdgeX = " << h.xMax() << endl;
    cout << "LowEdgeY = " << h.yMin() << " HighEdgeY = " << h.yMax() << endl;

    cout << "Sum of weights is " << h.sumW(true) << ", squared: " << h.sumW2(true) << endl;

    if (full) {
        cout << "Means: " << h.xMean(true) << " " << h.yMean(true) << endl;
        cout << "Variance: " << h.xVariance(true) << " " << h.yVariance(true) << endl;
        cout << "StdDevs: " << h.xStdDev(true) << " " << h.yStdDev(true) << endl;
    }
    cout << "-----------------------------" << endl;
}



int main() {

    cout << "-----------------------------" << endl;
    // Addition/subtraction:
    cout << "Creating histos to be added/subtracted/divided:" << endl;

    Histo2D first(10, 0, 100, 10, 0, 100);
    first.fill(1,1,1);
    first.fill(23,1,1);
    Histo2D second(10, 0, 100, 10, 0, 100);
    second.fill(1,1,1);

    cout << "Adding/Subtracting/Dividing" << endl;
    cout << "Testing the copy constructor:" << endl;
    Histo2D h(20, 0, 100, 20, 0, 100);
    Histo2D copyTest(h);
    cout << "Tested" << endl;

    Histo2D added(first+second);
    cout << "Addition? Copy constructor?" << endl;
    Histo2D subtracted(first-second);

    // cout << "Division crashes!" << endl;
    // Scatter3D divided(first/second);

    //const HistoBin2D& x =
    // static_cast<const Histo2D>(added).binAt(50,50);
    added.binAt(50,50);

    cout << "Done!" << endl;
    printStats(added);
    printStats(subtracted);

    // Write to stdout
    WriterYODA::write(cout, added);

    cout << "-----------------------------" << endl;
    return EXIT_SUCCESS;
}
