#ifndef density_profiler_H
#define density_profiler_H
#include <string>
#include <vector>
#include <map>
using namespace std;

/** Represents a training interval in a given training input file. Based on the functionality present
 * within this class, it appears that it exists to provide average data values on regions selected by
 * a human observer so that the model can make inferences about that data.
 */
class gap_interval{
public:
	double start, stop;
	vector<double> X;
	vector<double> Y;
	gap_interval();
	gap_interval(double, double);
};

double get_table_mean_var(string, string,double, double, double);



#endif
