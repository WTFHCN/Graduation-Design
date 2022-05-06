#include <random>
#include <time.h>
#define LOWBIT(x) (x & -x)
#define DEBUG(x) std::cout << x << "\n"
namespace util
{
    std::default_random_engine random(time(NULL));
    std::uniform_real_distribution<double> randOut(0.0, 1.0);
    double GetLap(double u, double b)
    {
        double x = randOut(random);
        if (x < 0.5)
            return b * std::log(2 * x) + u;
        else
            return u - b * std::log(2 - 2 * x);
    }
}
