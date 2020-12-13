#pragma once
#include <ctime>
inline double second()
{
    return (double) std::clock()/(double) CLOCKS_PER_SEC;
}
