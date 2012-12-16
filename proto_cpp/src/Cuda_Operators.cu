#include "Cuda_Operators.h"

std::ostream& operator<<(std::ostream& os, const cuComplex& z)
{
    os << "(" << z.x << " , " << z.y << ")";
    return os;
}