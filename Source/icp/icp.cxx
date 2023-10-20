#include "icp.hpp"

namespace icp
{
template class ICP<double, double, 3>;
template class ICP<double, double, 2>;
template class ICP<float, float, 3>;
template class ICP<float, float, 2>;
}
