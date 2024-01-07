#include "icp.h"

namespace icp
{
#if defined(ICP_DOUBLE_SUPPORT)
template class ICP<double, double, 3>;
template class ICP<double, double, 2>;
#endif
template class ICP<float, float, 3>;
template class ICP<float, float, 2>;
}  // namespace icp
