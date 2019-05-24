#ifndef BSPLINE_HPP
#define BSPLINE_HPP

#include <vector>
#include "stddef.h"

namespace cie
{
namespace splinekernel
{

/*! Evaluates one B-Spline basis function.
 *  @param t The parametric coordinate
 *  @param i The basis function index
 *  @param p The polynomial degree
 *  @param knotVector
 *  @return The function value
 */
double evaluateBSplineBasis( double t, size_t i, size_t p, const std::vector<double>& knotVector );

} // namespace splinekernel
} // namespace cie

#endif // BSPLINE_HPP
