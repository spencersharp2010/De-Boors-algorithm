#ifndef CURVE_HPP
#define CURVE_HPP

#include <vector>
#include <array>

namespace cie
{
namespace splinekernel
{

/*! Evaluate B-Spline curve by summing up basis functions times control points.
 *  @param tCoordinates The parametric coordinates at which the curve shall be evaluated
 *  @param xCoordinates The x coordinates of the control points
 *  @param yCoordinates The y coordinates of the control points
 *  @return A vector of x and a vector of y coordinates with one value for each parametric
 *          coordinate tCoordinates
 */
std::array<std::vector<double>, 2> evaluate2DCurve( const std::vector<double>& tCoordinates,
                                                    const std::vector<double>& xCoordinates,
                                                    const std::vector<double>& yCoordinates,
                                                    const std::vector<double>& knotVector );

//! Identical to evaluate2DCurve, but using De Boor's algorithm.
std::array<std::vector<double>, 2> evaluate2DCurveDeBoor( const std::vector<double>& tCoordinates,
                                                          const std::vector<double>& xCoordinates,
                                                          const std::vector<double>& yCoordinates,
                                                          const std::vector<double>& knotVector );

/*! De Boor's algorithm for evaluating (x, y) at one parametric coordinate t. The parameter
 *  recursionLevel has a default value of 1, which will be used if no argument is passed. */
std::array<double, 2> deBoor( double t,
                              size_t knotSpanIndex,
                              size_t polynomialDegree,
                              const std::vector<double>& knotVector,
                              const std::vector<double>& xCoordinates,
                              const std::vector<double>& yCoordinates,
                              size_t recursionLevel = 1 );

//! Same as above but without recursion.
std::array<double, 2> deBoorOptimized( double t,
                                       size_t i,
                                       size_t p,
                                       const std::vector<double>& knotVector,
                                       const std::vector<double>& xCoordinates,
                                       const std::vector<double>& yCoordinates );

//! Determines the knot span of the parametric coordinate t.
size_t findKnotSpan( double t,
                     size_t numberOfControlPoints,
                     const std::vector<double>& knotVector );

} // namespace splinekernel
} // namespace cie

#endif // CURVE_HPP
