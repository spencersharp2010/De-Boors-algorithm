#include "basisfunctions.hpp"

#include <string>
#include <cmath>
#include <stdexcept>

namespace cie
{
namespace splinekernel
{

double evaluateBSplineBasis( double t, size_t i, size_t p, const std::vector<double>& knotVector )
{
  double tolerance = 1e-12; // Numerically zero.
  size_t m = knotVector.size() - 1;

  // check if i is in interval 0 <= i <= n. (with n = m-p-1) and if t is in the interval of the knot vector
  if( i > ( m - p - 1 ) )
  {
    throw std::range_error( "Index " + std::to_string( i ) + " out of range!");
  }

  if( t < knotVector.front() || t > knotVector.back( ) )
  {
    throw std::range_error("t is not with the interval!");
  }

  if( p == 0 )
  {
    return ( ( t >= knotVector[i] ) && ( t < knotVector[i + 1] ) ) || // t is in [t_i, t_{i+1}]
           ( std::abs( knotVector.back( ) - knotVector[i + 1] ) < tolerance && // corner case of t = t_n
             std::abs( knotVector.back( ) - t             ) < tolerance ); 
  }
  else
  {
    double result = 0.0;

    double a = t - knotVector[i];
    double b = knotVector[i + p] - knotVector[i];

    if( std::abs( b ) > tolerance )
    {
      result += a / b * evaluateBSplineBasis( t, i, p - 1, knotVector );
    }

    a = knotVector[i + p + 1] - t;
    b = knotVector[i + p + 1] - knotVector[i + 1];

    if( std::abs( b ) > tolerance )
    {
      result += a / b * evaluateBSplineBasis( t, i + 1, p - 1, knotVector );
    }
    
    return result;
  }
}

} // namespace splinekernel
} // namespace cie
