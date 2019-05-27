#include "basisfunctions.hpp"
#include "curve.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace cie
{
namespace splinekernel
{

std::array<std::vector<double>, 2> evaluate2DCurve( const std::vector<double>& tCoordinates, 
                                                    const std::vector<double>& xCoordinates,
                                                    const std::vector<double>& yCoordinates, 
                                                    const std::vector<double>& knotVector )
{
    size_t numberOfSamples = tCoordinates.size( );
    size_t numberOfPoints = xCoordinates.size( );
    size_t m = knotVector.size( );
    size_t p = m - numberOfPoints - 1;

    if( yCoordinates.size( ) != numberOfPoints )
    {
        throw std::runtime_error( "Inconsistent size in evaluate2DCurve." );
    }

    std::vector<double> curveX( numberOfSamples, 0.0 );
    std::vector<double> curveY( numberOfSamples, 0.0 );
    
    for( size_t i = 0; i < numberOfSamples; ++i )
    {
        double t = tCoordinates[i];

        for( size_t j = 0; j < numberOfPoints; ++j )
        {
            double N = evaluateBSplineBasis( t, j, p, knotVector );

            curveX[i] += N * xCoordinates[j];
            curveY[i] += N * yCoordinates[j];
        }
    }

    return { curveX, curveY };
}

std::array<double, 2> deBoor( double t,
                              size_t knotSpanIndex,
                              size_t polynomialDegree,
                              const std::vector<double>& knotVector,
                              const std::vector<double>& xCoordinates,
                              const std::vector<double>& yCoordinates,
                              size_t recursionLevel )
{
    return { 0.0, 0.0 };
}

size_t findKnotSpan( double t,
                     size_t numberOfControlPoints,
                     const std::vector<double>& knotVector )
{
    return 0.0;
}

std::array<std::vector<double>, 2> evaluate2DCurveDeBoor( const std::vector<double>& tCoordinates,
                                                          const std::vector<double>& xCoordinates,
                                                          const std::vector<double>& yCoordinates,
                                                          const std::vector<double>& knotVector )
{
    return { std::vector<double>{ }, std::vector<double>{ } };
}

} // namespace splinekernel
} // namespace cie
