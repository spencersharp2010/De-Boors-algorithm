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
    if(recursionLevel == polynomialDegree + 1)
    {
        return {xCoordinates[knotSpanIndex], yCoordinates[knotSpanIndex]};
    }

    double gamma = (t - knotVector[knotSpanIndex])/(knotVector[knotSpanIndex+recursionLevel] - knotVector[knotSpanIndex]);
    std::array<double, 2> dummy1 = deBoor(t, knotSpanIndex - 1, polynomialDegree, knotVector, xCoordinates, yCoordinates, recursionLevel + 1);
    std::array<double, 2> dummy2 = deBoor(t, knotSpanIndex, polynomialDegree, knotVector, xCoordinates, yCoordinates, recursionLevel + 1);
    double resultX = (1 - gamma) * dummy1[0] + gamma * dummy2[0];
    double resultY = (1 - gamma) * dummy1[1] + gamma * dummy2[1];

    return { resultX, resultY };
}

size_t findKnotSpan( double t,
                     size_t numberOfControlPoints,
                     const std::vector<double>& knotVector )
{
    size_t index = 0;
    int size_of_knot_vector = knotVector.size();

    if (t > knotVector[size_of_knot_vector-1])
    {
        throw std::runtime_error("invalid input");
    }
    while(t >= knotVector[index])
    {
        if(index >= numberOfControlPoints)
        {
            break;
        }
        index++;
    }


    return index-1;
}

std::array<std::vector<double>, 2> evaluate2DCurveDeBoor( const std::vector<double>& tCoordinates,
                                                          const std::vector<double>& xCoordinates,
                                                          const std::vector<double>& yCoordinates,
                                                          const std::vector<double>& knotVector )
{
    size_t numberOfSamples = tCoordinates.size( );
    size_t numberOfPoints = xCoordinates.size( );
    size_t m = knotVector.size( );
    size_t p = m - numberOfPoints - 1;
    std::array<std::vector<double>, 2> returnValues;

    returnValues[0].reserve(numberOfSamples);
    returnValues[1].reserve(numberOfSamples);

    if( yCoordinates.size( ) != numberOfPoints )
    {
        throw std::runtime_error( "Inconsistent size in evaluate2DCurve." );
    }

    for(size_t i=0; i < numberOfSamples; i++)
    {
        size_t knotSpanIndex_x = findKnotSpan(tCoordinates[i],numberOfPoints,knotVector);

        std::array<double, 2> dummy = deBoor(tCoordinates[i],knotSpanIndex_x,p,knotVector,xCoordinates,yCoordinates);
        returnValues[0].push_back(dummy[0]);
        returnValues[1].push_back(dummy[1]);


    }

    return {returnValues};
    //return { std::vector<double>{ }, std::vector<double>{ } };
}

} // namespace splinekernel
} // namespace cie
