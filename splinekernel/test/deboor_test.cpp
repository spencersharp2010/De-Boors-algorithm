#include "catch.hpp"
#include "curve.hpp"

#include <array>
#include <vector>

namespace cie
{
namespace splinekernel
{

TEST_CASE("FindKnotSpan_test")
{
    std::vector<double> knotVector{ 0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 9.0, 9.0, 9.0, 9.0 };

    std::vector<double> x{ 2.0, 3.0, 0.5, 0.2, 0.5, 3.0 };
    std::vector<double> y{ 1.0, 3.0, 3.0, 4.0, 5.0, 6.0 };

    // So p = 3 but thats irrelevant here;
    size_t n = x.size( );

    CHECK( findKnotSpan( 0.0, n, knotVector ) == 3 );
    CHECK( findKnotSpan( 0.2, n, knotVector ) == 3 );
    CHECK( findKnotSpan( 1.1, n, knotVector ) == 4 );
    CHECK( findKnotSpan( 9.0, n, knotVector ) == 5 );

    CHECK_THROWS(findKnotSpan(9.1, n, knotVector));
}

TEST_CASE("DeBoor_test")
{
    std::vector<double> knotVector{ 0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 9.0, 9.0, 9.0, 9.0 };

    std::vector<double> x{ 0.5, 10.0, 9.0, 4.5, 1.5, 1.0 };
    std::vector<double> y{ 0.5,  1.0, 4.0, 7.5, 6.0, 1.0 };

    std::array<double, 3> t{ 0.0, 1.1, 9.0 };

    std::array<double, 2> P;

    size_t n = x.size( );
    size_t m = knotVector.size();
    size_t p = m - n - 1;
    size_t s = 77;

    // For t0 (Beginning of curve)
    REQUIRE_NOTHROW( s = findKnotSpan( t[0], n, knotVector ) );
    REQUIRE_NOTHROW( P = deBoor( t[0], s, p, knotVector, x, y ) );

    // x-coordinates of curve
    CHECK( P[0] == Approx( x[0] ) );
    CHECK( P[1] == Approx( y[0] ) );

    // For t1 in the middle of the curve (see skript)
    REQUIRE_NOTHROW( s = findKnotSpan( t[1], n, knotVector ) );
    REQUIRE_NOTHROW( P = deBoor(t[1], s, p, knotVector, x, y ) );

    CHECK( P[0] == Approx(9.3419010416666683 ) );
    CHECK( P[1] == Approx(2.6049366319444447 ) );

    // For t3, so at the end of the curve
    REQUIRE_NOTHROW( s = findKnotSpan(t[2], n, knotVector ) );
    REQUIRE_NOTHROW( P = deBoor(t[2], s, p, knotVector, x, y ) );

    CHECK( P[0] == Approx( x.back( ) ) );
    CHECK( P[1] == Approx( y.back( ) ) );
}

TEST_CASE( "DeBoorCurve_test" )
{
    // This is the same test as in the curve test because given the same setup,
    // DeBoor should deliver the same points on a curve as using the classic way: N * P
    std::vector<double> knotVector{ 0.0, 0.0, 0.5, 1.0, 1.0 };
    std::vector<double> x{ 2.0, 3.0, 0.5 };
    std::vector<double> y{ 1.0, 3.0, 3.0 };

    std::vector<double> t{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    std::array<std::vector<double>, 2> C;

    REQUIRE_NOTHROW( C = evaluate2DCurveDeBoor( t, x, y, knotVector ) );

    REQUIRE( C[0].size() == t.size( ) );
    REQUIRE( C[1].size() == t.size( ) );

    // x-coordinates of curve
    CHECK(C[0][0] == Approx( 2.0 ) );
    CHECK(C[0][1] == Approx( 2.2 ) );
    CHECK(C[0][2] == Approx( 2.4 ) );
    CHECK(C[0][3] == Approx( 2.6 ) );
    CHECK(C[0][4] == Approx( 2.8 ) );
    CHECK(C[0][5] == Approx( 3.0 ) );
    CHECK(C[0][6] == Approx( 2.5 ) );
    CHECK(C[0][7] == Approx( 2.0 ) );
    CHECK(C[0][8] == Approx( 1.5 ) );
    CHECK(C[0][9] == Approx( 1.0 ) );
    CHECK(C[0][10] == Approx( 0.5 ) );

    // y-coordinates of curve
    CHECK(C[1][0] == Approx( 1.0 ) );
    CHECK(C[1][1] == Approx( 1.4 ) );
    CHECK(C[1][2] == Approx( 1.8 ) );
    CHECK(C[1][3] == Approx( 2.2 ) );
    CHECK(C[1][4] == Approx( 2.6 ) );
    CHECK(C[1][5] == Approx( 3.0 ) );
    CHECK(C[1][6] == Approx( 3.0 ) );
    CHECK(C[1][7] == Approx( 3.0 ) );
    CHECK(C[1][8] == Approx( 3.0 ) );
    CHECK(C[1][9] == Approx( 3.0 ) );
    CHECK(C[1][10] == Approx( 3.0 ) );
}

TEST_CASE( "DeBoorCurveScript_test" )
{
    std::vector<double> knotVector { 0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 9.0, 9.0, 9.0, 9.0 };
    std::vector<double> x { 0.0, 10.0, 9.0, 4.5, 1.5, 1.0};
    std::vector<double> y { 0.0,  1.0, 4.0, 7.5, 6.0, 1.0};
    std::vector<double> t { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

    std::array<std::vector<double>, 2> C;

    REQUIRE_NOTHROW( C = evaluate2DCurveDeBoor( t, x, y, knotVector ) );

    REQUIRE( C[0].size( ) == t.size( ) );
    REQUIRE( C[1].size( ) == t.size( ) );

    // x-coordinates of curve
    CHECK( C[0][0] == Approx( 0.0 ) );
    CHECK( C[0][1] == Approx( 9.4375 ) );
    CHECK( C[0][2] == Approx( 8.3385 ) );
    CHECK( C[0][3] == Approx( 7.0208 ) );
    CHECK( C[0][4] == Approx( 5.6406 ) );
    CHECK( C[0][5] == Approx( 4.336 ) );
    CHECK( C[0][6] == Approx( 3.1724 ) );
    CHECK( C[0][7] == Approx( 2.197 ) );
    CHECK( C[0][8] == Approx( 1.45713 ) );
    CHECK( C[0][9] == Approx( 1.0 ) );

    CHECK( C[1][0] == Approx( 0.0 ) );
    CHECK( C[1][1] == Approx( 2.40972 ) );
    CHECK( C[1][2] == Approx( 4.12413 ) );
    CHECK( C[1][3] == Approx( 5.33333 ) );
    CHECK( C[1][4] == Approx( 6.07378 ) );
    CHECK( C[1][5] == Approx( 6.35778 ) );
    CHECK( C[1][6] == Approx( 6.10094 ) );
    CHECK( C[1][7] == Approx( 5.19472 ) );
    CHECK( C[1][8] == Approx( 3.53059 ) );
    CHECK( C[1][9] == Approx( 1.0 ) );
}

} //namespace splinekernel
} // namespace cie
