#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <rational/rational.hpp>

TEST_CASE("Rational", "[rational]")
{
    using namespace rational;
    Rational r;

    r = Rational(0.5);
    CHECK(r.numerator() == 1);
    CHECK(r.denominator() == 2);
    CHECK(r == 0.5);
    CHECK(double(r) == 0.5);

    r = Rational("1", "2");
    CHECK(r.numerator() == 1);
    CHECK(r.denominator() == 2);
    CHECK(r == 0.5);
    CHECK(double(r) == 0.5);

    r = Rational("1", "3");
    CHECK(r.numerator() == 1);
    CHECK(r.denominator() == 3);
    CHECK(r != 1.0 / 3.0);         // 1/3 is not exactly representable in binary
    CHECK(double(r) == 1.0 / 3.0); // match when converted to the closest double

    r = Rational("-1670137647214829", "36893488147419103232");
    CHECK(std::string(r) == "-1670137647214829/36893488147419103232");
    CHECK(r == -1670137647214829.0 / 36893488147419103232.0);

    r = Rational("1", "0");
    CHECK(r.numerator() == 1);
    CHECK(r.denominator() == 0);
    CHECK(std::string(r) == "1/0");
}