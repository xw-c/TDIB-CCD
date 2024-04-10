// Modified version of rational-gmp.cpp from CVC3
// (https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/rational-gmp_8cpp-source.html)
// by Sergey Berezin. Licensed under the following license:
// https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/LICENSE.html

#include "rational.hpp"

#include <cassert>
#include <iostream>
#include <cmath>

namespace rational {
Rational::Rational()
{
    mpq_init(value);
    mpq_set_d(value, 0);
}

Rational::Rational(double d)
{
    assert(std::isfinite(d));
    mpq_init(value);
    mpq_set_d(value, d);
    canonicalize();
}

Rational::Rational(float d)
{
    assert(std::isfinite(d));
    mpq_init(value);
    double ddouble = d; // convert (float)d to double
    mpq_set_d(value, ddouble);
    canonicalize();
}

Rational::Rational(int i)
{
    mpq_init(value);
    mpq_set_si(value, i, 1U);
    canonicalize();
}

Rational::Rational(long i)
{
    mpq_init(value);
    mpq_set_si(value, i, 1U);
    canonicalize();
}

Rational::Rational(const mpq_t& v_)
{
    mpq_init(value);
    mpq_set(value, v_);
    // canonicalize();
}

Rational::Rational(const Rational& other)
{
    mpq_init(value);
    mpq_set(value, other.value);
}

Rational::Rational(const std::string& numerator, const std::string& denominator)
    : Rational(numerator + "/" + denominator)
{
}

Rational::Rational(const std::string& str)
{
    mpq_init(value);
    // The return value is 0 if the entire string is a valid number.
    const int r = mpq_set_str(value, str.c_str(), 10);
    if (r != 0) {
        throw std::runtime_error("Invalid rational: {}" + str);
    }
}

Rational::~Rational() { mpq_clear(value); }

Rational operator-(const Rational& v)
{
    Rational r_out;
    mpq_neg(r_out.value, v.value);
    return r_out;
}

Rational operator+(const Rational& x, const Rational& y)
{
    Rational r_out;
    mpq_add(r_out.value, x.value, y.value);
    return r_out;
}

Rational operator-(const Rational& x, const Rational& y)
{
    Rational r_out;
    mpq_sub(r_out.value, x.value, y.value);
    return r_out;
}

Rational operator*(const Rational& x, const Rational& y)
{
    Rational r_out;
    mpq_mul(r_out.value, x.value, y.value);
    return r_out;
}

Rational operator/(const Rational& x, const Rational& y)
{
    Rational r_out;
    mpq_div(r_out.value, x.value, y.value);
    return r_out;
}

Rational& Rational::operator=(const Rational& x)
{
    if (this == &x)
        return *this;
    mpq_set(value, x.value);
    return *this;
}

Rational& Rational::operator=(const double x)
{
    mpq_set_d(value, x);
    // canonicalize();
    return *this;
}

Rational& Rational::operator=(const float x)
{
    double xd = x;
    mpq_set_d(value, xd);
    // canonicalize();
    return *this;
}

Rational& Rational::operator=(const int x)
{
    mpq_set_si(value, x, 1U);
    // canonicalize();
    return *this;
}

Rational& Rational::operator=(const long x)
{
    mpq_set_si(value, x, 1U);
    // canonicalize();
    return *this;
}

template <typename T> bool Rational::operator<(T r1)
{
    if constexpr (std::is_same<T, Rational>::value)
        return mpq_cmp(value, r1.value) < 0;
    else
        return *this < Rational(r1);
}

template <typename T> bool Rational::operator>(T r1)
{
    if constexpr (std::is_same<T, Rational>::value)
        return mpq_cmp(value, r1.value) > 0;
    else
        return *this > Rational(r1);
}

template <typename T> bool Rational::operator<=(T r1)
{
    if constexpr (std::is_same<T, Rational>::value)
        return mpq_cmp(value, r1.value) <= 0;
    else
        return *this <= Rational(r1);
}

template <typename T> bool Rational::operator>=(T r1)
{
    if constexpr (std::is_same<T, Rational>::value)
        return mpq_cmp(value, r1.value) >= 0;
    else
        return *this >= Rational(r1);
}

template <typename T> bool Rational::operator==(T r1)
{
    if constexpr (std::is_same<T, Rational>::value)
        return mpq_equal(value, r1.value);
    else
        return *this == Rational(r1);
}

template <typename T> bool Rational::operator!=(T r1)
{
    if constexpr (std::is_same<T, Rational>::value)
        return !mpq_equal(value, r1.value);
    else
        return *this != Rational(r1);
}

void Rational::canonicalize() { mpq_canonicalize(value); }

int Rational::sign() const { return mpq_sgn(value); }

// conversion operators

Rational::operator double() const { return mpq_get_d(value); }

Rational::operator std::string() const
{
    return std::string(mpq_get_str(NULL, 10, value));
}

// output operator

std::ostream& operator<<(std::ostream& os, const Rational& r)
{
    return (os << std::string(r));
}

// numerator and denominator

long long Rational::numerator() const { return Integer::from_numerator(*this); }

long long Rational::denominator() const
{
    return Integer::from_denominator(*this);
}

std::string Rational::numerator_str() const
{
    return Integer::from_numerator(*this);
}

std::string Rational::denominator_str() const
{
    return Integer::from_denominator(*this);
}

void Rational::print_numerator() const
{
    Integer::from_numerator(*this).print();
}

void Rational::print_denominator() const
{
    Integer::from_denominator(*this).print();
}

// --------------------- Integer ---------------------

Rational::Integer::Integer() { mpz_init(value); }

Rational::Integer::~Integer() { mpz_clear(value); }

Rational::Integer Rational::Integer::from_numerator(const Rational& r)
{
    Integer i;
    mpq_get_num(i.value, r.value);
    return i;
}

Rational::Integer Rational::Integer::from_denominator(const Rational& r)
{
    Integer i;
    mpq_get_den(i.value, r.value);
    return i;
}

void Rational::Integer::print() const { mpz_out_str(NULL, 10, value); }

Rational::Integer::operator long long() const
{
    assert(mpz_fits_slong_p(value));
    return mpz_get_si(value);
}

Rational::Integer::operator std::string() const
{
    return std::string(mpz_get_str(NULL, 10, value));
}

// ==================== Template instantiations ====================
template bool Rational::operator< <double>(double r1);
template bool Rational::operator< <float>(float r1);
template bool Rational::operator< <int>(int r1);
template bool Rational::operator< <long>(long r1);
template bool Rational::operator< <const mpq_t&>(const mpq_t& r1);
template bool Rational::operator< <const std::string&>(const std::string& r1);
template bool Rational::operator< <const Rational&>(const Rational& r1);

template bool Rational::operator><double>(double r1);
template bool Rational::operator><float>(float r1);
template bool Rational::operator><int>(int r1);
template bool Rational::operator><long>(long r1);
template bool Rational::operator><const mpq_t&>(const mpq_t& r1);
template bool Rational::operator><const std::string&>(const std::string& r1);
template bool Rational::operator><const Rational&>(const Rational& r1);

template bool Rational::operator<=<double>(double r1);
template bool Rational::operator<=<float>(float r1);
template bool Rational::operator<=<int>(int r1);
template bool Rational::operator<=<long>(long r1);
template bool Rational::operator<=<const mpq_t&>(const mpq_t& r1);
template bool Rational::operator<=<const std::string&>(const std::string& r1);
template bool Rational::operator<=<const Rational&>(const Rational& r1);

template bool Rational::operator>=<double>(double r1);
template bool Rational::operator>=<float>(float r1);
template bool Rational::operator>=<int>(int r1);
template bool Rational::operator>=<long>(long r1);
template bool Rational::operator>=<const mpq_t&>(const mpq_t& r1);
template bool Rational::operator>=<const std::string&>(const std::string& r1);
template bool Rational::operator>=<const Rational&>(const Rational& r1);

template bool Rational::operator==<double>(double r1);
template bool Rational::operator==<float>(float r1);
template bool Rational::operator==<int>(int r1);
template bool Rational::operator==<long>(long r1);
template bool Rational::operator==<const mpq_t&>(const mpq_t& r1);
template bool Rational::operator==<const std::string&>(const std::string& r1);
template bool Rational::operator==<const Rational&>(const Rational& r1);

template bool Rational::operator!=<double>(double r1);
template bool Rational::operator!=<float>(float r1);
template bool Rational::operator!=<int>(int r1);
template bool Rational::operator!=<long>(long r1);
template bool Rational::operator!=<const mpq_t&>(const mpq_t& r1);
template bool Rational::operator!=<const std::string&>(const std::string& r1);
template bool Rational::operator!=<const Rational&>(const Rational& r1);

} // namespace rational
