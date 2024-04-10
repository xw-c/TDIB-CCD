// Modified version of rational-gmp.cpp from CVC3
// (https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/rational-gmp_8cpp-source.html)
// by Sergey Berezin. Licensed under the following license:
// https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/LICENSE.html

#pragma once

#include <gmp.h>

#include <string>
#include <ostream>

namespace rational {

/// @brief Helper class that wraps GMP's mpq_t (rational) type
class Rational {
public:
    Rational();
    Rational(float d);
    Rational(double d);
    Rational(int i);
    Rational(long i);
    Rational(const mpq_t& v_);
    Rational(const Rational& other);
    Rational(const std::string& numerator, const std::string& denominator);
    Rational(const std::string& str);

    ~Rational();

    friend Rational operator-(const Rational& v);
    friend Rational operator+(const Rational& x, const Rational& y);
    friend Rational operator-(const Rational& x, const Rational& y);
    friend Rational operator*(const Rational& x, const Rational& y);
    friend Rational operator/(const Rational& x, const Rational& y);

    Rational& operator=(const Rational& x);
    Rational& operator=(const double x);
    Rational& operator=(const float x);
    Rational& operator=(const int x);
    Rational& operator=(const long x);

    template <typename T> bool operator<(T r1);
    template <typename T> bool operator>(T r1);
    template <typename T> bool operator<=(T r1);
    template <typename T> bool operator>=(T r1);
    template <typename T> bool operator==(T r1);
    template <typename T> bool operator!=(T r1);

    operator double() const;
    operator std::string() const;

    friend std::ostream& operator<<(std::ostream& os, const Rational& r);

    void canonicalize();

    int sign() const;

    long long numerator() const;
    long long denominator() const;

    std::string numerator_str() const;
    std::string denominator_str() const;

    void print_numerator() const;
    void print_denominator() const;

    /// @brief The actual GMP rational value
    mpq_t value;

protected:
    /// @brief Helper class that wraps GMP's mpz_t (integer) type
    class Integer {
    public:
        Integer();

        ~Integer();

        static Integer from_numerator(const Rational& r);
        static Integer from_denominator(const Rational& r);

        void print() const;

        operator long long() const;
        operator std::string() const;

        /// @brief The actual GMP integer value
        mpz_t value;
    };
};

} // namespace rational
