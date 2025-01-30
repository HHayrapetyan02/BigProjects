#include "rational.h"
#include <iostream>
#include <cmath>

Rational::Rational() : a_(0), b_(1) {
}

Rational::Rational(int a) : a_(a), b_(1) {
}

Rational::Rational(int a, int b) : a_(a), b_(b) {
  NormalView();
}

int Rational::Gcd(int a, int b) {
  if (a == 0 || b == 0) {
    return a + b;
  }
  if (a > b) {
    return Gcd(a - b, b);
  }
  return Gcd(a, b - a);
}

void Rational::NormalView() {
  if (b_ == 0) {
    throw RationalDivisionByZero();
  }
  if (b_ < 0) {
    a_ = -a_;
    b_ = -b_;
  }
  int nod = Gcd(std::abs(a_), b_);
  a_ /= nod;
  b_ /= nod;
}

int Rational::GetNumerator() const {
  return a_;
}

int Rational::GetDenominator() const {
  return b_;
}

void Rational::SetNumerator(int a) {
  a_ = a;
  NormalView();
}

void Rational::SetDenominator(int b) {
  b_ = b;
  NormalView();
}

Rational& Rational::operator=(Rational rational) {
  rational.NormalView();
  a_ = rational.a_;
  b_ = rational.b_;
  return *this;
}

Rational operator+(const Rational& first, const Rational& second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return {Rational(a_first * b_second + a_second * b_first, b_first * b_second)};
}

Rational operator-(const Rational& first, const Rational& second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return {Rational(a_first * b_second - a_second * b_first, b_first * b_second)};
}

Rational operator*(const Rational& first, const Rational& second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return {Rational(a_first * a_second, b_first * b_second)};
}

Rational operator/(const Rational& first, const Rational& second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return {Rational(a_first * b_second, b_first * a_second)};
}

Rational& operator+=(Rational& first, const Rational& second) {
  first = first + second;
  return first;
}

Rational& operator-=(Rational& first, const Rational& second) {
  first = first - second;
  return first;
}

Rational& operator*=(Rational& first, const Rational& second) {
  first = first * second;
  return first;
}

Rational& operator/=(Rational& first, const Rational& second) {
  first = first / second;
  return first;
}

Rational operator+(const Rational& rational) {
  int a = rational.GetNumerator();
  int b = rational.GetDenominator();
  return {Rational(a, b)};
}

Rational operator-(const Rational& rational) {
  int a = -rational.GetNumerator();
  int b = rational.GetDenominator();
  return {Rational(a, b)};
}

Rational& operator++(Rational& rational) {
  int a = rational.GetNumerator();
  int b = rational.GetDenominator();
  rational.SetNumerator(a + b);
  rational.NormalView();
  return rational;
}

Rational& operator--(Rational& rational) {
  int a = rational.GetNumerator();
  int b = rational.GetDenominator();
  rational.SetNumerator(a - b);
  rational.NormalView();
  return rational;
}

Rational operator++(Rational& rational, int) {
  Rational old = rational;
  int a = rational.GetNumerator();
  int b = rational.GetDenominator();
  rational.SetNumerator(a + b);
  return old;
}

Rational operator--(Rational& rational, int) {
  Rational old = rational;
  int a = rational.GetNumerator();
  int b = rational.GetDenominator();
  rational.SetNumerator(a - b);
  return old;
}

bool operator<=(const Rational first, const Rational second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return a_first * b_second <= a_second * b_first;
}

bool operator==(const Rational first, const Rational second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return a_first * b_second == a_second * b_first;
}

bool operator>=(const Rational first, const Rational second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return a_first * b_second >= a_second * b_first;
}

bool operator>(const Rational first, const Rational second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return a_first * b_second > a_second * b_first;
}

bool operator<(const Rational first, const Rational second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return a_first * b_second < a_second * b_first;
}

bool operator!=(const Rational first, const Rational second) {
  int a_first = first.GetNumerator();
  int b_first = first.GetDenominator();
  int a_second = second.GetNumerator();
  int b_second = second.GetDenominator();
  return a_first * b_second != a_second * b_first;
}

std::istream& operator>>(std::istream& is, Rational& rational) {
  int32_t a = 0;
  int32_t b = 0;
  char c = ' ';
  is >> a;
  if (is.peek() == '/') {
    is >> c;
    is >> b;
    rational = Rational(a, b);
  } else {
    rational = Rational(a);
  }
  return is;
}

std::ostream& operator<<(std::ostream& os, const Rational& rational) {
  if (rational.GetDenominator() != 1) {
    os << rational.GetNumerator() << '/' << rational.GetDenominator();
  } else {
    os << rational.GetNumerator();
  }
  return os;
}
