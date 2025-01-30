#include <stdexcept>

#ifndef RATIONAL
#define RATIONAL

class RationalDivisionByZero : public std::runtime_error {
 public:
  RationalDivisionByZero() : std::runtime_error("RationalDivisionByZero") {
  }
};

class Rational {
  int a_;
  int b_;
  int Gcd(int a, int b);

 public:
  void NormalView();
  Rational(const Rational&) = default;
  Rational();
  Rational(int a);  // NOLINT
  Rational(int a, int b);
  int GetNumerator() const;
  int GetDenominator() const;
  void SetNumerator(int a);
  void SetDenominator(int b);
  Rational& operator=(Rational rational);
};

Rational operator+(const Rational&, const Rational&);
Rational operator-(const Rational&, const Rational&);
Rational operator*(const Rational&, const Rational&);
Rational operator/(const Rational&, const Rational&);
Rational& operator+=(Rational&, const Rational&);
Rational& operator-=(Rational&, const Rational&);
Rational& operator*=(Rational&, const Rational&);
Rational& operator/=(Rational&, const Rational&);
Rational operator+(const Rational&);
Rational operator-(const Rational&);
Rational& operator++(Rational&);
Rational& operator--(Rational&);
Rational operator++(Rational&, int);
Rational operator--(Rational&, int);
bool operator<=(const Rational, const Rational);
bool operator>=(const Rational, const Rational);
bool operator<(const Rational, const Rational);
bool operator>(const Rational, const Rational);
bool operator==(const Rational, const Rational);
bool operator!=(const Rational, const Rational);

std::istream& operator>>(std::istream&, Rational&);
std::ostream& operator<<(std::ostream&, const Rational&);

#endif
