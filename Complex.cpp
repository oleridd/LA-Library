#include "Complex.h"

#include <iostream>


Comp::Comp() : Re(0), Im(0) {}
Comp::Comp(double Re, double Im) : Re(Re), Im(Im) {}

bool Comp::isZero() {
    return (Im == 0 && Re == 0);
}

Comp Comp::getConjugate() {
    return Comp(Re, -Im);
}

ostream& operator<<(ostream& os, Comp& c) {
    if(c.Im >= 0) {
        os << c.Re << "+" << c.Im << "i";
    } else {
        os << c.Re << "-" << c.Im * -1 << "i";
    }
    return os;
}

Comp& Comp::operator+=(Comp& c) {
    Re += c.Re;
    Im += c.Im;
    return *this;
}
Comp& Comp::operator+(Comp& c) {
    Comp b = *this;
    return b += c;
}
Comp& Comp::operator-=(Comp& c) {
    Re -= c.Re;
    Im -= c.Im;
    return *this;
}
Comp& Comp::operator-(Comp& c) {
    Comp b = *this;
    return b -= c;
}
Comp& Comp::operator*=(Comp& c) {
    double oRe = Re;
    double oIm = Im;
    Re = oRe*c.Re - oIm*c.Im;
    Im = oRe*c.Im + oIm*c.Re;
    return *this;
}
Comp& Comp::operator*(Comp& c) {
    return *this *= c;
}
Comp& Comp::operator/(Comp& c) {
    if(c.isZero()) throw runtime_error("Denominator can not be 0");
    double oRe = Re;
    double oIm = Im;
    Re = (oRe*c.Re + oIm*c.Im)/(pow(c.Re, 2) + pow(c.Im, 2));
    Im = (oIm*c.Re - oRe*c.Im)/(pow(c.Re, 2) + pow(c.Im, 2));
    return *this;
}