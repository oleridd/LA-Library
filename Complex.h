#pragma once
#include <iostream>

using namespace std;

class Comp {
    private:
        double Re;
        double Im;
    public:
        Comp();
        Comp(double Re, double Im);

        bool isZero();
        Comp getConjugate();

        friend ostream& operator<<(ostream& os, Comp& c);
        Comp& operator+=(Comp& c);
        Comp& operator+(Comp& c);
        Comp& operator-=(Comp& c);
        Comp& operator-(Comp& c);
        Comp& operator*=(Comp& c);
        Comp& operator*(Comp& c);
        Comp& operator/(Comp& c);
};