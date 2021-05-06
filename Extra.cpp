#include "Extra.h"

double getVectorLength(vector<double> v) {
    double sum {0};
    for(double i : v) sum += i*i;
    return sqrt(sum);
}

vector<double> unify(vector<double> v) {
    double length = getVectorLength(v);
    if(!length) return v;
    vector<double> newVec {};
    for(double i : v) newVec.push_back(i/length);
    return newVec;
}