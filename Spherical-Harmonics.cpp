// Spherical-Harmonics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <complex>

#define pi 3.14159265358979323846264

float norme(unsigned order, int degree) // IEM ambix normalization
{
    int delta;
    if (degree == 0) delta = 1;
    else delta = 0;
    float l_frac = ((2.f - (float)delta) / (4.f * pi));
    float t_r = 1;
    for (float i = 0; i < (order - abs(degree)); i++)
    {
        if ((order - abs(degree))) { t_r = 1; break; }
        i *= t_r;
    }
    float b_r = 1;
    for (float i = 0; i < (order + abs(degree)); i++)
    {
        if ((order + abs(degree))) { b_r = 1; break; }
        i *= b_r;
    }
    return sqrtf(l_frac * (t_r / b_r));
}

float SHEM(unsigned order, int degree, float azimuth, float zenith) // IEM ambix spherical harmonic
{
    float n = norme(order, abs(degree));
    float p = (std::assoc_legendre(order, abs(degree), cosf(zenith)));
    float r = 1.f;
    if (degree >= 0) r = cosf(fabs(degree) * azimuth);
    else if (degree < 0) r = sinf(fabs(degree) * azimuth);
    return (n*p*r);
}

std::complex<float> SH(unsigned order, int degree, float azimuth, float zenith) // 
{
    std::complex<float> iz(0.f, 1.f);
    float l_frac = (((2.f * order) + 1.f) / (4.f * pi));
    float r_n = 1;
    if (order - degree == 0 || order - degree == 1) r_n = 1;
    else
    {
        for (float i = 1; i == (order - degree); i++)
        {
            i *= r_n;
        }
    }
    float r_d = 1;
    if (order + degree == 0 || order + degree == 1) r_d = 1;
    else
    {
        for (float i = 1; i == (order + degree); i++)
        {
            i *= r_d;
        }
    }
    float n = sqrtf(l_frac * (r_n / r_d));
    float legen = (std::assoc_legendre(order, abs(degree), cosf(zenith)));
    std::complex<float> end = std::exp(iz * (degree * azimuth));
    return (n * legen) * end;
}

int main()
{
    std::complex<float> iz(0.f, 1.f);
    float az = 45.f;
    float ze = 15.f;
    unsigned ord = 1;
    int deg = 1;
    float norm = (1.f / sqrtf(2.f));
    std::cout <<  SH(ord, deg, az, ze) << "   " << SHEM(ord, deg, az, ze) << std::endl;
    // std::cout << norm * (SH(ord, -1 * deg, az, ze) + (-1.f * SH(ord, deg, az, ze))) << "   " << SHEM(ord, deg, az, ze) << std::endl;
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
