// Online C++ compiler to run C++ program online
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

#define pi 3.14159265358979323846264

float rad2degree(float radians)
{
    return (radians * 180.f / pi);
}

float degree2rad(float degrees)
{
    return (degrees * pi / 180.f);
}

float norme(unsigned order, int degree)
{

    float l_frac = ((2.f * order - 1) / (4.f * pi));
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

std::vector<float> SH(unsigned order, const float azimuth_, const float zenith_)
{
    float azimuth = degree2rad(azimuth_);
    float zenith = degree2rad(zenith_);
    std::vector<float> result = std::vector<float>((int)pow((order + 1), 2), 0); // instantiate and reserve a vector that is the size of the results that shall be returned
    result[0] = 1.f;
    for (int i = -(int)order; i <= (int)order; i++)
    {
        float n = norme(order, i);
        float p = (std::assoc_legendref(order, abs(i), cosf((pi / 2.f) - zenith))); // subtract 90 degrees to rotate
        float r = 0.f;
        if (i < 0) r = sinf(i * ((pi / 2.f) - azimuth)); // subtract 90 degrees to rotate
        else if (i >= 0) r = cosf(i * ((pi / 2.f) - azimuth)); // subtract 90 degrees to rotate
        result[(i + order) + 1] = (n * p * r); // place inside vector so it is ordered as Y^0_0, Y^1_-1, Y^1_0, Y^1_1
    }
    return result;
}

int main() {
    float azi = 270.f;
    float zeni = 0.f;
    unsigned order = 1;
    std::vector<float> spheric = SH(order, azi, zeni);
    for (int i = 0; i < spheric.capacity(); i++)
    {
        std::cout << spheric[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}