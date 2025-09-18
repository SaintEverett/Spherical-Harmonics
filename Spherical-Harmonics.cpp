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

constexpr int factorial(int n)
{
    return n <= 1 ? 1 : (n * factorial(n - 1));
}

float norme(unsigned order, int degree)
{
    float l_frac = ((2.f * order + 1) / (4.f * pi));
    float t_r = factorial(order - abs(degree));
    float b_r = factorial(order + abs(degree));
    return sqrtf(l_frac * (t_r / b_r));
}

float inversenorme(unsigned order, int degree)
{
    float t_r = factorial(order - abs(degree));
    float b_r = factorial(order + abs(degree));
    return sqrt(4*pi) * sqrtf(b_r / t_r);
}

std::vector<float> SH(unsigned order_, const float azimuth_, const float zenith_)
{
    float azimuth = degree2rad(azimuth_);
    float zenith = degree2rad(zenith_);
    std::vector<float> result = std::vector<float>((int)pow((order_ + 1), 2), 0); // instantiate and reserve a vector that is the size of the results that shall be returned
    for (int order = 0; order <= order_; order++)
    {
        if (order == 0) result[0] = norme(order, 0);
        for (int i = -order; i <= order; i++)
        {
            float n = norme(order, i);
            float p = (std::assoc_legendref(order, abs(i), cosf((pi / 2.f) - zenith))); // subtract 90 degrees to rotate
            float r = 0.f;
            if (i < 0) r = sinf(i * ((pi / 2.f) - azimuth)); // subtract 90 degrees to rotate
            else if (i >= 0) r = cosf(i * ((pi / 2.f) - azimuth)); // subtract 90 degrees to rotate
            std::cout << "Order: " << order << " Degree: " << i << " Complex component: " << r << " Legendre: " << p << " Normalization term: " << n << " ACN: " << (pow(order, 2) + order + i) << '\n';
            result[pow(order, 2) + order + i] = (n * p * r); // place inside vector so it is ordered as Y^0_0, Y^1_-1, Y^1_0, Y^1_1
        }
    }
    return result;
}

std::vector<float> BFormSH(std::vector<float> SH)
{
    std::vector<float> result = std::vector<float>(SH.capacity(), 0);
    unsigned order = 0;
    int degree;
    for (int ACN = 0; ACN < SH.capacity(); ACN++)
    {
        order = floor(sqrtf(ACN));
        degree = ACN - pow(order, 2) - order;
        std::cout << "Order: " << order << " Degree: " << degree << " Encode factor: " << inversenorme(order, degree) << '\n';
        result[ACN] = SH[ACN] * inversenorme(order, degree);
    }
    return result;
}


int main() 
{
    float azi = 45.f;
    float zeni = 26.f;
    unsigned order = 4;
    std::vector<float> spheric = SH(order, azi, zeni);
    for (int i = 0; i < spheric.capacity(); i++)
    {
        std::cout << "Raw: " << spheric[i] << '\n';
    }
    std::cout << std::endl;
    std::vector<float> bformat = BFormSH(spheric);
    for (int i = 0; i < bformat.capacity(); i++)
    {
        std::cout << "B-Format: " << bformat[i] << '\n';
    }
    return 0;
}