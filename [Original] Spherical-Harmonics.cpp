#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

#define pi 3.14159265358979323846264
#define fourpi 12.566370614359172953850

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

float danielNorm(unsigned order, int degree)
{
    int d = (degree == 0) ? 1 : 0; // Kronecker delta
    float ratio = static_cast<float>(factorial(order - abs(degree))) / factorial(order + abs(degree));
    return sqrtf((2.f - d) * ratio);
}

std::vector<float> SH(unsigned order_, const float azimuth_, const float zenith_)
{
    float azimuth = degree2rad(azimuth_);
    float zenith = degree2rad(zenith_);
    std::vector<float> result = std::vector<float>(pow(((int)order_ + 1), 2), 0); // instantiate and reserve a vector that is the size of the results that shall be returned
    for (int order = 0; order <= (int)order_; order++)
    {
        if (order == 0)
            result[0] = danielNorm(order, 0);
        for (int i = -order; i <= order; i++)
        {
            float n = danielNorm(order, i);
            float p = (std::assoc_legendref(order, abs(i), cosf(zenith - pi / 2)));
            float r = 0.f;
            if (i < 0)
                r = sinf(abs(i) * (azimuth - pi / 2));
            else if (i >= 0)
                r = cosf(i * (azimuth - pi / 2));
            // std::cout << "Order: " << order << " Degree: " << i << " Complex component: " << r << " Legendre: " << p << " Normalization term: " << n << " ACN: " << (pow(order, 2) + order + i) << '\n';
            result[pow(order, 2) + order + i] = (n * p) * r; // place inside vector so it is ordered as Y^0_0, Y^1_-1, Y^1_0, Y^1_1
        }
    }
    return result;
}

int main()
{
    float azi = 0.f;
    float zeni = 0.f;
    int order = 12;
    int degree;
    std::vector<float> spheric = SH(order, azi, zeni);
    int m_order = 0;
    int m_degree = 0;
    for (int i = 0; i < spheric.capacity(); i++)
    {
        m_order = (int)(floor(sqrt(i)));
        m_degree = (int)(i - (m_order * m_order) - m_order);
        std::cout << "ACN: " << i << " Order: " << m_order << " Degree: " << m_degree << " Value: " << spheric[i] << '\n';
        // std::cout << "Raw: " << spheric[i] << '\n';
    }
    std::cout << std::endl;
    return 0;
}