
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

#define pi 3.14159265358979323846264
#define fourpi 12.566370614359172953850
#define sqrt4pi  3.544907701811032 // iem

constexpr float decodeCorrection(const int N) { return sqrt4pi / (N + 1) / (N + 1); } // iem

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
    float lf = (2.f * (int)order) + 1;
    float t = factorial((int)order - abs(degree));
    float b = factorial((int)order + abs(degree));
    float res = (lf / fourpi) * (t / b);
    //std::cout << '\n' << "\tLeft fraction: " << lf << " Top factorial: " << t << " Bottom factorial: " << b << " Result: " << (sqrtf(res)) << '\n';
    return (sqrtf(res));
}

float danielNorm(unsigned order, int degree)
{
    int d = (degree == 0) ? 1 : 0;  // Kronecker delta
    float ratio = static_cast<float>(factorial(order - abs(degree))) / factorial(order + abs(degree));
    return sqrtf((2.f - d) * ratio);
}

float norme(int order, int degree, bool n3d)
{
    int d = (degree == 0) ? 1 : 0;  // Kronecker delta
    float ratio = static_cast<float>(factorial(order - abs(degree))) / factorial(order + abs(degree));
    return n3d ? sqrtf(2 * order + 1) * sqrtf((2.f - d) / (4.f * pi) * ratio) : sqrtf((2.f - d) / (4.f * pi) * ratio); // n3d = TRUE? then N3D (equal weighting),,,, n3d = FALSE? then SN3D
}

float inversenorme(unsigned order, int degree)
{
    float lf = (2.f * (int)order) + 1;
    float t = factorial((int)order - abs(degree));
    float b = factorial((int)order + abs(degree));
    float res = (lf / fourpi) * (t / b);
    //std::cout << "Order: " << order << " Degree: " << degree << " Inverse Norm: " << (sqrtf(res)) << std::endl;
    return (1.f / (sqrtf(res)));
}

std::vector<float> SH(unsigned order_, const float azimuth_, const float zenith_)
{
    float azimuth = degree2rad(azimuth_);
    float zenith = degree2rad(zenith_);
    std::vector<float> result = std::vector<float>(pow(((int)order_ + 1), 2), 0); // instantiate and reserve a vector that is the size of the results that shall be returned
    for (int order = 0; order <= (int)order_; order++)
    {
        if (order == 0) result[0] = danielNorm(order, 0);
        for (int i = -order; i <= order; i++)
        {
            float n = danielNorm(order, i);
            float p = (std::assoc_legendref(abs(i), order, sinf(zenith)));
            float r = 0.f;
            if (i < 0) r = sinf(abs(i) * (azimuth)); 
            else if (i >= 0) r = cosf(i * (azimuth)); 
            // std::cout << "Order: " << order << " Degree: " << i << " Complex component: " << r << " Legendre: " << p << " Normalization term: " << n << " ACN: " << (pow(order, 2) + order + i) << '\n';
            result[pow(order, 2) + order + i] = sqrtf((2 * order) + 1) * (n * p) * r; // place inside vector so it is ordered as Y^0_0, Y^1_-1, Y^1_0, Y^1_1
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
        int order = 0;
        while ((order + 1) * (order + 1) <= ACN) 
        {
            order++;
        }
        int degree = ACN - (order * order) - order;
        //std::cout << "Order: " << order << " Degree: " << degree << " Encode factor: " << inversenorme(order, degree) << '\n';
        result[ACN] = SH[ACN];
    }
    return result;
}


int main() 
{
    float azi = 90.f;
    float zeni = 0.f;
    int order = 2;
    int degree;
    /*
    for (order = 0; order <= 3; order++)
    {
        for(degree = -order; degree <= order; degree++)
        {
            std::cout << "Order: " << order << " Degree: " << degree << " Norme: " << norme(order, degree) << '\n';
        }
    }
    */
    std::vector<float> spheric = SH(order, azi, zeni);
    int m_order = 0;
    int m_degree = 0;
    for (int i = 0; i < spheric.capacity(); i++)
    {
        m_order = (int)(floor(sqrt(i)));
        m_degree = (int)(i - (m_order * m_order) - m_order);
        std::cout << "ACN: " << i << " Order: " << m_order << " Degree: " << m_degree << " Normalization: " << danielNorm(m_order, m_degree) << '\n';
        std::cout << "Raw: " << spheric[i] << '\n';
    }
    std::cout << std::endl;
    /**
    std::vector<float> bformat = BFormSH(spheric);
    for (int i = 0; i < bformat.capacity(); i++)
    {
        int order = 0;
        while ((order + 1) * (order + 1) <= i) { order++; }
        int degree = i - (order * order) - order;
        std::cout << "ACN: " << i << " B-Format: " << bformat[i] << '\n';
    }
    */
    return 0;
}