//-----------------------------------------------------------------------------------------------------------------------------//
//                                  Everett M. Carpenter's spherical harmonic library for C++                                  //
//                                                                                                                             //
//                              ~{-------------------------------------------------------------}~                              //
//                               |           __..--''``---....___   _..._    __                |                               //
//                               |    /// //_.-'    .-/";  `        ``<._  ``.''_ `. / // /    |                               //
//                               |    ///_.-' _..--.'_    \                    `( ) ) // //    |                               //
//                               |    / (_..-' // (< _     ;_..__               ; `' / ///     |                               //
//                               |    / // // //  `-._,_)' // / ``--...____..-' /// / //       |                               //
//                              ~{-------------------------------------------------------------}~                              //
//                                                                                                                             //
//-----------------------------------------------------------------------------------------------------------------------------//

#include <iostream>
#include <cmath>
#include <vector>
#include <array>

#define pi 3.14159265358979323846264
#define fourpi 12.566370614359172953850
#define MAX_ORDER 9
#define MAX_SIZE (MAX_ORDER + MAX_ORDER)

static const float degree2rad = pi / 180.0f;
static const float rad2deg = 180.0f / pi;

// constexpr function to generate factorial table at compile time
constexpr auto create_factorial_lut()
{
    std::array<unsigned long long, MAX_SIZE> lut{}; // create look up table of huge values (unsigned long long)
    lut[0] = 1;                                     // 0! = 1
    for (size_t i = 1; i < MAX_SIZE; i++)
    {
        lut[i] = lut[i - 1] * i; // clever way of calculating without having to do 1*2*3*4....*n
    }
    return lut;
}

constexpr auto FACTORIAL_LOOKUP = create_factorial_lut(); // actually make the lookup table

unsigned long long factorial(size_t n) // access this lookup table
{
    if (n > MAX_SIZE)
    {
        return 0;
        std::cout << "ERROR: FACTORIAL INDEX NOT AVAILABLE \n";
    }
    else
        return FACTORIAL_LOOKUP[n]; // since lut[0] = 0!, lut[1] = 1!, lut[2] = 2!, your n and index are interchangable
}

float SN3D(unsigned order, int degree) // SN3D normalization, returns [-1,1] normalized values
{
    int d = (degree == 0) ? 1 : 0;                                                                     // Kronecker delta
    float ratio = static_cast<float>(factorial(order - abs(degree))) / factorial(order + abs(degree)); // ratio of factorials
    return sqrtf((2.f - d) * ratio);                                                                   // final result, 1/4pi omitted
}

float N3D(unsigned order, int degree) // N3D normalization
{
    int d = (degree == 0) ? 1 : 0;                                                                     // Kronecker delta
    float ratio = static_cast<float>(factorial(order - abs(degree))) / factorial(order + abs(degree)); // ratio of factorials
    return sqrtf((2 * order) + 1) * sqrtf((2.f - d) * ratio);                                          // final result, N3D factor included
}

std::vector<float> SH(unsigned order_, const float azimuth_, const float zenith_, bool n3d)
{
    float azimuth_shift = (azimuth_ - 90.f) * degree2rad;      // shift "perspective" so that azi = 0 and zeni = 0 is a unity vector facing outwards from the listener (vector pointing from roughly the nose forward)
    float zenith_shift = (zenith_ - 90.f) * degree2rad;        // same here
    float coszeni = cosf(zenith_shift);                        // pre calculate cos(zenith)
    int size = (order_ + 1) * (order_ + 1);                    // pre-compute size of vector to be returned
    std::vector<float> result = std::vector<float>(size, 0.f); // instantiate vector that is the size of the results that shall be returned
    for (int order = 0; order <= (int)order_; order++)         // all orders from 0 - desired order
    {
        if (order == 0)
            result[0] = SN3D(order, 0);                      // Y^0_0 is omnidirectional
        for (int degree = -order; degree <= order; degree++) // all degrees of current order
        {
            float n = n3d ? N3D(SN3D(order, degree), order) : SN3D(order, degree);                         // normalization term if n3d bool = TRUE, return N3D else SN3D
            float p = (std::assoc_legendref(order, abs(degree), coszeni));                                 // legendre NOTE: degree of legendre is current ambisonic order & order of legendre is current ambisonic degree (very frustrating)
            float r = (degree < 0) ? sinf(abs(degree) * (azimuth_shift)) : cosf(degree * (azimuth_shift)); // degree positive? Re(exp(i*azimuth*degree)) degree negative? Im(exp(i*azimuth*degree))
            result[(order * order) + order + degree] = n * p * r;                                          // place inside vector so it is ordered as Y^0_0, Y^1_-1, Y^1_0, Y^1_1
        }
    }
    return result;
}

void print(std::vector<float> shs) // print function
{
    int m_order = 0;
    int m_degree = 0;
    for (int i = 0; i < shs.capacity(); i++)
    {
        m_order = (int)(floor(sqrt(i)));
        m_degree = (int)(i - (m_order * m_order) - m_order);
        // std::cout << "ACN: " << i << " Order: " << m_order << " Degree: " << m_degree << " Value: " << shs[i] << '\n';
    }
}
/*
int main()
{
    float azi = 0.f;                                      // test instance azimuth
    float zeni = 0.f;                                     // test instance zenith
    int order = 5;                                        // test instance order
    std::vector<float> spheric = SH(order, azi, zeni, 0); // as simple as that!
    print(spheric);
    return 0;
}
*/
