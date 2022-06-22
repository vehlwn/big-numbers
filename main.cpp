#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <iterator>
#include <vector>

using BigInt_t = boost::multiprecision::cpp_int;
using BigRational_t = boost::multiprecision::cpp_rational;
using BigCplx_t = boost::multiprecision::
    cpp_complex<101, boost::multiprecision::digit_base_10, std::allocator<void>>;
using BigFloat_t =
    boost::multiprecision::number<boost::multiprecision::cpp_bin_float<
        101,
        boost::multiprecision::digit_base_10,
        std::allocator<void>>>;

std::vector<BigCplx_t>
    slowDft(std::vector<BigCplx_t> input, const bool inverse = false)
{
    const std::size_t N = input.size();
    std::vector<BigCplx_t> ret(N);
    for(std::size_t k = 0; k < N; k++)
        for(std::size_t n = 0; n < N; n++)
            ret[k] += input[n]
                * exp(BigCplx_t{
                    0,
                    (inverse ? +1 : -1) * 2
                        * boost::math::constants::pi<BigCplx_t::value_type>() * k * n
                        / N})
                / sqrt(BigCplx_t::value_type{N});
    return ret;
}

BigCplx_t::value_type
    LinfMetric(const std::vector<BigCplx_t>& a, const std::vector<BigCplx_t>& b)
{
    std::vector<BigCplx_t::value_type> tmp;
    std::transform(
        a.begin(),
        a.end(),
        b.begin(),
        std::back_inserter(tmp),
        [](auto a, auto b) { return abs(a - b); });
    return *std::max_element(tmp.begin(), tmp.end());
}

int main()
{
    BigInt_t n = 1;
    for(int i = 1; i <= 1000; i++)
        n *= i;
    std::cout << n << std::endl;
    BigRational_t r = n;
    r /= n * 1001;
    std::cout << r << std::endl;
    std::cout.precision(101);
    std::cout << boost::math::constants::pi<BigFloat_t>() << std::endl;
    std::vector<BigCplx_t> v{0, 1, 2, 3, 4, 5};
    std::vector<BigCplx_t> v2 = slowDft(v, false);
    std::vector<BigCplx_t> v3 = slowDft(v2, true);
    std::cout << LinfMetric(v, v3) << std::endl;
}
