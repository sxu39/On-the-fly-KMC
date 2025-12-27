#ifndef KMC_UTILS_H
#define KMC_UTILS_H

#include <vector>
#include <iostream>

namespace KMC
{
    template <class PREC_T>
    static std::vector<PREC_T> operator+(const std::vector<PREC_T> &vec_a, const std::vector<PREC_T> &vec_b)
    {
        try
        {
            if (vec_a.size()!= vec_b.size())
                throw std::invalid_argument("Cannot add vectors of different sizes");
            std::vector<PREC_T> vec_res(vec_a.size());
            for (size_t i = 0 ; i < vec_a.size(); ++i)
            {
                vec_res[i] = vec_a[i] + vec_b[i];
            }
            return vec_res;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << "\n";
        }
    }

    template <class PREC_T>
    static std::vector<PREC_T> operator-(const std::vector<PREC_T> &vec_a, const std::vector<PREC_T> &vec_b)
    {
        try
        {
            if (vec_a.size()!= vec_b.size())
                throw std::invalid_argument("Cannot subtract vectors of different sizes");
            std::vector<PREC_T> vec_res(vec_a.size());
            for (size_t i = 0 ; i < vec_a.size(); ++i)
            {
                vec_res[i] = vec_a[i] - vec_b[i];
            }
            return vec_res;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << "\n";
        }
    }

    template <class PREC_T>
    static std::vector<PREC_T> fabs(const std::vector<PREC_T> &vec)
    {
        std::vector<PREC_T> vec_fabs(vec.size());
        for (size_t i = 0 ; i < vec.size() ; ++i)
        vec_fabs[i] = std::fabs(vec[i]);
        return vec_fabs;
    }

    template <class PREC_T>
    static std::vector<PREC_T> &&fabs(std::vector<PREC_T> &&vec)
    {
        for (auto &ele : vec)
            ele = std::fabs(ele);
        return std::move(vec);
    }

    template <class PREC_T>
    static PREC_T max_vec(const std::vector<PREC_T> &vec)
    {
        PREC_T maximum = INT64_MIN;
        for (int i = 0 ; i < vec.size() ; ++i)
            maximum = std::max(maximum, std::fabs(vec[i]));
        return maximum;
    }

    template <class PREC_T>
    static PREC_T min_vec(const std::vector<PREC_T> &vec)
    {
        PREC_T minimum = INT64_MAX;
        for (int i = 0 ; i < vec.size() ; ++i)
            minimum = std::min(minimum, std::fabs(vec[i]));
        return minimum;
    }
}

#endif