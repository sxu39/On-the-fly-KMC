#ifndef KMC_TIMER_H
#define KMC_TIMER_H

#include <time.h>
#include <string>
#include <unordered_map>

namespace KMC
{
    class timer
    {
        private:
            tm *time_p;
            std::unordered_map<int, std::string> month = {{0, "Jan"}, {1, "Feb"}, {2, "Mär"}, {3, "Apr"}, {4, "Mai"}, {5, "Jun"}, {6, "Jnl"}, {7, "Aug"},
            {8, "Sep"}, {9, "Okt"}, {10, "Nov"}, {11, "Dez"}};
        private:
            static inline const std::string two_digits(int);
        public:
            timer() = default;
            ~timer() = default;
            inline const std::string current_time();
    };
}

const std::string KMC::timer::two_digits(int time)
{
    return time < 10 ? "0" + std::to_string(time) : std::to_string(time);
}

const std::string KMC::timer::current_time()
{
    time_t t_curr = time(0);
    time_p = localtime(&t_curr);
    std::string t_str = month[time_p->tm_mon] + " " + std::to_string(time_p->tm_mday) + "/" + two_digits(time_p->tm_hour) + ":" + two_digits(time_p->tm_min)
     + ":" + two_digits(time_p->tm_sec) + " ";
    return t_str;
}

#endif