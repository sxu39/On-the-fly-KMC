#ifndef KMC_POTENTIAL_H
#define KMC_POTENTIAL_H

#include "cell.h"
#include <string>
#include "deepmd/DeepPot.h"

namespace KMC
{
    class potential
    {
        public:
            virtual prec_t infer(cell *) = 0;
            virtual std::vector<prec_t> infer(std::vector<cell> &) = 0;
            potential() = default;
            virtual ~potential() = default;
    };
}

namespace KMC
{
    class ml_potential : public potential
    {
        public:
            deepmd::DeepPot dp;
        public:
            ml_potential();
            explicit ml_potential(const std::string &);
            ~ml_potential() = default;
            void init(const std::string &);
            prec_t infer(cell *) override;
            std::vector<prec_t> infer(std::vector<cell> &) override;
    };
}

#endif
