//
// Created by Eric Moss on 6/26/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_INDEX_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_INDEX_H

#include <complex>
#include <sstream>
#include <string>
#include <iostream>
#include "Auxiliary.h"

class Index {
public:
    Index(int a, int b);

    [[nodiscard]] std::complex<double> getComplex(int d) const;
    [[nodiscard]] double getAbs(int d) const;
    [[nodiscard]] double getAngle(int d) const;

    [[nodiscard]] Index rotate(int d) const;
    [[nodiscard]] Index conj(int d) const;
    std::string to_string() const;

    Index mul(const Index& index, int d) const;

    [[nodiscard]] int getA() const;
    [[nodiscard]] int getB() const;

    friend std::ostream& operator<<(std::ostream&, const Index&);
    friend bool operator==(const Index& index1, const Index& index2);
    friend bool operator!=(const Index& index1, const Index& index2);
    friend Index operator-(const Index& lhs, const Index& rhs);
    friend bool operator<(const Index& index1, const Index& index2);

private:
    int a;
    int b;
};

template<>
struct std::hash<Index>
{
    std::size_t operator()(Index const& n) const
    {
        std::size_t hashA = std::hash<int>{}(n.getA());
        std::size_t hashB = std::hash<int>{}(n.getB());

        // Use prime numbers to mix the hash values
        // The choice of prime numbers can affect the distribution quality
        const std::size_t prime1 = 31;
        const std::size_t prime3 = 41;

        return hashA ^ (hashB * prime1);
    }
};

#endif //EUCLIDEAN_BIANCHI_MAASS_FORMS_INDEX_H
