//
// Created by Eric Moss on 6/26/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_INDEX_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_INDEX_H

#include <complex>
#include <sstream>
#include <string>
#include <iostream>
#include "Auxilliary.h"

class Index {
public:
    Index(int a, int b);

    std::complex<double> getComplex(int d) const;
    double getAbs(int d) const;
    double getAngle(int d) const;

    Index rotate(int d) const;
    Index conj(int d) const;
    std::string to_string() const;

    int getA() const;
    int getB() const;

private:
    int a;
    int b;

    friend std::ostream& operator<<(std::ostream&, const Index&);
    friend bool operator==(const Index& index1, const Index& index2);
    friend bool operator!=(const Index& index1, const Index& index2);
    friend Index operator-(const Index& lhs, const Index& rhs);
    friend bool operator<(const Index& index1, const Index& index2);
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
