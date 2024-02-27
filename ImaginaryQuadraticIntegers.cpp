//
// Created by Eric Moss on 2/3/24.
//

#include "ImaginaryQuadraticIntegers.h"
#include "Auxiliary.h"
#include <algorithm>

ImaginaryQuadraticIntegers::ImaginaryQuadraticIntegers(int d) {
    this->d = d;

    if (Auxiliary::mod(-d, 4) == 1) {
        //theta = 1/2 + I*sqrt(d)/2
        theta = std::complex<double> {1.0/2, sqrt(d)/2};
    } else if (Auxiliary::mod(-d, 4) == 0) {
        throw std::invalid_argument("d is incorrect");
    } else {
        //theta = I*sqrt(d)
        theta = std::complex<double> {0, sqrt(d)};
    }

    A = theta.imag();

    //These results are hardcoded output from John Cremona's bianchi-progs
    if (d == 3) {
        Y0 = sqrt(2.0/3);
    } else if (d == 19) {
        Y0 = sqrt(2.0/19);
    } else if (d == 43) {
        Y0 = sqrt(2.0/43);
    } else if (d == 67) {
        Y0 = sqrt(2.0/67);
    } else if (d == 163) {
        Y0 = sqrt(2.0/163);
    } else {
        //Y0 = sqrt(1- (1/2)^2 - (theta.imag()/2)^2)
        //=sqrt(3/4 - theta.imag()^2/4)
        Y0 = .75 - pow(this->getTheta().imag(),2)/4;
        Y0 = sqrt(Y0);
    }
}

ImaginaryQuadraticIntegers::ImaginaryQuadraticIntegers() {
    d = 0;
    A = 0;
    theta = {0,0};
    Y0 = 0;
}

vector<Index> ImaginaryQuadraticIntegers::indicesUpToM(const double M) {
    vector<Index> answer;
    // Computed bounds using Lagrange multipliers
    // aub = ceil(maxN * thetaModulus / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));
    int aUpperBound = ceil(M * abs(theta)/ sqrt(pow(abs(theta),2) - pow(theta.real(),2)));
    int aLowerBound = -aUpperBound;

    //bub = ceil(maxN / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));
    int bUpperBound = ceil(M /sqrt(pow(abs(theta),2) - pow(theta.real(),2)));
    int bLowerBound = -bUpperBound;

    for (int a = aLowerBound; a <= aUpperBound; a++) {
        for (int b = bLowerBound; b <= bUpperBound; b++) {
            if (a == 0 && b == 0) {
                continue;
            }
            Index index = Index(a, b);
            if (index.getAbs(d) <= M) {
                answer.push_back(index);
            }
        }
    }

    //Sort using lambda so that I don't have to store d in every Index object... waste of memory!
    auto indexComparator = [this](const Index& index1, const Index& index2) -> bool {
        if (index1.getAbs(d) < index2.getAbs(d)) {
            return true;
        } else if (index1.getAbs(d) == index2.getAbs(d) && index1.getAngle(d) < index2.getAngle(d)) {
            return true;
        } else {
            return false;
        }
    };

    /*Sort the indices by absolute value then by angle with positive real axis.*/
    std::sort(answer.begin(), answer.end(), indexComparator);
    return answer;
}

pair<vector<Index>, map<Index, vector<pair<Index, int>>>>
ImaginaryQuadraticIntegers::indexOrbitQuotientData(vector<Index> indices, const char symClass) {
    vector<Index> indexTransversal;
    map<Index, vector<pair<Index,int>>> orbitDataModSign;

    int rotationCoeff = symClass == 'D' || symClass == 'G' ? 1 : -1;
    int conjCoeff = symClass == 'D' || symClass == 'C' ? 1 : -1;


    while (!indices.empty()) {
        Index index = indices[0];
        vector<pair<Index, int>> orbit = {pair<Index, int>(index, 1)};

        Index tempIndex = index.rotate(d);
        int tempCoeff = 1*rotationCoeff;
        pair<Index, int> tempPair = {tempIndex, tempCoeff};
        while (tempIndex != index) {
            orbit.push_back(tempPair);

            tempIndex = tempIndex.rotate(d);
            tempCoeff = tempCoeff*rotationCoeff;
            tempPair = {tempIndex, tempCoeff};
        }

        Index conjIndex = index.conj(d);
        bool conjIsInRotations = false;
        for (const auto& tup : orbit) {
            if (conjIndex == tup.first) {
                conjIsInRotations = true;
                break;
            }
        }

        if (!conjIsInRotations) {
            /*add in the rotations of the conjugate*/
            tempPair = {conjIndex, conjCoeff};
            orbit.push_back(tempPair);

            tempIndex = conjIndex.rotate(d);
            tempCoeff = conjCoeff*rotationCoeff;
            tempPair = {tempIndex, tempCoeff};
            while (tempIndex != conjIndex) {
                orbit.push_back(tempPair);

                tempIndex = tempIndex.rotate(d);
                tempCoeff = tempCoeff*rotationCoeff;
                tempPair = {tempIndex, tempCoeff};
            }

        }
        /*Save our starting index to the transversal.*/
        indexTransversal.push_back(index);


        /*Compute the orbit mod +-1*/
        vector<pair<Index,int>> orbitModSign;
        vector<Index> alreadyGotten;
        for (const auto& tup : orbit) {
            Index l = tup.first;
            int pmOne = tup.second;
            if (std::find(alreadyGotten.begin(), alreadyGotten.end(),l) == alreadyGotten.end()) {
                //add it
                pair<Index,int> classModMinusOne = {l,pmOne};
                alreadyGotten.push_back(l);
                alreadyGotten.push_back(Index(-l.getA(), -l.getB()));
                orbitModSign.push_back(classModMinusOne);
            }
        }
        orbitDataModSign[index] = orbitModSign;

        /*Delete the indicesM0 in the orbit from the copied list of indicesM0.*/
        for (auto itr1 : orbit) {
            Index toDelete = itr1.first;
            auto toDeleteItr = std::find(indices.begin(), indices.end(), toDelete);
            indices.erase(toDeleteItr);
        }
    }

    return {indexTransversal, orbitDataModSign};
}



