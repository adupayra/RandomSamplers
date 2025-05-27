#pragma once

#include <stdint.h>
#include <vector>
#include <random>

/**
 * This file lists a collection of random sampling algorithms.
 * Note: The correctness and random properties of these samplers have not been tested thoroughly against errors as these solutions have eventualy been dropped, though they seem to fill their
 * purpose of sampling K elements from N with equal chances in O(K) time complexity
*/

class RandomSamplers
{
public:
    RandomSamplers(uint64_t N, uint64_t K);

    /**
     * Vitters' algorithms from Vitter, J. S. (1984). Faster methods for random sampling. Communications of the ACM, 27(7), 703-718.
    */
    void VitterD(); 
    std::vector<uint64_t> VitterA(uint64_t StartingIndex, uint64_t K, uint64_t N, bool shuffle = false);

    void Floyd(); // From Jon Bentley and Bob Floyd. 1987. Programming pearls: a sample of brilliance. Commun. ACM 30, 9 (Sept. 1987), 754�757. https://doi.org/10.1145/30401.315746
    void HiddenShuffle(); // from Shekelyan, M.& amp; Cormode, G.. (2021).Sequential Random Sampling Revisited : Hidden Shuffle Method . Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, in Proceedings of Machine Learning Research 130 : 3628 - 3636 Available from https ://proceedings.mlr.press/v130/shekelyan21a.html.
    
    void Reset();

    /**
     * Note: The PRNG used is the C++ default one, though it is definitely not the most robust in terms of randomness.
     * I intended to implement a xoshiro256++/** before I dropped these sampling algorithm solutions (or xoshiro256+ for algorithm only producing floats in [0, 1])
    */
    std::mt19937_64 gen{ std::random_device()() };
    std::uniform_real_distribution<double> dis{ 0.0, 1.0 };
    std::uniform_int_distribution<unsigned long long> disInt;

    std::vector<uint64_t> result;

private:
    uint64_t samplesSize;
    uint64_t recordSize;
};

