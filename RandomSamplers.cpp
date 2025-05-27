#include "RandomSamplers.h"
#include <unordered_set>

RandomSamplers::RandomSamplers(uint64_t N, uint64_t K)
{
    samplesSize = K;
    recordSize = N;
    result.reserve(K);
}

void RandomSamplers::HiddenShuffle()
{
    uint64_t N = this->recordSize;
    uint64_t n = this->samplesSize;

    uint64_t H = 0;
    uint64_t i = 0;

    if (N > n)
    {
        H = n;

        while (i < n)
        {
            double q = 1 - static_cast<double>(N - n) / (N - i);

            i = i + static_cast<uint64_t>(log(dis(gen)) / log(1 - q));
            double p_i = 1 - static_cast<double>(N - n) / (N - i);
            if (i < n && dis(gen) < p_i / q)
            {
                H--;
            }
            i++;
        }
    }
    uint64_t L = n - H;
    double a = 1.f;

    while (H > 0)
    {
        uint64_t S_old = static_cast<uint64_t>(n + a * (N - n));
        a = a * pow(dis(gen), 1.f / H);
        uint64_t S = static_cast<uint64_t>(n + a * (N - n));
        if (S < S_old)
        {
            result.push_back(N - 1 - S);
        }
        else
        {
            L++;
        }
        H--;
    }

    while (L > 0)
    {
        double u = dis(gen);
        uint64_t s = 0;
        double F = static_cast<double>(L) / n;

        while (F < u && s < (n - L))
        {
            F = 1 - (1 - static_cast<double>(L) / (n - s)) * (1 - F);
            s++;
        }
        L--;
        n = n - s - 1;
        result.push_back(N - 1 - n);
    }
    std::shuffle(std::begin(result), std::end(result), gen);
}

void RandomSamplers::Reset()
{
    result.clear();
    result.reserve(this->samplesSize);
}

void RandomSamplers::VitterD()
{
    long double nreal, ninv, Nreal, nmin1inv, Vprime, U, X, y1, y2, top, bottom, negSreal, qu1real;
    uint64_t qu1, threshold, S, n, N, limit;
    int negalphainv;
    uint64_t currentRecordIndex = 0;
    n = samplesSize;
    N = recordSize;
    nreal = static_cast<long double>(samplesSize);
    ninv = 1 / nreal;
    Nreal = static_cast<long double>(recordSize);
    Vprime = exp(log(dis(gen)) * ninv);
    qu1 = 1 - n + N;
    qu1real = Nreal - nreal + 1;
    negalphainv = -13;
    threshold = -negalphainv * n;
    while (n > 1 && threshold < N)
    {
        nmin1inv = 1 / (nreal - 1);
        while (true)
        {
            while (true)
            {
                X = Nreal * (-Vprime + 1);
                S = static_cast<uint64_t>(X);
                if (S < qu1)
                {
                    break;
                }
                Vprime = exp(log(dis(gen)) * ninv);
            }
            U = dis(gen);
            negSreal = -(long double)S;
            y1 = exp(log(U * Nreal / qu1real) * nmin1inv);
            Vprime = y1 * (-X / Nreal + 1) * (qu1real / (negSreal + qu1real));
            if (Vprime <= 1)
            {
                break;
            }
            y2 = 1;
            top = Nreal - 1;
            if (n - 1 > S)
            {
                bottom = Nreal - nreal;
                limit = N - S;
            }
            else
            {
                bottom = Nreal + negSreal - 1;
                limit = qu1;
            }
            for (size_t t = N - 1; t >= limit; t--)
            {
                y2 *= top / bottom;
                top -= 1;
                bottom -= 1;
            }
            if (Nreal / (Nreal - X) >= y1 * exp(log(y2) * nmin1inv))
            {
                Vprime = exp(log(dis(gen)) * nmin1inv);
                break;
            }
            Vprime = exp(log(dis(gen)) * ninv);

        }
        currentRecordIndex += S + 1;
        result.push_back(currentRecordIndex);
        N = N - S - 1;
        Nreal = negSreal + Nreal - 1;
        n -= 1;
        nreal -= 1;
        ninv = nmin1inv;
        qu1 -= S;
        qu1real += negSreal;
        threshold += negalphainv;
    }
    if (n > 1)
    {
        result.append_range(VitterA(currentRecordIndex, n, N));
    }
    else
    {
        S = static_cast<uint64_t>(N * Vprime);
        currentRecordIndex += S + 1;
        result.push_back(currentRecordIndex);
    }
    std::shuffle(std::begin(result), std::end(result), gen);

}

std::vector<uint64_t> RandomSamplers::VitterA(uint64_t startingIndex, uint64_t n, uint64_t N, bool shuffle)
{
    long double quot, Nreal, top;
    long double V;
    top = static_cast<long double>(N - n);
    Nreal = static_cast<long double>(N);
    uint64_t S = 0;
    uint64_t currentRecordIndex = startingIndex;
    std::vector<uint64_t> randomSamples;
    while (n >= 2)
    {
        V = dis(gen);
        S = 0;
        quot = top / Nreal;
        while (quot > V)
        {
            S += 1;
            top -= 1;
            Nreal -= 1;
            quot *= top / Nreal;
        }
        currentRecordIndex += S + 1;
        randomSamples.push_back(currentRecordIndex);
        Nreal -= 1;
        n -= 1;
    }

    S = static_cast<uint64_t>(Nreal * dis(gen));
    currentRecordIndex += S + 1;
    randomSamples.push_back(currentRecordIndex);
    if (shuffle)
    {
        std::shuffle(std::begin(randomSamples), std::end(randomSamples), gen);
    }
    return randomSamples;
}

void RandomSamplers::Floyd()
{
    /**
     * Floyd algorithm is the only one that requires the generation of pseudo random numbers that are 64 bits long.
     * Such numbers can represent some difficulties to be generated by pseudo random generator and the generator used should be random tested on the generation of such values.
     * However, this algorithm being the slowest of the one implemented in O(K) and the only one generating 64bits numbers, I did not bother doing these tests.
     * In conclusion, this algorithm might not work properly with big N, some tests should be conducted on disInt to make sure it works.
    */
    uint64_t N = recordSize;
    uint64_t K = samplesSize;
    std::unordered_set<uint64_t> randomSamples;
    randomSamples.reserve(K);
    for (size_t i = N - K; i < N; i++)
    {
        disInt = std::uniform_int_distribution<uint64_t>(0, i);
        uint64_t newValue = disInt(gen);

        uint64_t valueToAdd;
        auto itr = randomSamples.find(newValue);
        auto itr2 = randomSamples.end();

        if (itr == itr2)
        {
            valueToAdd = newValue;
        }
        else
        {
            valueToAdd = i + 1;
        }
        randomSamples.insert(valueToAdd);

    }

    for (auto itr : randomSamples)
    {
        result.push_back(itr);
    }
    std::shuffle(std::begin(result), std::end(result), gen);
}