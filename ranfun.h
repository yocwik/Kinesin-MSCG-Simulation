//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file contains functions for generating random numbers
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef RANFUN_INCLUDED
#define RANFUN_INCLUDED

#include <ctime>
#include <cmath>
#include <cstdlib>

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long idum)
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    if (idum < 0 || iff == 0)
    {
        iff=1;
        mj=labs(MSEED-labs(idum));
        mj %= MBIG;
        ma[55]=mj;
        mk=1;

        for (i=1;i<=54;i++) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;

            if (mk < MZ)
            {
                mk += MBIG;
            }

            mj=ma[ii];
        }

        for (k=1;k<=4;k++)
        {
            for(i=1;i<=55;i++)
            {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ)
                {
                    ma[i] += MBIG;
                }
            }
        }
        inext=0;
        inextp=31;
        idum=1;
    }

if (++inext == 56) inext=1;
if (++inextp == 56) inextp=1;
mj=ma[inext]-ma[inextp];
if (mj < MZ) mj += MBIG;
ma[inext]=mj;
return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

using namespace std;

float nrand(long& init)
{
    static int count = 0;
    static float nextGaussianVal;
    float firstGaussianVal, v1, v2, s;

    if (count == 0) {
       do {
           v1 = 2 * ran3(init) - 1;   // between -1.0 and 1.0
           v2 = 2 * ran3(init) - 1;   // between -1.0 and 1.0
           s = v1 * v1 + v2 * v2;
        } while (s > 1 || s == 0);
        float multiplier = sqrt(-2 * log(s)/s);
        nextGaussianVal = v2 * multiplier;
        firstGaussianVal = v1 * multiplier;
        count = 1;
        return firstGaussianVal;
    }

    count = 0;
    return nextGaussianVal;
}

#endif
