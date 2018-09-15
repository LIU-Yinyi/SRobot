#include "SBase.h"
#include "SType.h"
#include "SMath.h"
#include "SAsmMath.h"

namespace SRobot
{
    const SReal SMath::POS_INFINITY = std::numeric_limits<SReal>::infinity();
    const SReal SMath::NEG_INFINITY = -std::numeric_limits<SReal>::infinity();
    const SReal SMath::PI = SReal(4.0 * atan(1.0));
    const SReal SMath::TWO_PI = SReal(2.0 * PI);
    const SReal SMath::HALF_PI = SReal(0.5 * PI);
    const SReal SMath::fDeg2Rad = PI / SReal(180.0);
    const SReal SMath::fRad2Deg = SReal(180.0) / PI;
    const SReal SMath::LOG2 = log(SReal(2.0));

    int SMath::mTrigTableSize;

    SReal SMath::mTrigTableFactor;
    SReal *SMath::mSinTable = NULL;
    SReal *SMath::mTanTable = NULL;

    //-----------------------------------------------------------------------
    SMath::SMath(unsigned int trigTableSize)
    {
        mTrigTableSize = trigTableSize;
        mTrigTableFactor = mTrigTableSize / SMath::TWO_PI;

        mSinTable = new SReal[mTrigTableSize];
        mTanTable = new SReal[mTrigTableSize];

        buildTrigTables();
    }

    //-----------------------------------------------------------------------
    SMath::~SMath()
    {
        delete[] mSinTable;
        delete[] mTanTable;
    }

    //-----------------------------------------------------------------------
    void SMath::buildTrigTables(void)
    {
        SReal angle;
        for (int i = 0; i < mTrigTableSize; ++i)
        {
            angle = SMath::TWO_PI * i / mTrigTableSize;
            mSinTable[i] = sin(angle);
            mTanTable[i] = tan(angle);
        }
    }
    //-----------------------------------------------------------------------
    SReal SMath::SinTable(SReal fValue)
    {
        int idx;
        if (fValue >= 0)
        {
            idx = int(fValue * mTrigTableFactor) % mTrigTableSize;
        }
        else
        {
            idx = mTrigTableSize - (int(-fValue * mTrigTableFactor) % mTrigTableSize) - 1;
        }

        return mSinTable[idx];
    }
    //-----------------------------------------------------------------------
    SReal SMath::TanTable(SReal fValue)
    {
        int idx = int(fValue *= mTrigTableFactor) % mTrigTableSize;
        return mTanTable[idx];
    }
    //-----------------------------------------------------------------------
    int SMath::ISign(int iValue)
    {
        return (iValue > 0 ? +1 : (iValue < 0 ? -1 : 0));
    }
    //-----------------------------------------------------------------------
    SRadian SMath::ACos(SReal fValue)
    {
        if (-1.0 < fValue)
        {
            if (fValue < 1.0)
                return SRadian(acos(fValue));
            else
                return SRadian(0.0);
        }
        else
        {
            return SRadian(PI);
        }
    }
    //-----------------------------------------------------------------------
    SRadian SMath::ASin(SReal fValue)
    {
        if (-1.0 < fValue)
        {
            if (fValue < 1.0)
                return SRadian(asin(fValue));
            else
                return SRadian(HALF_PI);
        }
        else
        {
            return SRadian(-HALF_PI);
        }
    }
    //-----------------------------------------------------------------------
    SReal SMath::Sign(SReal fValue)
    {
        if (fValue > 0.0)
            return 1.0;

        if (fValue < 0.0)
            return -1.0;

        return 0.0;
    }
    //-----------------------------------------------------------------------
    SReal SMath::InvSqrt(SReal fValue)
    {
        return SReal(asm_rsq(fValue));
    }
    //-----------------------------------------------------------------------
    SReal SMath::UnitRandom()
    {
        return asm_rand() / asm_rand_max();
    }

    //-----------------------------------------------------------------------
    SReal SMath::RangeRandom(SReal fLow, SReal fHigh)
    {
        return (fHigh - fLow) * UnitRandom() + fLow;
    }

    //-----------------------------------------------------------------------
    SReal SMath::SymmetricRandom()
    {
        return 2.0f * UnitRandom() - 1.0f;
    }

    //-----------------------------------------------------------------------
    bool SMath::pointInTri2D(const SVector2 &p, const SVector2 &a, const SVector2 &b, const SVector2 &c)
    {

        SVector2 v1, v2;
        SReal dot[3];
        bool zeroDot[3];

        v1 = b - a;
        v2 = p - a;

        dot[0] = v1.crossProduct(v2);
        zeroDot[0] = SMath::RealEqual(dot[0], (SReal)0.0f, (SReal)1e-3);

        v1 = c - b;
        v2 = p - b;

        dot[1] = v1.crossProduct(v2);
        zeroDot[1] = SMath::RealEqual(dot[1], (SReal)0.0f, (SReal)1e-3);

        if (!zeroDot[0] && !zeroDot[1] && SMath::Sign(dot[0]) != SMath::Sign(dot[1]))
        {
            return false;
        }

        v1 = a - c;
        v2 = p - c;

        dot[2] = v1.crossProduct(v2);
        zeroDot[2] = SMath::RealEqual(dot[2], (SReal)0.0f, (SReal)1e-3);

        if ((!zeroDot[0] && !zeroDot[2] && SMath::Sign(dot[0]) != SMath::Sign(dot[2])) ||
            (!zeroDot[1] && !zeroDot[2] && SMath::Sign(dot[1]) != SMath::Sign(dot[2])))
        {
            return false;
        }

        return true;
    }
    //-----------------------------------------------------------------------
    bool SMath::pointInTri3D(const SVector3 &p, const SVector3 &a, const SVector3 &b, const SVector3 &c, const SVector3 &normal)
    {

        SVector3 v1, v2;
        SReal dot[3];
        bool zeroDot[3];

        v1 = b - a;
        v2 = p - a;

        dot[0] = v1.crossProduct(v2).dotProduct(normal);
        zeroDot[0] = SMath::RealEqual(dot[0], (SReal)0.0f, (SReal)1e-3);

        v1 = c - b;
        v2 = p - b;

        dot[1] = v1.crossProduct(v2).dotProduct(normal);
        zeroDot[1] = SMath::RealEqual(dot[1], (SReal)0.0f, (SReal)1e-3);

        if (!zeroDot[0] && !zeroDot[1] && SMath::Sign(dot[0]) != SMath::Sign(dot[1]))
        {
            return false;
        }

        v1 = a - c;
        v2 = p - c;

        dot[2] = v1.crossProduct(v2).dotProduct(normal);
        zeroDot[2] = SMath::RealEqual(dot[2], (SReal)0.0f, (SReal)1e-3);

        if ((!zeroDot[0] && !zeroDot[2] && SMath::Sign(dot[0]) != SMath::Sign(dot[2])) ||
            (!zeroDot[1] && !zeroDot[2] && SMath::Sign(dot[1]) != SMath::Sign(dot[2])))
        {
            return false;
        }

        return true;
    }
    //-----------------------------------------------------------------------
    bool SMath::RealEqual(SReal a, SReal b, SReal tolerance)
    {
        if (fabs(b - a) <= tolerance)
            return true;
        else
            return false;
    }
}
