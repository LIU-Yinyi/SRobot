#ifndef _SMATH_H_
#define _SMATH_H_

#include "SBase.h"
#include "SType.h"
#include "SAsmMath.h"
namespace SRobot
{
    class SMath
    {
      public:
        class RandomValueProvider
        {
          public:
            virtual ~RandomValueProvider() {}
            virtual SReal getRandomUnit() = 0;
        };

      protected:
        static int mTrigTableSize;

        static SReal mTrigTableFactor;
        static SReal *mSinTable;
        static SReal *mTanTable;

        static RandomValueProvider *mRandProvider;

        void buildTrigTables();  

        static SReal SinTable(SReal fValue);
        static SReal TanTable(SReal fValue);

      public:
        SMath(unsigned int trigTableSize = 4096);
        ~SMath();

        static inline int IAbs(int iValue) { return (iValue >= 0 ? iValue : -iValue); }
        static inline int ICeil(float fValue) { return int(ceil(fValue)); }
        static inline int IFloor(float fValue) { return int(floor(fValue)); }
        static int ISign(int iValue);

        static inline SReal Abs(SReal fValue) { return SReal(fabs(fValue)); }
        static inline SDegree Abs(const SDegree &dValue) { return SDegree(fabs(dValue.valueDegrees())); }
        static inline SRadian Abs(const SRadian &rValue) { return SRadian(fabs(rValue.valueRadians())); }

        static SRadian ACos(SReal fValue);
        static SRadian ASin(SReal fValue);
        static inline SRadian ATan(SReal fValue) { return SRadian(atan(fValue)); }
        static inline SRadian ATan2(SReal fY, SReal fX) { return SRadian(atan2(fY, fX)); }

        static inline SReal Ceil(SReal fValue) { return SReal(ceil(fValue)); }

        static inline bool isNaN(SReal f) { return f != f; }

        static inline SReal Cos(const SRadian &fValue, bool useTables = false)
        {
            return (!useTables) ? SReal(cos(fValue.valueRadians())) : SinTable(fValue.valueRadians() + HALF_PI);
        }

        static inline SReal Cos(SReal fValue, bool useTables = false)
        {
            return (!useTables) ? SReal(cos(fValue)) : SinTable(fValue + HALF_PI);
        }

        static inline SReal Exp(SReal fValue) { return SReal(exp(fValue)); }

        static inline SReal Floor(SReal fValue) { return SReal(floor(fValue)); }

        static inline SReal Log(SReal fValue) { return SReal(log(fValue)); }

        static const SReal LOG2;

        static inline SReal Log2(SReal fValue) { return SReal(log(fValue) / LOG2); }

        static inline SReal LogN(SReal base, SReal fValue) { return SReal(log(fValue) / log(base)); }

        static inline SReal Pow(SReal fBase, SReal fExponent) { return SReal(pow(fBase, fExponent)); }

        static SReal Sign(SReal fValue);
        static inline SRadian Sign(const SRadian &rValue) { return SRadian(Sign(rValue.valueRadians())); }
        static inline SDegree Sign(const SDegree &dValue) { return SDegree(Sign(dValue.valueDegrees())); }

        static inline float saturate(float t) { return (t < 0) ? 0 : ((t > 1) ? 1 : t); }
        static inline double saturate(double t) { return (t < 0) ? 0 : ((t > 1) ? 1 : t); }

        template <typename V, typename T>
        static V lerp(const V &v0, const V &v1, const T &t)
        {
            return v0 * (1 - t) + v1 * t;
        }

        static inline SReal Sin(const SRadian &fValue, bool useTables = false)
        {
#if CALCULATE_MODE == ASM
        return asm_sin(fValue.valueRadians());
#else
         return (!useTables) ? SReal(sin(fValue.valueRadians())) : SinTable(fValue.valueRadians());
#endif
        }

        static inline SReal Sin(SReal fValue, bool useTables = false)
        {
#if CALCULATE_MODE == ASM
            return asm_sin(fValue);
#else
            return (!useTables) ? SReal(sin(fValue)) : SinTable(fValue);
#endif
        }

        static inline SReal Sqr(SReal fValue) { return fValue * fValue; }
        static inline SReal Sqrt(SReal fValue) { return SReal(sqrt(fValue)); }
        static inline SRadian Sqrt(const SRadian &fValue) { return SRadian(Sqr(fValue.valueRadians())); }
        static inline SDegree Sqrt(const SDegree &fValue) { return SDegree(Sqrt(fValue.valueDegrees())); }

        static SReal InvSqrt(SReal fValue);

        static SReal UnitRandom();
        static SReal RangeRandom(SReal fLow, SReal fHigh);
        static SReal SymmetricRandom();

        static void SetRandomValueProvider(RandomValueProvider *provider);

        static inline SReal Tan(const SRadian &fValue, bool useTables = false)
        {
#if CALCULATE_MODE == ASM
            return asm_tan(fValue.valueRadians());
#else
            return (!useTables) ? SReal(tan(fValue.valueRadians())) : TanTable(fValue.valueRadians());
#endif
        }

        static inline SReal Tan(SReal fValue, bool useTables = false)
        {
#if CALCULATE_MODE == ASM
            return asm_tan(fValue);
#else
            return (!useTables) ? SReal(tan(fValue)) : TanTable(fValue);
#endif
        }

        static inline SReal DegreesToRadians(SReal degrees) { return degrees * fDeg2Rad; }
        static inline SReal RadiansToDegrees(SReal radians) { return radians * fRad2Deg; }

        static bool pointInTri2D(const SVector2 &p, const SVector2 &a,
                                 const SVector2 &b, const SVector2 &c);

        static bool pointInTri3D(const SVector3 &p, const SVector3 &a,
                                 const SVector3 &b, const SVector3 &c, const SVector3 &normal);

        static bool RealEqual(SReal a, SReal b,
                              SReal tolerance = std::numeric_limits<SReal>::epsilon());

        static SReal gaussianDistribution(SReal x, SReal offset = 0.0f, SReal scale = 1.0f);

        template <typename T>
        static T Clamp(T val, T minval, T maxval)
        {
            assert(minval <= maxval && "Invalid clamp range");
            return (std::max)((std::min)(val, maxval), minval);
        }

		static SUint64 Hash(SString String)
		{
			return CityHash64(String.c_str(),String.length());
		}

        static const SReal POS_INFINITY;
        static const SReal NEG_INFINITY;
        static const SReal PI;
        static const SReal TWO_PI;
        static const SReal HALF_PI;
        static const SReal fDeg2Rad;
        static const SReal fRad2Deg;
    };
}

#endif
