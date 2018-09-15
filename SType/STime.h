#ifndef _STIME_H_
#define _STIME_H_

#include "STypePreDefine.h"

namespace SRobot
{
    //时间
    class STime
    {
      public:
        SReal Time;

      public:
        STime()
        {
        }

        STime(const SReal fTime)
            : Time(fTime)
        {
        }
        STime &operator=(const STime &rkSTime)
        {
            Time = rkSTime.Time;
            return *this;
        }

        STime &operator=(const SReal fTime)
        {
            Time = fTime;
            return *this;
        }

        bool operator==(const STime &rkSTime) const { return (Time == rkSTime.Time); }
        bool operator!=(const STime &rkSTime) const { return (Time != rkSTime.Time); }
        STime operator+(const STime &rkSTime) const { return STime(Time + rkSTime.Time); }
        STime operator-(const STime &rkSTime) const { return STime(Time - rkSTime.Time); }
        STime operator*(const SReal fScalar) const { return STime(Time * fScalar); }
        STime operator*(const STime &rkSTime) const { return STime(Time * rkSTime); }

        STime operator/(const SReal fScalar) const
        {
            assert(fScalar != 0.0);
            return STime(Time / fScalar);
        }

        STime operator/(const STime &rkSTime) const
        {
            assert(rkSTime.Time != 0.0);
            return STime(Time / rkSTime.Time);
        }

        const STime &operator+() const { return *this; }
        STime operator-() const { return STime(-Time); }
        friend STime operator*(const SReal fScalar, const STime &rkSTime) { return STime(fScalar * rkSTime.Time); }
        friend STime operator/(const SReal fScalar, const STime &rkSTime)
        {
            assert(rkSTime.Time != 0.0);
            return STime(fScalar / rkSTime.Time);
        }

        friend STime operator+(const STime &rkSTime, const SReal fTime) { return STime(rkSTime.Time + fTime); }
        friend STime operator+(const SReal fTime, const STime &rkSTime) { return STime(fTime + rkSTime.Time); }
        friend STime operator-(const STime &rkSTime, const SReal fTime) { return STime(rkSTime.Time - fTime); }
        friend STime operator-(const SReal fTime, const STime &rkSTime) { return STime(fTime - rkSTime.Time); }

        STime &operator+=(const STime &rkSTime)
        {
            Time += rkSTime.Time;
            return *this;
        }

        STime &operator+=(const SReal fTime)
        {
            Time += fTime;
            return *this;
        }

        STime &operator-=(const STime &rkSTime)
        {
            Time -= rkSTime.Time;
            return *this;
        }

        STime &operator-=(const SReal fTime)
        {
            Time -= fTime;
            return *this;
        }

        STime &operator*=(const SReal fScalar)
        {
            Time *= fScalar;
            return *this;
        }

        STime &operator*=(const STime &rkSTime)
        {
            Time *= rkSTime.Time;
            return *this;
        }

        STime &operator/=(const SReal fScalar)
        {
            assert(fScalar != 0.0);
            Time /= fScalar;
            return *this;
        }

        STime &operator/=(const STime &rkSTime)
        {
            Time /= rkSTime.Time;
            return *this;
        }
    };
}

#endif