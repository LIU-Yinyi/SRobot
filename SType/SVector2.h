#ifndef _SVECTOR2_H_
#define _SVECTOR2_H_

#include "STypePreDefine.h"
#include "SAngle.h"

namespace SRobot
{
    class SVector2
    {
      public:
        SReal x;
        SReal y;

      public:
        SVector2(const SReal fX=0, const SReal fY=0);

        explicit SVector2(const SReal scaler);

        explicit SVector2(const SReal afCoordinate[2]);

        explicit SVector2(const int afCoordinate[2]);

        SVector2(SReal *const r);

        void swap(SVector2 &other);

        SReal operator[](const size_t i) const;

        SReal &operator[](const size_t i);

        SReal *ptr();

        const SReal *ptr() const;

        SVector2 &operator=(const SVector2 &rkSVector);
        SVector2 &operator=(const SReal fScalar);
        bool operator==(const SVector2 &rkSVector) const;
        bool operator!=(const SVector2 &rkSVector) const;
        const SVector2 &operator+() const;
        SVector2 operator+(const SVector2 &rkSVector) const;
        SVector2 operator-() const;
        SVector2 operator-(const SVector2 &rkSVector) const;
        SVector2 operator*(const SReal fScalar) const;
        SVector2 operator*(const SVector2 &rhs) const;
        SVector2 operator/(const SReal fScalar) const;
        SVector2 operator/(const SVector2 &rhs) const;

        friend SVector2 operator+(const SVector2 &lhs, const SReal rhs);
        friend SVector2 operator+(const SReal lhs, const SVector2 &rhs);
        friend SVector2 operator-(const SVector2 &lhs, const SReal rhs);
        friend SVector2 operator-(const SReal lhs, const SVector2 &rhs);
        friend SVector2 operator*(const SReal fScalar, const SVector2 &rkSVector);
        friend SVector2 operator/(const SReal fScalar, const SVector2 &rkSVector);


        SVector2 &operator+=(const SVector2 &rkSVector);
        SVector2 &operator+=(const SReal fScaler);
        SVector2 &operator-=(const SVector2 &rkSVector);
        SVector2 &operator-=(const SReal fScaler);
        SVector2 &operator*=(const SReal fScalar);
        SVector2 &operator*=(const SVector2 &rkSVector);
        SVector2 &operator/=(const SReal fScalar);
        SVector2 &operator/=(const SVector2 &rkSVector);

        SReal length() const;

        SReal squaredLength() const;

        SReal distance(const SVector2 &rhs) const;

        SReal squaredDistance(const SVector2 &rhs) const;

        SReal dotProduct(const SVector2 &vec) const;

        /**
         * 将该向量转化为单位向量
         * @return 
         * 返回向量在转换之前的长度
         */
        SReal normalise();

        SVector2 midPoint(const SVector2 &vec) const;

        bool operator<(const SVector2 &rhs) const;

        bool operator>(const SVector2 &rhs) const;

        void makeFloor(const SVector2 &cmp);

        void makeCeil(const SVector2 &cmp);

        SVector2 perpendicular(void) const;

        SReal crossProduct(const SVector2 &rkSVector) const;

        SVector2 randomDeviant(SRadian angle) const;

        bool isZeroLength(void) const;

        SVector2 normalisedCopy(void) const;

        SVector2 reflect(const SVector2 &normal) const;

        bool isNaN() const;

        SRadian angleBetween(const SVector2 &other) const;

        SRadian angleTo(const SVector2 &other) const;

        static const SVector2 ZERO;
        static const SVector2 UNIT_X;
        static const SVector2 UNIT_Y;
        static const SVector2 NEGATIVE_UNIT_X;
        static const SVector2 NEGATIVE_UNIT_Y;
        static const SVector2 UNIT_SCALE;

        friend std::ostream &operator<<(std::ostream &o, const SVector2 &v);
    };
}

#endif
