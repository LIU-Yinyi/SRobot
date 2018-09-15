#ifndef _SVECTOR3_H_
#define _SVECTOR3_H_

#include "STypePreDefine.h"
#include "SAngle.h"
#include "SQuaternion.h"
#include "SVector2.h"

namespace SRobot
{
    class SVector3
    {
      public:
        SReal x, y, z;

      public:
        explicit SVector3(const SReal fX=0, const SReal fY=0, const SReal fZ=0);

        explicit SVector3(const SReal afCoordinate[3]);

        explicit SVector3(const int afCoordinate[3]);

        explicit SVector3(SReal *const r);

        explicit SVector3(const SReal scaler);

        void swap(SVector3 &other);

        SReal operator[](const size_t i) const;

        SReal &operator[](const size_t i);

        SReal *ptr();

        const SReal *ptr() const;

        SVector3 &operator=(const SVector3 &rkSVector);
        SVector3 &operator=(const SReal fScaler);
        bool operator==(const SVector3 &rkSVector) const;
        bool operator!=(const SVector3 &rkSVector) const;

        const SVector3 &operator+() const;
        SVector3 operator-() const;        
        SVector3 operator+(const SVector3 &rkSVector) const;
        SVector3 operator-(const SVector3 &rkSVector) const;
        SVector3 operator*(const SReal fScalar) const;
        SVector3 operator*(const SVector3 &rhs) const;
        SVector3 operator/(const SReal fScalar) const;
        SVector3 operator/(const SVector3 &rhs) const;

        friend SVector3 operator+(const SVector3 &lhs, const SReal rhs);
        friend SVector3 operator+(const SReal lhs, const SVector3 &rhs);
        friend SVector3 operator-(const SVector3 &lhs, const SReal rhs);
        friend SVector3 operator-(const SReal lhs, const SVector3 &rhs);
        friend SVector3 operator*(const SReal fScalar, const SVector3 &rkSVector);
        friend SVector3 operator/(const SReal fScalar, const SVector3 &rkSVector);
        
        SVector3 &operator+=(const SVector3 &rkSVector);
        SVector3 &operator+=(const SReal fScalar);
        SVector3 &operator-=(const SVector3 &rkSVector);
        SVector3 &operator-=(const SReal fScalar);
        SVector3 &operator*=(const SReal fScalar);
        SVector3 &operator*=(const SVector3 &rkSVector);
        SVector3 &operator/=(const SReal fScalar);
        SVector3 &operator/=(const SVector3 &rkSVector);

        bool operator<(const SVector3 &rhs) const;
        bool operator>(const SVector3 &rhs) const;

        SReal length() const;

        SReal squaredLength() const;

        SReal distance(const SVector3 &rhs) const;

        SReal squaredDistance(const SVector3 &rhs) const;

        SReal dotProduct(const SVector3 &vec) const;

        SReal absDotProduct(const SVector3 &vec) const;

        SReal normalise();

        SVector3 crossProduct(const SVector3 &rkSVector) const;

        SVector3 midPoint(const SVector3 &vec) const;


        void makeFloor(const SVector3 &cmp);

        void makeCeil(const SVector3 &cmp);

        SVector3 perpendicular(void) const;

        SVector3 randomDeviant(const SRadian &angle, const SVector3 &up = SVector3::ZERO) const;

        SRadian angleBetween(const SVector3 &dest) const;

        SQuaternion getRotationTo(const SVector3 &dest, const SVector3 &fallbackAxis = SVector3::ZERO) const;

        bool isZeroLength(void) const;

        SVector3 normalisedCopy(void) const;

        SVector3 reflect(const SVector3 &normal) const;

        bool positionEquals(const SVector3 &rhs, SReal tolerance = 1e-03) const;

        bool positionCloses(const SVector3 &rhs, SReal tolerance = 1e-03f) const;

        bool directionEquals(const SVector3 &rhs, const SRadian &tolerance) const;

        bool isNaN() const;

        SVector3 primaryAxis() const;

        /**
         * @return
         * 将三位变量转换为去除z分量二维向量
         */
		SVector2 compressZ() const;

        static const SVector3 ZERO;
        static const SVector3 UNIT_X;
        static const SVector3 UNIT_Y;
        static const SVector3 UNIT_Z;
        static const SVector3 NEGATIVE_UNIT_X;
        static const SVector3 NEGATIVE_UNIT_Y;
        static const SVector3 NEGATIVE_UNIT_Z;
        static const SVector3 UNIT_SCALE;

        friend std::ostream &operator<<(std::ostream &o, const SVector3 &v);
    };
}

#endif
