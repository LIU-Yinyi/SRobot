#ifndef _SQUATERNION_H_
#define _SQUATERNION_H_

#include "STypePreDefine.h"

namespace SRobot
{
    //四元数
    class SQuaternion
    {
      public:
        SQuaternion()
            : w(1), x(0), y(0), z(0)
        {
        }

        SQuaternion(SReal fW, SReal fX, SReal fY, SReal fZ)
            : w(fW), x(fX), y(fY), z(fZ)
        {
        }

        SQuaternion(const SMatrix3 &rot)
        {
            this->FromRotationMatrix(rot);
        }

        SQuaternion(const SRadian &rfAngle, const SVector3 &rkAxis)
        {
            this->FromAngleAxis(rfAngle, rkAxis);
        }

        SQuaternion(const SVector3 &xaxis, const SVector3 &yaxis, const SVector3 &zaxis)
        {
            this->FromAxes(xaxis, yaxis, zaxis);
        }

        SQuaternion(const SVector3 *akAxis)
        {
            this->FromAxes(akAxis);
        }

        SQuaternion(SReal *valptr)
        {
            memcpy(&w, valptr, sizeof(SReal) * 4);
        }

        void swap(SQuaternion &other)
        {
            std::swap(w, other.w);
            std::swap(x, other.x);
            std::swap(y, other.y);
            std::swap(z, other.z);
        }

        SReal operator[](const size_t i) const
        {
            assert(i < 4);

            return *(&w + i);
        }

        SReal &operator[](const size_t i)
        {
            assert(i < 4);

            return *(&w + i);
        }

        SReal *ptr()
        {
            return &w;
        }

        const SReal *ptr() const
        {
            return &w;
        }

        void FromRotationMatrix(const SMatrix3 &kRot);
        void ToRotationMatrix(SMatrix3 &kRot) const;

        void FromAngleAxis(const SRadian &rfAngle, const SVector3 &rkAxis);
        void ToAngleAxis(SRadian &rfAngle, SVector3 &rkAxis) const;
        void ToAngleAxis(SDegree &dAngle, SVector3 &rkAxis) const
        {
            SRadian rAngle;
            ToAngleAxis(rAngle, rkAxis);
            dAngle = rAngle;
        }

        void FromAxes(const SVector3 *akAxis);
        void FromAxes(const SVector3 &xAxis, const SVector3 &yAxis, const SVector3 &zAxis);

        void ToAxes(SVector3 *akAxis) const;
        void ToAxes(SVector3 &xAxis, SVector3 &yAxis, SVector3 &zAxis) const;

        SVector3 xAxis(void) const;

        SVector3 yAxis(void) const;

        SVector3 zAxis(void) const;

        SQuaternion &operator=(const SQuaternion &rkQ)
        {
            w = rkQ.w;
            x = rkQ.x;
            y = rkQ.y;
            z = rkQ.z;
            return *this;
        }
        SQuaternion operator+(const SQuaternion &rkQ) const;
        SQuaternion operator-(const SQuaternion &rkQ) const;
        SQuaternion operator*(const SQuaternion &rkQ) const;
        SQuaternion operator*(SReal fScalar) const;
        friend SQuaternion operator*(SReal fScalar,
                                     const SQuaternion &rkQ);
        SQuaternion operator-() const;
        bool operator==(const SQuaternion &rhs) const
        {
            return (rhs.x == x) && (rhs.y == y) &&
                   (rhs.z == z) && (rhs.w == w);
        }
        bool operator!=(const SQuaternion &rhs) const
        {
            return !operator==(rhs);
        }

        SReal Dot(const SQuaternion &rkQ) const;

        SReal Norm() const;

        SReal normalise(void);
        SQuaternion Inverse() const;
        SQuaternion UnitInverse() const;
        SQuaternion Exp() const;
        SQuaternion Log() const;

        SVector3 operator*(const SVector3 &v) const;

        SRadian getRoll(bool reprojectAxis = true) const;

        SRadian getPitch(bool reprojectAxis = true) const;

        SRadian getYaw(bool reprojectAxis = true) const;

        bool equals(const SQuaternion &rhs, const SRadian &tolerance) const;

        static SQuaternion Slerp(SReal fT, const SQuaternion &rkP,
                                 const SQuaternion &rkQ, bool shortestPath = false);

        static SQuaternion SlerpExtraSpins(SReal fT,
                                           const SQuaternion &rkP, const SQuaternion &rkQ,
                                           int iExtraSpins);

        static void Intermediate(const SQuaternion &rkQ0,
                                 const SQuaternion &rkQ1, const SQuaternion &rkQ2,
                                 SQuaternion &rka, SQuaternion &rkB);

        static SQuaternion Squad(SReal fT, const SQuaternion &rkP,
                                 const SQuaternion &rkA, const SQuaternion &rkB,
                                 const SQuaternion &rkQ, bool shortestPath = false);

        static SQuaternion nlerp(SReal fT, const SQuaternion &rkP,
                                 const SQuaternion &rkQ, bool shortestPath = false);

        static const SReal msEpsilon;

        static const SQuaternion ZERO;
        static const SQuaternion IDENTITY;

        SReal w, x, y, z;

        bool isNaN() const;

        friend std::ostream &operator<<(std::ostream &o, const SQuaternion &q)
        {
            o << "SQuaternion(" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << ")";
            return o;
        }
    };
}

#endif