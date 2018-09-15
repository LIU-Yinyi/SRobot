#include "../SBase.h"
#include "../SType.h"
#include "../SMath.h"

namespace SRobot
{
    const SVector3 SVector3::ZERO(0, 0, 0);

    const SVector3 SVector3::UNIT_X(1, 0, 0);
    const SVector3 SVector3::UNIT_Y(0, 1, 0);
    const SVector3 SVector3::UNIT_Z(0, 0, 1);
    const SVector3 SVector3::NEGATIVE_UNIT_X(-1, 0, 0);
    const SVector3 SVector3::NEGATIVE_UNIT_Y(0, -1, 0);
    const SVector3 SVector3::NEGATIVE_UNIT_Z(0, 0, -1);
    const SVector3 SVector3::UNIT_SCALE(1, 1, 1);

    SVector3::SVector3(const SReal fX, const SReal fY, const SReal fZ)
        : x(fX), y(fY), z(fZ)
    {
    }

    SVector3::SVector3(const SReal afCoordinate[3])
        : x(afCoordinate[0]),
          y(afCoordinate[1]),
          z(afCoordinate[2])
    {
    }

    SVector3::SVector3(const int afCoordinate[3])
    {
        x = (SReal)afCoordinate[0];
        y = (SReal)afCoordinate[1];
        z = (SReal)afCoordinate[2];
    }

    SVector3::SVector3(SReal *const r)
        : x(r[0]), y(r[1]), z(r[2])
    {
    }

    SVector3::SVector3(const SReal scaler)
        : x(scaler), y(scaler), z(scaler)
    {
    }

    void SVector3::swap(SVector3 &other)
    {
        std::swap(x, other.x);
        std::swap(y, other.y);
        std::swap(z, other.z);
    }

    SReal SVector3::operator[](const size_t i) const
    {
        assert(i < 3);

        return *(&x + i);
    }

    SReal &SVector3::operator[](const size_t i)
    {
        assert(i < 3);

        return *(&x + i);
    }

    SReal *SVector3::ptr()
    {
        return &x;
    }

    const SReal *SVector3::ptr() const
    {
        return &x;
    }

    SVector3 &SVector3::operator=(const SVector3 &rkSVector)
    {
        x = rkSVector.x;
        y = rkSVector.y;
        z = rkSVector.z;

        return *this;
    }

    SVector3 &SVector3::operator=(const SReal fScaler)
    {
        x = fScaler;
        y = fScaler;
        z = fScaler;

        return *this;
    }

    bool SVector3::operator==(const SVector3 &rkSVector) const
    {
        return (x == rkSVector.x && y == rkSVector.y && z == rkSVector.z);
    }

    bool SVector3::operator!=(const SVector3 &rkSVector) const
    {
        return (x != rkSVector.x || y != rkSVector.y || z != rkSVector.z);
    }

    SVector3 SVector3::operator+(const SVector3 &rkSVector) const
    {
        return SVector3(x + rkSVector.x, y + rkSVector.y, z + rkSVector.z);
    }

    SVector3 SVector3::operator-(const SVector3 &rkSVector) const
    {
        return SVector3(x - rkSVector.x, y - rkSVector.y, z - rkSVector.z);
    }

    SVector3 SVector3::operator*(const SReal fScalar) const
    {
        return SVector3(x * fScalar, y * fScalar, z * fScalar);
    }

    SVector3 SVector3::operator*(const SVector3 &rhs) const
    {
        return SVector3(x * rhs.x, y * rhs.y, z * rhs.z);
    }

    SVector3 SVector3::operator/(const SReal fScalar) const
    {
        assert(fScalar != 0.0);
        SReal fInv = 1.0f / fScalar;
        return SVector3(x * fInv, y * fInv, z * fInv);
    }

    SVector3 SVector3::operator/(const SVector3 &rhs) const
    {
        return SVector3(x / rhs.x, y / rhs.y, z / rhs.z);
    }

    const SVector3 &SVector3::operator+() const
    {
        return *this;
    }

    SVector3 SVector3::operator-() const
    {
        return SVector3(-x, -y, -z);
    }

    SVector3 operator*(const SReal fScalar, const SVector3 &rkSVector)
    {
        return SVector3(fScalar * rkSVector.x, fScalar * rkSVector.y, fScalar * rkSVector.z);
    }

    SVector3 operator/(const SReal fScalar, const SVector3 &rkSVector)
    {
        return SVector3(fScalar / rkSVector.x, fScalar / rkSVector.y, fScalar / rkSVector.z);
    }

    SVector3 operator+(const SVector3 &lhs, const SReal rhs)
    {
        return SVector3(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
    }

    SVector3 operator+(const SReal lhs, const SVector3 &rhs)
    {
        return SVector3(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z);
    }

    SVector3 operator-(const SVector3 &lhs, const SReal rhs)
    {
        return SVector3(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs);
    }

    SVector3 operator-(const SReal lhs, const SVector3 &rhs)
    {
        return SVector3(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z);
    }

    SVector3 &SVector3::operator+=(const SVector3 &rkSVector)
    {
        x += rkSVector.x;
        y += rkSVector.y;
        z += rkSVector.z;

        return *this;
    }

    SVector3 &SVector3::operator+=(const SReal fScalar)
    {
        x += fScalar;
        y += fScalar;
        z += fScalar;
        return *this;
    }

    SVector3 &SVector3::operator-=(const SVector3 &rkSVector)
    {
        x -= rkSVector.x;
        y -= rkSVector.y;
        z -= rkSVector.z;

        return *this;
    }

    SVector3 &SVector3::operator-=(const SReal fScalar)
    {
        x -= fScalar;
        y -= fScalar;
        z -= fScalar;
        return *this;
    }

    SVector3 &SVector3::operator*=(const SReal fScalar)
    {
        x *= fScalar;
        y *= fScalar;
        z *= fScalar;
        return *this;
    }

    SVector3 &SVector3::operator*=(const SVector3 &rkSVector)
    {
        x *= rkSVector.x;
        y *= rkSVector.y;
        z *= rkSVector.z;

        return *this;
    }

    SVector3 &SVector3::operator/=(const SReal fScalar)
    {
        assert(fScalar != 0.0);

        SReal fInv = 1.0f / fScalar;

        x *= fInv;
        y *= fInv;
        z *= fInv;

        return *this;
    }

    SVector3 &SVector3::operator/=(const SVector3 &rkSVector)
    {
        x /= rkSVector.x;
        y /= rkSVector.y;
        z /= rkSVector.z;

        return *this;
    }

    SReal SVector3::length() const
    {
        return SMath::Sqrt(x * x + y * y + z * z);
    }

    SReal SVector3::squaredLength() const
    {
        return x * x + y * y + z * z;
    }

    SReal SVector3::distance(const SVector3 &rhs) const
    {
        return (*this - rhs).length();
    }

    SReal SVector3::squaredDistance(const SVector3 &rhs) const
    {
        return (*this - rhs).squaredLength();
    }

    SReal SVector3::dotProduct(const SVector3 &vec) const
    {
        return x * vec.x + y * vec.y + z * vec.z;
    }

    SReal SVector3::absDotProduct(const SVector3 &vec) const
    {
        return SMath::Abs(x * vec.x) + SMath::Abs(y * vec.y) + SMath::Abs(z * vec.z);
    }

    SReal SVector3::normalise()
    {
        SReal fLength = SMath::Sqrt(x * x + y * y + z * z);

        if (fLength > SReal(0.0f))
        {
            SReal fInvLength = 1.0f / fLength;
            x *= fInvLength;
            y *= fInvLength;
            z *= fInvLength;
        }

        return fLength;
    }

    SVector3 SVector3::crossProduct(const SVector3 &rkSVector) const
    {
        return SVector3(
            y * rkSVector.z - z * rkSVector.y,
            z * rkSVector.x - x * rkSVector.z,
            x * rkSVector.y - y * rkSVector.x);
    }

    SVector3 SVector3::midPoint(const SVector3 &vec) const
    {
        return SVector3(
            (x + vec.x) * 0.5f,
            (y + vec.y) * 0.5f,
            (z + vec.z) * 0.5f);
    }

    bool SVector3::operator<(const SVector3 &rhs) const
    {
        if (x < rhs.x && y < rhs.y && z < rhs.z)
            return true;
        return false;
    }

    bool SVector3::operator>(const SVector3 &rhs) const
    {
        if (x > rhs.x && y > rhs.y && z > rhs.z)
            return true;
        return false;
    }

    void SVector3::makeFloor(const SVector3 &cmp)
    {
        if (cmp.x < x) x = cmp.x;
        if (cmp.y < y) y = cmp.y;
        if (cmp.z < z) z = cmp.z;
    }

    void SVector3::makeCeil(const SVector3 &cmp)
    {
        if (cmp.x > x) x = cmp.x;
        if (cmp.y > y) y = cmp.y;
        if (cmp.z > z) z = cmp.z;
    }

    SVector3 SVector3::perpendicular(void) const
    {
        static const SReal fSquareZero = (SReal)(1e-06 * 1e-06);

        SVector3 perp = this->crossProduct(SVector3::UNIT_X);

        if (perp.squaredLength() < fSquareZero)
        {
            perp = this->crossProduct(SVector3::UNIT_Y);
        }
        perp.normalise();

        return perp;
    }

    SVector3 SVector3::randomDeviant(const SRadian &angle, const SVector3 &up) const
    {
        SVector3 newUp;

        if (up == SVector3::ZERO)
        {
            newUp = this->perpendicular();
        }
        else
        {
            newUp = up;
        }

        SQuaternion q;
        q.FromAngleAxis(SRadian(SMath::UnitRandom() * SMath::TWO_PI), *this);
        newUp = q * newUp;

        q.FromAngleAxis(angle, newUp);
        return q * (*this);
    }

    SRadian SVector3::angleBetween(const SVector3 &dest) const
    {
        SReal lenProduct = length() * dest.length();

        if (lenProduct < 1e-6f)
            lenProduct = 1e-6f;

        SReal f = dotProduct(dest) / lenProduct;

        f = SMath::Clamp(f, (SReal)-1.0, (SReal)1.0);
        return SMath::ACos(f);
    }

    SQuaternion SVector3::getRotationTo(const SVector3 &dest, const SVector3 &fallbackAxis) const
    {
        SQuaternion q;

        SVector3 v0 = *this;
        SVector3 v1 = dest;
        v0.normalise();
        v1.normalise();

        SReal d = v0.dotProduct(v1);

        if (d >= 1.0f)
        {
            return SQuaternion::IDENTITY;
        }
        if (d < (1e-6f - 1.0f))
        {
            if (fallbackAxis != SVector3::ZERO)
            {
                // rotate 180 degrees about the fallback axis
                q.FromAngleAxis(SRadian(SMath::PI), fallbackAxis);
            }
            else
            {
                // Generate an axis
                SVector3 axis = SVector3::UNIT_X.crossProduct(*this);
                if (axis.isZeroLength()) // pick another if colinear
                    axis = SVector3::UNIT_Y.crossProduct(*this);
                axis.normalise();
                q.FromAngleAxis(SRadian(SMath::PI), axis);
            }
        }
        else
        {
            SReal s = SMath::Sqrt((1 + d) * 2);
            SReal invs = 1 / s;

            SVector3 c = v0.crossProduct(v1);

            q.x = c.x * invs;
            q.y = c.y * invs;
            q.z = c.z * invs;
            q.w = s * 0.5f;
            q.normalise();
        }
        return q;
    }

    bool SVector3::isZeroLength(void) const
    {
        SReal sqlen = (x * x) + (y * y) + (z * z);
        return (sqlen < (1e-06 * 1e-06));
    }

    SVector3 SVector3::normalisedCopy(void) const
    {
        SVector3 ret = *this;
        ret.normalise();
        return ret;
    }

    SVector3 SVector3::reflect(const SVector3 &normal) const
    {
        return SVector3(*this - (2 * this->dotProduct(normal) * normal));
    }

    bool SVector3::positionEquals(const SVector3 &rhs, SReal tolerance) const
    {
        return SMath::RealEqual(x, rhs.x, tolerance) &&
               SMath::RealEqual(y, rhs.y, tolerance) &&
               SMath::RealEqual(z, rhs.z, tolerance);
    }

    bool SVector3::positionCloses(const SVector3 &rhs, SReal tolerance) const
    {
        return squaredDistance(rhs) <=
               (squaredLength() + rhs.squaredLength()) * tolerance;
    }

    bool SVector3::directionEquals(const SVector3 &rhs,
                                   const SRadian &tolerance) const
    {
        SReal dot = dotProduct(rhs);
        SRadian angle = SMath::ACos(dot);

        return SMath::Abs(angle.valueRadians()) <= tolerance.valueRadians();
    }

    bool SVector3::isNaN() const
    {
        return SMath::isNaN(x) || SMath::isNaN(y) || SMath::isNaN(z);
    }

    SVector3 SVector3::primaryAxis() const
    {
        SReal absx = SMath::Abs(x);
        SReal absy = SMath::Abs(y);
        SReal absz = SMath::Abs(z);
        if (absx > absy)
            if (absx > absz)
                return x > 0 ? SVector3::UNIT_X : SVector3::NEGATIVE_UNIT_X;
            else
                return z > 0 ? SVector3::UNIT_Z : SVector3::NEGATIVE_UNIT_Z;
        else // absx <= absy
            if (absy > absz)
            return y > 0 ? SVector3::UNIT_Y : SVector3::NEGATIVE_UNIT_Y;
        else
            return z > 0 ? SVector3::UNIT_Z : SVector3::NEGATIVE_UNIT_Z;
    }

	SVector2 SVector3::compressZ() const
	{
		return SVector2(x, y);
	}

    std::ostream &operator<<(std::ostream &o, const SVector3 &v)
    {
        o << "SVector3(" << v.x << ", " << v.y << ", " << v.z << ")";
        return o;
    }
}
