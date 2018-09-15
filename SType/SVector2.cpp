#include "../SBase.h"
#include "../SType.h"
#include "../SMath.h"

namespace SRobot
{
	const SVector2 SVector2::ZERO(0, 0);

	const SVector2 SVector2::UNIT_X(1, 0);
	const SVector2 SVector2::UNIT_Y(0, 1);
	const SVector2 SVector2::NEGATIVE_UNIT_X(-1, 0);
	const SVector2 SVector2::NEGATIVE_UNIT_Y(0, -1);
	const SVector2 SVector2::UNIT_SCALE(1, 1);

	SVector2::SVector2(const SReal fX, const SReal fY)
		: x(fX), y(fY)
	{
	}

	SVector2::SVector2(const SReal scaler)
		: x(scaler), y(scaler)
	{
	}

	SVector2::SVector2(const SReal afCoordinate[2])
		: x(afCoordinate[0]),
		  y(afCoordinate[1])
	{
	}

	SVector2::SVector2(const int afCoordinate[2])
	{
		x = (SReal)afCoordinate[0];
		y = (SReal)afCoordinate[1];
	}

	SVector2::SVector2(SReal *const r)
		: x(r[0]), y(r[1])
	{
	}

	void SVector2::swap(SVector2 &other)
	{
		std::swap(x, other.x);
		std::swap(y, other.y);
	}

	SReal SVector2::operator[](const size_t i) const
	{
		assert(i < 2);

		return *(&x + i);
	}

	SReal &SVector2::operator[](const size_t i)
	{
		assert(i < 2);

		return *(&x + i);
	}

	SReal *SVector2::ptr()
	{
		return &x;
	}

	const SReal *SVector2::ptr() const
	{
		return &x;
	}

	SVector2 &SVector2::operator=(const SVector2 &rkSVector)
	{
		x = rkSVector.x;
		y = rkSVector.y;

		return *this;
	}

	SVector2 &SVector2::operator=(const SReal fScalar)
	{
		x = fScalar;
		y = fScalar;

		return *this;
	}

	bool SVector2::operator==(const SVector2 &rkSVector) const
	{
		return (x == rkSVector.x && y == rkSVector.y);
	}

	bool SVector2::operator!=(const SVector2 &rkSVector) const
	{
		return (x != rkSVector.x || y != rkSVector.y);
	}

	SVector2 SVector2::operator+(const SVector2 &rkSVector) const
	{
		return SVector2(x + rkSVector.x, y + rkSVector.y);
	}

	SVector2 SVector2::operator-(const SVector2 &rkSVector) const
	{
		return SVector2(x - rkSVector.x, y - rkSVector.y);
	}

	SVector2 SVector2::operator*(const SReal fScalar) const
	{
		return SVector2(x * fScalar, y * fScalar);
	}

	SVector2 SVector2::operator*(const SVector2 &rhs) const
	{
		return SVector2(x * rhs.x, y * rhs.y);
	}

	SVector2 SVector2::operator/(const SReal fScalar) const
	{
		assert(fScalar != 0.0);
		SReal fInv = 1.0f / fScalar;
		return SVector2(x * fInv, y * fInv);
	}

	SVector2 SVector2::operator/(const SVector2 &rhs) const
	{
		return SVector2(x / rhs.x, y / rhs.y);
	}

	const SVector2 &SVector2::operator+() const
	{
		return *this;
	}

	SVector2 SVector2::operator-() const
	{
		return SVector2(-x, -y);
	}

	SVector2 operator*(const SReal fScalar, const SVector2 &rkSVector)
	{
		return SVector2(fScalar * rkSVector.x, fScalar * rkSVector.y);
	}

	SVector2 operator/(const SReal fScalar, const SVector2 &rkSVector)
	{
		return SVector2(fScalar / rkSVector.x, fScalar / rkSVector.y);
	}

	SVector2 operator+(const SVector2 &lhs, const SReal rhs)
	{
		return SVector2(lhs.x + rhs, lhs.y + rhs);
	}

	SVector2 operator+(const SReal lhs, const SVector2 &rhs)
	{
		return SVector2(lhs + rhs.x, lhs + rhs.y);
	}

	SVector2 operator-(const SVector2 &lhs, const SReal rhs)
	{
		return SVector2(lhs.x - rhs, lhs.y - rhs);
	}

	SVector2 operator-(const SReal lhs, const SVector2 &rhs)
	{
		return SVector2(lhs - rhs.x, lhs - rhs.y);
	}

	SVector2 &SVector2::operator+=(const SVector2 &rkSVector)
	{
		x += rkSVector.x;
		y += rkSVector.y;

		return *this;
	}

	SVector2 &SVector2::operator+=(const SReal fScaler)
	{
		x += fScaler;
		y += fScaler;

		return *this;
	}

	SVector2 &SVector2::operator-=(const SVector2 &rkSVector)
	{
		x -= rkSVector.x;
		y -= rkSVector.y;

		return *this;
	}

	SVector2 &SVector2::operator-=(const SReal fScaler)
	{
		x -= fScaler;
		y -= fScaler;

		return *this;
	}

	SVector2 &SVector2::operator*=(const SReal fScalar)
	{
		x *= fScalar;
		y *= fScalar;

		return *this;
	}

	SVector2 &SVector2::operator*=(const SVector2 &rkSVector)
	{
		x *= rkSVector.x;
		y *= rkSVector.y;

		return *this;
	}

	SVector2 &SVector2::operator/=(const SReal fScalar)
	{
		assert(fScalar != 0.0);
		SReal fInv = 1.0f / fScalar;
		x *= fInv;
		y *= fInv;
		return *this;
	}

	SVector2 &SVector2::operator/=(const SVector2 &rkSVector)
	{
		x /= rkSVector.x;
		y /= rkSVector.y;
		return *this;
	}

	SReal SVector2::length() const
	{
		return SMath::Sqrt(x * x + y * y);
	}

	SReal SVector2::squaredLength() const
	{
		return x * x + y * y;
	}

	SReal SVector2::distance(const SVector2 &rhs) const
	{
		return (*this - rhs).length();
	}

	SReal SVector2::squaredDistance(const SVector2 &rhs) const
	{
		return (*this - rhs).squaredLength();
	}

	SReal SVector2::dotProduct(const SVector2 &vec) const
	{
		return x * vec.x + y * vec.y;
	}

	SReal SVector2::normalise()
	{
		SReal fLength = SMath::Sqrt(x * x + y * y);

		if (fLength > SReal(0.0f))
		{
			SReal fInvLength = 1.0f / fLength;
			x *= fInvLength;
			y *= fInvLength;
		}

		return fLength;
	}

	SVector2 SVector2::midPoint(const SVector2 &vec) const
	{
		return SVector2((x + vec.x) * 0.5f, (y + vec.y) * 0.5f);
	}

	bool SVector2::operator<(const SVector2 &rhs) const
	{
		if (x < rhs.x && y < rhs.y)
			return true;
		return false;
	}

	bool SVector2::operator>(const SVector2 &rhs) const
	{
		if (x > rhs.x && y > rhs.y)
			return true;
		return false;
	}

	void SVector2::makeFloor(const SVector2 &cmp)
	{
		if (cmp.x < x) x = cmp.x;
		if (cmp.y < y) y = cmp.y;
	}

	void SVector2::makeCeil(const SVector2 &cmp)
	{
		if (cmp.x > x) x = cmp.x;
		if (cmp.y > y) y = cmp.y;
	}

	SVector2 SVector2::perpendicular(void) const
	{
		return SVector2(-y, x);
	}

	SReal SVector2::crossProduct(const SVector2 &rkSVector) const
	{
		return x * rkSVector.y - y * rkSVector.x;
	}

	SVector2 SVector2::randomDeviant(SRadian angle) const
	{
		angle *= SMath::RangeRandom(-1, 1);
		SReal cosa = SMath::Cos(angle);
		SReal sina = SMath::Sin(angle);
		return SVector2(cosa * x - sina * y,
						sina * x + cosa * y);
	}

	bool SVector2::isZeroLength(void) const
	{
		SReal sqlen = (x * x) + (y * y);
		return (sqlen < (1e-06 * 1e-06));
	}

	SVector2 SVector2::normalisedCopy(void) const
	{
		SVector2 ret = *this;
		ret.normalise();
		return ret;
	}

	SVector2 SVector2::reflect(const SVector2 &normal) const
	{
		return SVector2(*this - (2 * this->dotProduct(normal) * normal));
	}

	bool SVector2::isNaN() const
	{
		return SMath::isNaN(x) || SMath::isNaN(y);
	}

	SRadian SVector2::angleBetween(const SVector2 &other) const
	{
		SReal lenProduct = length() * other.length();

		if (lenProduct < 1e-6f)
			lenProduct = 1e-6f;

		SReal f = dotProduct(other) / lenProduct;

		f = SMath::Clamp(f, (SReal)-1.0, (SReal)1.0);
		return SMath::ACos(f);
	}

	SRadian SVector2::angleTo(const SVector2 &other) const
	{
		SRadian angle = angleBetween(other);

		if (crossProduct(other) < 0)
			angle = (SRadian)SMath::TWO_PI - angle;

		return angle;
	}

	std::ostream &operator<<(std::ostream &o, const SVector2 &v)
	{
		o << "SVector2(" << v.x << ", " << v.y << ")";
		return o;
	}
}