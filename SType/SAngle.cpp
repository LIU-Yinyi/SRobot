#include "../SBase.h"
#include "../SType.h"
#include "../SMath.h"

namespace SRobot
{
	SRadian::SRadian(const SDegree &d)
		: mRad(d.valueRadians())
	{
	}
	SRadian &SRadian::operator=(const SDegree &d)
	{
		mRad = d.valueRadians();
		return *this;
	}
	SRadian SRadian::operator+(const SDegree &d) const
	{
		return SRadian(mRad + d.valueRadians());
	}
	SRadian &SRadian::operator+=(const SDegree &d)
	{
		mRad += d.valueRadians();
		return *this;
	}
	SRadian SRadian::operator-(const SDegree &d) const
	{
		return SRadian(mRad - d.valueRadians());
	}
	SRadian &SRadian::operator-=(const SDegree &d)
	{
		mRad -= d.valueRadians();
		return *this;
	}

	SReal SRadian::valueDegrees() const
	{
		return SMath::RadiansToDegrees(mRad);
	}

	SReal SDegree::valueRadians() const
	{
		return SMath::DegreesToRadians(mDeg);
	}
}
