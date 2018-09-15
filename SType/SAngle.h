#ifndef _SANGLE_H_
#define _SANGLE_H_

#include "STypePreDefine.h"

namespace SRobot
{
	//弧度
	class SRadian
	{
	  private:
		SReal mRad;

	  public:
		explicit SRadian(SReal r = 0)
			: mRad(r) {}
		SRadian(const SDegree &d);

		operator SReal()
		{
			return valueRadians();
		}

		SRadian &operator=(const SReal &f)
		{
			mRad = f;
			return *this;
		}
		SRadian &operator=(const SRadian &r)
		{
			mRad = r.mRad;
			return *this;
		}
		SRadian &operator=(const SDegree &d);

		SReal valueDegrees() const;
		SReal valueRadians() const { return mRad; }

		const SRadian &operator+() const { return *this; }
		SRadian operator+(const SRadian &r) const { return SRadian(mRad + r.mRad); }
		SRadian operator+(const SDegree &d) const;
		SRadian &operator+=(const SRadian &r)
		{
			mRad += r.mRad;
			return *this;
		}
		SRadian &operator+=(const SDegree &d);

		SRadian operator-() const { return SRadian(-mRad); }
		SRadian operator-(const SRadian &r) const { return SRadian(mRad - r.mRad); }
		SRadian operator-(const SDegree &d) const;
		SRadian &operator-=(const SRadian &r)
		{
			mRad -= r.mRad;
			return *this;
		}
		SRadian &operator-=(const SDegree &d);

		SRadian operator*(SReal f) const { return SRadian(mRad * f); }
		SRadian operator*(const SRadian &r) const { return SRadian(mRad * r.mRad); }
		SRadian &operator*=(SReal f)
		{
			mRad *= f;
			return *this;
		}

		SRadian operator/(SReal f) const { return SRadian(mRad / f); }
		SRadian &operator/=(SReal f)
		{
			mRad /= f;
			return *this;
		}

		bool operator<(const SRadian &r) const { return mRad < r.mRad; }
		bool operator<=(const SRadian &r) const { return mRad <= r.mRad; }
		bool operator==(const SRadian &r) const { return mRad == r.mRad; }
		bool operator!=(const SRadian &r) const { return mRad != r.mRad; }
		bool operator>=(const SRadian &r) const { return mRad >= r.mRad; }
		bool operator>(const SRadian &r) const { return mRad > r.mRad; }

		friend SRadian operator*(SReal f, const SRadian &r)
		{
			return SRadian(f * r.mRad);
		}

		friend std::ostream &operator<<(std::ostream &o, const SRadian &v)
		{
			o << "Radian(" << v.valueRadians() << ")";
			return o;
		}
	};

	//角度
	class SDegree
	{
	  private:
		SReal mDeg;

	  public:
		explicit SDegree(SReal d = 0)
			: mDeg(d) {}
		SDegree(const SRadian &r)
			: mDeg(r.valueDegrees()) {}
		SDegree &operator=(const SReal &f)
		{
			mDeg = f;
			return *this;
		}
		SDegree &operator=(const SDegree &d)
		{
			mDeg = d.mDeg;
			return *this;
		}
		SDegree &operator=(const SRadian &r)
		{
			mDeg = r.valueDegrees();
			return *this;
		}


		SReal valueDegrees() const { return mDeg; }
		SReal valueRadians() const;

		const SDegree &operator+() const { return *this; }
		SDegree operator+(const SDegree &d) const { return SDegree(mDeg + d.mDeg); }
		SDegree operator+(const SRadian &r) const { return SDegree(mDeg + r.valueDegrees()); }
		SDegree &operator+=(const SDegree &d)
		{
			mDeg += d.mDeg;
			return *this;
		}
		SDegree &operator+=(const SRadian &r)
		{
			mDeg += r.valueDegrees();
			return *this;
		}

		SDegree operator-() const { return SDegree(-mDeg); }
		SDegree operator-(const SDegree &d) const { return SDegree(mDeg - d.mDeg); }
		SDegree operator-(const SRadian &r) const { return SDegree(mDeg - r.valueDegrees()); }
		SDegree &operator-=(const SDegree &d)
		{
			mDeg -= d.mDeg;
			return *this;
		}
		SDegree &operator-=(const SRadian &r)
		{
			mDeg -= r.valueDegrees();
			return *this;
		}

		SDegree operator*(SReal f) const { return SDegree(mDeg * f); }
		SDegree operator*(const SDegree &f) const { return SDegree(mDeg * f.mDeg); }
		SDegree &operator*=(SReal f)
		{
			mDeg *= f;
			return *this;
		}

		SDegree operator/(SReal f) const { return SDegree(mDeg / f); }
		SDegree &operator/=(SReal f)
		{
			mDeg /= f;
			return *this;
		}

		bool operator<(const SDegree &d) const { return mDeg < d.mDeg; }
		bool operator<=(const SDegree &d) const { return mDeg <= d.mDeg; }
		bool operator==(const SDegree &d) const { return mDeg == d.mDeg; }
		bool operator!=(const SDegree &d) const { return mDeg != d.mDeg; }
		bool operator>=(const SDegree &d) const { return mDeg >= d.mDeg; }
		bool operator>(const SDegree &d) const { return mDeg > d.mDeg; }

		friend SDegree operator*(SReal f, const SDegree &d)
		{
			return SDegree(f * d.mDeg);
		}

		friend std::ostream &operator<<(std::ostream &o, const SDegree &v)
		{
			o << "Degree(" << v.valueDegrees() << ")";
			return o;
		}
	};

	//角(弧度)
	class SAngle : public SRadian
	{
	  public:
		SAngle(SReal angle)
			: SRadian(angle) {}

		SAngle() : SRadian(0){}
		SAngle(SRadian angle)
			: SRadian(angle) {}


	};

	class SEulerAngle
	{
		public:

		SAngle roll;
		SAngle pitch;
		SAngle yaw;

		SEulerAngle(SReal roll=0, SReal pitch=0, SReal yaw=0)
		{
			this->roll = roll;
			this->pitch = pitch;
			this->yaw = yaw;
		}


		friend std::ostream &operator<<(std::ostream &o, const SEulerAngle &e)
		{
			o << "EulerAngle(roll:" << e.roll.valueDegrees() << ",pitch:" << e.pitch.valueDegrees() << ",yaw:" << e.yaw.valueDegrees() <<")";
			return o;
		}

	};
}

#endif
