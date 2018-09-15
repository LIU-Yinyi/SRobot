#include "../SType.h"
#include "../SMath.h"

namespace SRobot
{
	const SReal SQuaternion::msEpsilon = (SReal)1e-03;
	const SQuaternion SQuaternion::ZERO(0, 0, 0, 0);
	const SQuaternion SQuaternion::IDENTITY(1, 0, 0, 0);

	//-----------------------------------------------------------------------
	void SQuaternion::FromRotationMatrix(const SMatrix3 &kRot)
	{
		SReal fTrace = kRot[0][0] + kRot[1][1] + kRot[2][2];
		SReal fRoot;

		if (fTrace > 0.0)
		{
			fRoot = SMath::Sqrt(fTrace + 1.0f);
			w = 0.5f * fRoot;
			fRoot = 0.5f / fRoot;
			x = (kRot[2][1] - kRot[1][2]) * fRoot;
			y = (kRot[0][2] - kRot[2][0]) * fRoot;
			z = (kRot[1][0] - kRot[0][1]) * fRoot;
		}
		else
		{
			static size_t s_iNext[3] = {1, 2, 0};
			size_t i = 0;
			if (kRot[1][1] > kRot[0][0])
				i = 1;
			if (kRot[2][2] > kRot[i][i])
				i = 2;
			size_t j = s_iNext[i];
			size_t k = s_iNext[j];

			fRoot = SMath::Sqrt(kRot[i][i] - kRot[j][j] - kRot[k][k] + 1.0f);
			SReal *apkQuat[3] = {&x, &y, &z};
			*apkQuat[i] = 0.5f * fRoot;
			fRoot = 0.5f / fRoot;
			w = (kRot[k][j] - kRot[j][k]) * fRoot;
			*apkQuat[j] = (kRot[j][i] + kRot[i][j]) * fRoot;
			*apkQuat[k] = (kRot[k][i] + kRot[i][k]) * fRoot;
		}
	}
	//-----------------------------------------------------------------------
	void SQuaternion::ToRotationMatrix(SMatrix3 &kRot) const
	{
		SReal fTx = x + x;
		SReal fTy = y + y;
		SReal fTz = z + z;
		SReal fTwx = fTx * w;
		SReal fTwy = fTy * w;
		SReal fTwz = fTz * w;
		SReal fTxx = fTx * x;
		SReal fTxy = fTy * x;
		SReal fTxz = fTz * x;
		SReal fTyy = fTy * y;
		SReal fTyz = fTz * y;
		SReal fTzz = fTz * z;

		kRot[0][0] = 1.0f - (fTyy + fTzz);
		kRot[0][1] = fTxy - fTwz;
		kRot[0][2] = fTxz + fTwy;
		kRot[1][0] = fTxy + fTwz;
		kRot[1][1] = 1.0f - (fTxx + fTzz);
		kRot[1][2] = fTyz - fTwx;
		kRot[2][0] = fTxz - fTwy;
		kRot[2][1] = fTyz + fTwx;
		kRot[2][2] = 1.0f - (fTxx + fTyy);
	}
	//-----------------------------------------------------------------------
	void SQuaternion::FromAngleAxis(const SRadian &rfAngle, const SVector3 &rkAxis)
	{
		SRadian fHalfAngle(0.5 * rfAngle);
		SReal fSin = SMath::Sin(fHalfAngle);
		w = SMath::Cos(fHalfAngle);
		x = fSin * rkAxis.x;
		y = fSin * rkAxis.y;
		z = fSin * rkAxis.z;
	}
	//-----------------------------------------------------------------------
	void SQuaternion::ToAngleAxis(SRadian &rfAngle, SVector3 &rkAxis) const
	{
		SReal fSqrLength = x * x + y * y + z * z;
		if (fSqrLength > 0.0)
		{
			rfAngle = SMath::ACos(w) * 2.0;
			SReal fInvLength = SMath::InvSqrt(fSqrLength);
			rkAxis.x = x * fInvLength;
			rkAxis.y = y * fInvLength;
			rkAxis.z = z * fInvLength;
		}
		else
		{
			rfAngle = SRadian(0.0);
			rkAxis.x = 1.0;
			rkAxis.y = 0.0;
			rkAxis.z = 0.0;
		}
	}
	//-----------------------------------------------------------------------
	void SQuaternion::FromAxes(const SVector3 *akAxis)
	{
		SMatrix3 kRot;

		for (size_t iCol = 0; iCol < 3; iCol++)
		{
			kRot[0][iCol] = akAxis[iCol].x;
			kRot[1][iCol] = akAxis[iCol].y;
			kRot[2][iCol] = akAxis[iCol].z;
		}

		FromRotationMatrix(kRot);
	}
	//-----------------------------------------------------------------------
	void SQuaternion::FromAxes(const SVector3 &xaxis, const SVector3 &yaxis, const SVector3 &zaxis)
	{
		SMatrix3 kRot;

		kRot[0][0] = xaxis.x;
		kRot[1][0] = xaxis.y;
		kRot[2][0] = xaxis.z;

		kRot[0][1] = yaxis.x;
		kRot[1][1] = yaxis.y;
		kRot[2][1] = yaxis.z;

		kRot[0][2] = zaxis.x;
		kRot[1][2] = zaxis.y;
		kRot[2][2] = zaxis.z;

		FromRotationMatrix(kRot);
	}
	//-----------------------------------------------------------------------
	void SQuaternion::ToAxes(SVector3 *akAxis) const
	{
		SMatrix3 kRot;

		ToRotationMatrix(kRot);

		for (size_t iCol = 0; iCol < 3; iCol++)
		{
			akAxis[iCol].x = kRot[0][iCol];
			akAxis[iCol].y = kRot[1][iCol];
			akAxis[iCol].z = kRot[2][iCol];
		}
	}
	//-----------------------------------------------------------------------
	SVector3 SQuaternion::xAxis(void) const
	{
		//SReal fTx  = 2.0*x;
		SReal fTy = 2.0f * y;
		SReal fTz = 2.0f * z;
		SReal fTwy = fTy * w;
		SReal fTwz = fTz * w;
		SReal fTxy = fTy * x;
		SReal fTxz = fTz * x;
		SReal fTyy = fTy * y;
		SReal fTzz = fTz * z;

		return SVector3(1.0f - (fTyy + fTzz), fTxy + fTwz, fTxz - fTwy);
	}
	//-----------------------------------------------------------------------
	SVector3 SQuaternion::yAxis(void) const
	{
		SReal fTx = 2.0f * x;
		SReal fTy = 2.0f * y;
		SReal fTz = 2.0f * z;
		SReal fTwx = fTx * w;
		SReal fTwz = fTz * w;
		SReal fTxx = fTx * x;
		SReal fTxy = fTy * x;
		SReal fTyz = fTz * y;
		SReal fTzz = fTz * z;

		return SVector3(fTxy - fTwz, 1.0f - (fTxx + fTzz), fTyz + fTwx);
	}
	//-----------------------------------------------------------------------
	SVector3 SQuaternion::zAxis(void) const
	{
		SReal fTx = 2.0f * x;
		SReal fTy = 2.0f * y;
		SReal fTz = 2.0f * z;
		SReal fTwx = fTx * w;
		SReal fTwy = fTy * w;
		SReal fTxx = fTx * x;
		SReal fTxz = fTz * x;
		SReal fTyy = fTy * y;
		SReal fTyz = fTz * y;

		return SVector3(fTxz + fTwy, fTyz - fTwx, 1.0f - (fTxx + fTyy));
	}
	//-----------------------------------------------------------------------
	void SQuaternion::ToAxes(SVector3 &xaxis, SVector3 &yaxis, SVector3 &zaxis) const
	{
		SMatrix3 kRot;

		ToRotationMatrix(kRot);

		xaxis.x = kRot[0][0];
		xaxis.y = kRot[1][0];
		xaxis.z = kRot[2][0];

		yaxis.x = kRot[0][1];
		yaxis.y = kRot[1][1];
		yaxis.z = kRot[2][1];

		zaxis.x = kRot[0][2];
		zaxis.y = kRot[1][2];
		zaxis.z = kRot[2][2];
	}

	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::operator+(const SQuaternion &rkQ) const
	{
		return SQuaternion(w + rkQ.w, x + rkQ.x, y + rkQ.y, z + rkQ.z);
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::operator-(const SQuaternion &rkQ) const
	{
		return SQuaternion(w - rkQ.w, x - rkQ.x, y - rkQ.y, z - rkQ.z);
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::operator*(const SQuaternion &rkQ) const
	{
		return SQuaternion(
			w * rkQ.w - x * rkQ.x - y * rkQ.y - z * rkQ.z,
			w * rkQ.x + x * rkQ.w + y * rkQ.z - z * rkQ.y,
			w * rkQ.y + y * rkQ.w + z * rkQ.x - x * rkQ.z,
			w * rkQ.z + z * rkQ.w + x * rkQ.y - y * rkQ.x);
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::operator*(SReal fScalar) const
	{
		return SQuaternion(fScalar * w, fScalar * x, fScalar * y, fScalar * z);
	}
	//-----------------------------------------------------------------------
	SQuaternion operator*(SReal fScalar, const SQuaternion &rkQ)
	{
		return SQuaternion(fScalar * rkQ.w, fScalar * rkQ.x, fScalar * rkQ.y,
						   fScalar * rkQ.z);
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::operator-() const
	{
		return SQuaternion(-w, -x, -y, -z);
	}
	//-----------------------------------------------------------------------
	SReal SQuaternion::Dot(const SQuaternion &rkQ) const
	{
		return w * rkQ.w + x * rkQ.x + y * rkQ.y + z * rkQ.z;
	}
	//-----------------------------------------------------------------------
	SReal SQuaternion::Norm() const
	{
		return w * w + x * x + y * y + z * z;
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::Inverse() const
	{
		SReal fNorm = w * w + x * x + y * y + z * z;
		if (fNorm > 0.0)
		{
			SReal fInvNorm = 1.0f / fNorm;
			return SQuaternion(w * fInvNorm, -x * fInvNorm, -y * fInvNorm, -z * fInvNorm);
		}
		else
		{
			return ZERO;
		}
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::UnitInverse() const
	{
		return SQuaternion(w, -x, -y, -z);
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::Exp() const
	{
		SRadian fAngle(SMath::Sqrt(x * x + y * y + z * z));
		SReal fSin = SMath::Sin(fAngle);

		SQuaternion kResult;
		kResult.w = SMath::Cos(fAngle);

		if (SMath::Abs(fSin) >= msEpsilon)
		{
			SReal fCoeff = fSin / (fAngle.valueRadians());
			kResult.x = fCoeff * x;
			kResult.y = fCoeff * y;
			kResult.z = fCoeff * z;
		}
		else
		{
			kResult.x = x;
			kResult.y = y;
			kResult.z = z;
		}

		return kResult;
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::Log() const
	{
		SQuaternion kResult;
		kResult.w = 0.0;

		if (SMath::Abs(w) < 1.0)
		{
			SRadian fAngle(SMath::ACos(w));
			SReal fSin = SMath::Sin(fAngle);
			if (SMath::Abs(fSin) >= msEpsilon)
			{
				SReal fCoeff = fAngle.valueRadians() / fSin;
				kResult.x = fCoeff * x;
				kResult.y = fCoeff * y;
				kResult.z = fCoeff * z;
				return kResult;
			}
		}

		kResult.x = x;
		kResult.y = y;
		kResult.z = z;

		return kResult;
	}
	//-----------------------------------------------------------------------
	SVector3 SQuaternion::operator*(const SVector3 &v) const
	{
		SVector3 uv, uuv;
		SVector3 qvec(x, y, z);
		uv = qvec.crossProduct(v);
		uuv = qvec.crossProduct(uv);
		uv *= (2.0f * w);
		uuv *= 2.0f;

		return v + uv + uuv;
	}
	//-----------------------------------------------------------------------
	bool SQuaternion::equals(const SQuaternion &rhs, const SRadian &tolerance) const
	{
		SReal fCos = Dot(rhs);
		SRadian angle = SMath::ACos(fCos);

		return (SMath::Abs(angle.valueRadians()) <= tolerance.valueRadians()) || SMath::RealEqual(angle.valueRadians(), SMath::PI, tolerance.valueRadians());
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::Slerp(SReal fT, const SQuaternion &rkP,
								   const SQuaternion &rkQ, bool shortestPath)
	{
		SReal fCos = rkP.Dot(rkQ);
		SQuaternion rkT;

		if (fCos < 0.0f && shortestPath)
		{
			fCos = -fCos;
			rkT = -rkQ;
		}
		else
		{
			rkT = rkQ;
		}

		if (SMath::Abs(fCos) < 1 - msEpsilon)
		{
			SReal fSin = SMath::Sqrt(1 - SMath::Sqr(fCos));
			SRadian fAngle = SMath::ATan2(fSin, fCos);
			SReal fInvSin = 1.0f / fSin;
			SReal fCoeff0 = SMath::Sin((1.0f - fT) * fAngle) * fInvSin;
			SReal fCoeff1 = SMath::Sin(fT * fAngle) * fInvSin;
			return fCoeff0 * rkP + fCoeff1 * rkT;
		}
		else
		{
			SQuaternion t = (1.0f - fT) * rkP + fT * rkT;

			t.normalise();
			return t;
		}
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::SlerpExtraSpins(SReal fT,
											 const SQuaternion &rkP, const SQuaternion &rkQ, int iExtraSpins)
	{
		SReal fCos = rkP.Dot(rkQ);
		SRadian fAngle(SMath::ACos(fCos));

		if (SMath::Abs(fAngle.valueRadians()) < msEpsilon)
			return rkP;

		SReal fSin = SMath::Sin(fAngle);
		SRadian fPhase(SMath::PI * iExtraSpins * fT);
		SReal fInvSin = 1.0f / fSin;
		SReal fCoeff0 = SMath::Sin((1.0f - fT) * fAngle - fPhase) * fInvSin;
		SReal fCoeff1 = SMath::Sin(fT * fAngle + fPhase) * fInvSin;
		return fCoeff0 * rkP + fCoeff1 * rkQ;
	}
	//-----------------------------------------------------------------------
	void SQuaternion::Intermediate(const SQuaternion &rkQ0,
								   const SQuaternion &rkQ1, const SQuaternion &rkQ2,
								   SQuaternion &rkA, SQuaternion &rkB)
	{

		SQuaternion kQ0inv = rkQ0.UnitInverse();
		SQuaternion kQ1inv = rkQ1.UnitInverse();
		SQuaternion rkP0 = kQ0inv * rkQ1;
		SQuaternion rkP1 = kQ1inv * rkQ2;
		SQuaternion kArg = 0.25 * (rkP0.Log() - rkP1.Log());
		SQuaternion kMinusArg = -kArg;

		rkA = rkQ1 * kArg.Exp();
		rkB = rkQ1 * kMinusArg.Exp();
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::Squad(SReal fT,
								   const SQuaternion &rkP, const SQuaternion &rkA,
								   const SQuaternion &rkB, const SQuaternion &rkQ, bool shortestPath)
	{
		SReal fSlerpT = 2.0f * fT * (1.0f - fT);
		SQuaternion kSlerpP = Slerp(fT, rkP, rkQ, shortestPath);
		SQuaternion kSlerpQ = Slerp(fT, rkA, rkB);
		return Slerp(fSlerpT, kSlerpP, kSlerpQ);
	}
	//-----------------------------------------------------------------------
	SReal SQuaternion::normalise(void)
	{
		SReal len = Norm();
		SReal factor = 1.0f / SMath::Sqrt(len);
		*this = *this *factor;
		return len;
	}
	//-----------------------------------------------------------------------
	SRadian SQuaternion::getRoll(bool reprojectAxis) const
	{
		if (reprojectAxis)
		{
			SReal fTy = 2.0f * y;
			SReal fTz = 2.0f * z;
			SReal fTwz = fTz * w;
			SReal fTxy = fTy * x;
			SReal fTyy = fTy * y;
			SReal fTzz = fTz * z;

			return SRadian(SMath::ATan2(fTxy + fTwz, 1.0f - (fTyy + fTzz)));
		}
		else
		{
			return SRadian(SMath::ATan2(2 * (x * y + w * z), w * w + x * x - y * y - z * z));
		}
	}
	//-----------------------------------------------------------------------
	SRadian SQuaternion::getPitch(bool reprojectAxis) const
	{
		if (reprojectAxis)
		{

			SReal fTx = 2.0f * x;

			SReal fTz = 2.0f * z;
			SReal fTwx = fTx * w;
			SReal fTxx = fTx * x;
			SReal fTyz = fTz * y;
			SReal fTzz = fTz * z;

			return SRadian(SMath::ATan2(fTyz + fTwx, 1.0f - (fTxx + fTzz)));
		}
		else
		{
			return SRadian(SMath::ATan2(2 * (y * z + w * x), w * w - x * x - y * y + z * z));
		}
	}
	//-----------------------------------------------------------------------
	SRadian SQuaternion::getYaw(bool reprojectAxis) const
	{
		if (reprojectAxis)
		{
			SReal fTx = 2.0f * x;
			SReal fTy = 2.0f * y;
			SReal fTz = 2.0f * z;
			SReal fTwy = fTy * w;
			SReal fTxx = fTx * x;
			SReal fTxz = fTz * x;
			SReal fTyy = fTy * y;

			return SRadian(SMath::ATan2(fTxz + fTwy, 1.0f - (fTxx + fTyy)));
		}
		else
		{
			return SRadian(SMath::ASin(-2 * (x * z - w * y)));
		}
	}
	//-----------------------------------------------------------------------
	SQuaternion SQuaternion::nlerp(SReal fT, const SQuaternion &rkP,
								   const SQuaternion &rkQ, bool shortestPath)
	{
		SQuaternion result;
		SReal fCos = rkP.Dot(rkQ);
		if (fCos < 0.0f && shortestPath)
		{
			result = rkP + fT * ((-rkQ) - rkP);
		}
		else
		{
			result = rkP + fT * (rkQ - rkP);
		}
		result.normalise();
		return result;
	}

	bool SQuaternion::isNaN() const
	{
		return SMath::isNaN(x) || SMath::isNaN(y) || SMath::isNaN(z) || SMath::isNaN(w);
	}
}
