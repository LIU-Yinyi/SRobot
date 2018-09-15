#include "../SBase.h"
#include "../SType.h"
#include "../SMath.h"

namespace SRobot
{
    const SReal SMatrix3::EPSILON = (SReal)1e-06;
    const SMatrix3 SMatrix3::ZERO(0, 0, 0, 0, 0, 0, 0, 0, 0);
    const SMatrix3 SMatrix3::IDENTITY(1, 0, 0, 0, 1, 0, 0, 0, 1);
    const SReal SMatrix3::msSvdEpsilon = (SReal)1e-04;
    const unsigned int SMatrix3::msSvdMaxIterations = 32;

    //-----------------------------------------------------------------------
    SVector3 SMatrix3::GetColumn(size_t iCol) const
    {
        assert(iCol < 3);
        return SVector3(m[0][iCol], m[1][iCol], m[2][iCol]);
    }
    //-----------------------------------------------------------------------
    void SMatrix3::SetColumn(size_t iCol, const SVector3 &vec)
    {
        assert(iCol < 3);
        m[0][iCol] = vec.x;
        m[1][iCol] = vec.y;
        m[2][iCol] = vec.z;
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromAxes(const SVector3 &xAxis, const SVector3 &yAxis, const SVector3 &zAxis)
    {
        SetColumn(0, xAxis);
        SetColumn(1, yAxis);
        SetColumn(2, zAxis);
    }

    //-----------------------------------------------------------------------
    bool SMatrix3::operator==(const SMatrix3 &rkSMatrix) const
    {
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
            {
                if (m[iRow][iCol] != rkSMatrix.m[iRow][iCol])
                    return false;
            }
        }

        return true;
    }
    //-----------------------------------------------------------------------
    SMatrix3 SMatrix3::operator+(const SMatrix3 &rkSMatrix) const
    {
        SMatrix3 kSum;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
            {
                kSum.m[iRow][iCol] = m[iRow][iCol] +
                                     rkSMatrix.m[iRow][iCol];
            }
        }
        return kSum;
    }
    //-----------------------------------------------------------------------
    SMatrix3 SMatrix3::operator-(const SMatrix3 &rkSMatrix) const
    {
        SMatrix3 kDiff;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
            {
                kDiff.m[iRow][iCol] = m[iRow][iCol] - rkSMatrix.m[iRow][iCol];
            }
        }
        return kDiff;
    }
    //-----------------------------------------------------------------------
    SMatrix3 SMatrix3::operator*(const SMatrix3 &rkSMatrix) const
    {
        SMatrix3 kProd;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
            {
                kProd.m[iRow][iCol] =
                    m[iRow][0] * rkSMatrix.m[0][iCol] +
                    m[iRow][1] * rkSMatrix.m[1][iCol] +
                    m[iRow][2] * rkSMatrix.m[2][iCol];
            }
        }
        return kProd;
    }
    //-----------------------------------------------------------------------
    SVector3 SMatrix3::operator*(const SVector3 &rkPoint) const
    {
        SVector3 kProd;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            kProd[iRow] =
                m[iRow][0] * rkPoint[0] +
                m[iRow][1] * rkPoint[1] +
                m[iRow][2] * rkPoint[2];
        }
        return kProd;
    }
    //-----------------------------------------------------------------------
    SVector3 operator*(const SVector3 &rkPoint, const SMatrix3 &rkSMatrix)
    {
        SVector3 kProd;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            kProd[iRow] =
                rkPoint[0] * rkSMatrix.m[0][iRow] +
                rkPoint[1] * rkSMatrix.m[1][iRow] +
                rkPoint[2] * rkSMatrix.m[2][iRow];
        }
        return kProd;
    }
    //-----------------------------------------------------------------------
    SMatrix3 SMatrix3::operator-() const
    {
        SMatrix3 kNeg;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
                kNeg[iRow][iCol] = -m[iRow][iCol];
        }
        return kNeg;
    }
    //-----------------------------------------------------------------------
    SMatrix3 SMatrix3::operator*(SReal fScalar) const
    {
        SMatrix3 kProd;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
                kProd[iRow][iCol] = fScalar * m[iRow][iCol];
        }
        return kProd;
    }
    //-----------------------------------------------------------------------
    SMatrix3 operator*(SReal fScalar, const SMatrix3 &rkSMatrix)
    {
        SMatrix3 kProd;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
                kProd[iRow][iCol] = fScalar * rkSMatrix.m[iRow][iCol];
        }
        return kProd;
    }
    //-----------------------------------------------------------------------
    SMatrix3 SMatrix3::Transpose() const
    {
        SMatrix3 kTranspose;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
                kTranspose[iRow][iCol] = m[iCol][iRow];
        }
        return kTranspose;
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::Inverse(SMatrix3 &rkInverse, SReal fTolerance) const
    {
        rkInverse[0][0] = m[1][1] * m[2][2] -
                          m[1][2] * m[2][1];
        rkInverse[0][1] = m[0][2] * m[2][1] -
                          m[0][1] * m[2][2];
        rkInverse[0][2] = m[0][1] * m[1][2] -
                          m[0][2] * m[1][1];
        rkInverse[1][0] = m[1][2] * m[2][0] -
                          m[1][0] * m[2][2];
        rkInverse[1][1] = m[0][0] * m[2][2] -
                          m[0][2] * m[2][0];
        rkInverse[1][2] = m[0][2] * m[1][0] -
                          m[0][0] * m[1][2];
        rkInverse[2][0] = m[1][0] * m[2][1] -
                          m[1][1] * m[2][0];
        rkInverse[2][1] = m[0][1] * m[2][0] -
                          m[0][0] * m[2][1];
        rkInverse[2][2] = m[0][0] * m[1][1] -
                          m[0][1] * m[1][0];

        SReal fDet = m[0][0] * rkInverse[0][0] +
                     m[0][1] * rkInverse[1][0] +
                     m[0][2] * rkInverse[2][0];

        if (SMath::Abs(fDet) <= fTolerance)
            return false;

        SReal fInvDet = 1.0f / fDet;
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
                rkInverse[iRow][iCol] *= fInvDet;
        }

        return true;
    }
    //-----------------------------------------------------------------------
    SMatrix3 SMatrix3::Inverse(SReal fTolerance) const
    {
        SMatrix3 kInverse = SMatrix3::ZERO;
        Inverse(kInverse, fTolerance);
        return kInverse;
    }
    //-----------------------------------------------------------------------
    SReal SMatrix3::Determinant() const
    {
        SReal fCofactor00 = m[1][1] * m[2][2] -
                            m[1][2] * m[2][1];
        SReal fCofactor10 = m[1][2] * m[2][0] -
                            m[1][0] * m[2][2];
        SReal fCofactor20 = m[1][0] * m[2][1] -
                            m[1][1] * m[2][0];

        SReal fDet =
            m[0][0] * fCofactor00 +
            m[0][1] * fCofactor10 +
            m[0][2] * fCofactor20;

        return fDet;
    }
    //-----------------------------------------------------------------------
    void SMatrix3::Bidiagonalize(SMatrix3 &kA, SMatrix3 &kL,
                                 SMatrix3 &kR)
    {
        SReal afV[3], afW[3];
        SReal fLength, fSign, fT1, fInvT1, fT2;
        bool bIdentity;

        fLength = SMath::Sqrt(kA[0][0] * kA[0][0] +
                              kA[1][0] * kA[1][0] +
                              kA[2][0] * kA[2][0]);
        if (fLength > 0.0)
        {
            fSign = (kA[0][0] > 0.0f ? 1.0f : -1.0f);
            fT1 = kA[0][0] + fSign * fLength;
            fInvT1 = 1.0f / fT1;
            afV[1] = kA[1][0] * fInvT1;
            afV[2] = kA[2][0] * fInvT1;

            fT2 = -2.0f / (1.0f + afV[1] * afV[1] + afV[2] * afV[2]);
            afW[0] = fT2 * (kA[0][0] + kA[1][0] * afV[1] + kA[2][0] * afV[2]);
            afW[1] = fT2 * (kA[0][1] + kA[1][1] * afV[1] + kA[2][1] * afV[2]);
            afW[2] = fT2 * (kA[0][2] + kA[1][2] * afV[1] + kA[2][2] * afV[2]);
            kA[0][0] += afW[0];
            kA[0][1] += afW[1];
            kA[0][2] += afW[2];
            kA[1][1] += afV[1] * afW[1];
            kA[1][2] += afV[1] * afW[2];
            kA[2][1] += afV[2] * afW[1];
            kA[2][2] += afV[2] * afW[2];

            kL[0][0] = 1.0f + fT2;
            kL[0][1] = kL[1][0] = fT2 * afV[1];
            kL[0][2] = kL[2][0] = fT2 * afV[2];
            kL[1][1] = 1.0f + fT2 * afV[1] * afV[1];
            kL[1][2] = kL[2][1] = fT2 * afV[1] * afV[2];
            kL[2][2] = 1.0f + fT2 * afV[2] * afV[2];
            bIdentity = false;
        }
        else
        {
            kL = SMatrix3::IDENTITY;
            bIdentity = true;
        }

        fLength = SMath::Sqrt(kA[0][1] * kA[0][1] + kA[0][2] * kA[0][2]);
        if (fLength > 0.0)
        {
            fSign = (kA[0][1] > 0.0f ? 1.0f : -1.0f);
            fT1 = kA[0][1] + fSign * fLength;
            afV[2] = kA[0][2] / fT1;

            fT2 = -2.0f / (1.0f + afV[2] * afV[2]);
            afW[0] = fT2 * (kA[0][1] + kA[0][2] * afV[2]);
            afW[1] = fT2 * (kA[1][1] + kA[1][2] * afV[2]);
            afW[2] = fT2 * (kA[2][1] + kA[2][2] * afV[2]);
            kA[0][1] += afW[0];
            kA[1][1] += afW[1];
            kA[1][2] += afW[1] * afV[2];
            kA[2][1] += afW[2];
            kA[2][2] += afW[2] * afV[2];

            kR[0][0] = 1.0;
            kR[0][1] = kR[1][0] = 0.0;
            kR[0][2] = kR[2][0] = 0.0;
            kR[1][1] = 1.0f + fT2;
            kR[1][2] = kR[2][1] = fT2 * afV[2];
            kR[2][2] = 1.0f + fT2 * afV[2] * afV[2];
        }
        else
        {
            kR = SMatrix3::IDENTITY;
        }

        fLength = SMath::Sqrt(kA[1][1] * kA[1][1] + kA[2][1] * kA[2][1]);
        if (fLength > 0.0)
        {
            fSign = (kA[1][1] > 0.0f ? 1.0f : -1.0f);
            fT1 = kA[1][1] + fSign * fLength;
            afV[2] = kA[2][1] / fT1;

            fT2 = -2.0f / (1.0f + afV[2] * afV[2]);
            afW[1] = fT2 * (kA[1][1] + kA[2][1] * afV[2]);
            afW[2] = fT2 * (kA[1][2] + kA[2][2] * afV[2]);
            kA[1][1] += afW[1];
            kA[1][2] += afW[2];
            kA[2][2] += afV[2] * afW[2];

            SReal fA = 1.0f + fT2;
            SReal fB = fT2 * afV[2];
            SReal fC = 1.0f + fB * afV[2];

            if (bIdentity)
            {
                kL[0][0] = 1.0;
                kL[0][1] = kL[1][0] = 0.0;
                kL[0][2] = kL[2][0] = 0.0;
                kL[1][1] = fA;
                kL[1][2] = kL[2][1] = fB;
                kL[2][2] = fC;
            }
            else
            {
                for (int iRow = 0; iRow < 3; iRow++)
                {
                    SReal fTmp0 = kL[iRow][1];
                    SReal fTmp1 = kL[iRow][2];
                    kL[iRow][1] = fA * fTmp0 + fB * fTmp1;
                    kL[iRow][2] = fB * fTmp0 + fC * fTmp1;
                }
            }
        }
    }
    //-----------------------------------------------------------------------
    void SMatrix3::GolubKahanStep(SMatrix3 &kA, SMatrix3 &kL,
                                  SMatrix3 &kR)
    {
        SReal fT11 = kA[0][1] * kA[0][1] + kA[1][1] * kA[1][1];
        SReal fT22 = kA[1][2] * kA[1][2] + kA[2][2] * kA[2][2];
        SReal fT12 = kA[1][1] * kA[1][2];
        SReal fTrace = fT11 + fT22;
        SReal fDiff = fT11 - fT22;
        SReal fDiscr = SMath::Sqrt(fDiff * fDiff + 4.0f * fT12 * fT12);
        SReal fRoot1 = 0.5f * (fTrace + fDiscr);
        SReal fRoot2 = 0.5f * (fTrace - fDiscr);

        SReal fY = kA[0][0] - (SMath::Abs(fRoot1 - fT22) <=
                                       SMath::Abs(fRoot2 - fT22)
                                   ? fRoot1
                                   : fRoot2);
        SReal fZ = kA[0][1];
        SReal fInvLength = SMath::InvSqrt(fY * fY + fZ * fZ);
        SReal fSin = fZ * fInvLength;
        SReal fCos = -fY * fInvLength;

        SReal fTmp0 = kA[0][0];
        SReal fTmp1 = kA[0][1];
        kA[0][0] = fCos * fTmp0 - fSin * fTmp1;
        kA[0][1] = fSin * fTmp0 + fCos * fTmp1;
        kA[1][0] = -fSin * kA[1][1];
        kA[1][1] *= fCos;

        size_t iRow;
        for (iRow = 0; iRow < 3; iRow++)
        {
            fTmp0 = kR[0][iRow];
            fTmp1 = kR[1][iRow];
            kR[0][iRow] = fCos * fTmp0 - fSin * fTmp1;
            kR[1][iRow] = fSin * fTmp0 + fCos * fTmp1;
        }

        fY = kA[0][0];
        fZ = kA[1][0];
        fInvLength = SMath::InvSqrt(fY * fY + fZ * fZ);
        fSin = fZ * fInvLength;
        fCos = -fY * fInvLength;

        kA[0][0] = fCos * kA[0][0] - fSin * kA[1][0];
        fTmp0 = kA[0][1];
        fTmp1 = kA[1][1];
        kA[0][1] = fCos * fTmp0 - fSin * fTmp1;
        kA[1][1] = fSin * fTmp0 + fCos * fTmp1;
        kA[0][2] = -fSin * kA[1][2];
        kA[1][2] *= fCos;

        size_t iCol;
        for (iCol = 0; iCol < 3; iCol++)
        {
            fTmp0 = kL[iCol][0];
            fTmp1 = kL[iCol][1];
            kL[iCol][0] = fCos * fTmp0 - fSin * fTmp1;
            kL[iCol][1] = fSin * fTmp0 + fCos * fTmp1;
        }

        fY = kA[0][1];
        fZ = kA[0][2];
        fInvLength = SMath::InvSqrt(fY * fY + fZ * fZ);
        fSin = fZ * fInvLength;
        fCos = -fY * fInvLength;

        kA[0][1] = fCos * kA[0][1] - fSin * kA[0][2];
        fTmp0 = kA[1][1];
        fTmp1 = kA[1][2];
        kA[1][1] = fCos * fTmp0 - fSin * fTmp1;
        kA[1][2] = fSin * fTmp0 + fCos * fTmp1;
        kA[2][1] = -fSin * kA[2][2];
        kA[2][2] *= fCos;

        for (iRow = 0; iRow < 3; iRow++)
        {
            fTmp0 = kR[1][iRow];
            fTmp1 = kR[2][iRow];
            kR[1][iRow] = fCos * fTmp0 - fSin * fTmp1;
            kR[2][iRow] = fSin * fTmp0 + fCos * fTmp1;
        }

        fY = kA[1][1];
        fZ = kA[2][1];
        fInvLength = SMath::InvSqrt(fY * fY + fZ * fZ);
        fSin = fZ * fInvLength;
        fCos = -fY * fInvLength;

        kA[1][1] = fCos * kA[1][1] - fSin * kA[2][1];
        fTmp0 = kA[1][2];
        fTmp1 = kA[2][2];
        kA[1][2] = fCos * fTmp0 - fSin * fTmp1;
        kA[2][2] = fSin * fTmp0 + fCos * fTmp1;

        for (iCol = 0; iCol < 3; iCol++)
        {
            fTmp0 = kL[iCol][1];
            fTmp1 = kL[iCol][2];
            kL[iCol][1] = fCos * fTmp0 - fSin * fTmp1;
            kL[iCol][2] = fSin * fTmp0 + fCos * fTmp1;
        }
    }
    //-----------------------------------------------------------------------
    void SMatrix3::SingularValueDecomposition(SMatrix3 &kL, SVector3 &kS,
                                              SMatrix3 &kR) const
    {
        size_t iRow, iCol;

        SMatrix3 kA = *this;
        Bidiagonalize(kA, kL, kR);

        for (unsigned int i = 0; i < msSvdMaxIterations; i++)
        {
            SReal fTmp, fTmp0, fTmp1;
            SReal fSin0, fCos0, fTan0;
            SReal fSin1, fCos1, fTan1;

            bool bTest1 = (SMath::Abs(kA[0][1]) <=
                           msSvdEpsilon * (SMath::Abs(kA[0][0]) + SMath::Abs(kA[1][1])));
            bool bTest2 = (SMath::Abs(kA[1][2]) <=
                           msSvdEpsilon * (SMath::Abs(kA[1][1]) + SMath::Abs(kA[2][2])));
            if (bTest1)
            {
                if (bTest2)
                {
                    kS[0] = kA[0][0];
                    kS[1] = kA[1][1];
                    kS[2] = kA[2][2];
                    break;
                }
                else
                {
                    fTmp = (kA[1][1] * kA[1][1] - kA[2][2] * kA[2][2] +
                            kA[1][2] * kA[1][2]) /
                           (kA[1][2] * kA[2][2]);
                    fTan0 = 0.5f * (fTmp + SMath::Sqrt(fTmp * fTmp + 4.0f));
                    fCos0 = SMath::InvSqrt(1.0f + fTan0 * fTan0);
                    fSin0 = fTan0 * fCos0;

                    for (iCol = 0; iCol < 3; iCol++)
                    {
                        fTmp0 = kL[iCol][1];
                        fTmp1 = kL[iCol][2];
                        kL[iCol][1] = fCos0 * fTmp0 - fSin0 * fTmp1;
                        kL[iCol][2] = fSin0 * fTmp0 + fCos0 * fTmp1;
                    }

                    fTan1 = (kA[1][2] - kA[2][2] * fTan0) / kA[1][1];
                    fCos1 = SMath::InvSqrt(1.0f + fTan1 * fTan1);
                    fSin1 = -fTan1 * fCos1;

                    for (iRow = 0; iRow < 3; iRow++)
                    {
                        fTmp0 = kR[1][iRow];
                        fTmp1 = kR[2][iRow];
                        kR[1][iRow] = fCos1 * fTmp0 - fSin1 * fTmp1;
                        kR[2][iRow] = fSin1 * fTmp0 + fCos1 * fTmp1;
                    }

                    kS[0] = kA[0][0];
                    kS[1] = fCos0 * fCos1 * kA[1][1] -
                            fSin1 * (fCos0 * kA[1][2] - fSin0 * kA[2][2]);
                    kS[2] = fSin0 * fSin1 * kA[1][1] +
                            fCos1 * (fSin0 * kA[1][2] + fCos0 * kA[2][2]);
                    break;
                }
            }
            else
            {
                if (bTest2)
                {
                    fTmp = (kA[0][0] * kA[0][0] + kA[1][1] * kA[1][1] -
                            kA[0][1] * kA[0][1]) /
                           (kA[0][1] * kA[1][1]);
                    fTan0 = 0.5f * (-fTmp + SMath::Sqrt(fTmp * fTmp + 4.0f));
                    fCos0 = SMath::InvSqrt(1.0f + fTan0 * fTan0);
                    fSin0 = fTan0 * fCos0;

                    for (iCol = 0; iCol < 3; iCol++)
                    {
                        fTmp0 = kL[iCol][0];
                        fTmp1 = kL[iCol][1];
                        kL[iCol][0] = fCos0 * fTmp0 - fSin0 * fTmp1;
                        kL[iCol][1] = fSin0 * fTmp0 + fCos0 * fTmp1;
                    }

                    fTan1 = (kA[0][1] - kA[1][1] * fTan0) / kA[0][0];
                    fCos1 = SMath::InvSqrt(1.0f + fTan1 * fTan1);
                    fSin1 = -fTan1 * fCos1;

                    for (iRow = 0; iRow < 3; iRow++)
                    {
                        fTmp0 = kR[0][iRow];
                        fTmp1 = kR[1][iRow];
                        kR[0][iRow] = fCos1 * fTmp0 - fSin1 * fTmp1;
                        kR[1][iRow] = fSin1 * fTmp0 + fCos1 * fTmp1;
                    }

                    kS[0] = fCos0 * fCos1 * kA[0][0] -
                            fSin1 * (fCos0 * kA[0][1] - fSin0 * kA[1][1]);
                    kS[1] = fSin0 * fSin1 * kA[0][0] +
                            fCos1 * (fSin0 * kA[0][1] + fCos0 * kA[1][1]);
                    kS[2] = kA[2][2];
                    break;
                }
                else
                {
                    GolubKahanStep(kA, kL, kR);
                }
            }
        }

        for (iRow = 0; iRow < 3; iRow++)
        {
            if (kS[iRow] < 0.0)
            {
                kS[iRow] = -kS[iRow];
                for (iCol = 0; iCol < 3; iCol++)
                    kR[iRow][iCol] = -kR[iRow][iCol];
            }
        }
    }
    //-----------------------------------------------------------------------
    void SMatrix3::SingularValueComposition(const SMatrix3 &kL,
                                            const SVector3 &kS, const SMatrix3 &kR)
    {
        size_t iRow, iCol;
        SMatrix3 kTmp;

        for (iRow = 0; iRow < 3; iRow++)
        {
            for (iCol = 0; iCol < 3; iCol++)
                kTmp[iRow][iCol] = kS[iRow] * kR[iRow][iCol];
        }

        for (iRow = 0; iRow < 3; iRow++)
        {
            for (iCol = 0; iCol < 3; iCol++)
            {
                m[iRow][iCol] = 0.0;
                for (int iMid = 0; iMid < 3; iMid++)
                    m[iRow][iCol] += kL[iRow][iMid] * kTmp[iMid][iCol];
            }
        }
    }
    //-----------------------------------------------------------------------
    void SMatrix3::Orthonormalize()
    {
        SReal fInvLength = SMath::InvSqrt(m[0][0] * m[0][0] + m[1][0] * m[1][0] +
                                          m[2][0] * m[2][0]);

        m[0][0] *= fInvLength;
        m[1][0] *= fInvLength;
        m[2][0] *= fInvLength;

        SReal fDot0 =
            m[0][0] * m[0][1] +
            m[1][0] * m[1][1] +
            m[2][0] * m[2][1];

        m[0][1] -= fDot0 * m[0][0];
        m[1][1] -= fDot0 * m[1][0];
        m[2][1] -= fDot0 * m[2][0];

        fInvLength = SMath::InvSqrt(m[0][1] * m[0][1] +
                                    m[1][1] * m[1][1] +
                                    m[2][1] * m[2][1]);

        m[0][1] *= fInvLength;
        m[1][1] *= fInvLength;
        m[2][1] *= fInvLength;

        SReal fDot1 =
            m[0][1] * m[0][2] +
            m[1][1] * m[1][2] +
            m[2][1] * m[2][2];

        fDot0 =
            m[0][0] * m[0][2] +
            m[1][0] * m[1][2] +
            m[2][0] * m[2][2];

        m[0][2] -= fDot0 * m[0][0] + fDot1 * m[0][1];
        m[1][2] -= fDot0 * m[1][0] + fDot1 * m[1][1];
        m[2][2] -= fDot0 * m[2][0] + fDot1 * m[2][1];

        fInvLength = SMath::InvSqrt(m[0][2] * m[0][2] +
                                    m[1][2] * m[1][2] +
                                    m[2][2] * m[2][2]);

        m[0][2] *= fInvLength;
        m[1][2] *= fInvLength;
        m[2][2] *= fInvLength;
    }
    //-----------------------------------------------------------------------
    void SMatrix3::QDUDecomposition(SMatrix3 &kQ,
                                    SVector3 &kD, SVector3 &kU) const
    {

        SReal fInvLength = m[0][0] * m[0][0] + m[1][0] * m[1][0] + m[2][0] * m[2][0];
        if (SMath::RealEqual(fInvLength, 0)) fInvLength = SMath::InvSqrt(fInvLength);

        kQ[0][0] = m[0][0] * fInvLength;
        kQ[1][0] = m[1][0] * fInvLength;
        kQ[2][0] = m[2][0] * fInvLength;

        SReal fDot = kQ[0][0] * m[0][1] + kQ[1][0] * m[1][1] +
                     kQ[2][0] * m[2][1];
        kQ[0][1] = m[0][1] - fDot * kQ[0][0];
        kQ[1][1] = m[1][1] - fDot * kQ[1][0];
        kQ[2][1] = m[2][1] - fDot * kQ[2][0];
        fInvLength = kQ[0][1] * kQ[0][1] + kQ[1][1] * kQ[1][1] + kQ[2][1] * kQ[2][1];
        if (SMath::RealEqual(fInvLength, 0)) fInvLength = SMath::InvSqrt(fInvLength);

        kQ[0][1] *= fInvLength;
        kQ[1][1] *= fInvLength;
        kQ[2][1] *= fInvLength;

        fDot = kQ[0][0] * m[0][2] + kQ[1][0] * m[1][2] +
               kQ[2][0] * m[2][2];
        kQ[0][2] = m[0][2] - fDot * kQ[0][0];
        kQ[1][2] = m[1][2] - fDot * kQ[1][0];
        kQ[2][2] = m[2][2] - fDot * kQ[2][0];
        fDot = kQ[0][1] * m[0][2] + kQ[1][1] * m[1][2] +
               kQ[2][1] * m[2][2];
        kQ[0][2] -= fDot * kQ[0][1];
        kQ[1][2] -= fDot * kQ[1][1];
        kQ[2][2] -= fDot * kQ[2][1];
        fInvLength = kQ[0][2] * kQ[0][2] + kQ[1][2] * kQ[1][2] + kQ[2][2] * kQ[2][2];
        if (SMath::RealEqual(fInvLength, 0)) fInvLength = SMath::InvSqrt(fInvLength);

        kQ[0][2] *= fInvLength;
        kQ[1][2] *= fInvLength;
        kQ[2][2] *= fInvLength;

        SReal fDet = kQ[0][0] * kQ[1][1] * kQ[2][2] + kQ[0][1] * kQ[1][2] * kQ[2][0] +
                     kQ[0][2] * kQ[1][0] * kQ[2][1] - kQ[0][2] * kQ[1][1] * kQ[2][0] -
                     kQ[0][1] * kQ[1][0] * kQ[2][2] - kQ[0][0] * kQ[1][2] * kQ[2][1];

        if (fDet < 0.0)
        {
            for (size_t iRow = 0; iRow < 3; iRow++)
                for (size_t iCol = 0; iCol < 3; iCol++)
                    kQ[iRow][iCol] = -kQ[iRow][iCol];
        }

        SMatrix3 kR;
        kR[0][0] = kQ[0][0] * m[0][0] + kQ[1][0] * m[1][0] +
                   kQ[2][0] * m[2][0];
        kR[0][1] = kQ[0][0] * m[0][1] + kQ[1][0] * m[1][1] +
                   kQ[2][0] * m[2][1];
        kR[1][1] = kQ[0][1] * m[0][1] + kQ[1][1] * m[1][1] +
                   kQ[2][1] * m[2][1];
        kR[0][2] = kQ[0][0] * m[0][2] + kQ[1][0] * m[1][2] +
                   kQ[2][0] * m[2][2];
        kR[1][2] = kQ[0][1] * m[0][2] + kQ[1][1] * m[1][2] +
                   kQ[2][1] * m[2][2];
        kR[2][2] = kQ[0][2] * m[0][2] + kQ[1][2] * m[1][2] +
                   kQ[2][2] * m[2][2];

        kD[0] = kR[0][0];
        kD[1] = kR[1][1];
        kD[2] = kR[2][2];

        SReal fInvD0 = 1.0f / kD[0];
        kU[0] = kR[0][1] * fInvD0;
        kU[1] = kR[0][2] * fInvD0;
        kU[2] = kR[1][2] / kD[1];
    }
    //-----------------------------------------------------------------------
    SReal SMatrix3::MaxCubicRoot(SReal afCoeff[3])
    {
        const SReal fOneThird = (SReal)1.0 / (SReal)3.0;
        const SReal fEpsilon = (SReal)1e-06;
        SReal fDiscr = afCoeff[2] * afCoeff[2] - 3.0f * afCoeff[1];
        if (fDiscr <= fEpsilon)
            return -fOneThird * afCoeff[2];

        SReal fX = 1.0;
        SReal fPoly = afCoeff[0] + fX * (afCoeff[1] + fX * (afCoeff[2] + fX));
        if (fPoly < 0.0)
        {
            fX = SMath::Abs(afCoeff[0]);
            SReal fTmp = 1.0f + SMath::Abs(afCoeff[1]);
            if (fTmp > fX)
                fX = fTmp;
            fTmp = 1.0f + SMath::Abs(afCoeff[2]);
            if (fTmp > fX)
                fX = fTmp;
        }

        SReal fTwoC2 = 2.0f * afCoeff[2];
        for (int i = 0; i < 16; i++)
        {
            fPoly = afCoeff[0] + fX * (afCoeff[1] + fX * (afCoeff[2] + fX));
            if (SMath::Abs(fPoly) <= fEpsilon)
                return fX;

            SReal fDeriv = afCoeff[1] + fX * (fTwoC2 + 3.0f * fX);
            fX -= fPoly / fDeriv;
        }

        return fX;
    }
    //-----------------------------------------------------------------------
    SReal SMatrix3::SpectralNorm() const
    {
        SMatrix3 kP;
        size_t iRow, iCol;
        SReal fPmax = 0.0;
        for (iRow = 0; iRow < 3; iRow++)
        {
            for (iCol = 0; iCol < 3; iCol++)
            {
                kP[iRow][iCol] = 0.0;
                for (int iMid = 0; iMid < 3; iMid++)
                {
                    kP[iRow][iCol] +=
                        m[iMid][iRow] * m[iMid][iCol];
                }
                if (kP[iRow][iCol] > fPmax)
                    fPmax = kP[iRow][iCol];
            }
        }

        SReal fInvPmax = 1.0f / fPmax;
        for (iRow = 0; iRow < 3; iRow++)
        {
            for (iCol = 0; iCol < 3; iCol++)
                kP[iRow][iCol] *= fInvPmax;
        }

        SReal afCoeff[3];
        afCoeff[0] = -(kP[0][0] * (kP[1][1] * kP[2][2] - kP[1][2] * kP[2][1]) +
                       kP[0][1] * (kP[2][0] * kP[1][2] - kP[1][0] * kP[2][2]) +
                       kP[0][2] * (kP[1][0] * kP[2][1] - kP[2][0] * kP[1][1]));
        afCoeff[1] = kP[0][0] * kP[1][1] - kP[0][1] * kP[1][0] +
                     kP[0][0] * kP[2][2] - kP[0][2] * kP[2][0] +
                     kP[1][1] * kP[2][2] - kP[1][2] * kP[2][1];
        afCoeff[2] = -(kP[0][0] + kP[1][1] + kP[2][2]);

        SReal fRoot = MaxCubicRoot(afCoeff);
        SReal fNorm = SMath::Sqrt(fPmax * fRoot);
        return fNorm;
    }
    //-----------------------------------------------------------------------
    void SMatrix3::ToAngleAxis(SVector3 &rkAxis, SRadian &rfSRadians) const
    {
        SReal fTrace = m[0][0] + m[1][1] + m[2][2];
        SReal fCos = 0.5f * (fTrace - 1.0f);
        rfSRadians = SMath::ACos(fCos); // in [0,PI]

        if (rfSRadians > SRadian(0.0))
        {
            if (rfSRadians < SRadian(SMath::PI))
            {
                rkAxis.x = m[2][1] - m[1][2];
                rkAxis.y = m[0][2] - m[2][0];
                rkAxis.z = m[1][0] - m[0][1];
                rkAxis.normalise();
            }
            else
            {
                float fHalfInverse;
                if (m[0][0] >= m[1][1])
                {
                    if (m[0][0] >= m[2][2])
                    {
                        rkAxis.x = 0.5f * SMath::Sqrt(m[0][0] -
                                                      m[1][1] - m[2][2] + 1.0f);
                        fHalfInverse = 0.5f / rkAxis.x;
                        rkAxis.y = fHalfInverse * m[0][1];
                        rkAxis.z = fHalfInverse * m[0][2];
                    }
                    else
                    {
                        rkAxis.z = 0.5f * SMath::Sqrt(m[2][2] -
                                                      m[0][0] - m[1][1] + 1.0f);
                        fHalfInverse = 0.5f / rkAxis.z;
                        rkAxis.x = fHalfInverse * m[0][2];
                        rkAxis.y = fHalfInverse * m[1][2];
                    }
                }
                else
                {
                    if (m[1][1] >= m[2][2])
                    {
                        rkAxis.y = 0.5f * SMath::Sqrt(m[1][1] -
                                                      m[0][0] - m[2][2] + 1.0f);
                        fHalfInverse = 0.5f / rkAxis.y;
                        rkAxis.x = fHalfInverse * m[0][1];
                        rkAxis.z = fHalfInverse * m[1][2];
                    }
                    else
                    {
                        rkAxis.z = 0.5f * SMath::Sqrt(m[2][2] -
                                                      m[0][0] - m[1][1] + 1.0f);
                        fHalfInverse = 0.5f / rkAxis.z;
                        rkAxis.x = fHalfInverse * m[0][2];
                        rkAxis.y = fHalfInverse * m[1][2];
                    }
                }
            }
        }
        else
        {
            rkAxis.x = 1.0;
            rkAxis.y = 0.0;
            rkAxis.z = 0.0;
        }
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromAngleAxis(const SVector3 &rkAxis, const SRadian &fSRadians)
    {
        SReal fCos = SMath::Cos(fSRadians);
        SReal fSin = SMath::Sin(fSRadians);
        SReal fOneMinusCos = 1.0f - fCos;
        SReal fX2 = rkAxis.x * rkAxis.x;
        SReal fY2 = rkAxis.y * rkAxis.y;
        SReal fZ2 = rkAxis.z * rkAxis.z;
        SReal fXYM = rkAxis.x * rkAxis.y * fOneMinusCos;
        SReal fXZM = rkAxis.x * rkAxis.z * fOneMinusCos;
        SReal fYZM = rkAxis.y * rkAxis.z * fOneMinusCos;
        SReal fXSin = rkAxis.x * fSin;
        SReal fYSin = rkAxis.y * fSin;
        SReal fZSin = rkAxis.z * fSin;

        m[0][0] = fX2 * fOneMinusCos + fCos;
        m[0][1] = fXYM - fZSin;
        m[0][2] = fXZM + fYSin;
        m[1][0] = fXYM + fZSin;
        m[1][1] = fY2 * fOneMinusCos + fCos;
        m[1][2] = fYZM - fXSin;
        m[2][0] = fXZM - fYSin;
        m[2][1] = fYZM + fXSin;
        m[2][2] = fZ2 * fOneMinusCos + fCos;
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::ToEulerAnglesXYZ(SRadian &rfYAngle, SRadian &rfPAngle,
                                    SRadian &rfRAngle) const
    {
        rfPAngle = SRadian(SMath::ASin(m[0][2]));
        if (rfPAngle < SRadian(SMath::HALF_PI))
        {
            if (rfPAngle > SRadian(-SMath::HALF_PI))
            {
                rfYAngle = SMath::ATan2(-m[1][2], m[2][2]);
                rfRAngle = SMath::ATan2(-m[0][1], m[0][0]);
                return true;
            }
            else
            {
                SRadian fRmY = SMath::ATan2(m[1][0], m[1][1]);
                rfRAngle = SRadian(0.0);
                rfYAngle = rfRAngle - fRmY;
                return false;
            }
        }
        else
        {
            SRadian fRpY = SMath::ATan2(m[1][0], m[1][1]);
            rfRAngle = SRadian(0.0);
            rfYAngle = fRpY - rfRAngle;
            return false;
        }
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::ToEulerAnglesXZY(SRadian &rfYAngle, SRadian &rfPAngle,
                                    SRadian &rfRAngle) const
    {
        rfPAngle = SMath::ASin(-m[0][1]);
        if (rfPAngle < SRadian(SMath::HALF_PI))
        {
            if (rfPAngle > SRadian(-SMath::HALF_PI))
            {
                rfYAngle = SMath::ATan2(m[2][1], m[1][1]);
                rfRAngle = SMath::ATan2(m[0][2], m[0][0]);
                return true;
            }
            else
            {
                SRadian fRmY = SMath::ATan2(-m[2][0], m[2][2]);
                rfRAngle = SRadian(0.0);
                rfYAngle = rfRAngle - fRmY;
                return false;
            }
        }
        else
        {
            SRadian fRpY = SMath::ATan2(-m[2][0], m[2][2]);
            rfRAngle = SRadian(0.0);
            rfYAngle = fRpY - rfRAngle;
            return false;
        }
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::ToEulerAnglesYXZ(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const
    {

        rfPAngle = SMath::ASin(-m[1][2]);
        if (rfPAngle < SRadian(SMath::HALF_PI))
        {
            if (rfPAngle > SRadian(-SMath::HALF_PI))
            {
                rfYAngle = SMath::ATan2(m[0][2], m[2][2]);
                rfRAngle = SMath::ATan2(m[1][0], m[1][1]);
                return true;
            }
            else
            {
                SRadian fRmY = SMath::ATan2(-m[0][1], m[0][0]);
                rfRAngle = SRadian(0.0);
                rfYAngle = rfRAngle - fRmY;
                return false;
            }
        }
        else
        {
            SRadian fRpY = SMath::ATan2(-m[0][1], m[0][0]);
            rfRAngle = SRadian(0.0);
            rfYAngle = fRpY - rfRAngle;
            return false;
        }
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::ToEulerAnglesYZX(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const
    {
        rfPAngle = SMath::ASin(m[1][0]);
        if (rfPAngle < SRadian(SMath::HALF_PI))
        {
            if (rfPAngle > SRadian(-SMath::HALF_PI))
            {
                rfYAngle = SMath::ATan2(-m[2][0], m[0][0]);
                rfRAngle = SMath::ATan2(-m[1][2], m[1][1]);
                return true;
            }
            else
            {
                SRadian fRmY = SMath::ATan2(m[2][1], m[2][2]);
                rfRAngle = SRadian(0.0);
                rfYAngle = rfRAngle - fRmY;
                return false;
            }
        }
        else
        {
            SRadian fRpY = SMath::ATan2(m[2][1], m[2][2]);
            rfRAngle = SRadian(0.0);
            rfYAngle = fRpY - rfRAngle;
            return false;
        }
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::ToEulerAnglesZXY(SRadian &rfYAngle, SRadian &rfPAngle,
                                    SRadian &rfRAngle) const
    {
        rfPAngle = SMath::ASin(m[2][1]);
        if (rfPAngle < SRadian(SMath::HALF_PI))
        {
            if (rfPAngle > SRadian(-SMath::HALF_PI))
            {
                rfYAngle = SMath::ATan2(-m[0][1], m[1][1]);
                rfRAngle = SMath::ATan2(-m[2][0], m[2][2]);
                return true;
            }
            else
            {
                SRadian fRmY = SMath::ATan2(m[0][2], m[0][0]);
                rfRAngle = SRadian(0.0);
                rfYAngle = rfRAngle - fRmY;
                return false;
            }
        }
        else
        {
            SRadian fRpY = SMath::ATan2(m[0][2], m[0][0]);
            rfRAngle = SRadian(0.0);
            rfYAngle = fRpY - rfRAngle;
            return false;
        }
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::ToEulerAnglesZYX(SRadian &rfYAngle, SRadian &rfPAngle,
                                    SRadian &rfRAngle) const
    {
        rfPAngle = SMath::ASin(-m[2][0]);
        if (rfPAngle < SRadian(SMath::HALF_PI))
        {
            if (rfPAngle > SRadian(-SMath::HALF_PI))
            {
                rfYAngle = SMath::ATan2(m[1][0], m[0][0]);
                rfRAngle = SMath::ATan2(m[2][1], m[2][2]);
                return true;
            }
            else
            {
                SRadian fRmY = SMath::ATan2(-m[0][1], m[0][2]);
                rfRAngle = SRadian(0.0);
                rfYAngle = rfRAngle - fRmY;
                return false;
            }
        }
        else
        {
            SRadian fRpY = SMath::ATan2(-m[0][1], m[0][2]);
            rfRAngle = SRadian(0.0);
            rfYAngle = fRpY - rfRAngle;
            return false;
        }
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromEulerAnglesXYZ(const SRadian &fYAngle, const SRadian &fPAngle,
                                      const SRadian &fRAngle)
    {
        SReal fCos, fSin;

        fCos = SMath::Cos(fYAngle);
        fSin = SMath::Sin(fYAngle);
        SMatrix3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

        fCos = SMath::Cos(fPAngle);
        fSin = SMath::Sin(fPAngle);
        SMatrix3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

        fCos = SMath::Cos(fRAngle);
        fSin = SMath::Sin(fRAngle);
        SMatrix3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

        *this = kXMat *(kYMat * kZMat);
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromEulerAnglesXZY(const SRadian &fYAngle, const SRadian &fPAngle,
                                      const SRadian &fRAngle)
    {
        SReal fCos, fSin;

        fCos = SMath::Cos(fYAngle);
        fSin = SMath::Sin(fYAngle);
        SMatrix3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

        fCos = SMath::Cos(fPAngle);
        fSin = SMath::Sin(fPAngle);
        SMatrix3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

        fCos = SMath::Cos(fRAngle);
        fSin = SMath::Sin(fRAngle);
        SMatrix3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

        *this = kXMat *(kZMat * kYMat);
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromEulerAnglesYXZ(const SRadian &fYAngle, const SRadian &fPAngle,
                                      const SRadian &fRAngle)
    {
        SReal fCos, fSin;

        fCos = SMath::Cos(fYAngle);
        fSin = SMath::Sin(fYAngle);
        SMatrix3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

        fCos = SMath::Cos(fPAngle);
        fSin = SMath::Sin(fPAngle);
        SMatrix3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

        fCos = SMath::Cos(fRAngle);
        fSin = SMath::Sin(fRAngle);
        SMatrix3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

        *this = kYMat *(kXMat * kZMat);
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromEulerAnglesYZX(const SRadian &fYAngle, const SRadian &fPAngle,
                                      const SRadian &fRAngle)
    {
        SReal fCos, fSin;

        fCos = SMath::Cos(fYAngle);
        fSin = SMath::Sin(fYAngle);
        SMatrix3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

        fCos = SMath::Cos(fPAngle);
        fSin = SMath::Sin(fPAngle);
        SMatrix3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

        fCos = SMath::Cos(fRAngle);
        fSin = SMath::Sin(fRAngle);
        SMatrix3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

        *this = kYMat *(kZMat * kXMat);
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromEulerAnglesZXY(const SRadian &fYAngle, const SRadian &fPAngle,
                                      const SRadian &fRAngle)
    {
        SReal fCos, fSin;

        fCos = SMath::Cos(fYAngle);
        fSin = SMath::Sin(fYAngle);
        SMatrix3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

        fCos = SMath::Cos(fPAngle);
        fSin = SMath::Sin(fPAngle);
        SMatrix3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

        fCos = SMath::Cos(fRAngle);
        fSin = SMath::Sin(fRAngle);
        SMatrix3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

        *this = kZMat *(kXMat * kYMat);
    }
    //-----------------------------------------------------------------------
    void SMatrix3::FromEulerAnglesZYX(const SRadian &fYAngle, const SRadian &fPAngle,
                                      const SRadian &fRAngle)
    {
        SReal fCos, fSin;

        fCos = SMath::Cos(fYAngle);
        fSin = SMath::Sin(fYAngle);
        SMatrix3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

        fCos = SMath::Cos(fPAngle);
        fSin = SMath::Sin(fPAngle);
        SMatrix3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

        fCos = SMath::Cos(fRAngle);
        fSin = SMath::Sin(fRAngle);
        SMatrix3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

        *this = kZMat *(kYMat * kXMat);
    }
    //-----------------------------------------------------------------------
    void SMatrix3::Tridiagonal(SReal afDiag[3], SReal afSubDiag[3])
    {
        SReal fA = m[0][0];
        SReal fB = m[0][1];
        SReal fC = m[0][2];
        SReal fD = m[1][1];
        SReal fE = m[1][2];
        SReal fF = m[2][2];

        afDiag[0] = fA;
        afSubDiag[2] = 0.0;
        if (SMath::Abs(fC) >= EPSILON)
        {
            SReal fLength = SMath::Sqrt(fB * fB + fC * fC);
            SReal fInvLength = 1.0f / fLength;
            fB *= fInvLength;
            fC *= fInvLength;
            SReal fQ = 2.0f * fB * fE + fC * (fF - fD);
            afDiag[1] = fD + fC * fQ;
            afDiag[2] = fF - fC * fQ;
            afSubDiag[0] = fLength;
            afSubDiag[1] = fE - fB * fQ;
            m[0][0] = 1.0;
            m[0][1] = 0.0;
            m[0][2] = 0.0;
            m[1][0] = 0.0;
            m[1][1] = fB;
            m[1][2] = fC;
            m[2][0] = 0.0;
            m[2][1] = fC;
            m[2][2] = -fB;
        }
        else
        {
            afDiag[1] = fD;
            afDiag[2] = fF;
            afSubDiag[0] = fB;
            afSubDiag[1] = fE;
            m[0][0] = 1.0;
            m[0][1] = 0.0;
            m[0][2] = 0.0;
            m[1][0] = 0.0;
            m[1][1] = 1.0;
            m[1][2] = 0.0;
            m[2][0] = 0.0;
            m[2][1] = 0.0;
            m[2][2] = 1.0;
        }
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::QLAlgorithm(SReal afDiag[3], SReal afSubDiag[3])
    {
        for (int i0 = 0; i0 < 3; i0++)
        {
            const unsigned int iMaxIter = 32;
            unsigned int iIter;
            for (iIter = 0; iIter < iMaxIter; iIter++)
            {
                int i1;
                for (i1 = i0; i1 <= 1; i1++)
                {
                    SReal fSum = SMath::Abs(afDiag[i1]) +
                                 SMath::Abs(afDiag[i1 + 1]);
                    if (SMath::Abs(afSubDiag[i1]) + fSum == fSum)
                        break;
                }
                if (i1 == i0)
                    break;

                SReal fTmp0 = (afDiag[i0 + 1] - afDiag[i0]) / (2.0f * afSubDiag[i0]);
                SReal fTmp1 = SMath::Sqrt(fTmp0 * fTmp0 + 1.0f);
                if (fTmp0 < 0.0)
                    fTmp0 = afDiag[i1] - afDiag[i0] + afSubDiag[i0] / (fTmp0 - fTmp1);
                else
                    fTmp0 = afDiag[i1] - afDiag[i0] + afSubDiag[i0] / (fTmp0 + fTmp1);
                SReal fSin = 1.0;
                SReal fCos = 1.0;
                SReal fTmp2 = 0.0;
                for (int i2 = i1 - 1; i2 >= i0; i2--)
                {
                    SReal fTmp3 = fSin * afSubDiag[i2];
                    SReal fTmp4 = fCos * afSubDiag[i2];
                    if (SMath::Abs(fTmp3) >= SMath::Abs(fTmp0))
                    {
                        fCos = fTmp0 / fTmp3;
                        fTmp1 = SMath::Sqrt(fCos * fCos + 1.0f);
                        afSubDiag[i2 + 1] = fTmp3 * fTmp1;
                        fSin = 1.0f / fTmp1;
                        fCos *= fSin;
                    }
                    else
                    {
                        fSin = fTmp3 / fTmp0;
                        fTmp1 = SMath::Sqrt(fSin * fSin + 1.0f);
                        afSubDiag[i2 + 1] = fTmp0 * fTmp1;
                        fCos = 1.0f / fTmp1;
                        fSin *= fCos;
                    }
                    fTmp0 = afDiag[i2 + 1] - fTmp2;
                    fTmp1 = (afDiag[i2] - fTmp0) * fSin + 2.0f * fTmp4 * fCos;
                    fTmp2 = fSin * fTmp1;
                    afDiag[i2 + 1] = fTmp0 + fTmp2;
                    fTmp0 = fCos * fTmp1 - fTmp4;

                    for (int iRow = 0; iRow < 3; iRow++)
                    {
                        fTmp3 = m[iRow][i2 + 1];
                        m[iRow][i2 + 1] = fSin * m[iRow][i2] +
                                          fCos * fTmp3;
                        m[iRow][i2] = fCos * m[iRow][i2] -
                                      fSin * fTmp3;
                    }
                }
                afDiag[i0] -= fTmp2;
                afSubDiag[i0] = fTmp0;
                afSubDiag[i1] = 0.0;
            }

            if (iIter == iMaxIter)
            {
                return false;
            }
        }

        return true;
    }
    //-----------------------------------------------------------------------
    void SMatrix3::EigenSolveSymmetric(SReal afEigenvalue[3],
                                       SVector3 akEigenSVector[3]) const
    {
        SMatrix3 kSMatrix = *this;
        SReal afSubDiag[3];
        kSMatrix.Tridiagonal(afEigenvalue, afSubDiag);
        kSMatrix.QLAlgorithm(afEigenvalue, afSubDiag);

        for (size_t i = 0; i < 3; i++)
        {
            akEigenSVector[i][0] = kSMatrix[0][i];
            akEigenSVector[i][1] = kSMatrix[1][i];
            akEigenSVector[i][2] = kSMatrix[2][i];
        }

        SVector3 kCross = akEigenSVector[1].crossProduct(akEigenSVector[2]);
        SReal fDet = akEigenSVector[0].dotProduct(kCross);
        if (fDet < 0.0)
        {
            akEigenSVector[2][0] = -akEigenSVector[2][0];
            akEigenSVector[2][1] = -akEigenSVector[2][1];
            akEigenSVector[2][2] = -akEigenSVector[2][2];
        }
    }
    //-----------------------------------------------------------------------
    void SMatrix3::TensorProduct(const SVector3 &rkU, const SVector3 &rkV,
                                 SMatrix3 &rkProduct)
    {
        for (size_t iRow = 0; iRow < 3; iRow++)
        {
            for (size_t iCol = 0; iCol < 3; iCol++)
                rkProduct[iRow][iCol] = rkU[iRow] * rkV[iCol];
        }
    }
    //-----------------------------------------------------------------------
    bool SMatrix3::hasScale() const
    {
        SReal t = m[0][0] * m[0][0] + m[1][0] * m[1][0] + m[2][0] * m[2][0];
        if (!SMath::RealEqual(t, 1.0, (SReal)1e-04))
            return true;
        t = m[0][1] * m[0][1] + m[1][1] * m[1][1] + m[2][1] * m[2][1];
        if (!SMath::RealEqual(t, 1.0, (SReal)1e-04))
            return true;
        t = m[0][2] * m[0][2] + m[1][2] * m[1][2] + m[2][2] * m[2][2];
        if (!SMath::RealEqual(t, 1.0, (SReal)1e-04))
            return true;

        return false;
    }
}
