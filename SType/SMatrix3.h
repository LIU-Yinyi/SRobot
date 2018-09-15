#ifndef _SMATRIX3_H_
#define _SMATRIX3_H_

#include "STypePreDefine.h"
#include "SAngle.h"
#include "SVector3.h"

namespace SRobot
{
    class SMatrix3
    {
      public:
        SMatrix3() {}
        explicit SMatrix3(const SReal arr[3][3])
        {
            memcpy(m, arr, 9 * sizeof(SReal));
        }
        SMatrix3(const SMatrix3 &rkSMatrix)
        {
            memcpy(m, rkSMatrix.m, 9 * sizeof(SReal));
        }
        SMatrix3(SReal fEntry00, SReal fEntry01, SReal fEntry02,
                 SReal fEntry10, SReal fEntry11, SReal fEntry12,
                 SReal fEntry20, SReal fEntry21, SReal fEntry22)
        {
            m[0][0] = fEntry00;
            m[0][1] = fEntry01;
            m[0][2] = fEntry02;
            m[1][0] = fEntry10;
            m[1][1] = fEntry11;
            m[1][2] = fEntry12;
            m[2][0] = fEntry20;
            m[2][1] = fEntry21;
            m[2][2] = fEntry22;
        }

        void swap(SMatrix3 &other)
        {
            std::swap(m[0][0], other.m[0][0]);
            std::swap(m[0][1], other.m[0][1]);
            std::swap(m[0][2], other.m[0][2]);
            std::swap(m[1][0], other.m[1][0]);
            std::swap(m[1][1], other.m[1][1]);
            std::swap(m[1][2], other.m[1][2]);
            std::swap(m[2][0], other.m[2][0]);
            std::swap(m[2][1], other.m[2][1]);
            std::swap(m[2][2], other.m[2][2]);
        }

        const SReal *operator[](size_t iRow) const { return m[iRow]; }
        SReal *operator[](size_t iRow) { return m[iRow]; }

        SVector3 GetColumn(size_t iCol) const;
        void SetColumn(size_t iCol, const SVector3 &vec);

        void FromAxes(const SVector3 &xAxis, const SVector3 &yAxis, const SVector3 &zAxis);

        SMatrix3 &operator=(const SMatrix3 &rkSMatrix)
        {
            memcpy(m, rkSMatrix.m, 9 * sizeof(SReal));
            return *this;
        }

        bool operator==(const SMatrix3 &rkSMatrix) const;
        bool operator!=(const SMatrix3 &rkSMatrix) const { return !operator==(rkSMatrix); }

        SMatrix3 operator+(const SMatrix3 &rkSMatrix) const;

        SMatrix3 operator-() const;
        SMatrix3 operator-(const SMatrix3 &rkSMatrix) const;

        SMatrix3 operator*(SReal fScalar) const;
        SMatrix3 operator*(const SMatrix3 &rkSMatrix) const;
        SVector3 operator*(const SVector3 &rkSVector) const;
        friend SVector3 operator*(const SVector3 &rkSVector, const SMatrix3 &rkSMatrix);
        friend SMatrix3 operator*(SReal fScalar, const SMatrix3 &rkSMatrix);

        SMatrix3 Transpose() const;
        bool Inverse(SMatrix3 &rkInverse, SReal fTolerance = 1e-06) const;
        SMatrix3 Inverse(SReal fTolerance = 1e-06) const;
        SReal Determinant() const;

        void SingularValueDecomposition(SMatrix3 &rkL, SVector3 &rkS, SMatrix3 &rkR) const;
        void SingularValueComposition(const SMatrix3 &rkL, const SVector3 &rkS, const SMatrix3 &rkR);

        void Orthonormalize();

        void QDUDecomposition(SMatrix3 &rkQ, SVector3 &rkD, SVector3 &rkU) const;

        SReal SpectralNorm() const;

        void ToAngleAxis(SVector3 &rkAxis, SRadian &rfAngle) const;
        void ToAngleAxis(SVector3 &rkAxis, SDegree &rfAngle) const
        {
            SRadian r;
            ToAngleAxis(rkAxis, r);
            rfAngle = r;
        }
        void FromAngleAxis(const SVector3 &rkAxis, const SRadian &fSRadians);

        bool ToEulerAnglesXYZ(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const;
        bool ToEulerAnglesXZY(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const;
        bool ToEulerAnglesYXZ(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const;
        bool ToEulerAnglesYZX(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const;
        bool ToEulerAnglesZXY(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const;
        bool ToEulerAnglesZYX(SRadian &rfYAngle, SRadian &rfPAngle, SRadian &rfRAngle) const;

        void FromEulerAnglesXYZ(const SRadian &fYAngle, const SRadian &fPAngle, const SRadian &fRAngle);
        void FromEulerAnglesXZY(const SRadian &fYAngle, const SRadian &fPAngle, const SRadian &fRAngle);
        void FromEulerAnglesYXZ(const SRadian &fYAngle, const SRadian &fPAngle, const SRadian &fRAngle);
        void FromEulerAnglesYZX(const SRadian &fYAngle, const SRadian &fPAngle, const SRadian &fRAngle);
        void FromEulerAnglesZXY(const SRadian &fYAngle, const SRadian &fPAngle, const SRadian &fRAngle);
        void FromEulerAnglesZYX(const SRadian &fYAngle, const SRadian &fPAngle, const SRadian &fRAngle);

        void EigenSolveSymmetric(SReal afEigenvalue[3], SVector3 akEigenSVector[3]) const;

        static void TensorProduct(const SVector3 &rkU, const SVector3 &rkV, SMatrix3 &rkProduct);

        bool hasScale() const;

        friend std::ostream &operator<<(std::ostream &o, const SMatrix3 &mat)
        {
            o << "SMatrix3(" << mat[0][0] << ", " << mat[0][1] << ", " << mat[0][2] << ", "
              << mat[1][0] << ", " << mat[1][1] << ", " << mat[1][2] << ", "
              << mat[2][0] << ", " << mat[2][1] << ", " << mat[2][2] << ")";
            return o;
        }

        friend SStream Serialize(const SMatrix3 &Matrix3);
        //template <>
        //friend SMatrix3 UnSerialize <SMatrix3> (SStream Stream);

        static const SReal EPSILON;
        static const SMatrix3 ZERO;
        static const SMatrix3 IDENTITY;

      public:
        void Tridiagonal(SReal afDiag[3], SReal afSubDiag[3]);
        bool QLAlgorithm(SReal afDiag[3], SReal afSubDiag[3]);

        static const SReal msSvdEpsilon;
        static const unsigned int msSvdMaxIterations;
        static void Bidiagonalize(SMatrix3 &kA, SMatrix3 &kL, SMatrix3 &kR);
        static void GolubKahanStep(SMatrix3 &kA, SMatrix3 &kL, SMatrix3 &kR);

        static SReal MaxCubicRoot(SReal afCoeff[3]);

        SReal m[3][3];

        friend class SMatrix4;
    };
}

#endif