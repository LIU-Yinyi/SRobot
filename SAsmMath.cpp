#include "SAsmMath.h"

typedef long long __int64;

namespace SRobot
{

	const float pi = 4.0f * atan(1.0f);
	const float half_pi = 0.5f * pi;

	float asm_arccos(float r)
	{
		// return half_pi + arctan( r / -sqr( 1.f - r * r ) );
#if  CALCULATE_MODE == ASM
		float asm_one = 1.f;
		float asm_half_pi = half_pi;
		__asm
		{
			fld r // r0 = r
			fld r // r1 = r0, r0 = r
				fmul r // r0 = r0 * r
				fsubr asm_one // r0 = r0 - 1.f
				fsqrt // r0 = sqrtf( r0 )
				fchs // r0 = - r0
				fdiv // r0 = r1 / r0
				fld1 // {{ r0 = atan( r0 )
				fpatan // }}
				fadd asm_half_pi // r0 = r0 + pi / 2
		} // returns r0
#else
		return float(acos(r));
#endif
	}

	//-----------------------------------------------------------------------

	float asm_arcsin(float r)
	{
		// return arctan( r / sqr( 1.f - r * r ) );
#if  CALCULATE_MODE == ASM
		const float asm_one = 1.f;
		__asm
		{
			fld r // r0 = r
			fld r // r1 = r0, r0 = r
				fmul r // r0 = r0 * r
				fsubr asm_one // r0 = r0 - 1.f
				fsqrt // r0 = sqrtf( r0 )
				fdiv // r0 = r1 / r0
				fld1 // {{ r0 = atan( r0 )
				fpatan // }}
		} // returns r0
#else
		return float(asin(r));
#endif
	}

	//-----------------------------------------------------------------------

	float asm_arctan(float r)
	{
#if  CALCULATE_MODE == ASM
		__asm
		{
			fld r // r0 = r
			fld1 // {{ r0 = atan( r0 )
				fpatan // }}
		} // returns r0
#else
		return float(atan(r));
#endif
	}

	//-----------------------------------------------------------------------

	float asm_sin(float r)
	{
#if  CALCULATE_MODE == ASM
		__asm
		{
			fld r // r0 = r
			fsin // r0 = sinf( r0 )
		} // returns r0
#else
		return sin(r);
#endif
	}

	//-----------------------------------------------------------------------

	float asm_cos(float r)
	{
#if  CALCULATE_MODE == ASM
		__asm
		{
			fld r // r0 = r
			fcos // r0 = cosf( r0 )
		} // returns r0
#else
		return cos(r);
#endif
	}

	//-----------------------------------------------------------------------

	float asm_tan(float r)
	{
#if  CALCULATE_MODE == ASM
		// return sin( r ) / cos( r );
		__asm
		{
			fld r // r0 = r
			fsin // r0 = sinf( r0 )
				fld r // r1 = r0, r0 = r
				fcos // r0 = cosf( r0 )
				fdiv // r0 = r1 / r0
		} // returns r0
#else
		return tan(r);
#endif
	}

	//-----------------------------------------------------------------------

	float asm_sqrt(float r)
	{
		// returns a for a * a = r
#if  CALCULATE_MODE == ASM
		__asm
		{
			fld r // r0 = r
			fsqrt // r0 = sqrtf( r0 )
		} // returns r0
#else
		return sqrt(r);
#endif
	}

	//-----------------------------------------------------------------------

	// returns 1 / a for a * a = r
	// -- Use this for Vector normalisation!!!
	float asm_rsq(float r)
	{
#if  CALCULATE_MODE == ASM
		__asm
		{
			fld1 // r0 = 1.f
			fld r // r1 = r0, r0 = r
				fsqrt // r0 = sqrtf( r0 )
				fdiv // r0 = r1 / r0
		} // returns r0
#else
		return 1. / sqrt(r);
#endif
	}

	//-----------------------------------------------------------------------

	// returns 1 / a for a * a = r
	// Another version
	float apx_rsq(float r)
	{
#if  CALCULATE_MODE == ASM
		const float asm_dot5 = 0.5f;
		const float asm_1dot5 = 1.5f;

		__asm
		{
			fld r // r0 = r
			fmul asm_dot5 // r0 = r0 * .5f
				mov eax, r // eax = r
				shr eax, 0x1 // eax = eax >> 1
				neg eax // eax = -eax
				add eax, 0x5F400000 // eax = eax & MAGICAL NUMBER
				mov r, eax // r = eax
				fmul r // r0 = r0 * r
				fmul r // r0 = r0 * r
				fsubr asm_1dot5 // r0 = 1.5f - r0
				fmul r // r0 = r0 * r
		} // returns r0
#else
		return 1. / sqrt(r);
#endif
	}

	inline unsigned long long GetCycleCount()
	{
		//注: 以下语句由于g++无法编译被注释掉
		/*
		__asm _emit 0x0F
		__asm _emit 0x31
		*/
		return 0;
	}

	//-----------------------------------------------------------------------

	// returns a random number
	/*FORCEINLINE*/ float asm_rand()
	{
#if  CALCULATE_MODE == ASM
		static unsigned __int64 q = GetCycleCount();
		_asm
		{
			movq mm0, q

			// do the magic MMX thing
				pshufw mm1, mm0, 0x1E
				paddd mm0, mm1

				// move to integer memory location and free MMX
				movq q, mm0
				emms
		}
		return float(q);
#else
		return float(rand());
#endif

	}

	//-----------------------------------------------------------------------

	// returns the maximum random number
	/*FORCEINLINE*/ float asm_rand_max()
	{
#if  CALCULATE_MODE == ASM
		return (float)(std::numeric_limits<unsigned __int64>::max)();
		return 9223372036854775807.0f;
#else
		return float(RAND_MAX);
#endif
	}

	//-----------------------------------------------------------------------

	// returns log2( r ) / log2( e )
	float asm_ln(float r)
	{
#if  CALCULATE_MODE == ASM

		const float asm_1_div_log2_e = .693147180559f;
		const float asm_neg1_div_3 = -.33333333333333333333333333333f;
		const float asm_neg2_div_3 = -.66666666666666666666666666667f;
		const float asm_2 = 2.f;

		int log_2 = 0;

		__asm
		{
			// log_2 = ( ( r >> 0x17 ) & 0xFF ) - 0x80;
			mov eax, r
			sar eax, 0x17
				and eax, 0xFF
				sub eax, 0x80
				mov log_2, eax

				// r = ( r & 0x807fffff ) + 0x3f800000;
				mov ebx, r
				and ebx, 0x807FFFFF
				add ebx, 0x3F800000
				mov r, ebx

				// r = ( asm_neg1_div_3 * r + asm_2 ) * r + asm_neg2_div_3;   // (1)
				fld r
				fmul asm_neg1_div_3
				fadd asm_2
				fmul r
				fadd asm_neg2_div_3
				fild log_2
				fadd
				fmul asm_1_div_log2_e
		}
#else
		return log(r);
#endif
	}
}
