#ifndef _SASMMATH_H_
#define _SASMMATH_H_

#include "SBase.h"

#define ASM 0
#define STD 1

namespace SRobot
{

	float asm_arccos(float r);
	//-----------------------------------------------------------------------

	float asm_arcsin(float r);
	//-----------------------------------------------------------------------

	float asm_arctan(float r);
	//-----------------------------------------------------------------------

	float asm_sin(float r);
	//-----------------------------------------------------------------------

	float asm_cos(float r);
	//-----------------------------------------------------------------------

	float asm_tan(float r);
	//-----------------------------------------------------------------------

	float asm_sqrt(float r);
	//-----------------------------------------------------------------------

	   // returns 1 / a for a * a = r
	   // -- Use this for Vector normalisation!!!
	float asm_rsq(float r);
	//-----------------------------------------------------------------------

	   // returns 1 / a for a * a = r
	   // Another version
	float apx_rsq(float r);
	//-----------------------------------------------------------------------
	inline unsigned long long GetCycleCount();
	//-----------------------------------------------------------------------

	   // returns a random number
	/*FORCEINLINE*/ float asm_rand();
	//-----------------------------------------------------------------------

	   // returns the maximum random number
	/*FORCEINLINE*/ float asm_rand_max();
	//-----------------------------------------------------------------------

	   // returns log2( r ) / log2( e )
	float asm_ln(float r);

}
#endif
