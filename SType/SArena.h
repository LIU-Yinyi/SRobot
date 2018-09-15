#ifndef _SARENA_H_
#define _SARENA_H_

#include "SBase.h"
#include "STypePreDefine.h"
#include "SMath.h"

namespace SRobot
{
	class SArena
	{
	public:
		SReal XUpperBound;
		SReal XLowerBound;
		SReal YUpperBound;
		SReal YLowerBound;

		SLine2 GreenLine;
		SLine2 RedLine;
		SLine2 WhiteLineLeft;
		SLine2 WhiteLineRight;

		SReal MaxXAllowed;
		SReal MinXAllowed;
		SReal MaxYAllowed;
		SReal MinYAllowed;
		SReal MaxZAllowed;
		SReal MinZAllowed;
 	};
}
#endif
