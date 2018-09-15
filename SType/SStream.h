#ifndef _SSTREAM_H
#define _SSTREAM_H

#include "STypePreDefine.h"

namespace SRobot
{
	//数据流
	class SStream : public SString
	{
	  public:
		template <typename T>
		SStream &operator&=(const T &Object);
	};
}

#endif