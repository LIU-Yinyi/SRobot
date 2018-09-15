#include "SStream.h"

namespace SRobot
{
	template <typename T>
	SStream &SStream::operator&=(const T &Object)
	{
		(*this) += Serialize(Object);
		return (*this);
	}
}
