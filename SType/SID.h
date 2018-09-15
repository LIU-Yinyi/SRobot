#ifndef _SID_H
#define _SID_H

#include "STypePreDefine.h"

namespace SRobot
{
	//数据流
	class SID
	{
	  public:
		SUint64 ID;

		SID();
		SID(SUint64 UID);
		SID(SInt UID);
		SID &operator=(const SInt &UID);
		SID &operator=(const SUint64 &UID);
		SID &operator=(const SUint &UID);

		bool operator==(const SID &UID) const;
		bool operator!=(const SID &UID) const;
		bool operator>(const SID &UID) const;
		bool operator>=(const SID &UID) const;
		bool operator<(const SID &UID) const;
		bool operator<=(const SID &UID) const;

		friend std::istream &operator>>(std::istream &is, SID &ID);
		friend std::ostream &operator<<(std::ostream &os, SID &ID);

        /*friend SUint64 operator=(const SID &UID)
        {
        	return UID.ID;
        }*/

		static const SID Anyone;
		static const SID Anonymous;
	};
}

#endif