#include "SID.h"

namespace SRobot
{
	//Êý¾ÝÁ÷

	SID::SID()
	{
		ID = 0;
	}
	SID::SID(SInt UID)
	{
		ID = (SUint64)UID;
	}

	SID::SID(SUint64 UID)
	{
		ID = UID;
	}

	SID &SID::operator=(const SInt &UID)
	{
		ID = (SUint64)UID;
		return *this;
	}

	SID &SID::operator=(const SUint64 &UID)
	{
		ID = UID;
		return *this;
	}

	SID &SID::operator=(const SUint &UID)
	{
		ID = (SUint64)UID;
		return *this;
	}

	bool SID::operator==(const SID &UID) const
	{
		return this->ID == UID.ID;
	}
	bool SID::operator!=(const SID &UID) const
	{
		return this->ID != UID.ID;
	}
	bool SID::operator>(const SID &UID) const
	{
		return this->ID > UID.ID;
	}
	bool SID::operator>=(const SID &UID) const
	{
		return this->ID >= UID.ID;
	}
	bool SID::operator<(const SID &UID) const
	{
		return this->ID < UID.ID;
	}
	bool SID::operator<=(const SID &UID) const
	{
		return this->ID <= UID.ID;
	}

	std::istream &operator>>(std::istream &is, SID &ID)
	{
		is >> ID.ID;
		return is;
	}

	std::ostream &operator<<(std::ostream &os, SID &ID)
	{
		os << ID.ID;
		return os;
	}

	/*friend SUint64 operator=(const SID &UID)
	{
	return UID.ID;
	}*/

	const SID SID::Anyone = 0;
	const SID SID::Anonymous = 1;
}