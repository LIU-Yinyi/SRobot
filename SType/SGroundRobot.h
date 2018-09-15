#ifndef _SGROUNDROBOT_H_
#define _SGROUNDROBOT_H_

#include "SBase.h"
#include "STypePreDefine.h"
#include "SMath.h"

namespace SRobot
{
	class SScene;

	class SGroundRobot
	{
	  public:
		SPosition Position;
		SDisplacement DisplacementFromLastFrame;
		SVelocity Velocity;

		SPosition LastFramePosition;
		SVelocity LastFrameVelocity;

		SReal MaxVelocity;

		STime TravelTime;
		STime MaxTravelTime;

	  public:
		SGroundRobot();
		~SGroundRobot();

		void Run(STime Interval, SScene &Scene);

		SPosition GetPosition() const;
		SDisplacement GetDisplacementFromLastFrame() const;
		SVelocity GetVelocity() const;

		SPosition GetLastFramePosition() const;
		SVelocity GetLastFrameVelocity() const;

		SReal GetMaxVelocity() const;

		STime GetTravelTime() const;
		STime GetMaxTravelTime() const;

		SPosition SetPosition(SPosition NewPosition);
		SDisplacement SetDisplacementFromLastFrame(SDisplacement NewDisplacementFromLastFrame);
		SVelocity SetVelocity(SVelocity NewVelocity);

		SPosition SetLastFramePosition(SPosition NewLastFramePosition);
		SVelocity SetLastFrameVelocity(SVelocity NewLastFrameCelocity);

		SReal SetMaxVelocity(SReal NewMaxVelocity);

		STime SetTravelTime(STime NewTravelTime);
		STime SetMaxTravelTime(STime NewMaxTravelTime);
	};
}

#endif
