#include "SType.h"
#include "SGroundRobot.h"

namespace SRobot
{
	SGroundRobot::SGroundRobot()
	{
		Position = SPosition(0.0f, 0.0f);
		DisplacementFromLastFrame = SDisplacement(0.0f, 0.0f);
		Velocity = SVelocity(0.0f, 0.0f);

		LastFramePosition = SPosition(0.0f, 0.0f);
		LastFrameVelocity = SVelocity(0.0f, 0.0f);

		MaxVelocity = 0.0f;
	}

	SGroundRobot::~SGroundRobot()
	{
	}

	SPosition SGroundRobot::GetPosition() const
	{
		return Position;
	}
	SDisplacement SGroundRobot::GetDisplacementFromLastFrame() const
	{
		return DisplacementFromLastFrame;
	}
	SVelocity SGroundRobot::GetVelocity() const
	{
		return Velocity;
	}

	SPosition SGroundRobot::GetLastFramePosition() const
	{
		return LastFramePosition;
	}
	SVelocity SGroundRobot::GetLastFrameVelocity() const
	{
		return LastFrameVelocity;
	}

	SReal SGroundRobot::GetMaxVelocity() const
	{
		return MaxVelocity;
	}

	STime SGroundRobot::GetTravelTime() const
	{
		return TravelTime;
	}
	STime SGroundRobot::GetMaxTravelTime() const
	{
		return MaxTravelTime;
	}

	SPosition SGroundRobot::SetPosition(SPosition NewPosition)
	{
		return Position = NewPosition;
	}
	SDisplacement SGroundRobot::SetDisplacementFromLastFrame(SDisplacement NewDisplacementFromLastFrame)
	{
		return DisplacementFromLastFrame = NewDisplacementFromLastFrame;
	}
	SVelocity SGroundRobot::SetVelocity(SVelocity NewVelocity)
	{
		return Velocity = NewVelocity;
	}

	SPosition SGroundRobot::SetLastFramePosition(SPosition NewLastFramePosition)
	{
		return LastFramePosition = NewLastFramePosition;
	}
	SVelocity SGroundRobot::SetLastFrameVelocity(SVelocity NewLastFrameVelocity)
	{
		return LastFrameVelocity = NewLastFrameVelocity;
	}

	SReal SGroundRobot::SetMaxVelocity(SReal NewMaxVelocity)
	{
		return MaxVelocity = NewMaxVelocity;
	}

	STime SGroundRobot::SetTravelTime(STime NewTravelTime)
	{
		return TravelTime = NewTravelTime;
	}
	STime SGroundRobot::SetMaxTravelTime(STime NewMaxTravelTime)
	{
		return MaxTravelTime = NewMaxTravelTime;
	}
}
