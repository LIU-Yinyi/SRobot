#include "SType.h"
#include "SAerialRobot.h"

namespace SRobot
{
	SAerialRobot::SAerialRobot()
	{
		Position = SPosition(0.0f,0.0f);	
		DisplacementFromLastFrame = SDisplacement(0.0f,0.0f);
		Velocity = SVelocity(0.0f,0.0f);
		Acceleration = SAcceleration(0.0f,0.0f);

		LastFramePosition = SPosition(0.0f,0.0f);
		LastFrameVelocity = SVelocity(0.0f,0.0f);
		LastFrameAcceleration = SAcceleration(0.0f,0.0f);

		MaxHorizontalVelocity = 0.0f;
		MaxVerticalVelocity = 0.0f;

		MaxHorizontalAcceleration = 0.0f;
		MaxVerticalAcceleration = 0.0f;
	}

	SAerialRobot::~SAerialRobot()
	{
	}

	SPosition		SAerialRobot::GetPosition() const { return Position; }
	SDisplacement	SAerialRobot::GetDisplacementFromLastFrame() const {return DisplacementFromLastFrame;}
	SVelocity		SAerialRobot::GetVelocity() const { return Velocity; }
	SAcceleration	SAerialRobot::GetAcceleration() const { return Acceleration; }

	SPosition		SAerialRobot::GetLastFramePosition() const { return LastFramePosition; }
	SVelocity		SAerialRobot::GetLastFrameVelocity() const { return LastFrameVelocity; }
	SAcceleration	SAerialRobot::GetLastFrameAcceleration() const { return LastFrameAcceleration; }

	SReal			SAerialRobot::GetMaxHorizontalVelocity() const { return MaxHorizontalVelocity; }
	SReal			SAerialRobot::GetMaxVerticalVelocity() const { return MaxVerticalVelocity; }

	SReal			SAerialRobot::GetMaxHorizontalAcceleration() const { return MaxHorizontalAcceleration; }
	SReal			SAerialRobot::GetMaxVerticalAcceleration() const { return MaxVerticalAcceleration; }

	SPosition		SAerialRobot::SetPosition(SPosition NewPosition) { return Position = NewPosition; }
	SDisplacement	SAerialRobot::SetDisplacementFromLastFrame(SDisplacement NewDisplacementFromLastFrame) { return DisplacementFromLastFrame = NewDisplacementFromLastFrame; }
	SVelocity		SAerialRobot::SetVelocity(SVelocity NewVelocity) { return Velocity = NewVelocity; }
	SAcceleration	SAerialRobot::SetAcceleration(SAcceleration NewAcceleration) { return Acceleration = NewAcceleration; }

	SPosition		SAerialRobot::SetLastFramePosition(SPosition NewLastFramePosition)	{ return LastFramePosition = NewLastFramePosition; }
	SVelocity		SAerialRobot::SetLastFrameVelocity(SVelocity NewLastFrameVelocity) { return LastFrameVelocity = NewLastFrameVelocity; }
	SAcceleration	SAerialRobot::SetLastFrameAcceleration(SAcceleration NewLastFrameAcceleration) { return LastFrameAcceleration = NewLastFrameAcceleration; }

	SReal			SAerialRobot::SetMaxHorizontalVelocity(SReal NewMaxHorizontalVelocity) { return MaxHorizontalVelocity = NewMaxHorizontalVelocity; }
	SReal			SAerialRobot::SetMaxVerticalVelocity(SReal NewMaxVerticalVelocity) { return MaxVerticalVelocity = NewMaxVerticalVelocity; }

	SReal			SAerialRobot::SetMaxHorizontalAcceleration(SReal NewMaxHorizontalAcceleration) { return MaxHorizontalAcceleration = NewMaxHorizontalAcceleration; }
	SReal			SAerialRobot::SetMaxVerticalAcceleration(SReal NewMaxVerticalAcceleration) { return MaxVerticalAcceleration = NewMaxVerticalAcceleration; }
}
