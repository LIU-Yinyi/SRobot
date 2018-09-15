#ifndef _SAERIALROBOT_H_
#define _SAERIALROBOT_H_

#include "SBase.h"
#include "STypePreDefine.h"
#include "SMath.h"


namespace SRobot
{
	class SScene;

	class SAerialRobot
	{
	public:
		SPosition		Position;	
		SDisplacement	DisplacementFromLastFrame;
		SVelocity		Velocity;
		SAcceleration	Acceleration;

		SPosition		LastFramePosition;	
		SVelocity		LastFrameVelocity;
		SAcceleration	LastFrameAcceleration;

		SReal			MaxHorizontalVelocity;
		SReal			MaxVerticalVelocity;

		SReal			MaxHorizontalAcceleration;
		SReal			MaxVerticalAcceleration;
	public:
		SAerialRobot();
		~SAerialRobot();

		void			Run(STime Interval,SScene& Scene);

		SPosition		GetPosition() const;
		SDisplacement	GetDisplacementFromLastFrame() const;
		SVelocity		GetVelocity() const;
		SAcceleration	GetAcceleration() const;

		SPosition		GetLastFramePosition() const;	
		SVelocity		GetLastFrameVelocity() const;
		SAcceleration	GetLastFrameAcceleration() const;

		SReal			GetMaxHorizontalVelocity() const; 
		SReal			GetMaxVerticalVelocity() const;

		SReal			GetMaxHorizontalAcceleration() const; 
		SReal			GetMaxVerticalAcceleration() const;

		SPosition		SetPosition(SPosition NewPosition);
		SDisplacement	SetDisplacementFromLastFrame(SDisplacement NewDisplacementFromLastFrame);
		SVelocity		SetVelocity(SVelocity NewVelocity);
		SAcceleration	SetAcceleration(SAcceleration NewAcceleration);

		SPosition		SetLastFramePosition(SPosition NewLastFramePosition);	
		SVelocity		SetLastFrameVelocity(SVelocity NewLastFrameCelocity);
		SAcceleration	SetLastFrameAcceleration(SAcceleration NewLastFrameAcceleration);

		SReal			SetMaxHorizontalVelocity(SReal NewMaxHorizontalVelocity); 
		SReal			SetMaxVerticalVelocity(SReal NewMaxVerticalVelocity);

		SReal			SetMaxHorizontalAcceleration(SReal NewMaxHorizontalAcceleration); 
		SReal			SetMaxVerticalAcceleration(SReal NewMaxVerticalAcceleration);
	};
}

#endif