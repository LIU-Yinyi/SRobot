#ifndef _STYPEPREDEFINE_H_
#define _STYPEPREDEFINE_H_

#include "../SBase.h"

namespace SRobot
{
	typedef bool SBool;					//布尔型
	typedef double SReal;				//浮点数
	typedef unsigned char SUchar; 		//字符

	typedef unsigned int SUint;			//无符号整形
	typedef unsigned char SUint8;		//无符号8位整形
	typedef unsigned short Suint16;		//无符号16位整形
	typedef unsigned int SUint32;		//无符号32位整形
	typedef unsigned long long SUint64;	//无符号64位整形

	typedef int SInt;					//整形
	typedef signed char SInt8;			//8位整形
	typedef short SInt16;				//16位整形
	typedef int SInt32;					//32位整形
	typedef long long SInt64;			//64位整形

	typedef std::string SString;						//字符串
	typedef std::ostringstream SStringWriteStream;
	typedef std::istringstream SStringReadStream;

	typedef unsigned int SSize;	//大小

	class SID;

	class SStream;						//数据流

	class SRadian;						//弧度
	class SDegree;						//角度
	class SAngle;						//角
	class SQuaternion;					//四元数
	class SMatrix3;						//3*3矩阵
	class SMatrix4;						//4*4矩阵
	class SVector2;						//二维向量
	class SVector3;						//三维向量
	class SVector4;						//四维向量
	class SPoint4;						//四维向量
	class SLine2;						//平面线段
	class SLine3;						//空间线段
	class SLine4;						//四维空间线段
	class SRectangle2;					//平面矩形
	class STime;						//时间

	class SMath;						//数学库

	class SInputInterface;				//数据输入接口定义
	class SOutputInterface;				//数据输出接口定义

	//class SMessage;						//消息

	class SAerialRobot;
	//class SArena;
	class SGroundRobot;
	//class SScene;
	typedef SVector2 SPoint2;			//平面点
	typedef SVector3 SPoint3;			//空间点
	typedef SVector3 SPosition;			//位置
	typedef SVector3 SDisplacement;		//位移
	typedef SVector3 SVelocity;			//速度
	typedef SVector3 SAcceleration;		//加速度

}

#endif
