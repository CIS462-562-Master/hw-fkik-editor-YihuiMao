#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	
	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	//vec3 tempPos= m_Guide.getGlobalTranslation();
	AJoint *root = m_pSkeleton->getRootNode();
	vec3 guidePos = m_Guide.getGlobalTranslation();
	mat3 guiderot= m_Guide.getGlobalRotation();
	vec3 rootPos = root->getGlobalTranslation();
    vec3 updated = guidePos + guiderot * rootPos;


	guideTargetPos[1] = 0;
	
	// 2.	Set the y component of the guide position to 0
	updated[1] = 0.f;
	m_Guide.setGlobalTranslation(updated);
	
	// 3.	Set the global rotation of the guide joint towards the guideTarget


	vec3 z = vec3(0, 0, 1);
	vec3 l = (guideTargetPos - guidePos).Normalize();
	vec3 axis = z ^ l;
	double angle = acos(z * l);
	mat3 Rotation;
	Rotation.FromAxisAngle(axis, angle);
	m_Guide.setGlobalRotation(Rotation);




}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	AJoint *root = m_pSkeleton->getRootNode();
	vec3 pos = root->getGlobalTranslation();

	float height;




	pos[1] = pos[1] +  (leftHeight+rightHeight)/2;
	root->setLocalTranslation(pos);

	m_pSkeleton->update();


	
	// 2.	Update the character with Limb-based IK 

	
	
	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal

		ATarget target;
		vec3 lFoot = leftFoot->getGlobalTranslation();
		lFoot[1] = leftHeight;
		target.setGlobalTranslation(lFoot);
		m_IKController->IKSolver_Limb(leftFoot->getID(), target);

		vec3 y_axis = (0.f, 1.f, 0.f);
		vec3 axis = y_axis ^ leftNormal;
		vec3 unit_leftNormal = unit_leftNormal.Normalize();
		double angle = acos(y_axis*unit_leftNormal);

		
		mat3 guideRotation;
		guideRotation.FromAxisAngle(axis, angle);
		leftFoot->setLocalRotation(guideRotation);
		
		
	}
	if (rotateRight)
	{


		ATarget target;	
		vec3 rFoot = rightFoot->getGlobalTranslation();
		rFoot[1] = rightHeight;
		target.setGlobalTranslation(rFoot);
		m_IKController->IKSolver_Limb(rightFoot->getID(), target);



		vec3 y_axis = (0.f, 1.f, 0.f);
		vec3 axis = y_axis ^ rightNormal;
		vec3 unit_rightNormal = unit_rightNormal.Normalize();
		double angle = acos(y_axis * unit_rightNormal);



		mat3 guideRotation;
		guideRotation.FromAxisAngle(axis, angle);
	    rightFoot->setLocalRotation(guideRotation);
		
		
	}
	m_pSkeleton->update();
}
