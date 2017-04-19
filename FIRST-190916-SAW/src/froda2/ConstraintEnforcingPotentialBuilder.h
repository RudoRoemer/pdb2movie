#ifndef CONSTRAINTENFORCINGPOTENTIALBUILDER_H_
#define CONSTRAINTENFORCINGPOTENTIALBUILDER_H_

class RigidUnitSystem;
class Energy;
class Gradient;
class Mismatch;
class Repulsion;
class VL;
class RopesAssignment;
class ConstraintEnforcingPotential;
class InequalityAngleConstraints;

class ConstraintEnforcingPotentialBuilder
{
public:
	ConstraintEnforcingPotentialBuilder( RigidUnitSystem *rigidUnitSystem );
	virtual ~ConstraintEnforcingPotentialBuilder();
  void setSharedPointsEnergy();
  void setRopesEnergy( const RopesAssignment* ropesAssignment );
  void setOverlapEnergy( Repulsion *repulsion );
  void setInequalityAngleConstraints( InequalityAngleConstraints *con );
  ConstraintEnforcingPotential *getConstraintEnforcingPotential();
protected:
  RigidUnitSystem *rigidUnitSystem;
  Energy *energy;
  Gradient *gradient;
  Mismatch *mismatch;
private:
  ConstraintEnforcingPotential *cep;
  VL *verletList;
};

#endif /*CONSTRAINTENFORCINGPOTENTIALBUILDER_H_*/
