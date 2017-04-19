#ifndef FRODA2BUILDER_FROMMOLFRAMEWORK_H_
#define FRODA2BUILDER_FROMMOLFRAMEWORK_H_

class MolFramework;
class RigidUnitSystem;
class ConstraintEnforcingPotential;
class PDB;
class Parameters;
class SymmetryMatrices;
class AmberPrmtop;

class Froda2Builder_FromMolFramework
{
public:
	Froda2Builder_FromMolFramework(
				       MolFramework &structure, Parameters &parameters, SymmetryMatrices *symmetryMatrices );
	virtual ~Froda2Builder_FromMolFramework();

  RigidUnitSystem *getRigidUnitSystem();
  ConstraintEnforcingPotential *getConstraintEnforcingPotential();
  PDB *getPDB();
  AmberPrmtop *getAmberPrmtop();

private:
  RigidUnitSystem *rigidUnitSystem;
  ConstraintEnforcingPotential *cep;
  PDB *pdb;
  AmberPrmtop *prmtop;
};

#endif /*FRODA2BUILDER_FROMMOLFRAMEWORK_H_*/
