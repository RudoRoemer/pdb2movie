#ifndef RANDOMCENTERPERTURBER_H_
#define RANDOMCENTERPERTURBER_H_

class RigidUnitSystem;

class RandomCenterPerturber
{
public:
  RandomCenterPerturber( RigidUnitSystem *rigidUnitSystem_,
                         double size_ );   
	virtual ~RandomCenterPerturber() {}
  
  void perturb();
private:
  RigidUnitSystem *rigidUnitSystem;
  double size;  
};

#endif /*RANDOMCENTERPERTURBER_H_*/
