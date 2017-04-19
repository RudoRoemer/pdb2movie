#ifndef RANDOMROTORPERTURBER_H_
#define RANDOMROTORPERTURBER_H_

class RigidUnitSystem;

class RandomRotorPerturber
{
public:
  RandomRotorPerturber( RigidUnitSystem *rigidUnitSystem_,
                        double size_ );   
  virtual ~RandomRotorPerturber() {}

  void perturb();
private:
  RigidUnitSystem *rigidUnitSystem;
  double size;  
};

#endif /*RANDOMROTORPERTURBER_H_*/
