#ifndef OBSERVABLE_H_
#define OBSERVABLE_H_

#include <list>

class Observable;
class Observer {
public:
  Observer() {}
  virtual ~Observer() {}
  virtual void receiveNotification( Observable *observable )=0;
};

class Observable {
public:
  Observable() {}
  virtual ~Observable() {}
  void registerObserver( Observer *obs ) {
    observers.push_back( obs );
  }
  void notifyObservers() {
    for ( std::list<Observer*>::iterator it = observers.begin();
          it != observers.end();
          it++ )
    {
      (*it)->receiveNotification( this );
    }
  }
private:
  std::list<Observer*> observers;
};

#endif /*OBSERVABLE_H_*/
