#ifndef COMMANDMEMBERFUNCTION_
#define COMMANDMEMBERFUNCTION_

#include "Command.h"

template <class T> 
class CommandMemberFunction : public Command
{
public:
  CommandMemberFunction( T *t_, void (T::*f_)() ) :
    t(t_),
    f(f_) {}
  ~CommandMemberFunction() {}
  void operator()() { (t->*f)(); }
private:
  T *t;
  void (T::*f)();
};

#endif /*COMMANDMEMBERFUNCTION_*/
