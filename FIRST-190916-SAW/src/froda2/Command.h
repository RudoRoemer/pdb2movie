#ifndef COMMAND_H_
#define COMMAND_H_

class Command
{
public:
  virtual ~Command() {}
  virtual void operator()() = 0;
};

#endif /*COMMAND_H_*/
