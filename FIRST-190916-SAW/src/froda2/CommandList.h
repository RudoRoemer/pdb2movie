#ifndef COMMANDLIST_H_
#define COMMANDLIST_H_

#include "Command.h"
#include <list>

class CommandList : public Command
{
public:
  CommandList() {}
  virtual ~CommandList() {}
  void addCommand( Command* command ) {
    commands.push_back( command );
  }
  void operator()() {
    for ( std::list<Command*>::iterator it = commands.begin();
          it != commands.end();
          it++ ) (**it)();
  }
private:
  std::list<Command*> commands;
};

#endif /*COMMANDLIST_H_*/
