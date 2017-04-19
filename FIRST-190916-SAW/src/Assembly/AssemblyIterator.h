#ifndef _ITERATOR_
#define _ITERATOR_

using namespace std;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
class AssemblyIterator {

public:
  AssemblyIterator()
  {};
  ~AssemblyIterator()
  {};

public:
  bool hasNext();
  R next();
};

#endif
