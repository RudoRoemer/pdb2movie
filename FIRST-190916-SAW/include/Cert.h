/* Copyright (c) 2006 Stanford University and An Nguyen.
 *
 * Permission is hereby granted, free of charge, to any person  obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be  included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND  NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN  ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN  CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * Certificate to be used in KDS
 * Each certificate must define 
 *   getCBound(): conservative bound on when this certificate may fail
 *                getCBound must return a value less than current time
 *                   when the certificate has already failed
 *   update: what to do when the certificate fails
 */

#ifndef _CERT_H
#define _CERT_H

#include <iostream>
#include <assert.h>

// #define myassert(e, m, id1, id2) ;
#define myassert(e, m, id1, id2) if (!(e)) \
   fprintf(stderr, "Assert failed, %s:%d, %s - %x %x\n", __FILE__, __LINE__, \
           m, id1, id2)

using namespace std;

class Cert {
public:
  virtual float getCBound(float now) const { 
    assert(false);
    return 0; 
  }
  virtual void update() { assert(false); }
  virtual void print(ostream &os) const { assert(false); }
  friend ostream &operator<<(ostream &os, const Cert& cert) {
    cert.print(os); return os;
  }
  virtual ~Cert() {}
};

template <class N>
class ParentCert : public Cert {
  N *node;
public:
  ParentCert(N *node_) : node(node_) {}
  float getCBound(float now) const;
  void update();
  void print(ostream &os) const;
  bool isOn(const N* node_) const { return node == node_; }
  virtual ~ParentCert() {}
};

template <class N>
class EdgeCert : public Cert {
  N *node1, *node2;
public:
  EdgeCert(N *node1_, N *node2_) : node1(node1_), node2(node2_) {}
  float getCBound(float now) const;
  void update();
  void print(ostream &os) const;
  bool isOn(const N* node_) const { return node1 == node_ || node2 == node_; }
  virtual ~EdgeCert() {}
};

template <class N>
class NEdgeCert : public Cert {
  N *node1, *node2;
public:
  NEdgeCert(N *node1_, N *node2_) : node1(node1_), node2(node2_) {}
  float getCBound(float now) const;
  void update();
  void print(ostream &os) const;
  bool isOn(const N* node_) const { return node1 == node_ || node2==node_; }
  virtual ~NEdgeCert() {}
};

template <class N>
void ParentCert<N>::print(ostream &os) const {
  os << "parentcert " << *node;
}

template <class N>
void EdgeCert<N>::print(ostream &os) const {
  os << "edgecert " << *node1 << " " << *node2;
}

template <class N>
void NEdgeCert<N>::print(ostream &os) const {
  os << "nedgecert " << *node1 << " " << *node2;
}

template <class N>
float ParentCert<N>::getCBound(float now) const {
  return node->firstTimeFurther(
    now, *node->getParent(), node->getParent()->getRadius2());
}

template <class N>
float EdgeCert<N>::getCBound(float now) const {
  return min(
    node1->firstTimeFurther(
      now, *node2, node1->getHierarchy()->getC2() * node1->getRadius2()),
    node1->firstTimeNearer(now, *node2, node1->getRadius2()));
}

template <class N>
float NEdgeCert<N>::getCBound(float now) const {
  return node1->firstTimeNearer(
    now, *node2, node1->getHierarchy()->getC2() * node1->getRadius2());
}

template <class N>
void ParentCert<N>::update() {
  N* oldparent = node->getParent();

  // After the call to node->removeParent(), this certificate is destroyed. 
  N* savednode = node;
  node->removeParent();
  savednode->findNewParent(oldparent);
}

template <class N>
void EdgeCert<N>::update() {
  if (node1->dist2(*node2) < node1->getRadius2()) {
    // edge too short, demote nodes with lower level
    if (node1->getParent()->getPoint() == node1->getPoint()) {
      // cerr << "    demote " << node2->getId() << endl;
      node2->demote(*node1);
    } else {
      // cerr << "    demote " << node1->getId() << endl;
      node1->demote(*node2);
    }
  } else if (node1->dist2(*node2) >
      node1->getHierarchy()->getC2() * node1->getRadius2()) {
    // edge is too long but is still a potential edge
    // cerr << "    removeedge " << endl;

    // After the call to fullRemoveNeighbor(), this certificate is destroyed. 
    N* savednode1 = node1;
    N* savednode2 = node2;
    savednode1->fullRemoveNeighbor(*savednode2);
    savednode1->addPotentialNeighbor(savednode2);
  } else {
    // this should not happen
    myassert(true, "invalid edgecert update", -1, -1);
  }
}

template <class N>
void NEdgeCert<N>::update() {
  myassert((node1->dist2(*node2) <
      node1->getHierarchy()->getC2() * node1->getRadius2()), 
    "invalid nedgecert update", -1, -1);

  // After call to removePotentialNeighbor(), this certificate is destroyed. 
  N* savednode1 = node1;
  N* savednode2 = node2;
  savednode1->removePotentialNeighbor(savednode2);
  savednode1->fullAddNeighbor(*savednode2);
}

#endif
