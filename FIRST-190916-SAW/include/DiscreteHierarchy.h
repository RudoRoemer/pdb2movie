// Created by An Nguyen
//
// Modifications by Dan Farrell
//   June 16, 2006
//
// Here is the paper that describes the Discrete Centers Heirarchy
// and its parameters:
//
// Jie Gao, Leonidas J. Guibas, An Nguyen.  "Deformable Spanners
// and Applications."  Computational Geometry, volume 35, issues 1-2,
// August 2006, pp. 2-19.
//
// Also see www.cs.sunysb.edu/~jgao, the website of the paper's first author
//
// The original dch code was obtained from
//   simtk.org (project name for the code is "dch")
// At that site, the "dch" original code comes with a helpful README.txt file
// and a very basic sample program that uses the dch.
//
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

#include <iterator>
#include <iostream>
#include <vector>
#include <queue>
#include <math.h>
#include <string>
#include <algorithm>
#include <climits>
#include <cstdio>

#include "Cert.h"
#include "PQueue.h"

#define USEMAP

#ifdef USEMAP
#include <map>
#elif defined __GNUC__
#include <ext/hash_map>
using namespace __gnu_cxx;
#elif defined _MSC_VER
#include <hash_map>
using namespace stdext;
#else
#include <hash_map>
#endif

using namespace std;

template <class MP>
class DiscreteHierarchy;

template <class S, class T>
void full_remove(S& s, const T& t) {
  s.erase(remove(s.begin(), s.end(), t), s.end());
}

#ifndef USEMAP
#ifdef __GNUC__
template <class T>
struct PtrHash {
  typedef const T *TPtr;
  hash<int> H;
  size_t operator() (const TPtr &t) const {
    return H((int) t);
  }
};

template <class T>
struct PtrEqual {
  bool operator()(const T* s1, const T* s2) const {
    return s1 == s2;
  }
};
#endif
#endif

template <class MP> class Node;
template <class MP> class DiscreteHierarchy;
template <class MP> ostream& operator<< (ostream& os, const Node<MP>& node);
template <class MP> ostream& operator<< (ostream& os, 
                                         const DiscreteHierarchy<MP>& node);

template <class MP>
class Node {
public:
  typedef Node<MP> NodeMP;
  typedef DiscreteHierarchy<MP> DH;
  typedef ParentCert< NodeMP > ParentCertMP;
  typedef EdgeCert< NodeMP > EdgeCertMP;
  typedef NEdgeCert< NodeMP > NEdgeCertMP;
  typedef PriorityQNode<Cert> *PQLink;
  typedef PriorityQ<Cert> PQueue;

private:
  static int idcount;
  int id;
  DH *dh;			// the DH this node is a part of
  const MP* point;		// pointer to the point asso. with this node 
  float radius2;		// square radius of this node
                                //   all children must be within this radius
  Node* parent;			// parent in the hierarchy, self if root
  vector< Node* > children;	// list of all children nodes
  vector< Node* > neighbors;	// list of all neighbor nodes

  // certificates for parent, children, and neighbors
  PQLink p_cert;		// ParentCerts
  vector<PQLink> e_certs;	// list of EdgeCerts, one for each neighbor
  vector<PQLink> ne_certs;	// list of NECerts, one for each potential nb
    
public:
  Node(const MP* point_, Node* parent_, DH *dh_, float radius2_) : 
      dh(dh_), point(point_), radius2(radius2_), parent(parent_), p_cert(NULL) {
    id = idcount++; 
    // cerr << "  node created " << id << endl;
  }
  ~Node() {
    myassert(id >= 0, "bad call to deconstructor", id, -1);
    // cerr << "  node removed " << id << endl;
    id = -1;
  }
  
  static void resetIDcount() { idcount=0; }
  /*
   * Wrapper functions dist2(), firstTimeFurther(), and firstTimeNearer()
   * Class MP must have these functions defined.
   */
  float dist2(const Node& other) const { return point->dist2(*other.point); }
  float firstTimeFurther(float currtime, const Node& other, float dist2) const {
    return point->firstTimeFurther(currtime, *other.point, dist2); 
  }
  float firstTimeNearer(float currtime, const Node& other, float dist2) const {
    return point->firstTimeNearer(currtime, *other.point, dist2); 
  }

  int getId() const { return id; }
  const MP* getPoint() const { return point; }
  float getRadius2() const { return radius2; }
  void setRadius2(float r2) { radius2 = r2; }
  const vector< Node* >& getChildren() const { return children; }
  const vector< Node* >& getNeighbors() const { return neighbors; }
  const Node* getParent() const { return parent; }
  Node* getParent() { return parent; }
  const DH* getHierarchy() const { return dh; }

  // the constant 1.1 is used to avoid numerical difficulty
  // perhaps the level number should be maintain exclusively instead
  //   of implicitly through the radius value
  bool isBottom() const { return radius2 < dh->getMinRadius2() * 1.1 ; }
  bool isNextToBottom() const {
    return radius2 < dh->getMinRadius2() * dh->getBeta2() * 1.1 && !isBottom(); 
  }

  bool isNeighbor(Node* nb) const {
    return find(neighbors.begin(), neighbors.end(), nb) != neighbors.end();
  }

  // find various certificates
  PQLink findParentCert() const { return p_cert; }
  PQLink findEdgeCert(Node* neighbor) const;
  PQLink findNEdgeCert(Node* neighbor) const;

  void addNeighbor(Node* neighbor);		// simply add a neighbor
  void addPotentialNeighbor(Node* neighbor);	// simply add a potential nb
  void removeNeighbor(Node* neighbor);		// simply remove a neighbor
  void removePotentialNeighbor(Node* neighbor);	// simply remove a potential nb

  // Setting the parent of this node.  If findNeighbors is true, 
  // the neighbors of the this node is automatically computed.
  // During maintenance of the hierarchy, findNeighbors should
  //   be false as neighbors have been computed already
  void setParent(Node* newparent, bool findNeighbors);

  // set the current parent to NULL, remove the potential neighbors
  // as the result of this unlink
  void removeParent();		

  // Compute a new parent for the node and compute all potential neighbors
  // associated with this new parent.  The parent of the
  // node must be NULL before the call.  If newnode is true,
  // the neighbors of the node is also computed.
  void findNewParent(Node* relaxedparent, bool newnode = false);

  // Add/remove a neighbor.  Also update the potential neighbors
  // affected as the result of this add/remove neighbor
  void fullRemoveNeighbor(Node& other);
  void fullAddNeighbor(Node& other);

  // demote the node by one level.
  void demote(Node& oldneighbor);

  friend ostream& operator<< <MP> (ostream& os, const Node<MP>& node);
};

template <class MP>
class DiscreteHierarchy {
public:
  typedef Node<MP> NodeMP;
  typedef ParentCert< NodeMP> ParentCertMP;
  typedef EdgeCert< NodeMP> EdgeCertMP;
  typedef NEdgeCert< NodeMP> NEdgeCertMP;
  typedef PriorityQNode<Cert> *PQLink;
  typedef PriorityQ<Cert> PQueue;

#ifdef USEMAP
  typedef map<const MP*, NodeMP*> NodeHash;
#elif defined __GNUC__
  typedef hash_map<const MP*, NodeMP*, PtrHash<MP>, PtrEqual<MP> > NodeHash;
#elif defined _MSC_VER
  typedef hash_map<const MP*, NodeMP*> NodeHash;
#else
  typedef hash_map<const MP*, NodeMP*> NodeHash;
#endif

private:
  float c2;			// square of c 
  float beta2;			// square of beta
  float now;			// current time
  float minradius2;		// minimum square radius at the bottom level
  float beta; //square root of beta2
				//   radius2 of nodes at bottom <= minradius2
  float bottomLevelMaxDistance;
  int countUpdates;
  NodeMP* root;			// top most node
  PQueue certs;			// queue of ALL certificates
  NodeHash bottom;		// store all bottom level nodes, i.e. all points
public:
  //DWF changed ordering of initializers
  DiscreteHierarchy(float beta, float c, float minradius) :
       c2(c*c), beta2(beta*beta), now(0), minradius2(minradius*minradius),
       beta(sqrt(beta*beta)), countUpdates(0), root(NULL)  {
    bottomLevelMaxDistance = c*minradius;
  }

  ~DiscreteHierarchy() {
    queue< NodeMP* > candidates;
    candidates.push(root);
    while (!candidates.empty()) {
      NodeMP* curr = candidates.front();
      candidates.pop();

      for (unsigned int i = 0; i < curr->getChildren().size(); i++) {
        candidates.push(curr->getChildren()[i]);
      }
      delete curr;
    }
  }

  float getC2() const { return c2; }
  float getBeta2() const { return beta2; }
  float getMinRadius2() const { return minradius2; }

  bool exists(const MP* p) { return (bottom.find(p) != bottom.end()); }
  void insert(const MP* p);
  void remove(const MP* p);
  void update(float time);
  void verify() const;
  void printstat() const;
  int numBottomEdges() const;

  //Dan Farrell added member function getBottomLevelNeighbors.
  ////NOTE: this member function avoids double counting by only
  ////outputting the neighbors j of a node i that satisfy i < j.
  //NOTE SAW has removed this feature 11 Sept. 06
  //so now returns all neighbors
  //Also, dch must be up-to-date before using this function,
  //meaning that if points have moved you should call 'update'
  //before using this function to query neighbors.
  void getBottomLevelNeighbors(const MP* p, float distance, vector<int>& outputPoints ) const;

private:
  float getTime() const { return now; }
  const PQueue& getCerts() const { return certs; }
  PQueue& getCerts() { return certs; }
  const NodeMP* getRoot() const { return root; }
  NodeMP* getRoot() { return root; }
  void setRoot(NodeMP* root_) { root = root_; }

  void extendroot();
  void collapseroot();

  friend class Node<MP>;
  friend ostream& operator<< <MP> (ostream& os, 
                                   const DiscreteHierarchy<MP>& hier);
};

template <class MP>
ostream& operator<< (ostream& os, const Node<MP>& node) {
  os << '(' << node.id << ':'
            << ((node.parent) ? node.parent->id : -1)
            << ',' << sqrt(node.radius2) << ')';
  for (unsigned int i = 0; i < node.neighbors.size(); i++)
    os << ' ' << node.neighbors[i]->id;
  os << " | ";
  for (unsigned int i = 0; i < node.children.size(); i++)
    os << ' ' << node.children[i]->id;
  return os;
}

template <class MP>
ostream& operator<< (ostream& os, const DiscreteHierarchy<MP>& hier) {
  os << "== Begin ==" << endl;
  os << "  minradius = " << sqrt(hier.minradius2)
     << ", c = " << sqrt(hier.c2) 
     << ", beta = " << sqrt(hier.beta2) << endl;
  queue< Node<MP>* > candidates;			// remain nodes
  candidates.push(hier.root);
  while (!candidates.empty()) {
    Node<MP>* curr = candidates.front();
    candidates.pop();
    os << *curr << endl;

    for (unsigned int i = 0; i < curr->getChildren().size(); i++) {
      candidates.push(curr->getChildren()[i]);
    }
  }

/*
  os << "  == Certs ==" << endl;
  for (unsigned int i = 0; i < hier.getCerts().size(); i++) {
    const Cert* cert = hier.getCerts().nth(i);
    os << "    " << cert << " * " << *cert << endl;
  }

  for (typename DiscreteHierarchy<MP>::NodeHash::const_iterator 
       iter = hier.bottom.begin(); iter != hier.bottom.end(); iter++) {
    os << "   bottom " << iter->second << endl;
  }
*/

  os << "== End ==" << endl;
  return os;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Searches all bottom-level neighbors of a given point to find
//   points that are closer than a given distance.  
//   The dch must be up-to-date before using this function,
//   meaning that if points have moved you should call 'update'
//   before using this function.
// Note:
//   The template parameter class MP must be a class that implements
//   "int getID() const"   
// Parameters:
//   p - the point whose neighbors you are requesting
//   targetDist - points that are closer than this distance will be output,
//      with some extra points that are actually farther than this distance
// Return Value List:
//   a vector containing the neighbors of p
////////////////////////////////////////////////////////////////////////////////
template <class MP>
void DiscreteHierarchy<MP>::getBottomLevelNeighbors(
    const MP* p, 
    float targetDist, 
    vector<int>& outputPoints ) const {
  
  if (targetDist > bottomLevelMaxDistance) {
    cerr << "Error, distance queried exceeds bottom level neighbor distance" << endl;
    cerr << targetDist << " > " << bottomLevelMaxDistance << endl;
    exit(0);
  }

  outputPoints.clear();
  
  typename NodeHash::const_iterator iter = bottom.find(p);
  myassert(iter != bottom.end(), "query of non-existing point", -1, -1);
  const Node<MP> *node = iter->second;
  const Node<MP> *parent = node->getParent();
  const Node<MP> *uncle;
  float targDistPlusBeta = targetDist + beta;
  float targDistPlusBeta2 = targDistPlusBeta*targDistPlusBeta;

  //add all bottom-level siblings to the output list,
  //since these siblings can be arbitrarily close.
  for (size_t j = 0; j < parent->getChildren().size(); j++) {
    outputPoints.push_back( parent->getChildren()[j]->getPoint()->getID() );
  }

  //for each uncle
  for (size_t j = 0; j < parent->getNeighbors().size(); j++) {
    uncle = parent->getNeighbors()[j];
    if ( node->getPoint()->dist2( *uncle->getPoint() )  > targDistPlusBeta2 ) continue;

    //for each cousin
    for (size_t k = 0; k < uncle->getChildren().size(); k++) {
      outputPoints.push_back( uncle->getChildren()[k]->getPoint()->getID() );
    }
  }
}

/*
 * Compute the number of edges at the bottom level
 * i.e. number of pairs of nodes that are within c * minradius
 * of each other.
 *
 * This function could be easily modified to list such pairs
 */
template <class MP>
int DiscreteHierarchy<MP>::numBottomEdges() const {
  int count = 0;
  float threshold2 = minradius2 * c2;

  for (typename NodeHash::const_iterator iter = bottom.begin(); 
       iter != bottom.end(); iter++) {
    Node<MP>* node = iter->second;

    // DWF changed ints to size_t
    // all pairs such that their parents are neighbors
    for (size_t j = 0; j < node->getParent()->getNeighbors().size(); j++) {
      Node<MP> *uncle = node->getParent()->getNeighbors()[j];
      for (size_t k = 0; k < uncle->getChildren().size(); k++) {
        Node<MP> *cousin = uncle->getChildren()[k];
        if (node->getId() > cousin->getId() && 
            node->dist2(*cousin) < threshold2) count++;
      }
    }

    // all pairs with the same parent
    for (size_t j = 0; j < node->getParent()->getChildren().size(); j++) {
      Node<MP> *sibling = node->getParent()->getChildren()[j];
      // we don't really have to check, as we know these pairs are neighbors
      if (node->getId() > sibling->getId() && 
          node->dist2(*sibling) < threshold2) count++;
    }
  }

  return count;
}

// extend the hierarchy upward one level
template <class MP>
void DiscreteHierarchy<MP>::extendroot() {
  Node<MP>* old = root;
  root = new Node<MP>(root->getPoint(), NULL, this, beta2*root->getRadius2());
  old->setParent(root, true);
}

// collapse the hierarchy upward one level if neccessary
template <class MP>
void DiscreteHierarchy<MP>::collapseroot() {
  if (root->getChildren().size() == 1) {
    Node<MP>* old = root;
    root = root->getChildren()[0];
    root->removeParent();
    delete old;
  }
}

template <class MP>
void DiscreteHierarchy<MP>::printstat() const {
  int count = 0;
  int totaledges = 0, edges = 0;
  vector< Node<MP>* > level;
  level.push_back(root);

  cout << "beta = " << sqrt(beta2) << endl;
  cout << "c = " << sqrt(c2) << endl;
  cout << "minradius = " << sqrt(minradius2) << endl; 

  while (level.size() != 0) {
    edges = 0;
    vector< Node<MP>* > nextlevel;

    if (!level[0]->isBottom()) {
      for (unsigned int i = 0; i < level.size(); i++) {
        // cout << "  " << *level[i] << endl;
        edges += (int) level[i]->getNeighbors().size();
        for (unsigned int j = 0; j < level[i]->getChildren().size(); j++) {
          nextlevel.push_back(level[i]->getChildren()[j]);
        }
      }
      edges /= 2;
    }

    cout << "Level " << count << ": nodes " << (int) level.size() 
          << ", edges: " << edges 
          << ", avg degree " << (1.0*edges/level.size()) 
          << ", radius: " << sqrt(level[0]->getRadius2()) 
          << endl;
    count++;
    level = nextlevel;
    totaledges += edges;
  }
  
  cout << "Total number of edges (without bottom edges) " << totaledges << endl;
  cout << "Number of certs: " << certs.size() << endl;
}

// only pointer to the point inserted is used, the point is NOT copied
template <class MP>
void DiscreteHierarchy<MP>::insert(const MP* p) {
  // cout << "insert " << *p << endl;
  Node<MP> *node = new Node<MP>(p, NULL, this, minradius2);

  // deal with special cases first
  if (root == NULL) {				// this is the first node
    root = new Node<MP>(p, NULL, this, minradius2 * beta2);
    node->setParent(root, false);
    bottom[p] = node;
    return;
  }

  // extends root up if neccessary, making sure root is a potential parent
  float initd2 = root->dist2(*node);
  while (initd2 > root->getRadius2()) extendroot();

  // find the true parent of the node
  Node<MP>* parent = root;

  // a node with radius 1 has all its descendant within a distance of
  // beta/(beta-1). exp_ratio2 is the square of this value
  float exp_ratio2 = beta2 / (beta2 + 1 - 2*sqrt(beta2));

  queue< Node<MP>* > candidates; // list of potential parents
  candidates.push(root);

  // find parent, which is the potential parent with the lowest level
  while (!candidates.empty()) {
    Node<MP>* curr = candidates.front();
    candidates.pop();
    float d2 = node->dist2(*curr);
    if (d2 < curr->getRadius2()) parent = curr;
    if (curr->isNextToBottom()) continue;
    if (d2 < exp_ratio2 * curr->getRadius2()) {
      for (unsigned int i = 0; i < curr->getChildren().size(); i++) {
        candidates.push(curr->getChildren()[i]);
      }
    }
  }

  // insert the node and find all neighbors
  node->setRadius2(parent->getRadius2() / beta2);
  node->setParent(parent, true);

  // set parent and find neighbors in lower levels
  while (!node->isBottom()) {
    Node<MP>* next = new Node<MP>(node->getPoint(), NULL, 
                          this, node->getRadius2()/beta2);
    next->setParent(node, true);
    node = next;
  }

  // all nodes must be in bottom
  bottom[p] = node;

  // cout << *this << endl;
  // verify();
}

template <class MP>
void DiscreteHierarchy<MP>::remove(const MP* p) {
  typename NodeHash::iterator iter = bottom.find(p);
  myassert(iter != bottom.end(), "remove non-existing point", -1, -1);

  // remove the point from the bottom
  Node<MP> *node = iter->second;
  // cout << "remove " << *node << endl;

  Node<MP> *pnode = node->getParent();
  bottom.erase(iter);
  node->removeParent();
  delete node;

  // now remove the point on higher levels
  // all nodes of the point p has been removed in previous level
  while (pnode->getPoint() == p) {
    node = pnode;

    // handle special cases, deleting root node
    if (node == root) {
      if (node->getChildren().size() == 0) {
        // cout << "remove root" << *p << endl;
        root = NULL;
        delete node;
        break;
      } else if (node->getChildren().size() == 1 && 
                 !node->getChildren()[0]->isBottom()) {
        // cout << "change root" << *p << endl;
        root = node->getChildren()[0];
        root->removeParent();
        delete node;
        break;
      } else {
        extendroot();
      }
    }

    // remove the node from the hierarchy
    pnode = node->getParent();
    node->removeParent();

    vector<Node<MP>*> savedneighbors = node->getNeighbors();
    vector<Node<MP>*> savedchildren = node->getChildren();

    for (int i = (int) savedchildren.size()-1; i >= 0; i--)
      savedchildren[i]->removeParent();
    for (int i = (int) savedneighbors.size()-1; i >= 0; i--)
      node->removeNeighbor(savedneighbors[i]);

    delete node;

    // now find new parents for former children
    for (unsigned int i = 0; i < savedchildren.size(); i++) {
      Node<MP>* child = savedchildren[i];
      for (unsigned int j = 0; j < savedneighbors.size(); j++) {
        if (savedneighbors[j]->dist2(*child) < savedneighbors[j]->getRadius2()){
          child->setParent(savedneighbors[j], false);
          break;
        }
      }

      if (child->getParent() == NULL) { // orphan
        Node<MP>* upchild = new Node<MP>(child->getPoint(), NULL, this, 
                                 child->getRadius2() * beta2);
        upchild->findNewParent(pnode, true);
        child->setParent(upchild, false);
        savedneighbors.push_back(upchild);
      }
    }
  }

  // cout << *this << endl;
  // verify();
}

template <class MP>
void DiscreteHierarchy<MP>::verify() const {
  queue< Node<MP>* > candidates;  // remain nodes to be verified
  candidates.push(root);
  int count = 0;

  while (!candidates.empty()) {
    Node<MP>* curr = candidates.front();
    candidates.pop();
    if (curr->isBottom()) count++;

    for (unsigned int i = 0; i < curr->getChildren().size(); i++) {
      candidates.push(curr->getChildren()[i]);
    }

    // verify structure, parent relation
    if (curr != root) {
      myassert(curr->getParent() != NULL, "parent not null", curr->getId(), -1);
      // make sure that the cert exists
      myassert(curr->findParentCert() != NULL, "parent certs must exist", 
        curr->getId(), curr->getParent()->getId());
      if (curr->findParentCert() != NULL) {
        myassert(
          ((ParentCertMP*) certs.inf(curr->findParentCert()))->isOn(curr), 
          "valid parent cert stored", curr->getId(), 
          curr->getParent()->getId());
      }
    } else {
      myassert(curr->getParent() == NULL, "parent is null", 
        curr->getId(), curr->getParent()->getId());
      myassert(curr->findParentCert() == NULL, "parent cert is null", 
        curr->getId(), -1);
    }

    // verify parent distance relation
    if (curr != root) {
      myassert( curr->dist2(*curr->getParent()) 
                  <= 1.0001 * curr->getParent()->getRadius2(),
        "parent doesn't cover node", curr->getId(), curr->getParent()->getId());
    }

    // verify structure, children relation
    for (unsigned int i = 0; i < curr->getChildren().size(); i++) {
      myassert(curr->getChildren()[i]->getParent() == curr,
        "parent-child dual", curr->getId(), curr->getChildren()[i]->getId());
    }

    // verify structure, neighbor relation
    for (unsigned int i = 0; i < curr->getNeighbors().size(); i++) {
      myassert(
        find(curr->getNeighbors()[i]->getNeighbors().begin(),
             curr->getNeighbors()[i]->getNeighbors().end(), curr) !=
          curr->getNeighbors()[i]->getNeighbors().end(), 
        "symmetric neighbor", curr->getId(), curr->getNeighbors()[i]->getId());
      myassert(curr->findEdgeCert(curr->getNeighbors()[i]) != NULL,
        "missing neighb cert", curr->getId(), curr->getNeighbors()[i]->getId());
    }

    // verify neighbor distance relation
    for (unsigned int i = 0; i < curr->getNeighbors().size(); i++) {
      if (curr->getId() < curr->getNeighbors()[i]->getId()) continue;
      myassert(curr->dist2(*curr->getNeighbors()[i]) < c2*curr->getRadius2(), 
           "neighbor near", curr->getId(), curr->getNeighbors()[i]->getId());
      if ( curr->dist2(*curr->getNeighbors()[i]) >= c2*curr->getRadius2() ) {
        cout << "dch time: " << now << endl;
        cout << "Difference: " << sqrt(curr->dist2(*curr->getNeighbors()[i])) - sqrt(c2*curr->getRadius2());
        cout << "  Distance: " << sqrt(c2*curr->getRadius2()) << endl;
        cout << "  EdgeCert key: " << curr->findEdgeCert( curr->getNeighbors()[i] )->key();
        cout << "  Root cert key: " << certs.firstKey();
        cout << "  Cert CBound: " << certs.inf(curr->findEdgeCert( curr->getNeighbors()[i] ))->getCBound( now ) << endl;
      }

      myassert(curr->dist2(*curr->getNeighbors()[i]) >= curr->getRadius2(), 
          "neighbor density", curr->getId(), curr->getNeighbors()[i]->getId());
      if ( curr->dist2(*curr->getNeighbors()[i]) < curr->getRadius2() ) {
        cout << "dch time: " << now << endl;
        cout << "Difference: " << sqrt(curr->getRadius2()) - sqrt(curr->dist2(*curr->getNeighbors()[i]));
        cout << "  Distance: " << sqrt(curr->getRadius2()) << endl;
        cout << "  EdgeCert key: " << curr->findEdgeCert( curr->getNeighbors()[i] )->key();
        cout << "  Root cert key: " << certs.firstKey();
        cout << "  Cert CBound: " << certs.inf(curr->findEdgeCert( curr->getNeighbors()[i] ))->getCBound( now ) << endl;
      }
    }

    // verify structure, neighbor and potential neighbor complimentary relation
    if (curr != root && !curr->isBottom()) {
      for (unsigned int i = 0; i < curr->getParent()->getNeighbors().size(); 
              i++) {
        Node<MP>* nb = curr->getParent()->getNeighbors()[i];
        for (unsigned int j = 0; j < nb->getChildren().size(); j++) {
          Node<MP>* cousin = nb->getChildren()[j];
          if (curr->findEdgeCert(cousin) == NULL) {
            myassert(!curr->isNeighbor(cousin),
              "nonneighbor associate with no cert", 
              curr->getId(), cousin->getId());
            myassert(curr->findNEdgeCert(cousin) != NULL, 
              "non-neighbor exists", curr->getId(), cousin->getId());
          } else {
            myassert(curr->isNeighbor(cousin),
              "neighbor associate with cert", curr->getId(), cousin->getId());
            myassert(curr->findNEdgeCert(cousin) == NULL, 
              "both neighbor and non-neighbor", curr->getId(), cousin->getId());
          }
        }
      }

      for (unsigned int i = 0; i < curr->getParent()->getChildren().size(); 
             i++) {
        Node<MP>* cousin = curr->getParent()->getChildren()[i];
        if (cousin == curr) continue;
        if (curr->findEdgeCert(cousin) == NULL) {
          myassert(!curr->isNeighbor(cousin),
            "nonneighbor associate with no cert", 
            curr->getId(), cousin->getId());
          myassert(curr->findNEdgeCert(cousin) != NULL, 
            "non-neighbor exists", curr->getId(), cousin->getId());
        } else {
          myassert(curr->isNeighbor(cousin),
            "neighbor associate with cert", curr->getId(), cousin->getId());
          myassert(curr->findNEdgeCert(cousin) == NULL, 
            "both neighbor and non-neighbor", curr->getId(), cousin->getId());
        }
      }
    } else {
      assert(curr->getNeighbors().size() == 0);
    }
  }
  //DWF cast bottom.size() as int
  myassert(count == (int)bottom.size(), "inconsistent bottom", count, bottom.size());
}

template <class MP>
void Node<MP>::addNeighbor(Node<MP>* neighbor) {
  // cerr << "  addneighbor " << id << " " << neighbor->id << endl;
  neighbors.push_back(neighbor);
  neighbor->neighbors.push_back(this);
  EdgeCertMP *cert = new EdgeCertMP(this, neighbor);
  PQLink pqlink = dh->getCerts().insert(cert->getCBound(dh->getTime()), cert);
  e_certs.push_back(pqlink);
  neighbor->e_certs.push_back(pqlink);
}

template <class MP>
void Node<MP>::addPotentialNeighbor(Node<MP>* neighbor) {
  NEdgeCertMP *cert = new NEdgeCertMP(this, neighbor);
  PQLink pqlink = dh->getCerts().insert(cert->getCBound(dh->getTime()), cert);
  ne_certs.push_back(pqlink);
  neighbor->ne_certs.push_back(pqlink);
}

template <class MP>
void Node<MP>::removeNeighbor(Node<MP>* neighbor) {
  // cerr << "  removeneighbor " << id << " " << neighbor->id << endl;

  full_remove(neighbors, neighbor);
  full_remove(neighbor->neighbors, this);
  PQLink pqlink = findEdgeCert(neighbor);
  dh->getCerts().erase(pqlink);
  full_remove(e_certs, pqlink);
  full_remove(neighbor->e_certs, pqlink);
}

template <class MP>
void Node<MP>::removePotentialNeighbor(Node<MP>* neighbor) {
  PQLink pqlink = findNEdgeCert(neighbor);
  dh->getCerts().erase(pqlink);
  full_remove(ne_certs, pqlink);
  full_remove(neighbor->ne_certs, pqlink);
}

template <class MP>
PriorityQNode< Cert >* 
Node<MP>::findEdgeCert(Node<MP>* neighbor) const {
  for (unsigned int i = 0; i < e_certs.size(); i++) {
    PQLink pqlink = e_certs[i];
    EdgeCertMP* cert = (EdgeCertMP*) dh->getCerts().inf(pqlink);
    if (cert->isOn(neighbor)) return pqlink;
  }
  return NULL;
}

template <class MP>
PriorityQNode< Cert >* 
Node<MP>::findNEdgeCert(Node<MP>* neighbor) const {
  for (unsigned int i = 0; i < ne_certs.size(); i++) {
    PQLink pqlink = ne_certs[i];
    NEdgeCertMP* cert = (NEdgeCertMP*) dh->getCerts().inf(pqlink);
    if (cert->isOn(neighbor)) return pqlink;
  }
  return NULL;
}

template <class MP>
void DiscreteHierarchy<MP>::update(float time) {
  now = time;
  if ( countUpdates == INT_MAX ) {
    countUpdates = 0;
  }
  else {
    countUpdates++;
  }
  while (true) {
    PQLink pqlink = certs.find_min();
    if (certs.prio(pqlink) >= now) break;

    float newbound = certs.inf(pqlink)->getCBound(now);
    
    //DWF Comments:
    //the newbound just calculated
    //may appear to be equal to "now" due to finite precision even if
    //the newbound may actually be less than "now",
    //(meaning that the certificate is invalid).  This
    //can happen whenever a certificate goes invalid,
    //if the amount by which the certificate is invalid is less than
    //epsilon*now.  For type "float", epsilon is about 1e-7.
    //If you desire to maintain a certain time resolution,
    //say "res", then make sure that time does not surpass res/eps
    //(If you add res to time, but res < eps*time,
    //it will be as if you added nothing.)
    if (newbound >= now) {
      certs.update(pqlink, newbound);
    } else {
      // cout << "Handling: " << certs.inf(pqlink) << " * "
      //                      << *certs.inf(pqlink) << endl;
      certs.inf(pqlink)->update();
      collapseroot();
    }
  }
}

// assumption: parent = NULL
template <class MP>
void Node<MP>::findNewParent(Node<MP>* relaxedparent, bool newnode) {
  myassert(parent == NULL, "findNewParent", id, -1);

  Node<MP> *newparent = NULL;
  for (unsigned int i = 0; i < relaxedparent->neighbors.size(); i++) {
    if (dist2(*relaxedparent->neighbors[i]) < relaxedparent->radius2) {
      newparent = relaxedparent->neighbors[i];
      break;
    }
  }

  if (dist2(*relaxedparent) < relaxedparent->radius2) {
    newparent = relaxedparent;
  }

  if (newparent != NULL) { 		// change parent
    //cout << "  new parent for " << id << " is " << newparent->id << endl;
    setParent(newparent, newnode);
    return;
  }

  // no parent, have to create one
  newparent = new Node<MP>(getPoint(), NULL, dh, dh->getBeta2()*radius2);
  // cout << "  create new node for parent of " 
       // << id  << ": " << newparent->id << endl;


  // make sure relaxedparent is not the root, extend if neccessary
  if (relaxedparent == dh->getRoot()) dh->extendroot();
  newparent->findNewParent(relaxedparent->parent, true);

  setParent(newparent, newnode);
}

template <class MP>
void Node<MP>::removeParent() {
  if (!isBottom()) {
    for (unsigned int i = 0; i < parent->neighbors.size(); i++) {
      Node<MP> *olduncle = parent->neighbors[i];
      for (unsigned int j = 0; j < olduncle->children.size(); j++) {
        Node<MP> *oldcousin = olduncle->children[j];
        if (find(neighbors.begin(), neighbors.end(), oldcousin) 
              == neighbors.end()) {	// oldcousin is only a potential neighb
          removePotentialNeighbor(oldcousin);
        }
      }
    }
    for (unsigned int i = 0; i < parent->children.size(); i++) {
      Node<MP> *oldcousin = parent->children[i];
      if (oldcousin == this) continue;
      if (find(neighbors.begin(), neighbors.end(), oldcousin) 
            == neighbors.end()) {	// oldcousin is only a potential neighb
        removePotentialNeighbor(oldcousin);
      }
    }
  }

  dh->getCerts().erase(p_cert);
  full_remove(parent->children, this);
  parent = NULL;
  p_cert = NULL;
}

template <class MP>
void Node<MP>::setParent(Node<MP>* newparent, bool findNeighbors) {
  myassert(parent == NULL, "setparent", id, parent->id);
  parent = newparent;
  parent->children.push_back(this);

  ParentCertMP *cert = new ParentCertMP(this);
  p_cert = dh->getCerts().insert(cert->getCBound(dh->getTime()), cert);

  if (isBottom()) return;
  
  for (unsigned int j = 0; j < parent->getChildren().size(); j++) {
    Node<MP>* other = parent->getChildren()[j];
    if (other == this) continue;

    if (findNeighbors) {
      if (dist2(*other) < dh->getC2()*radius2) {
        addNeighbor(other);
      } else {
        addPotentialNeighbor(other);
      }
    } else {
      if (find(neighbors.begin(), neighbors.end(), other) 
            == neighbors.end()) {	// newcousin is only a potential neighb
        addPotentialNeighbor(other);
      }
    }
  }

  for (unsigned int i = 0; i < parent->neighbors.size(); i++) {
    Node<MP>* nb = parent->neighbors[i];
    for (unsigned int j = 0; j < nb->getChildren().size(); j++) {
      Node<MP>* cousin = nb->getChildren()[j];
      if (findNeighbors) {
        if (dist2(*cousin) < dh->getC2()*radius2) {
          addNeighbor(cousin);
        } else {
          addPotentialNeighbor(cousin);
        }
      } else {
        if (find(neighbors.begin(), neighbors.end(), cousin) 
              == neighbors.end()) {	// newcousin is only a potential neighb
          addPotentialNeighbor(cousin);
        }
      }
    }
  }
}

template <class MP>
void Node<MP>::fullRemoveNeighbor(Node<MP>& other) {
  removeNeighbor(&other);
  if (isNextToBottom()) return;

  for (unsigned int i = 0; i < children.size(); i++) {
    Node<MP>* child = children[i];
    for (unsigned int j = 0; j < other.children.size(); j++) {
      Node<MP>* otherchild = other.children[j];
      if (find(child->neighbors.begin(), child->neighbors.end(), 
               otherchild) == child->neighbors.end()) {
        child->removePotentialNeighbor(otherchild);
      }
    }
  }
}

template <class MP>
void Node<MP>::fullAddNeighbor(Node<MP>& other) {
  addNeighbor(&other);
  if (isNextToBottom()) return;

  for (unsigned int i = 0; i < children.size(); i++) {
    Node<MP>* child = children[i];
    for (unsigned int j = 0; j < other.children.size(); j++) {
      Node<MP>* otherchild = other.children[j];
      if (find(child->neighbors.begin(), child->neighbors.end(), 
               otherchild) == child->neighbors.end()) {
        child->addPotentialNeighbor(otherchild);
      }
    }
  }
}

template <class MP>
void Node<MP>::demote(Node<MP>& oldneighbor) {
  // cout << "  demote node " << *this << endl;
  myassert(!isBottom(), "demote bottom node", id, parent->id);

  // what to do when the nodes are TOO close?
  vector< Node<MP>* > savedchildren = children;
  vector< Node<MP>* > savedneighbors = neighbors;

  removeParent();
  myassert(ne_certs.size() == 0, "no more potential edges left", id, -1);

/*
  if (ne_certs.size() != 0) {
    for (unsigned int i = 0; i < ne_certs.size(); i++) {
      const Cert* cert = dh->getCerts().inf(ne_certs[i]);
      cerr << "    " << cert << " * " << *cert << endl;
    }
  }
*/

  for (int i = (int) savedchildren.size() - 1; i >= 0; i--) {
    savedchildren[i]->removeParent();
  }
  for (int i = (int) savedneighbors.size() - 1; i >= 0; i--) {
    removeNeighbor(savedneighbors[i]);	// no need to use fullRemoveNeighbor
  }

  myassert(e_certs.size() == 0, "no more edges left", id, -1);

  for (unsigned int i = 0; i < savedchildren.size(); i++) {
    savedchildren[i]->findNewParent(&oldneighbor);
  }

  delete this;
}
