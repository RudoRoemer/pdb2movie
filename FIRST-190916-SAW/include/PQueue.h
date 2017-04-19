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
 * An adaptation of the standard priority queue code
 * to be used as a queue of certificates in Kinetic Data Structure (KDS).
 */

#ifndef _PQUEUE_H
#define _PQUEUE_H

using namespace std;

template <class _T> 
class PriorityQ;

template <class _T>
class PriorityQNode 
{
  friend class PriorityQ<_T>;

private:

  // If mode = 1, delete the data pointed to when this node is destroyed.
  // If mode = 0, don't delete the data.
  int mode;

  PriorityQNode<_T>* parent;
  PriorityQNode<_T>* left;
  PriorityQNode<_T>* rght;

  _T *data;
  float key_;

public:

  PriorityQNode(_T *d, const float &k, int m) {data = d; key_ = k; mode = m;}
  PriorityQNode() {}
  ~PriorityQNode() {
	  if (mode) delete data;
  }

  const float &key() const {return key_;}
  bool isleaf() const { return (left == NULL) && (rght == NULL); }
  bool isroot() const { return parent == NULL; }
};


template <class _T>
class PriorityQ {

public:

  typedef PriorityQNode<_T>* link_type;

private:

  unsigned int _size;  // number of nodes in queue
  unsigned int mask;  // 2^n, where n is max depth of tree (root is at depth 0)

  link_type root;

  //list< PriorityQNode<_T> > q;

  // The mode assigned to all PriorityQNodes in this PriorityQ.  If
  // and only if mode is TRUE, the data pointed to the nodes will be
  // deleted when the nodes are destroyed.  Nodes are destroyed when
  // the PriorityQ is cleared or destroyed, or e.g. when the
  // application deletes a PriorityQNode after removing from the
  // PriorityQ.  Mode is TRUE by default.
  int mode; 

  void promote(link_type const kid);


  void recursDelete(link_type node)
  {
    if (!node) return;
    recursDelete(node->left);
    recursDelete(node->rght);
    delete node;
  }

public:

  void clear() {if (mode) recursDelete(root); _size = mask = 0; root = NULL;}

  PriorityQ(int m = 1) {_size = mask = 0; root = NULL; mode = m;}

  ~PriorityQ() {if (mode) recursDelete(root);}

  int empty() const {return (_size == 0);}

  int length() const {return _size;}
  int size() const { return _size;}

  link_type insert(float key, _T* data = NULL);

  link_type removeFirst();

  void eraseFirst() {delete removeFirst();}

  void remove(link_type node);  // remove arbitrary node

  void erase(const link_type node) {  // remove & delete arbtrary node
    remove(node); 
	delete node;
  }

  const _T *first() const {return root->data;} // return 1st node w/o removal
        _T *first()       {return root->data;}

  const float &firstKey() const {return root->key_;}

  void update(const link_type node, const float &newkey);

  // return node in nth tree position, 0 <= n < _size
  const _T *nth(int n) const {
    unsigned int bit;
    link_type node;
    if (n < 0 || n >= _size) return NULL;
    n++;
    for (bit = mask; !(n & bit); bit >>= 1);
    for (bit >>= 1, node = root; bit; bit >>= 1)
      node = (n & bit) ? node->rght : node->left;
    return node->data;
  }

  link_type find_min() const { return root;}
  float prio(link_type a) const { return a->key_;}
  _T *inf(link_type a) const { return a->data;} 
};

template <class _T> 
void PriorityQ<_T>::promote(PriorityQNode<_T> *const kid)
{
  PriorityQNode<_T> *par, *tmp;

  par = kid->parent;
  if ( (kid->parent = par->parent) ) 
    if (par->parent->left == par) par->parent->left = kid; 
    else par->parent->rght = kid;
  else root = kid;
  par->parent = kid;

  if (par->left == kid) {
    if ( (par->left = kid->left) ) par->left->parent = par;
    kid->left = par;
    tmp = par->rght;
    if ( (par->rght = kid->rght) ) par->rght->parent = par;
    if ( (kid->rght = tmp) ) kid->rght->parent = kid;
  }
  else { 
    if ( (par->rght = kid->rght) ) par->rght->parent = par;
    kid->rght = par;
    tmp = par->left;
    if ( (par->left = kid->left) ) par->left->parent = par;
    if ( (kid->left = tmp) ) kid->left->parent = kid;
  }
}

template <class _T> 
PriorityQNode<_T>* 
PriorityQ<_T>::insert(float key, _T* data)
{
  PriorityQNode<_T> *par;
  unsigned int bit;
  PriorityQNode<_T> *node;
  
  node = new PriorityQNode<_T>(data, key, mode);
  node->left = node->rght = NULL;

  if (!(++_size & mask)) mask = _size;

  for (bit = mask >> 1, par = root; bit > 1; bit >>= 1)
    par = (_size & bit) ? par->rght : par->left;
  node->parent = par;
  if (par)
    if (_size & 1) par->rght = node; else par->left = node;
  else root = node;
    
  // percolate up
  while (node->parent && node->key_ < node->parent->key_) promote(node);

  return node;
}


// remove node from front of priority queue
template <class _T> 
PriorityQNode<_T> *PriorityQ<_T>::removeFirst()
{
  unsigned int bit;
  PriorityQNode<_T> *node, *minKid, *oldRoot;

  oldRoot = root;

  if (_size <= 1) {
    _size = mask = 0;
    root = NULL;
  }
  else {
    // locate last node in Q
    for (bit = mask >> 1, node = root; bit > 0; bit >>= 1)
      node = (_size & bit) ? node->rght : node->left;

    // pop bottom element up to root position ...
    if (_size & 1) node->parent->rght = NULL; else node->parent->left = NULL;
    if (node->left = root->left) root->left->parent = node;
    if (node->rght = root->rght) root->rght->parent = node;
    node->parent = NULL;
    root = node;
    if (!(--_size & mask)) mask >>= 1;

    // ... push down if necessary
    while (node->left) {
      minKid = (node->rght && node->rght->key_ < node->left->key_) ? 
	node->rght : node->left;
      if (node->key_ <= minKid->key_) break;
      promote(minKid);
    }
  }
  
  return oldRoot;
}


// Remove an arbitrary node from priority queue.  Obviously, the
// PriorityQ is not empty when this routine is called.
template <class _T> 
void 
PriorityQ<_T>::remove(typename PriorityQ<_T>::link_type node)
{
  unsigned int bit;
  PriorityQNode<_T> *last, *minKid;

  // locate last node in Q
  for (bit = mask >> 1, last = root; bit > 0; bit >>= 1)
    last = (_size & bit) ? last->rght : last->left;

  if (last == node) {
    // we're just removing the last node - easy
    if (last->parent)
      if (_size & 1) last->parent->rght = NULL; 
      else last->parent->left = NULL;
    else root = NULL;
    if (!(--_size & mask)) mask >>= 1;
    return;
  }

  // pop last node up to deleted node's position ...
  // (at this point, we know last node is not the root)
  if (_size & 1) last->parent->rght = NULL; else last->parent->left = NULL;
  if ( (last->left = node->left) ) last->left->parent = last;
  if ( (last->rght = node->rght) ) last->rght->parent = last;
  if ( (last->parent = node->parent) ) 
    if (node->parent->left == node) last->parent->left = last;
    else last->parent->rght = last;
  else root = last;
  if (!(--_size & mask)) mask >>= 1;

  // Li debugged: need to promote before push down
  //
  if(last->key_ < node->key_) {   
    //... percolate up
    while(last->parent) {
      if(last->key_< last->parent->key_) 
	promote(last);
      else return;
    }
  } else {
    // ... push down if necessary
    while (last->left) {
      minKid = (last->rght && last->rght->key_ < last->left->key_) ? 
	last->rght : last->left;
      if (last->key_ <= minKid->key_) break;
      promote(minKid);
    }
  }
}

//  Update a priority queue node with its new key.  Percolate node up
//  or push node down, as appropriate.
template <class _T> 
void PriorityQ<_T>::update(PriorityQNode<_T> *const node, 
			     const float &newkey)
{
  int increase;
  PriorityQNode<_T> *minKid;

  increase = node->key_ < newkey;
  node->key_ = newkey;
  if (increase)
    // push down
    while (node->left) {
      minKid = (node->rght && node->rght->key_ < node->left->key_) ? 
	node->rght : node->left;
      if (node->key_ <= minKid->key_) break;
      promote(minKid);
    }
  else
    // percolate up
    while (node->parent && node->key_ < node->parent->key_) promote(node);
}

#endif
