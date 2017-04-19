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
 * This class represents a single 3D point/vector
 */


#ifndef __Vect_H
#define __Vect_H

#include <iostream>

template <class T>
class Vect {
public:
  T x,y,z;
  Vect() : x(0), y(0), z(0) {}
  Vect(T *v) : x(v[0]), y(v[1]), z(v[2]) {}
  Vect(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}
  Vect(const Vect<T>& p) : x(p.x), y(p.y), z(p.z) {}

  ~Vect(){};

  Vect<T>& operator=(const Vect<T>& p) {
    x = p.x; y = p.y; z = p.z;
    return *this;
  }

  Vect<T> operator+(const Vect<T>& other) const {
    return Vect<T>(x + other.x, y + other.y, z + other.z);
  }
  Vect<T> operator-(const Vect<T>& other) const {
    return Vect<T>(x - other.x, y - other.y, z - other.z);
  }
  void operator+=(const Vect<T>& other) {
    x += other.x; y += other.y; z += other.z;
  }
  void operator-=(const Vect<T>& other) {
    x -= other.x; y -= other.y; z -= other.z;
  }
  void operator*=(T f) {
    x *= f; y *= f; z *= f;
  }
  void operator/=(T f) {
    x /= f; y /= f; z /= f;
  }
  Vect<T> operator*(T f) const {
    return Vect<T>(x*f, y*f, z*f);
  }

  T dot(const Vect<T> other) const {
    return x*other.x + y*other.y + z*other.z;
  }

  T norm2() const {
    return x*x + y*y + z*z;
  }

  // compute distance square
  T dist2(const Vect<T>&p) const {
    return (p.x-x)*(p.x-x)+(p.y-y)*(p.y-y)+(p.z-z)*(p.z-z);
  }

  friend std::ostream& operator<< (std::ostream& os, const Vect<T>& p) {
    os << '(' << p.x << ',' << p.y << ',' << p.z << ')';
    return os;
  }
};
#endif
