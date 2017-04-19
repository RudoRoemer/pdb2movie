#ifndef _MOVINGPOINT_H_
#define _MOVINGPOINT_H_

#include <cmath>
#include "Vect.h"
typedef Vect<double> Vec3;

class MovingPoint {
public:
  MovingPoint();
  MovingPoint(const MovingPoint& other);
  MovingPoint( const Vec3 *position_, int id_ );
  virtual ~MovingPoint();
  const MovingPoint& operator=(const MovingPoint& other);
  float distFromPastPosition2() const;
  float dist2(const MovingPoint& other) const;
  Vec3 delta( const MovingPoint& other ) const;
  float firstTimeNearer(float now, const MovingPoint& other, float threshold2) const;
  float firstTimeFurther(float now, const MovingPoint& other, float threshold2) const;
  void updatePastPosition();
  int getID() const {return id;}
  const Vec3 *getPosition() const {return position;}
  
protected:
  const Vec3 *position;
  Vec3 pastPosition;
  int id;
};

inline MovingPoint::MovingPoint() : position(NULL), pastPosition(), id(0) {}

inline MovingPoint::MovingPoint(const MovingPoint& other) : position(other.position),
                                                   pastPosition(other.pastPosition),
                                                   id(other.id) {}
                                                   
inline MovingPoint::MovingPoint( const Vec3 *position_, int id_ ) : 
  position(position_), pastPosition( *position_ ), id(id_) {}


inline const MovingPoint& MovingPoint::operator=(const MovingPoint& other) {
  position = other.position;
  pastPosition = other.pastPosition;
  id = other.id;
  return other;
}


inline MovingPoint::~MovingPoint() {
}

inline void MovingPoint::updatePastPosition() {
  if ( position != NULL )
  pastPosition = *position;
}

inline float MovingPoint::distFromPastPosition2() const {
  return (float) position->dist2( pastPosition );
}

inline float MovingPoint::dist2(const MovingPoint& other) const {
  return (float) position->dist2( *other.position );
}

inline float MovingPoint::firstTimeNearer(float now, const MovingPoint& other, float threshold2) const {
  return now + ( sqrt(dist2(other)) - sqrt(threshold2) ) / 2.0f;
}

inline float MovingPoint::firstTimeFurther(float now, const MovingPoint& other, float threshold2) const {
  return now + ( sqrt(threshold2) - sqrt(dist2(other)) ) / 2.0f;
}

inline Vec3 MovingPoint::delta( const MovingPoint& other ) const {
  return *other.position - *position;
}


#endif 
