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
 * Simple timer, measuring user and system time 
 * between successive start/stop calls
 */

#ifndef TIMER_H 
#define TIMER_H

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;

// time unit is millisecond
class Timer {
  long utime, stime;
  struct rusage startu, pauseu;
  bool running;
public:
  Timer() { reset(); }
  long gettime() { 
    if (!running) return utime; 
    stop(); start();
    return utime + stime;
  }
  long getutime() { 
    if (!running) return utime; 
    stop(); start();
    return utime;
  }
  long getstime() {
    if (!running) return stime; 
    stop(); start();
    return stime;
  }
  void start() {
    if (!running) {
      running = true;
      getrusage(RUSAGE_SELF, &startu);
    }
  }
  void stop() {
    if (running) {
      running = false;
      getrusage(RUSAGE_SELF, &pauseu);
      utime += ( (pauseu.ru_utime.tv_sec-startu.ru_utime.tv_sec) * 1000+ 
                 (pauseu.ru_utime.tv_usec-startu.ru_utime.tv_usec) / 1000 );
      stime += ( (pauseu.ru_stime.tv_sec-startu.ru_stime.tv_sec) * 1000+ 
                 (pauseu.ru_stime.tv_usec-startu.ru_stime.tv_usec) / 1000 );
    }
  }
  void reset() {
    running = false;
    utime = stime = 0;
  }
};

#endif
