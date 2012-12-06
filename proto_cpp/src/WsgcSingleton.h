/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>

     This file is part of WSGC. A Weak Signal transmission mode using Gold Codes

     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 51 Franklin Street, Boston, MA  02110-1301  USA

     Static not real time prototype in C++

     WsgcSingleton

     Implementation of the singleton pattern
*/

#ifndef __WSGC_SINGLETON_H__
#define __WSGC_SINGLETON_H__

template <typename T>
class WsgcSingleton
{
protected:
  // Constructeur/destructeur
  WsgcSingleton () { }
  ~WsgcSingleton () {}

public:
  // Interface publique
  static T *getInstance ()
  {
    if (0 == _singleton)
      {
        _singleton = new T;
      }

    return (static_cast<T*> (_singleton));
  }

  static void kill ()
  {
    if (0 != _singleton)
      {
        delete _singleton;
        _singleton = 0;
      }
  }

private:
  // Unique instance
  static T *_singleton;
};

template <typename T>
T *WsgcSingleton<T>::_singleton = 0;
#endif /* WSGCSINGLETON_H_ */
