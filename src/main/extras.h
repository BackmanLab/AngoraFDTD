/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EXTRAS_H
#define EXTRAS_H

//template declaration for the function that reads a material region from a file

#include "headers.h"

//base Angora exception class
#include "angora_excp.h"
// https://github.com/lteacy/maxsum-cpp/blob/master/include/maxsum/common.h

int * ind2sub(int matSize[], int index){
  int ndims = 3;
  static int sub[3];

  int cumProd[ndims];
  int nextProd = 1;
  for(int j = 0; j<ndims; j++){
    cumProd[j] = nextProd;
    nextProd *= matSize[j];
  }


  int rem = index;
  for(int k = ndims-1; k>=0; k--)
  {
     sub[k] = rem / cumProd[k];
     rem = rem % cumProd[k];
   }

  return sub;
}


int * ind2sub(int matSize[], int index, int order[]){
  int ndims = 3;
  static int sub[3];
  static int subReordered[3];

  int cumProd[ndims];
  int nextProd = 1;
  for(int j = 0; j<ndims; j++){
    cumProd[j] = nextProd;
    nextProd *= matSize[j];
  }

  int rem = index;
  for(int k = ndims-1; k>=0; k--)
  {
     sub[k] = rem / cumProd[k];
     rem = rem % cumProd[k];
   }

 for(int j = 0; j<ndims; j++){
   subReordered[j] = sub[order[j]];
 }
  return subReordered;
}

bool poolOnGLobalEdge(int currentLoc[], int gridSizes[]){
  static bool anyOnEdge = false;
  for(int i=0; i<3; i++){
    if(currentLoc[i] == gridSizes[i]){
      anyOnEdge = true;
      return anyOnEdge;
    }
  }
  return anyOnEdge;
}

bool * isDimOnGlobalEdge(int currentLoc[], int gridSizes[]){
  static bool onEdge[3];
  for(int i=0; i<3; i++){
    onEdge[i] = currentLoc[i] == gridSizes[i];
  }
  return onEdge;
}

bool * isDimNotOnGlobalEdge(int currentLoc[], int gridSizes[]){
  static bool onEdge[3];
  for(int i=0; i<3; i++){
    onEdge[i] = currentLoc[i] < gridSizes[i];
  }
  return onEdge;
}

bool isWithinBounds(int currentPos[], int stPos[], int edPos[]){
  static bool isIn = 1;
  for(int i=0; i<3; i++){
    if((currentPos[i] < stPos[i]) && (currentPos[i] > edPos[i])){
      isIn = 0;
    }
  }
  return isIn;
}


bool isWithinBounds(int currentPos, int stPos, int edPos){
  if((currentPos < stPos) || (currentPos > edPos)){
    return 0;
  }
  return 1;
}

/*
int maxSize(int stPos[], int endPos[], int npts){
  int sz = endPos[0] - stPos[0] + 1;
  for(int i=1; i<npts; i++){
    int newSz = endPos[i] - stPos[i] + 1;
    if(newSz > sz){
      sz = newSz;
    }
  }
  return sz;
}
*/

int maxSize(int *stPos, int *endPos, unsigned int n){
  static int sz = 0;
  for(int i=0; i<n; i++){
    int newSz = endPos[i] - stPos[i] + 1;
    if(newSz > sz){
      sz = newSz;
    }
    //stPos++;
    //endPos++;
  }
  return sz;
}


double find_max(double a[ ],int len){
   double max; /* Current max */

   max = a[0];
   for (int i=1;i<len;++i)
     if (a[i] > max)
        max = a[i];
   return max;
}


double find_min(double a[ ],int len){
   double min; /* Current max */
   min = a[0];
   for (int i=1;i<len;++i)
     if (a[i] < min)
        min = a[i];
   return min;
}

#endif // EXTRAS_H
