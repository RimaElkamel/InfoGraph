#include "Vector.h"
#include <vector>

 Vector operator+(const Vector& A, const Vector& B ){
     return Vector(A[0]+B[0], A[1]+B[1], A[2]+B[2]);
 }
 Vector operator-(const Vector& A, const Vector& B ){
     return Vector(A[0]-B[0], A[1]-B[1], A[2]-B[2]);
 }

 Vector operator*(double k, const Vector& B ){
     return Vector(k*B[0], k*B[1], k*B[2]);
 }
 Vector operator*(const Vector& B , double k){
     return Vector(k*B[0], k*B[1], k*B[2]);
 }
 Vector operator*(const Vector& B , const Vector& A){
     return Vector(A[0]*B[0], A[1]*B[1], A[2]*B[2]);
 }
 Vector operator/(const Vector& B , double k){
     return Vector(B[0]/k, B[1]/k, B[2]/k);
 }

 double dot(const Vector& A, const Vector& B){
     return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
 }

 Vector cr(const Vector& A, const Vector& B){
     return Vector(A[1]*B[2]-A[2]*B[1],A[2]*B[0]- A[0]*B[2],A[0]*B[1]-A[1]*B[0]);
 }
