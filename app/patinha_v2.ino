#include <cmath>
#include <vector>
#include <ArduinoEigen.h>

using namespace Eigen;

#define L1 82.0
#define L2 104.0
#define L3 98.0
#define L4 28.0
#define L5 50.0
#define OFFSET1 (M_PI/6)
#define OFFSET2 (M_PI/6)
#define OFFSET4 (-M_PI/2)

std::vector<float> o;
Matrix<float, 4, 4> pose;


void forward_kine(Matrix<float,4,4> &T, std::vector<float> o){
  T(0,0) = cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[0]+OFFSET1)*cos(o[4]) - sin(o[0]+OFFSET1)*sin(o[4]);
  T(0,1) = -cos(o[4])*sin(o[0]+OFFSET1) - cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[0]+OFFSET1)*sin(o[4]);
  T(0,2) = -sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[0]+OFFSET1);
  T(0,3) = cos(o[0]+OFFSET1)*(L3*cos(o[0]+OFFSET1+o[2]) + L2*cos(o[1]+OFFSET2) + L4*cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4) - L5*sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4));

  T(1,0) = cos(o[0]+OFFSET1)*sin(o[4]) + cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[4])*sin(o[0]+OFFSET1);
  T(1,1) = cos(o[0]+OFFSET1)*cos(o[4]) - cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*sin(o[0]+OFFSET1)*sin(o[4]);
  T(1,2) = -sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*sin(o[0]+OFFSET1);
  T(1,3) = sin(o[0]+OFFSET1)*(L3*cos(o[0]+OFFSET1+o[2]) + L2*cos(o[1]+OFFSET2) + L4*cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4) - L5*sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4));

  T(2,0) = sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[4]);
  T(2,1) = -sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*sin(o[4]);
  T(2,2) = cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4);
  T(2,3) = L1 + L3*sin(o[1]+OFFSET2+o[2]) + L2*sin(o[1]+OFFSET2) + L5*cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4) + L4*sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4);

  T(3,0) = 0;
  T(3,1) = 0;
  T(3,2) = 0;
  T(3,3) = 1;
}


std::vector<float> inverse_kine(Matrix<float,4,4> T){
  float o1, o2, o3, o4, o5, o234, C3, S3, num2, den2;
  std::vector<float> angles;
  
  o1 = atan2(T(1,3), T(0,3));
  while(o1 < 0)
    o1 = o1 + M_PI;

  o5 = -atan2(T(3,1),T(3,0));
  while(o5 < 0)
    o5 = o5 + M_PI;
  
  o234 = atan2(-(T(0,2)*cos(o1) + T(1,2)*sin(o1)), T(2,2));
  
  C3 = ( (pow(T(0,3)*cos(o1) + T(1,3)*sin(o1) - L4*cos(o234) + L5*sin(o234), 2) + pow(T(2,3)-L1-L5*cos(o234)-L4*sin(o234), 2) - pow(L2, 2) - pow(L3, 2)) / (2*L2*L3));
  S3 = sqrt(abs(1.0 - pow(C3, 2)));
  o3 = atan2(S3, C3);
  while(o3 < 0)
    o3 = o3 + M_PI;

  num2 = (T(2,3)-L1-L4*sin(o234)-L5*cos(o234))*(L3*cos(o3)+L2) - (L3*sin(o3))*(T(0,3)*cos(o1)+T(1,3)*sin(o1)-L4*cos(o234)+L5*sin(o234));
  den2 = (T(2,3)-L1-L4*sin(o234)-L5*cos(o234))*(L3*sin(o3)) + (L3*cos(o3)+L2)*(T(0,3)*cos(o1)+T(1,3)*sin(o1)-L4*cos(o234)+L5*sin(o234));
  o2 = atan2(num2,den2);
  while(o2 < 0)
    o2 = o2 + M_PI;  

  o4 = o234 - o2 - o3;
  while(o4 < 0)
    o4 = o4 + M_PI;   
  
  while(o1 > M_PI)
    o1 = o1 - M_PI;
  while(o2 > M_PI)
    o2 = o2 - M_PI;
  while(o3 > M_PI)
    o3 = o3 - M_PI;
  while(o4 > M_PI)
    o4 = o4 - M_PI;
  while(o5 > M_PI)
    o5 = o5 - M_PI;
  
  Serial.printf("\n\no1 = %f \n", o1);
  Serial.printf("\no2 = %f \n", o2);
  Serial.printf("\no3 = %f \n", o3);
  Serial.printf("\no4 = %f \n", o4);
  Serial.printf("\no5 = %f \n", o5);  
  
  return {o1,o2,o3,o4,o5}; 
}


void calculate_diff_operator(Matrix<float,4,4> &DELTA){
  float dx, dy, dz, drx, dry, drz;
 
  /*
  DELTA = {
    {0, -drz, dry, dx},
    {drz, 0, -drx, dy},
    {-dry, drx, 0, dz},
    {0,   0,   0,   0}
  };*/
  DELTA(0,0) = 0;
  DELTA(0,1) = -drz;
  DELTA(0,2) = dry;
  DELTA(0,3) = dx;
  DELTA(1,0) = drz;
  DELTA(1,1) = 0;
  DELTA(1,2) = -drx;
  DELTA(1,3) = dy;
  DELTA(2,0) = -dry;
  DELTA(2,1) = drx;
  DELTA(2,2) = 0;
  DELTA(2,3) = dz;
  DELTA(3,0) = 0;
  DELTA(3,1) = 0;
  DELTA(3,2) = 0;
  DELTA(3,3) = 0;
}


void calculate_new_pose(Matrix<float,4,4> T, Matrix<float,4,4> DELTA){
  Matrix<float,4,4> dT;
  dT = DELTA * T;
  T = dT + T;
}


void setup(){
  Serial.begin(115200);
  Serial.println("hemlo!");
  o = {0,0,0,0,0};
  forward_kine(pose, o);
  for(int i=0; i<4; i++){
    Serial.print("| ");
    for(int j=0; j<4; j++)
      Serial.printf("%f ", pose(i,j));
    Serial.print("| \n");
  }
  inverse_kine(pose);
}

void loop() {
  delay(10000);
  Serial.println("hemlo");
}
