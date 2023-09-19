#include <cmath>
#include <vector>
#include <ArduinoEigen.h>
#include <PS4Controller.h>

using namespace Eigen;

#define l1 82.0
#define l2 104.0
#define l3 98.0
#define l4 28.0
#define l5 50.0
#define OFFSET1 (M_PI/6)
#define OFFSET2 (M_PI/6)
#define OFFSET4 (-M_PI/2)

#define PWM_RESOLUTION 12
#define PWM_FREQUENCY  50
#define MG996R_MIN_POSITION  1*(4096*50/1000)   /* 1 milissecond pulse width */
#define MG996R_MID_POSITION  1.6*(4096*50/1000) /* 1.5 milissecond pulse width */
#define MG996R_MAX_POSITION  2.2*(4096*50/1000) /* 2 milisseconds pulse width */
#define MG996R_RAD_TO_PWM(a) (round((1.2*(4096*50/1000)/(2*M_PI/3))*a + (4096*50/1000)))
#define PS4_INT8_TO_DELTA(x) (1.0/63.0 * x)

#define SENSOR    32
#define GRIPPER_O 2
#define GRIPPER_S 15

/* mathematical model variables */
std::vector<float> o;
Matrix<float, 4, 4> pose;
Matrix<float, 4, 4> DELTA = Matrix4f::Zero();

/* ps4 controller variables*/
int8_t left_stick_x, left_stick_y, right_stick_y;
bool L1_bt, L2_bt, R1_bt, R2_bt;

void forward_kine(Matrix<float,4,4> &T, std::vector<float> o){
  T(0,0) = cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[0]+OFFSET1)*cos(o[4]) - sin(o[0]+OFFSET1)*sin(o[4]);
  T(0,1) = -cos(o[4])*sin(o[0]+OFFSET1) - cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[0]+OFFSET1)*sin(o[4]);
  T(0,2) = -sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[0]+OFFSET1);
  T(0,3) = cos(o[0]+OFFSET1)*(l3*cos(o[0]+OFFSET1+o[2]) + l2*cos(o[1]+OFFSET2) + l4*cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4) - l5*sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4));

  T(1,0) = cos(o[0]+OFFSET1)*sin(o[4]) + cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[4])*sin(o[0]+OFFSET1);
  T(1,1) = cos(o[0]+OFFSET1)*cos(o[4]) - cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*sin(o[0]+OFFSET1)*sin(o[4]);
  T(1,2) = -sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*sin(o[0]+OFFSET1);
  T(1,3) = sin(o[0]+OFFSET1)*(l3*cos(o[0]+OFFSET1+o[2]) + l2*cos(o[1]+OFFSET2) + l4*cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4) - l5*sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4));

  T(2,0) = sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*cos(o[4]);
  T(2,1) = -sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4)*sin(o[4]);
  T(2,2) = cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4);
  T(2,3) = l1 + l3*sin(o[1]+OFFSET2+o[2]) + l2*sin(o[1]+OFFSET2) + l5*cos(o[1]+OFFSET2+o[2]+o[3]+OFFSET4) + l4*sin(o[1]+OFFSET2+o[2]+o[3]+OFFSET4);

  T(3,0) = 0;
  T(3,1) = 0;
  T(3,2) = 0;
  T(3,3) = 1;
}


std::vector<float> inverse_kine(Matrix<float,4,4> T){
  float o1, o2, o3, o4, o5, o234, C3, S3, num2, den2;
  std::vector<float> old_values = o;
  
  o1 = atan2(T(1,3), T(0,3));
  while(o1 < 0)
    o1 = o1 + M_PI;

  o5 = -atan2(T(3,1),T(3,0));
  while(o5 < 0)
    o5 = o5 + M_PI;
  
  o234 = atan2(-(T(0,2)*cos(o1) + T(1,2)*sin(o1)), T(2,2));
  
  C3 = ( (pow(T(0,3)*cos(o1) + T(1,3)*sin(o1) - l4*cos(o234) + l5*sin(o234), 2) + pow(T(2,3)-l1-l5*cos(o234)-l4*sin(o234), 2) - pow(l2, 2) - pow(l3, 2)) / (2*l2*l3));
  S3 = sqrt(abs(1.0 - pow(C3, 2)));
  o3 = atan2(S3, C3);
  while(o3 < 0)
    o3 = o3 + M_PI;

  num2 = (T(2,3)-l1-l4*sin(o234)-l5*cos(o234))*(l3*cos(o3)+l2) - (l3*sin(o3))*(T(0,3)*cos(o1)+T(1,3)*sin(o1)-l4*cos(o234)+l5*sin(o234));
  den2 = (T(2,3)-l1-l4*sin(o234)-l5*cos(o234))*(l3*sin(o3)) + (l3*cos(o3)+l2)*(T(0,3)*cos(o1)+T(1,3)*sin(o1)-l4*cos(o234)+l5*sin(o234));
  o2 = atan2(num2,den2);
  while(o2 < 0)
    o2 = o2 + M_PI;  

  o4 = o234 - o2 - o3;
  while(o4 < 0)
    o4 = o4 + M_PI;   
  
  o1 = o1 - M_PI/6;
  o2 = o2 - M_PI/6;
  o4 = o4 + M_PI/2;

  while(o1 >= M_PI)
    o1 = o1 - M_PI;
  while(o2 >= M_PI)
    o2 = o2 - M_PI;
  while(o3 >= M_PI)
    o3 = o3 - M_PI;
  while(o4 >= M_PI)
    o4 = o4 - M_PI;
  while(o5 >= M_PI)
    o5 = o5 - M_PI;  

  if(o[0]<=2*M_PI/3 || o[1]<=2*M_PI/3 || o[2]<=2*M_PI/3 || o[3]<=2*M_PI/3 || o[4]<=2*M_PI/3 ||
     o[0]>=0 || o[1]>=0 || o[2]>=0 || o[3]>=0 || o[4]>=0)
    return {o1,o2,o3,o4,o5};

  return old_values;
}


void calculate_diff_operator(Matrix<float,4,4> &DELTA){
  float dx=0, dy=0, dz=0, drx=0, dry=0, drz=0;
  
  dz = PS4_INT8_TO_DELTA(left_stick_y);
  dy = PS4_INT8_TO_DELTA(left_stick_x);
  dx = PS4_INT8_TO_DELTA(right_stick_y);

  Serial.println(dx);
  Serial.println(dy);
  Serial.println(dz);
  
  DELTA << 0, -drz, dry, dx,
           drz, 0, -drx, dy,
           -dry, drx, 0, dz,
           0,   0,   0,   0;
}


void calculate_new_pose(Matrix<float,4,4> &T, Matrix<float,4,4> DELTA){
  Matrix<float,4,4> dT;
  dT = DELTA * T;
  T = T + dT;
}

bool is_pose_in_workspace(std::vector<float> joints){
  
}

void close_gripper(){
  digitalWrite(GRIPPER_O, LOW);
  delay(500);
  digitalWrite(GRIPPER_S, HIGH);
}


void open_gripper(){
  digitalWrite(GRIPPER_S, LOW);
  delay(500);
  digitalWrite(GRIPPER_O, HIGH);  
}


void notify(){
  L1_bt = PS4.L1();
  R1_bt = PS4.R1();
  L2_bt = PS4.L2();
  R2_bt = PS4.R2();
  left_stick_x = PS4.LStickX();
  left_stick_y = PS4.LStickY();
  right_stick_y = PS4.RStickY();
}


void onConnect(){
  Serial.println("Connected!.");
}


void onDisConnect(){
  Serial.println("Disconnected!.");    
}


void setup(){
  Serial.begin(115200);
  PS4.attach(notify);
  PS4.attachOnConnect(onConnect);
  PS4.attachOnDisconnect(onDisConnect);
  PS4.begin();
  
  pinMode(GRIPPER_O, OUTPUT);
  pinMode(GRIPPER_S, OUTPUT);
  ledcSetup(0, PWM_FREQUENCY, PWM_RESOLUTION);
  ledcSetup(1, PWM_FREQUENCY, PWM_RESOLUTION);
  ledcSetup(2, PWM_FREQUENCY, PWM_RESOLUTION);
  ledcSetup(3, PWM_FREQUENCY, PWM_RESOLUTION);
  ledcSetup(4, PWM_FREQUENCY, PWM_RESOLUTION);
  //ledcSetup(5, PWM_FREQUENCY, PWM_RESOLUTION);
  ledcAttachPin(13, 0);
  ledcAttachPin(12, 1);
  ledcAttachPin(14, 2);
  ledcAttachPin(27, 3);
  ledcAttachPin(26, 4);
  //ledcAttachPin(25, 5); 
  delay(100);

  o = {M_PI/3, M_PI/3, M_PI/2, M_PI/6, 0};
  forward_kine(pose, o);
  inverse_kine(pose);

  for(int i=0; i<5; i++){
    Serial.printf("%d ", (unsigned int) MG996R_RAD_TO_PWM(o[i]));
    ledcWrite(i, (unsigned int) MG996R_RAD_TO_PWM( o[i] ));
    }
}


void loop() { 
  if (PS4.isConnected()) {
  //  if (PS4.L1()) o[0] = o[0] - 0,1;  // Decremento base
  //  if (PS4.R1()) o[0] = o[0] + 0,1;  // Incremento base
  //  if (PS4.L2()) o[4] = o[4] - 0,1;  // Decremento pulso
  //  if (PS4.R2()) o[4] = o[4] + 0,1;  // Incremento pulso
  //  forward_kine(pose, o);
  
    /* joysticks linear moviments */
    forward_kine(pose, o);
    calculate_diff_operator(DELTA);
    calculate_new_pose(pose, DELTA);
    o = inverse_kine(pose);        
  }
  
  Serial.printf("\n\no1 = %f \n", o[0]);
  Serial.printf("\no2 = %f \n", o[1]);
  Serial.printf("\no3 = %f \n", o[2]);
  Serial.printf("\no4 = %f \n", o[3]);
  Serial.printf("\no5 = %f \n", o[4]);  
  for(int i=0; i<4; i++){
    ledcWrite(i,(unsigned int) MG996R_RAD_TO_PWM( o[i] ));
  }
  
  //analogRead(SENSOR)>600? close_gripper() : open_gripper();
  
  delay(50);
}
