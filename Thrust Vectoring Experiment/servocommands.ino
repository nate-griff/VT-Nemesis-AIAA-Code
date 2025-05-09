#include <Servo.h>
Servo tiltservo;
#define tiltPin 3 //servo pin number



int angle = 0;
void setup() {
  // put your setup code here, to run once:
  tiltservo.attach(tiltPin);

}

void loop() {
  // put your main code here, to run repeatedly:
  //reset servos to default position (90 degrees)
  tiltservo.write(0);


  // 45 Degree non instantaneous rotation 
  for (int angle = 90; angle <= 93; angle += 1) {
   
    tiltservo.write(angle);
    //Delay is in mili seconds, 1000 =1 second
    
    delay(1);
  }
delay (1000);

 for (int angle = 93; angle <= 96; angle += 1) {
   
    tiltservo.write(angle);
    //Delay is in mili seconds, 1000 =1 second
    
    delay(1);
  }
delay (1000);

 for (int angle = 96; angle <= 99; angle += 1) {
   
    tiltservo.write(angle);
    //Delay is in mili seconds, 1000 =1 second
    
    delay(1);
  }
delay (1000);

 for (int angle = 99; angle <= 102; angle += 1) {
   
    tiltservo.write(angle);
    //Delay is in mili seconds, 1000 =1 second
    
    delay(1);
  }
delay (1000);

 for (int angle = 102; angle <= 105; angle += 1) {
   
    tiltservo.write(angle);
    //Delay is in mili seconds, 1000 =1 second
    
    delay(1);
  }
delay (1000);

  //135 back to 90
  for (int angle = 105; angle >= 90; angle -= 1) {
    
    tiltservo.write(angle);
    delay(1);
  }
  delay(1000);
 for (int angle = 90; angle >= 87; angle -= 1) {
  
    tiltservo.write(angle);
    delay(1);
  }
delay(1000);

 for (int angle = 87; angle >= 84; angle -= 1) {
  
    tiltservo.write(angle);
    delay(1);
  }
delay(1000);

 for (int angle = 84; angle >= 81; angle -= 1) {
  
    tiltservo.write(angle);
    delay(1);
  }
delay(1000);

 for (int angle = 81; angle >= 78; angle -= 1) {
  
    tiltservo.write(angle);
    delay(1);
  }
delay(1000);

 for (int angle = 78; angle >= 75; angle -= 1) {
  
    tiltservo.write(angle);
    delay(1);
  }
delay(1000);

  for (int angle = 75; angle <= 90; angle += 1) {

    tiltservo.write(angle);
    //Delay is in mili seconds, 1000 =1 second
    
    delay(1);
  }
  delay(10000);
 

 
}
