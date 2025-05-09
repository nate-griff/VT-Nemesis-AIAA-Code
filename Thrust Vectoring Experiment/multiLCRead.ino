#include <HX711_ADC.h>

/*
 * Code for Reading multiple loadcells
 */

#include "HX711.h"

// load cell 0 is the thrust load cell
// load cells 1 & 2 are servo load cells


#define VCC3 4 //Set Pin 52 to be additional 5V source
#define GND3 11 //Set Pin 51 to be additional GND source

#define CQ0 8   // Load Cell 0 SCK0
#define DT0 7   // Load Cell 0 DOUT0

#define CQ1 9   // Load Cell 1 SCK1
#define DT1 6   // Load Cell 1 DOUT1

#define CQ2 10  // Load Cell 2 SCK2
#define DT2 5   // Load Cell 2 DOUT2

HX711 tr;       // this section tares the load cell
HX711 tq1;
HX711 tq2;


long clF0 = 14.12;     // Calibration value: Thrust
long clF1 = -4253;  // Calibration value: servoLC1
long clF2 = -1870;      // Calibration value: servoLC2

float trVal;
float tq1P;
float tq2P;
float tqVal;
float tqP;

void setup() {
  Serial.print("Start of Setup");
  Serial.begin(115200);

  pinMode(VCC3, OUTPUT); // Finish process of 52 5V 
  digitalWrite(VCC3, HIGH); // Finish process of 52 5V
  
  pinMode(GND3, OUTPUT); // Finish process of GND 
  digitalWrite(GND3, LOW); // Finish process of GND

  tr.begin(DT0, CQ0); //Thrust Load Cell
  tq1.begin(DT1, CQ1); //Servo Load Cells
  tq2.begin(DT2, CQ2); //Servo Load Cells

  tr.set_scale(clF0);
  tq1.set_scale(clF1);
  tq2.set_scale(clF2);  

  tr.tare();
  tq1.tare();
  tq2.tare();
  Serial.print("End Setup.");    

}

void loop() {
    trVal = tr.get_units(); //Force (N)
    tq1P = tq1.get_units(); //Torque (N*M)
    tq2P = tq2.get_units(); //Torque (N*m)

    tqP = (tq1P + tq2P)/2 ; // Torque
    // tqVal = tqP; // Torque in N.m

    if (trVal == 0 && tqP == 0 && tq1P == 0 && tq2P == 0){
      return ;
    }

    // Get the current time in milliseconds
    unsigned long t = millis();

    // Print values in CSV format with a timestamp
    Serial.print(trVal);
    Serial.print(",");
    Serial.print(tqP);
    Serial.print(",");
    Serial.print(tq1P);
    Serial.print(",");
    Serial.println(tq2P); // Use println to end the line
    //delay(300);
}
