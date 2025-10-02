#include <Servo.h>
#include <math.h>

// Arduino pin assignment
#define PIN_LED  9
#define PIN_TRIG 12
#define PIN_ECHO 13

#define PIN_SERVO 10

// configurable parameters
#define SND_VEL 346.0
#define INTERVAL 25
#define PULSE_DURATION 10
#define _DIST_MIN 100

#define TIMEOUT ((INTERVAL / 2) * 1000.0)
#define SCALE (0.001 * 0.5 * SND_VEL)

#define _MOVING_TIME 1000
#define _WAIT_FOR 2000
#define _MULTIPLIER 10000

// global variables
unsigned long last_sampling_time;
unsigned long last_sensed_time;

Servo myServo;
unsigned long moveStartTime;
int startAngle = 0;
int stopAngle  = 100;
boolean isOpen = false;

float sigmoid(float value, int gap)
{
  float x = (value / gap) * 20 - 10;
  return (1 / (1 + exp(-x))) * _MULTIPLIER;
}

float easeInOutCubic(float value, int gap)
{
  float x = value / gap;
  return (x < 0.5 ? 4 * x * x * x : 1 - pow(-2 * x + 2, 3) / 2) * _MULTIPLIER;
}

void setup() {
  Serial.begin(57600);
  
  // initialize GPIO pins
  pinMode(PIN_LED,OUTPUT);
  pinMode(PIN_TRIG,OUTPUT);
  pinMode(PIN_ECHO,INPUT);
  digitalWrite(PIN_TRIG, LOW);
  
  myServo.attach(PIN_SERVO);
  moveStartTime = millis(); // start moving

  myServo.write(startAngle); // Set position
}

void loop() {
  if (millis() >= last_sampling_time + INTERVAL)
  {
    float dist_raw = USS_measure(PIN_TRIG,PIN_ECHO);

    if (dist_raw > 20 && dist_raw < _DIST_MIN)
    {
      digitalWrite(PIN_LED, 0);
      last_sensed_time = millis();

      if (!isOpen)
      {
        startAngle = 0;
        stopAngle = 120;
        moveStartTime = millis(); // reset start time for next movement
        isOpen = true;
      }
    }
    else
    {
      digitalWrite(PIN_LED, 1);

      if (isOpen && millis() >= last_sensed_time + _WAIT_FOR)
      {
        startAngle = 120;
        stopAngle = 0;
        moveStartTime = millis(); // reset start time for next movement
        isOpen = false;
      }
    }

    last_sampling_time += INTERVAL;
  }
  
  unsigned long progress = millis() - moveStartTime;

  if (progress <= _MOVING_TIME) {
    // change below to use another function
    long angle = map(sigmoid(progress, _MOVING_TIME), 0, _MULTIPLIER,startAngle, stopAngle);
    myServo.write(angle);
    Serial.println(angle);
  }
}

float USS_measure(int TRIG, int ECHO)
{
  digitalWrite(TRIG, HIGH);
  delayMicroseconds(PULSE_DURATION);
  digitalWrite(TRIG, LOW);
  
  return pulseIn(ECHO, HIGH, TIMEOUT) * SCALE; // unit: mm
}
