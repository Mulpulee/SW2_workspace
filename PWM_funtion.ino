#define PIN_LED 7

int period = 100;
int onTime = 0;
int brightness = 0;
int timer = 0;
bool goingUp = true;

void setup()
{
  Serial.begin(115200);
  pinMode(PIN_LED, OUTPUT);

  set_period(100);
}

void loop()
{
    set_duty(brightness);
    digitalWrite(PIN_LED, 1);
    delayMicroseconds(onTime);
    digitalWrite(PIN_LED, 0);
    delayMicroseconds(period - onTime);

    timer += period;
    if (timer >= 5000)
    {
      if (goingUp) brightness += timer / 5000;
      else brightness -= timer / 5000;
      timer = 0;
      Serial.print("period : "); Serial.print(period);
      Serial.print(" / brightness : "); Serial.println(brightness);
    }

    if (brightness <= 0) goingUp = true;
    if (brightness >= 100) goingUp = false;
}

void set_period(int p)
{
  period = p;
}

void set_duty(int d)
{
  onTime = period / 100 * d;
}
