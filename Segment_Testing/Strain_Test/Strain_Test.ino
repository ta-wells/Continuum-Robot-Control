#define analogPin      A4
#define chargePin      7
#define dischargePin   8
#define resistorValue  6000000.0F   // actual resistor value in ohms
 
unsigned long startTime;
unsigned long elapsedTime;
float microFarads;
float picoFarads;
 
const int sampleCount = 100; // rolling window size
float readings[100];         // circular buffer
int indexPos = 0;            // current index in buffer
float total = 0;             // sum of current window
bool filled = false;         // tracks if buffer is full yet
 
void setup() {
  pinMode(chargePin, OUTPUT);
  digitalWrite(chargePin, LOW);
  Serial.begin(9600);
 
  // Initialize readings to 0
  for (int i = 0; i < sampleCount; i++) {
    readings[i] = 0;
  }
}
 
void loop() {
  // Charge the capacitor
  digitalWrite(chargePin, HIGH);
  startTime = micros();
  while (analogRead(analogPin) < 648) { }
  elapsedTime = micros() - startTime;
 
  // Convert to capacitance
  microFarads = ((float)elapsedTime / resistorValue);
  picoFarads = microFarads * 1000000.0;
  picoFarads = picoFarads - 54; // account for noise and offset
 // change the -54 value to calibrate it for each sensor
 
  microFarads = ((float)elapsedTime / resistorValue);
  Serial.print(elapsedTime);       // print the value to serial port
  Serial.print(" microS    ");
 
  // Update rolling sum: subtract old value, add new value
  total -= readings[indexPos];
  readings[indexPos] = picoFarads;
  total += picoFarads;
 
  // Advance circular buffer index
  indexPos++;
  if (indexPos >= sampleCount) {
    indexPos = 0;
    filled = true; // after first wrap, buffer is full
  }
 
  // Compute average
  float avgPico;
  if (filled) {
    avgPico = total / sampleCount;
  } else {
    avgPico = total / indexPos; // average over partial fill
  }
 
  // Print smoothed reading
 
  Serial.println(String(avgPico, 6) + " pF");
 
  // Discharge capacitor
  digitalWrite(chargePin, LOW);
  pinMode(dischargePin, OUTPUT);
  digitalWrite(dischargePin, LOW);
  while (analogRead(analogPin) > 0) { }
  pinMode(dischargePin, INPUT);
}