#include <Arduino.h>
//Biblioteker til BLE
#include <BLEDevice.h>
#include <BLEServer.h>
#include <BLEUtils.h>
#include <BLE2902.h>
#include <math.h> 
#include <stdio.h> 
#define LOWBYTE 0
#define HIGHBYTE 1
bool current_byte = LOWBYTE;
const int16_t arr_size = 3400; //Der læses 3400 samples ad gangen
const int16_t seg_size = arr_size / 2; //Hvert segment er 1700 samples langt
const int SKG_window_size = 300; //SKG-segmentet der bruges til alignment er 300 samples langt
const int EKG_window_size = 600; //EKG-segmentet der bruges til alignment er 600 samples langt
int16_t y_r = -32767; //Maksimumværdien for R-takken sættes til -32767 til at starte med
int16_t x_r = 0; //x-værdien til R-takkens maksimum
int16_t y_t = -32767; //Maksimumværdien for T-bølgen sættes til -32767 til at starte med
int16_t x_t = 0; //x-værdien til T-bølgens maksimum
int16_t y_ac_max = -32767; //Maksimumværdien for AC-amplituden sættes til -32767 til at starte med
int16_t x_ac_max = 0; //x-værdien til AC-amplitudens maksimum
float mean_arr_skg[SKG_window_size] = {}; //Array til gennemsnitligt SKG-signal
float mean_arr_ekg[EKG_window_size] = {}; //Array til gennemsnitligt EKG-signal
int16_t count = 0; //
int16_t ekg_data[seg_size] = {}; //EKG-data
int16_t skg_data[seg_size] = {}; //SKG-data
int16_t value = 0;
int loopCheck = 0;
const int loops = 25; //Antallet af signalsegmenter AC-amplituden skal bestemmes udfra
int16_t aligned_SKG[SKG_window_size][loops]; //Matrice med 300 rækker og 25 kolonner til de alignede SKG-segmenter
int16_t aligned_EKG[EKG_window_size][loops]; //Matrice med 600 rækker og 25 kolonner til de alignede EKG-segmenter
int current_loop = 0; 
int32_t final_mean_skg[SKG_window_size]; //Endeligt gennemsnit af SKG-signal
int32_t final_mean_ekg[EKG_window_size]; //Endeligt gennemsnit af SKG-signal
uint32_t sendToBLESKG[SKG_window_size]; //Det SKG-signal der sendes til GUI'en med BLE
uint32_t sendToBLEEKG[EKG_window_size]; //Det EKG-signal der sendes til GUI'en med BLE
int16_t final_max_ac; //AC-amplitudens maksimum
int16_t final_min_ac = 32767; //AC-amplitudens minimum
uint16_t final_ac_amp; //AC-amplituden
uint16_t y = 0;
uint32_t c = 0;
uint32_t defaultData = 10000;

// Opsætning af BLE

/* Definer UUID for Custom Service */
#define serviceID BLEUUID((uint16_t)0x1700)

/* Definer custom characteristic til EKG-data */
BLECharacteristic customCharacteristic(
  BLEUUID((uint16_t)0x1A00), // <-- UUID til characteristic for EKG
  BLECharacteristic::PROPERTY_READ
);

/* Definer custom characteristic til SKG-data */
BLECharacteristic customCharacteristic2(
  BLEUUID((uint16_t)0x2A00), // <-- nyt (SKG UUID)
  BLECharacteristic::PROPERTY_READ 
);

/* Definer custom characteristic til AC-data */
BLECharacteristic customCharacteristic3(
  BLEUUID((uint16_t)0x3A00), // <-- nyt (AC UUID)
  BLECharacteristic::PROPERTY_READ 
);

/* Håndterer server callbacks */
bool deviceConnected = false;
class ServerCallbacks: public BLEServerCallbacks {
    void onConnect(BLEServer* MyServer) {
      deviceConnected = true;
    };

    void onDisconnect(BLEServer* MyServer) {
      deviceConnected = false;
    }
};

void recieve_data() //Læser data fra den serielle port - Læser 3400 samples ad gangen
{
  while (count < arr_size)
  {
    if (Serial.available() > 0)
    {
      if (current_byte == LOWBYTE)
      {
        value = (int16_t)(Serial.read());
        current_byte = HIGHBYTE;
      }
      else
      {
        value += (int16_t)(Serial.read() << 8); // hent ny highbyte
        current_byte = LOWBYTE;

        if (count < seg_size) //De første 1700 samples gemmes i ekg_data[]
        {
          ekg_data[count] = value;
          count++;
        }

        else //De næste 1700 samples gemmes i skg_data[]
        {
          skg_data[count - seg_size] = value;
          count++;
        }
      }
    }
  }
}


//Til filtrering
#define FILTER_LENGTH 3
// koeficienter til EKG filter - Filterkoefficienter er fundet ved hjælp af Butter funktion i MATLAB
double a[FILTER_LENGTH + 1] = {1, -2.98869028151567, 2.97744442748529, -0.988753966159255};
double b[FILTER_LENGTH + 1] = {0.994361084395028, -2.98308325318508, 2.98308325318508, -0.994361084395028};

// array der kommer til at indeholde output af filtrering
double output_ekg[seg_size] = {};

// array der kommer til at indeholde det filtrerede data
int16_t ekg_filt[seg_size] = {};

// fikser afrundingsfejl;
int16_t floor_and_convert(double value)
{
  if (value > 0) // positiv
  {
    return (int16_t)(value + 0.5);
  }
  else // negativ
  {
    return (int16_t)(value - 0.5);
  }
}

// filtrering af EKG
void iir_ekg_filter(void) //3. ordens Butterworth højpasfilter med en knækfrekvens på 0,9 Hz.
{
  for (int i = 0; i < seg_size; i++)
  {
    if (i == 0) //Det første sample sammenlignes ikke med andre samples
    {
      output_ekg[i] = b[0] * ekg_data[i];
    }

    else if (i < FILTER_LENGTH)
    {
      output_ekg[i] = b[0] * ekg_data[i];

      for (int j = 1; j <= i; j++)
      {
        output_ekg[i] = output_ekg[i] + b[j] * ekg_data[i - j] - a[j] * output_ekg[i - j];
      }
    }
    else
    {
      output_ekg[i] = b[0] * ekg_data[i];
      for (int j = 1; j <= FILTER_LENGTH; j++)
      {
        output_ekg[i] = output_ekg[i] + b[j] * ekg_data[i - j] - a[j] * output_ekg[i - j];
      }
    }
  }
}

void HP_filt_ekg() //Filtrerer EKG-signalet med filteret ovenfor
{
  iir_ekg_filter();
  for (int i = 0; i < seg_size; i++)
  {
    ekg_filt[i] = floor_and_convert(output_ekg[i]);
  }
}

void find_r_t() //Finder R-takker, T-bølger og tilhørende x-værdier i EKG-dataet
{
  for (int i = 0; i < (seg_size - 500); i++) //Finder den maksimale værdi indenfor de første 1200 samples af EKG-signalet - Denne svarer til toppen af R-takken
  {
    if (ekg_filt[i] > y_r) //Hvis EKG på plads i er større end den tidligere maksimumværdi for R-takken
    {
      y_r = ekg_filt[i]; //Gemmer den nye maksimumværdi
      x_r = i; //Gemmer x-værdien for den fundne maksimale værdi
    }
  }

  for (int idx = (x_r + 33); idx < (x_r + 400); idx++) //Finder den maksimale værdi fra 33 samples til 400 samples efter R-takken - Denne svarer til toppen af T-bølgen
  {
    if (ekg_filt[idx] > y_t) //Hvis EKG på plads i er større end den tidligere maksimumværdi for T-bølgen
    {
      y_t = ekg_filt[idx]; //Gemmer den nye maksimumværdi
      x_t = idx; //Gemmer x-værdien for den fundne maksimale værdi
    }
  }
}

void find_ac_max() //Efter samme princip som find_r_t() finder denne funktion maksimum og tilhørende x-værdi i SKG-dataet
{ 

  for (int i = x_t; i < (x_t + 160); i++) //Der ledes efter maksimum fra T-bølgen og 160 samples frem
  {
    if (skg_data[i] > y_ac_max)
    {
      y_ac_max = skg_data[i];
      x_ac_max = i;
    }
  }
}

void align() //Aligner SKG- og EKG-signalerne ift. de funde maksimum-værdier
{
  for (int j = 0; j < SKG_window_size; j++) //Laver segmenter på 300 samples, fra 150 samples før maksimum til 150 samples efter maksimum
  {
    aligned_SKG[j][current_loop] = skg_data[x_ac_max - SKG_window_size / 2 + j]; //Segmenterne gemmes i en matrice, hvert segment i en ny kolonne
  }
 for (int j = 0; j < EKG_window_size; j++) 
 {
  aligned_EKG[j][current_loop] = ekg_data[x_r - (EKG_window_size / 2) + j]; 
 }
}

void mean_arrays() //Beregner gennemsnitssignaler ud fra de alignede signaler - Både for EKG og SKG
{
  for (int i = 0; i < loops; i++)
  {
    for (int k = 0; k < SKG_window_size; k++) //Summerer de alignede signaler
    {
      mean_arr_skg[k] = mean_arr_skg[k] + aligned_SKG[k][i];
    }
     for (int k = 0; k < EKG_window_size; k++)
    {
      mean_arr_ekg[k] = mean_arr_ekg[k] + aligned_EKG[k][i];
    }
  }

  for (int k = 0; k < SKG_window_size; k++)
  {
    mean_arr_skg[k] = mean_arr_skg[k] / loops; //Deler de summerede signaler med antallet af loops der er kørt igennem, for at få et gennemsnit
  }
   for (int k = 0; k < EKG_window_size; k++)
  {
    mean_arr_ekg[k] = mean_arr_ekg[k] / loops;
  }

  for (int k = 0; k < SKG_window_size; k++)
  {
    final_mean_skg[k] = mean_arr_skg[k]; //Laver datatypen om for gennemsnitssignalet - fra float til int32
    sendToBLESKG[k] = (final_mean_skg[k]+10000); //Lægger 10000 til før signalet skal sendes med BLE - desuden laves det til uint32

  }
   for (int k = 0; k < EKG_window_size; k++)
  {
    final_mean_ekg[k] = mean_arr_ekg[k];
    sendToBLEEKG[k] = (final_mean_ekg[k]+10000);
  }
}

void find_ac_amp() //Finder den endelige AC-amplitude, ud fra det gennemsnitlige SKG-signal
{

  final_max_ac = final_mean_skg[SKG_window_size / 2]; //På baggrund af måden signalerne er alignet på, vil maks-værdien for AC-amplituden findes på den midterste plads i gennemsnitssignalet

  for (int i = (SKG_window_size / 2); i > (SKG_window_size / 2) - 20; i--) //Fra maksværdien og 20 samples tilbage ledes der efter minimumsværdien
  {
    if (final_mean_skg[i] < final_min_ac)
    {
      final_min_ac = final_mean_skg[i];
    }
  }

  final_ac_amp = final_max_ac - final_min_ac; //Den endelige AC-amplitude beregnes
}

void reset() //Denne funktion bruges til at resette nogle variable når  loopet starter forfra
{
  count = 0;
  y_t = -32767;
  y_r = -32767;
  y_ac_max = -32767;
}

void setup() 
{
  Serial.begin(115200);
  pinMode(02, OUTPUT);

  BLEDevice::init("ESP BLE server enhed"); //Navn til BLE-enheden

  /* Opretter BLE Server */
  BLEServer *MyServer = BLEDevice::createServer();
  MyServer->setCallbacks(new ServerCallbacks()); // funktion der håndterer Server Callbacks

  /* Tilføj service til server */
  BLEService *customService = MyServer->createService(BLEUUID((uint16_t)0x1700)); //  Et tilfældigt ID vælges

  /* Tilføj characteristic til service */
  customService->addCharacteristic(&customCharacteristic);  
  customService->addCharacteristic(&customCharacteristic2); 
  customService->addCharacteristic(&customCharacteristic3); 

  customCharacteristic.setValue(defaultData);  // Sætter characteristic til at være 10000 indtil data kommer
  customCharacteristic2.setValue(defaultData); // Sætter characteristic til at være 10000 indtil data kommer
  customCharacteristic3.setValue(y); //Sætter characteristic til at være nul indtil data kommer
  
  /* Konfigurerer Advertising med Services der skal advertises */
  MyServer->getAdvertising()->addServiceUUID(serviceID);

  customService->start(); // Start service

  MyServer->getAdvertising()->start(); // Start Server/Advertising
}

void loop()
{
  current_loop = 0; 
  if (loopCheck == 0) { //Sørger for at der kun køres databehandling på 25 EKG- og SKG-segmenter

    for (int i = 0; i < loops; i++) //For-loop til de funktioner der skal køre alle 25 gange
    {
      current_loop = i;
      reset();
      recieve_data();
      HP_filt_ekg();
      find_r_t();
      find_ac_max();
      align();
    }
  mean_arrays(); 
  find_ac_amp();
  customCharacteristic3.setValue(final_ac_amp);

   for (int i = 0; i<5; i++) { //Blinker 5 gange med LED på ESP32 for at vise at databehandlingen er færdig
  digitalWrite(02, HIGH);
  delay(100);
  digitalWrite(02, LOW);
  delay(100);
  }

  loopCheck = 1; 
  }
  
if(deviceConnected){ //
    if (c < EKG_window_size-1){
      digitalWrite(02,HIGH);
       delay(500);
        for ( c = 0 ; c < EKG_window_size ; c++){
        customCharacteristic.setValue(sendToBLEEKG[c]); //Sender gennemsnitligt EKG-signal
        delay(60); 
        }
        for ( c = 0 ; c < SKG_window_size ; c++){
        customCharacteristic2.setValue(sendToBLESKG[c]); //Sender gennemsnitligt SKG-signal
        delay(80);
        digitalWrite(02,LOW);
        }
    }
    else {
      Serial.println("No device connected");
      delay(1000);
    } 
}
}
