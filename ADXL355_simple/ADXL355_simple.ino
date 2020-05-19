
/**
******************************************************************************
    @file     main.c
    @brief    Project main source file
    @version  V0.1
    @author   ADI
    @date     May 2016
   @par Revision History:
   - V0.1, May 2016: initial version.

*******************************************************************************
  Copyright 2016(c) Analog Devices, Inc.

  All rights reserved.

  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:
   - Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   - Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in
     the documentation and/or other materials provided with the
     distribution.
   - Neither the name of Analog Devices, Inc. nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.
   - The use of this software may or may not infringe the patent rights
     of one or more patent holders.  This license does not release you
     from the requirement that you obtain separate licenses from these
     patent holders to use this software.
   - Use of the software either in source or binary form, must be run
     on or directly connected to an Analog Devices Inc. component.

  THIS SOFTWARE IS PROVIDED BY ANALOG DEVICES "AS IS" AND ANY EXPRESS OR IMPLIED
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, NON-INFRINGEMENT, MERCHANTABILITY
  AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL ANALOG DEVICES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  INTELLECTUAL PROPERTY RIGHTS, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*******************************************************************************
**/

/***************************** Include Files **********************************/
#include <stdio.h>
#include <stdlib.h>
#include <SPI.h>

#include "ADXL355.h"

/*******************************************************************************
**************************** Internal types ************************************
********************************************************************************/

/* Write data mode */
typedef enum {
  SPI_WRITE_ONE_REG = 1,         /* Write 1 ACC register */
  SPI_WRITE_TWO_REG,             /* Write 2 ACC register */
} enWriteData;

typedef enum {
  SPI_READ_ONE_REG = 1,            /* Read one ACC register */
  SPI_READ_TWO_REG,                /* Read two ACC registers */
  SPI_READ_THREE_REG,              /* Read X,Y,Z ACC registers */
} enRegsNum;

/*******************************************************************************
**************************** Internal definitions ******************************
********************************************************************************/
#define int8_t    char
#define uint8_t   unsigned char
#define int32_t   int
#define uint32_t  unsigned int

/* Accelerometer write command */
#define ADXL355_WRITE         0x0

/* Accelerometer read command */
#define ADXL355_READ          0x1

#define CSACC_PIN           5   /* CSADXL355 - output */
#define INT1ACC_PIN         0   /* INT1 - input */
#define INT2ACC_PIN         0   /* INT2 - input */
#define DATARDYACC_PIN      16   /* DATA RDY - input */

#define ADXL_RANGE 2           // Define the accelerometer range to +/- 2g

const int ledPin = 2;             // ledPin refers to ESP32 GPIO 2

int32_t volatile i32SensorX;

int32_t volatile i32SensorX_init;

uint32_t volatile ui32SensorX;

float adxl355Scale;

/**************************** Function Definitions ****************************/

/**
   @brief SPI initialization

   @return none

**/
void SPI_Init(void)
{
  SPI.begin();
  //  SPI.setClockDivider(SPI_CLOCK_DIV128);
  SPI.setClockDivider(SPI_CLOCK_DIV8);

  //  SPI.begin();
  digitalWrite(CSACC_PIN, HIGH);         /* Deselect accelerometer */
}

/**
   @brief Writes a data, a command or a register to the LCD or to ACC via SPI.

   @param ui8address - ACC register address
   @param ui8Data - value to be written in 1 register write
   @param ui8Data2 - 2nd value to be written in 2 register write
   @enMode enWriteData - write mode

   @return none

**/
void SPI_Write(uint8_t ui8address, uint8_t ui8Data, uint8_t ui8Data2, enWriteData enMode)
{
  uint8_t ui8writeAddress;
  ui8writeAddress = ((ui8address << 1) | ADXL355_WRITE);

  if (enMode == SPI_WRITE_ONE_REG) {
    digitalWrite(CSACC_PIN, LOW);         /* Select accelerometer */

    SPI.transfer(ui8writeAddress);     /* Send register address */
    SPI.transfer(ui8Data);             /* Send value to be written */

    digitalWrite(CSACC_PIN, HIGH);         /* Deselect accelerometer */
  }

  if (enMode == SPI_WRITE_TWO_REG) {
    digitalWrite(CSACC_PIN, LOW);         /* Select accelerometer */

    SPI.transfer(ui8writeAddress);     /* Send register address */
    SPI.transfer(ui8Data);             /* Send 1st value to be written */
    SPI.transfer(ui8Data2);             /* Send 2nd value to be written */

    digitalWrite(CSACC_PIN, HIGH);         /* Deselect accelerometer */
  }
}

/**
   @brief Reads a specified register or two registers address in the accelerometer via SPI.

   @param ui8address - register address
   @param enRegs - register number

   @return reading result

**/
uint32_t SPI_Read(uint8_t ui8address, enRegsNum enRegs)
{
  uint32_t ui32Result = 0;

  uint32_t ui32valueL = 0; // LOW BYTE
  uint32_t ui32valueM = 0; // MEDIUM BYTE
  uint32_t ui32valueH = 0; // HIGH BYTE

  uint8_t ui8writeAddress;
  ui8writeAddress = ((ui8address << 1) | ADXL355_READ);

  digitalWrite(CSACC_PIN, LOW);         /* Select accelerometer */

  SPI.transfer(ui8writeAddress);       /* Send register address */

  if (enRegs == SPI_READ_THREE_REG) {          /* Only used for X,Y,Z axis data registers*/
    ui32valueH = SPI.transfer(0x00);             /* Read the register value */
    ui32valueM = SPI.transfer(0x00);
    ui32valueL = SPI.transfer(0x00);

    ui32Result = ((ui32valueH << 16) | (ui32valueM << 8) | ui32valueL); /* Set read result*/
  }

  digitalWrite(CSACC_PIN, HIGH);         /* Deselect accelerometer */

  return ui32Result;
}

/**
   @brief Initialization the accelerometer sensor

   @return none

**/
void ADXL355_Init(void)
{
  pinMode(CSACC_PIN, OUTPUT);                  /* Set CSACC pin as output */
  pinMode(DATARDYACC_PIN, INPUT);              /* Set DRDY pin as input */
}

/**
   @brief Turns on accelerometer measurement mode.

   @return none

**/
void ADXL355_Start_Sensor(void)
{
  uint8_t ui8temp;

  ui8temp = (uint8_t)SPI_Read(POWER_CTL, SPI_READ_ONE_REG);       /* Read POWER_CTL register, before modifying it */

  ui8temp = ui8temp & 0xFE;                                       /* Set measurement bit in POWER_CTL register */
  SPI_Write(POWER_CTL, ui8temp, 0x00, SPI_WRITE_ONE_REG);         /* Write the new value to POWER_CTL register */

}

/**
   @brief Puts the accelerometer into standby mode.

   @return none

**/
void ADXL355_Stop_Sensor(void)
{
  uint8_t ui8temp;

  ui8temp = (uint8_t)SPI_Read(POWER_CTL, SPI_READ_ONE_REG);        /*Read POWER_CTL register, before modifying it */
  ui8temp = ui8temp | 0x01;                                      /* Clear measurement bit in POWER_CTL register */
  SPI_Write(POWER_CTL, ui8temp, 0x00, SPI_WRITE_ONE_REG);                 /* Write the new value to POWER_CTL register */
}

/**
   @brief Reads the accelerometer data.

   @return none

**/
void ADXL355_Data_Scan(void)
{
  ui32SensorX = SPI_Read(XDATA3, SPI_READ_THREE_REG);

  i32SensorX = ADXL355_Acceleration_Data_Conversion(ui32SensorX);

}

/**
   @brief Convert the two's complement data in X,Y,Z registers to signed integers

   @param ui32SensorData - raw data from register

   @return int32_t - signed integer data

**/
int32_t ADXL355_Acceleration_Data_Conversion (uint32_t ui32SensorData)
{
  int32_t volatile i32Conversion = 0;

  ui32SensorData = (ui32SensorData  >> 4);
  ui32SensorData = (ui32SensorData & 0x000FFFFF);

  if ((ui32SensorData & 0x00080000)  == 0x00080000) {
    i32Conversion = (ui32SensorData | 0xFFF00000);
  } else {
    i32Conversion = ui32SensorData;
  }

  return i32Conversion;
}

void ADXL355_Data_Init()
{
  ADXL355_Data_Scan();
  i32SensorX_init = 0;
}

/**
   @brief The main application function

   @return the function contains infinite loop and never returns/

**/
void setup()
{

  /* Initialize UART */
  Serial.begin(230400);
  //Serial.println("***** ADXL355 Simple Test *****");

  pinMode(ledPin, OUTPUT);        // initialize digital pin ledPin as an output.
  for (int i = 0; i < 10; i++)
  {
    digitalWrite(ledPin, HIGH);
    delay(50);
    digitalWrite(ledPin, LOW);
    delay(50);
  }
  /* Initialize accelerometer */
  ADXL355_Init();

  /* Initialize SPI */
  SPI_Init();
  delay(1000);

#if ADXL_RANGE == 2
  SPI_Write(RANGE, 0x81, 0x00, SPI_WRITE_ONE_REG);          /* Set sensor range within RANGE register */
  adxl355Scale = 256.0f;
#endif
}

void loop()
{
  byte x[4];

  if ( digitalRead(DATARDYACC_PIN) == HIGH )
  {
    ADXL355_Data_Scan();

#if (0)


#else

    x[0] = i32SensorX & 0xFF; // 0x78 SEND LOW BYTE
    x[1] = (i32SensorX >> 8) & 0xFF; // 0x56
    x[2] = (i32SensorX >> 16) & 0xFF; // 0x34
    x[3] = (i32SensorX >> 24) & 0xFF; // 0x12 SEND HIGH BYTE


    Serial.write(x, 4);

#endif

  }

  serialEvent();
}

void serialEvent()
{
  volatile uint32_t ui32test = 0x00;                /* Read the ID register */
  volatile uint32_t ui32test2;                  /* Read the Manufacturer register */
  volatile uint32_t ui32test3;                  /* Read the Part ID register */
  while (Serial.available())
  {
    int inChar = Serial.read();       // get the new byte
    switch (inChar)
    {
      case 'S':
        digitalWrite(ledPin, HIGH);   // turn the LED on (HIGH is the voltage level)
        do {
          /* Start accelerometer measurement mode */
          ADXL355_Start_Sensor();

          ui32test = SPI_Read(DEVID_AD, SPI_READ_ONE_REG);                  /* Read the ID register */
          ui32test2 = SPI_Read(DEVID_MST, SPI_READ_ONE_REG);                  /* Read the Manufacturer register */
          ui32test3 = SPI_Read(PARTID, SPI_READ_ONE_REG);                  /* Read the Part ID register */

          delay(1000);
        } while ( ui32test == 0x00 );

        ADXL355_Data_Init();

        break;
      case 'E':
        digitalWrite(ledPin, LOW);    // turn the LED off by making the voltage LOW
        ADXL355_Stop_Sensor();
        break;
      default:
        break;
    }
  }
}

/* End Of File *///
