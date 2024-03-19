// GLOBAL VARIABLES
#ifndef MYGLOBALS_HH
#define MYGLOBALS_HH

// Set the number of sensors in the X-Y plane
#define NX_SENSORS 1
#define NY_SENSORS 2
#define N_TOTAL_SENSORS (NX_SENSORS*NY_SENSORS)

// Define the number of Plastic Scintillators
// The veto at the bottom is not considered
#define N_PL_SCINT_NO_VETO 1

// Set to 1 to enable Scintillation processes
// Scintillation is computationally expensive
#define OPTICAL_PROCESSES 0

// Generate the GDML file to export the geometry into a CAD
#define GENERATE_GDML 1

// Define the thickness of the Thin and the Thick detectors
// Thin: Delta E Detector
// Thick: E Detector
#define TK_THIN 100*um
#define TK_THICK 1*nm
#define TK_PLASTIC_VETO 1*cm

// Select the geometry
// Set to 1 to view the new geometry
// Set to 0 to view the geometry from the Master thesis
//      Here, the geometry is implemented with Ametek PIPS
#define NEW_GEOMETRY 1

#define CZT_DETECTOR 1


#endif