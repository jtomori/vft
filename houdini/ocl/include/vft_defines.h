/*
-----------------------------------------------------------------------------
This source file has been developed within the scope of the
Technical Director course at Filmakademie Baden-Wuerttemberg.
http://technicaldirector.de
    
Written by Juraj Tomori.
Copyright (c) 2019 Animationsinstitut of Filmakademie Baden-Wuerttemberg
-----------------------------------------------------------------------------
*/

#ifndef VFT_DEFINES
#define VFT_DEFINES

// constants

#define NULL                    0

#define LARGE_NUMBER            10e10
#define ORBITS_OFFSET           0.006f // is used to hide colored seams in sdf subtracting

#define ORBITS_ARRAY_LENGTH     9
#define ENABLE_DELTA_DE         0

// functions

#define SIN(x) native_sin(x)
#define COS(x) native_cos(x)
#define TAN(x) native_tan(x)
#define POWR(x, y) native_powr((x), (y))
#define SQRT(x) native_sqrt(x)
#define LOG(x) native_log(x)
#define EXP(x) native_exp(x)
#define DIV(x, y) native_divide((x), (y))
#define NORMALIZE(x) fast_normalize(x)
#define LENGTH(x) fast_length(x)

#endif
