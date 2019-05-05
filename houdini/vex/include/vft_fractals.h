/*
-----------------------------------------------------------------------------
This source file has been developed within the scope of the
Technical Director course at Filmakademie Baden-Wuerttemberg.
http://technicaldirector.de
    
Written by Juraj Tomori.
Copyright (c) 2019 Animationsinstitut of Filmakademie Baden-Wuerttemberg
-----------------------------------------------------------------------------
*/

#ifndef VFT_FRACTALS
#define VFT_FRACTALS

struct computeFogParms
{
    vector pos = {0.0, 0.0, 0.0};
    int max_iter = 20;
    float max_dist = 2;
    float power = 8;
    float orbit_volume_scale = 1.0;
    vector orbit_volume_offset = {0.0, 0.0, 0.0};
}

function float computeFog(computeFogParms parms; vector colors)
{
    vector z = parms.pos;
    int i = 0;

    float dr = 1.0;
    float r = 0.0;
    float orbit_volume = 10000000000000000000;
        
    for (; i < parms.max_iter; i++) {
        // convert to polar coordinates
        float theta = acos(z.z/r);
        float phi = atan(z.y, z.x);
        dr =  pow(r, parms.power-1.0) * parms.power * dr + 1.0;
        
        // scale and rotate the point
        float zr = pow(r, parms.power);
        theta = theta * parms.power;
        phi = phi * parms.power;
        
        // convert back to cartesian coordinates
        z = zr * set(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
        z += parms.pos;

        // check for exit
        r = length(z);
        orbit_volume = min(orbit_volume, volumesample(1, 0, (z + parms.orbit_volume_offset) * parms.orbit_volume_scale));

        if (r > parms.max_dist)
            break;
    }

    colors.x = orbit_volume;
    return i == parms.max_iter ? 1.0 : 0.0;
}

#endif
