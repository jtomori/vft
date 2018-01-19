#ifndef _VFT_FRACTALS
#define _VFT_FRACTALS



////////////// primitives


// sphere: position, radius, center
static float sphere( float3 P, float rad, float3 center )
{
    float dist = length(P - center) - rad;
    return dist;
}

// box: position, size
static float box( float3 P, float3 b )
{
  float3 d = fabs(P) - b;
  return min( max( d.x, max(d.y, d.z) ), 0.0f) + length( max(d, 0.0f) );
}

// round box: position, size, roundness
static float roundBox( float3 P, float3 b, float r )
{
  return length( max( fabs(P) - b, 0.0f ) )-r;
}

// torus: position, size (radius, thickness)
static float torus( float3 P, float2 t )
{
  float2 q = (float2)(length(P.xz)-t.x,P.y);
  return length(q)-t.y;
}

// infinite cone
static float cone( float3 P, float2 c )
{
    c = normalize(c);
    float q = length(P.xy);
    return dot( c, (float2)(q, P.z) );
}


////////////// fractals


static float mandelbulb( float3 P, float Power, float size )
{
    P /= size;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    
    const int Iterations = 60; // increase to remove banding
    const int Bailout = 4;
    //float Power = power; // animatable
    
    for (int i = 0; i < Iterations ; i++)
    {
            r = length(z);
            if (r > Bailout) break;
            
            // convert to polar coordinates
            float theta = acos(z.z/r);
            float phi = atan2(z.y, z.x);

            dr =  pow(r, Power-1) * Power * dr + 1;
            
            // scale and rotate the point
            float zr = pow(r, Power);
            theta *= Power;
            phi *= Power;
            
            // convert back to cartesian coordinates
            z = zr * (float3)(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
            z += P;
    }
    float out = 0.5 * log(r) * r/dr;
    return out * size;
}

static float mandelbox( float3 P, float scale, float size)
{
    P /= size;

    const int Iterations = 30;
    const int Bailout = 60;
    float3 P_orig = P;
    float DEfactor = scale;

    for (int i = 0; i < Iterations ; i++)
    {
        float fixedRadius = 1.0;
        float fR2 = fixedRadius * fixedRadius;
        float minRadius = 0.5;
        float mR2 = minRadius * minRadius;

        if (P.x > 1.0) P.x = 2.0 - P.x;
        else if (P.x < -1.0) P.x = -2.0 - P.x;

        if (P.y > 1.0) P.y = 2.0 - P.y;
        else if (P.y < -1.0) P.y = -2.0 - P.y;

        if (P.z > 1.0) P.z = 2.0 - P.z;
        else if (P.z < -1.0) P.z = -2.0 - P.z;

        float r2 = P.x*P.x + P.y*P.y + P.z*P.z;

        if (r2 < mR2)
        {
           P = P * fR2 / mR2;
           DEfactor = DEfactor * fR2 / mR2;
        }
        else if (r2 < fR2)
        {
           P = P * fR2 / r2;
           DEfactor *= fR2 / r2;
        }

        P = P * scale + P_orig;

        DEfactor *= scale;
    }

    return length(P)/fabs(DEfactor) * size;
}

static void mandelbulbPower2Iter(float3* Z, float* de, const bool julia, const float3* P_in) {
    float distance = length(*Z);

    *de = *de * 2.0f * distance;
    float x2 = Z->x * Z->x;
    float y2 = Z->y * Z->y;
    float z2 = Z->z * Z->z;
    float temp = 1.0 - z2 / (x2 + y2);
    float newx = (x2 - y2) * temp;
    float newy = 2.0 * Z->x * Z->y * temp;
    float newz = -2.0 * Z->z * sqrt(x2 + y2);
    Z->x = newx;
    Z->y = newy;
    Z->z = newz;

    *Z += *P_in;
}

// aux.r_dz -> dr (DE factor or something)
// aux.r -> r (length of Z)
// dr -> de
// r -> distance
// Bailout -> max_distance
// Iterations -> max_iterations

static float mandelbulbPower2(float3 P, float size)
{
    P /= size;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    const int Iterations = 60; // increase to remove banding
    const int Bailout = 6;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;

        dr = dr * 2.0 * r;
        float x2 = z.x * z.x;
        float y2 = z.y * z.y;
        float z2 = z.z * z.z;
        float temp = 1.0 - z2 / (x2 + y2);
        float newx = (x2 - y2) * temp;
        float newy = 2.0 * z.x * z.y * temp;
        float newz = -2.0 * z.z * sqrt(x2 + y2);
        z.x = newx;
        z.y = newy;
        z.z = newz;

        z += P;
    }

    float out = 0.5 * log(r) * r/dr;
    return out * size;
}

static void mengerSpongeIter(float3* Z, float* de, const bool julia, const float3* P_in) {
    float distance = length(*Z);

    Z->x = fabs(Z->x);
    Z->y = fabs(Z->y);
    Z->z = fabs(Z->z);

    if (Z->x - Z->y < 0.0f) Z->xy = Z->yx;
    if (Z->x - Z->z < 0.0f) Z->xz = Z->zx;
    if (Z->y - Z->z < 0.0f) Z->yz = Z->zy;

    *Z *= 3.0f;

    Z->x -= 2.0f;
    Z->y -= 2.0f;
    if (Z->z > 1.0f) Z->z -= 2.0f;

    *de *= 3.0;

    if (julia) *Z += *P_in; // julia mode
}

static float mengerSponge(float3 P, float size)
{
    P /= size;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    const int Iterations = 30;
    const int Bailout = 100;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;

        z.x = fabs(z.x);
        z.y = fabs(z.y);
        z.z = fabs(z.z);
    
        if (z.x - z.y < 0.0f) z.xy = z.yx;
        if (z.x - z.z < 0.0f) z.xz = z.zx;
        if (z.y - z.z < 0.0f) z.yz = z.zy;
    
        z *= 3.0f;
    
        z.x -= 2.0f;
        z.y -= 2.0f;
        if (z.z > 1.0f) z.z -= 2.0f;
    
        dr *= 3.0;
        //z += P; // julia mode

    }

    //float out = 0.5f * log(r) * r/dr;
    float out = r / dr;
    return out * size;
}

// bristorbrot
// shape matches well, but normals are noisy, using log de function helps
static float bristorbrot(float3 P, float size)
{
    P /= size;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    const int Iterations = 25;
    const int Bailout = 100;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;

        dr = dr * 2.0f * r;
        float newx = z.x * z.x - z.y * z.y - z.z * z.z;
        float newy = z.y * (2.0f * z.x - z.z);
        float newz = z.z * (2.0f * z.x + z.y);
        z.x = newx; // change to += for julia mode and comment out "z + =P;"
        z.y = newy;
        z.z = newz;

        z += P;
    }

    float out = 0.5f * log(r) * r/dr;
    //float out = r / dr;
    return out * size;
}

// xenodreambuie
static float xenodreambuie(float3 P, float power, float alpha, float beta, float size)
{
    P /= size;

    alpha = radians(alpha);
    beta = radians(beta);

    float pi = 3.14159265359;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    const int Iterations = 100;
    const int Bailout = 100;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;

        float rp = pow(r, power - 1.0f);
        dr = rp * dr * power + 1.0f;
        rp *= r;

        float th = atan2(z.y, z.x) + beta;
        float ph = acos(z.z / r) + alpha;

        if (fabs(ph) > 0.5f * pi) ph = sign(ph) * pi - ph;

        z.x = rp * cos(th * power) * sin(ph * power);
        z.y = rp * sin(th * power) * sin(ph * power);
        z.z = rp * cos(ph * power);

        z += P;
    }

    float out = 0.5f * log(r) * r/dr;
    //float out = r / dr;
    return out * size;
}

// Coastalbrot
// looks weird, but seems to work fine
static float coastalbrot(float3 P, float size)
{
    P /= size;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    const int Iterations = 30;
    const int Bailout = 10;
    float pi = 3.14159265359;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;
        
        float temp = r;
        temp = pow(temp, 7.7f);
        dr = temp * dr * 7.7f;
        temp *= r;

        z.x = sin(sin(sin(pi / 3.0f + z.x * pi)));
        z.y = sin(sin(sin(pi / 3.0f + z.y * pi)));
        z.z = sin(sin(sin(pi / 3.0f + z.z * pi)));

        z *= temp;

        z += P;
    }

    float out = 0.5f * log(r) * r/dr;
    //float out = r / dr;
    return out * size;
}

// sierpinski3d
static float sierpinski3d(float3 P, float scale, float3 offset, float3 rot, float size)
{
    P /= size;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    const int Iterations = 30;
    const int Bailout = 100;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;
        
        float3 temp = z;


        if (z.x - z.y < 0.0f) z.xy = z.yx;
        if (z.x - z.z < 0.0f) z.xz = z.zx;
        if (z.y - z.z < 0.0f) z.yz = z.zy;
        if (z.x + z.y < 0.0f)
        {
            temp.x = -z.y;
            z.y = -z.x;
            z.x = temp.x;
        }
        if (z.x + z.z < 0.0f)
        {
            temp.x = -z.z;
            z.z = -z.x;
            z.x = temp.x;
        }
        if (z.y + z.z < 0.0f)
        {
            temp.y = -z.z;
            z.z = -z.y;
            z.y = temp.y;
        }


        z = z * scale;
        dr *= scale;

        z -= offset;

        z = mtxPtMult( mtxRotate(rot) , z );

        //z += P; // enable for julia mode
    }

    //float out = 0.5f * log(r) * r/dr;
    float out = r / dr;
    return out * size;
}


#endif