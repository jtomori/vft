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
// aux.r_dz -> dr (DE factor or something)
// aux.r -> r (length of Z)
// dr -> de
// r -> distance
// Bailout -> max_distance
// Iterations -> max_iterations

static void mandelbulbIter(float3* Z, float* de, const bool julia, const float3* P_in, const float power) {
    float distance = length(*Z);

    // convert to polar coordinates
    float theta = acos(Z->z/distance);
    float phi = atan2(Z->y, Z->x);

    *de =  pow(distance, power-1) * power * (*de) + 1;
    
    // scale and rotate the point
    float zr = pow(distance, power);
    theta *= power;
    phi *= power;
    
    // convert back to cartesian coordinates
    float3 newZ = zr * (float3)( sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta) );

    // regular fractal
    //*Z = newZ + *P_in;

    // julia fractal
    *Z = newZ + (float3)(1,1,0);
}

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

static void mandelboxIter(float3* Z, float* de, const bool julia, const float3* P_in, const float scale) {
    float distance = length(*Z);

    float fixedRadius = 1.0;
    float fR2 = fixedRadius * fixedRadius;
    float minRadius = 0.5;
    float mR2 = minRadius * minRadius;

    if (Z->x > 1.0) Z->x = 2.0 - Z->x;
    else if (Z->x < -1.0) Z->x = -2.0 - Z->x;

    if (Z->y > 1.0) Z->y = 2.0 - Z->y;
    else if (Z->y < -1.0) Z->y = -2.0 - Z->y;

    if (Z->z > 1.0) Z->z = 2.0 - Z->z;
    else if (Z->z < -1.0) Z->z = -2.0 - Z->z;

    float r2 = Z->x*Z->x + Z->y*Z->y + Z->z*Z->z;

    if (r2 < mR2)
    {
        *Z = *Z * fR2 / mR2;
        *de = *de * fR2 / mR2;
    }
    else if (r2 < fR2)
    {
        *Z = *Z * fR2 / r2;
        *de *= fR2 / r2;
    }

    *de *= scale;

    // regular mode
    *Z = *Z * scale + *P_in;

    // julia mode
    //*Z = *Z * scale + (float3)(5,0,0);
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
    float3 new;
    new.x = (x2 - y2) * temp;
    new.y = 2.0 * Z->x * Z->y * temp;
    new.z = -2.0 * Z->z * sqrt(x2 + y2);

    // regular
    *Z = new + *P_in;

    // julia
    //*Z = new + (float3)(0,1,0);
}

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

    //if (julia) *Z += *P_in; // julia mode

    // regular
    // works without adding anything
    //*Z += *P_in;

    // julia
    *Z += (float3)(0,1,0);
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

static void bristorbrotIter(float3* Z, float* de, const bool julia, const float3* P_in) {
    float distance = length(*Z);

    *de = *de * 2.0f * distance;
    float3 new;
    new.x = Z->x * Z->x - Z->y * Z->y - Z->z * Z->z;
    new.y = Z->y * (2.0f * Z->x - Z->z);
    new.z = Z->z * (2.0f * Z->x + Z->y);
    
    *Z = new;

    // regular
    *Z += *P_in;

    // julia
    //*Z += (float3)(1,0,0);
}

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

static void xenodreambuieIter(float3* Z, float* de, const bool julia, const float3* P_in, const float power, float alpha, float beta) {
    alpha = radians(alpha);
    beta = radians(beta);

    float distance = length(*Z);

    float rp = pow(distance, power - 1.0f);
    *de = rp * (*de) * power + 1.0f;
    rp *= distance;

    float th = atan2(Z->y, Z->x) + beta;
    float ph = acos(Z->z / distance) + alpha;

    if (fabs(ph) > 0.5f * M_PI_F) ph = sign(ph) * M_PI_F - ph;

    Z->x = rp * cos(th * power) * sin(ph * power);
    Z->y = rp * sin(th * power) * sin(ph * power);
    Z->z = rp * cos(ph * power);

    // regular
    *Z += *P_in;

    // julia
    //*Z += (float3)(1,0,1);
}

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

static void coastalbrotIter(float3* Z, float* de, const bool julia, const float3* P_in) {
    float distance = length(*Z);

    float temp = distance;
    temp = pow(temp, 7.7f);
    *de = temp * (*de) * 7.7f;
    temp *= distance;

    Z->x = sin(sin(sin(M_PI_F / 3.0f + Z->x * M_PI_F)));
    Z->y = sin(sin(sin(M_PI_F / 3.0f + Z->y * M_PI_F)));
    Z->z = sin(sin(sin(M_PI_F / 3.0f + Z->z * M_PI_F)));

    *Z *= temp;

    // regular
    *Z += *P_in;

    // julia
    //*Z += (float3)(1,0,1);
}

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

static void sierpinski3dIter(float3* Z, float* de, const bool julia, const float3* P_in, const float scale, const float3 offset, const float3 rot) {
    float distance = length(*Z);

    float3 temp = *Z;

    if (Z->x - Z->y < 0.0f) Z->xy = Z->yx;
    if (Z->x - Z->z < 0.0f) Z->xz = Z->zx;
    if (Z->y - Z->z < 0.0f) Z->yz = Z->zy;
    if (Z->x + Z->y < 0.0f)
    {
        temp.x = -Z->y;
        Z->y = -Z->x;
        Z->x = temp.x;
    }
    if (Z->x + Z->z < 0.0f)
    {
        temp.x = -Z->z;
        Z->z = -Z->x;
        Z->x = temp.x;
    }
    if (Z->y + Z->z < 0.0f)
    {
        temp.y = -Z->z;
        Z->z = -Z->y;
        Z->y = temp.y;
    }

    *Z *= scale;
    *de *= scale;

    *Z -= offset;

    *Z = mtxPtMult( mtxRotate(rot) , *Z );

    // regular
    // works without adding anything
    //*Z += *P_in;

    // julia
    *Z += (float3)(0,0,1);
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