#ifndef VFT_FRACTALS
#define VFT_FRACTALS

// mapping of variables
// aux.r_dz -> de
// aux.r -> distance
// aux.de -> de
// dr -> de
// r -> distance
// Bailout -> max_distance
// Iterations -> max_iterations
// positive log_lin -> log, negative -> lin

// fractal formula sources:
// [M2] - Mandelbulber v2
// [M3D] - Mandelbulb 3D
// [WEB] - From the web

////////////// primitives
// [WEB] - http://iquilezles.org/www/articles/distfunctions/distfunctions.htm

// sphere: position, radius, center
static float sphere( float3 P, float rad )
{
    float dist = length(P) - rad;
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

// [WEB] - http://blog.hvidtfeldts.net/index.php/2011/09/
static void mandelbulbIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float power)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    // convert to polar coordinates
    float theta = acos((*Z).z/distance);
    float phi = atan2((*Z).y, (*Z).x);

    *de =  pow(distance, power-1) * power * (*de) + 1.0f;
    
    // scale and rotate the point
    float zr = pow(distance, power);
    theta *= power;
    phi *= power;
    
    // convert back to cartesian coordinates
    float3 newZ = zr * (float3)( sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta) );

    if (julia.x == 0.0f)
    {
        *Z = newZ + *P_in;
    }
    else {
        *Z = newZ + julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [WEB] - http://www.fractalforums.com/index.php?topic=2785.msg14893#msg14893
static void mandelboxIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float scale)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    float fixedRadius = 1.0f;
    float fR2 = fixedRadius * fixedRadius;
    float minRadius = 0.5f;
    float mR2 = minRadius * minRadius;

    if ((*Z).x > 1.0f) (*Z).x = 2.0f - (*Z).x;
    else if ((*Z).x < -1.0f) (*Z).x = -2.0f - (*Z).x;

    if ((*Z).y > 1.0f) (*Z).y = 2.0f - (*Z).y;
    else if ((*Z).y < -1.0f) (*Z).y = -2.0f - (*Z).y;

    if ((*Z).z > 1.0f) (*Z).z = 2.0f - (*Z).z;
    else if ((*Z).z < -1.0f) (*Z).z = -2.0f - (*Z).z;

    float r2 = (*Z).x*(*Z).x + (*Z).y*(*Z).y + (*Z).z*(*Z).z;

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

    if (julia.x == 0.0f)
    {
        *Z = *Z * scale + *P_in;
    }
    else
    {
        *Z = *Z * scale + julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Classic Mandelbulb Power 2 fractal - MandelbulbPower2Iteration
static void mandelbulbPower2Iter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia) 
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    *de = *de * 2.0f * distance;
    float x2 = (*Z).x * (*Z).x;
    float y2 = (*Z).y * (*Z).y;
    float z2 = (*Z).z * (*Z).z;
    float temp = 1.0f - z2 / (x2 + y2);
    float3 new;
    new.x = (x2 - y2) * temp;
    new.y = 2.0f * (*Z).x * (*Z).y * temp;
    new.z = -2.0f * (*Z).z * sqrt(x2 + y2);

    if (julia.x == 0.0f)
    {
        *Z = new + *P_in;
    }
    else
    {
        *Z = new + julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Menger Sponge formula created by Knighty - MengerSpongeIteration
static void mengerSpongeIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    (*Z).x = fabs((*Z).x);
    (*Z).y = fabs((*Z).y);
    (*Z).z = fabs((*Z).z);

    if ((*Z).x - (*Z).y < 0.0f) (*Z).xy = (*Z).yx;
    if ((*Z).x - (*Z).z < 0.0f) (*Z).xz = (*Z).zx;
    if ((*Z).y - (*Z).z < 0.0f) (*Z).yz = (*Z).zy;

    *Z *= 3.0f;

    (*Z).x -= 2.0f;
    (*Z).y -= 2.0f;
    if ((*Z).z > 1.0f) (*Z).z -= 2.0f;

    *de *= 3.0f;

    if (julia.x == 1.0f)
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Bristorbrot formula - BristorbrotIteration
static void bristorbrotIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    float3 new;
    new.x = (*Z).x * (*Z).x - (*Z).y * (*Z).y - (*Z).z * (*Z).z;
    new.y = (*Z).y * (2.0f * (*Z).x - (*Z).z);
    new.z = (*Z).z * (2.0f * (*Z).x + (*Z).y);
    
    *de = *de * 2.0f * distance;
    *Z = new;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Xenodreambuie - XenodreambuieIteration
static void xenodreambuieIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float power, float alpha, float beta)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    alpha = radians(alpha);
    beta = radians(beta);

    float rp = pow(distance, power - 1.0f);
    *de = rp * (*de) * power + 1.0f;
    rp *= distance;

    float th = atan2((*Z).y, (*Z).x) + beta;
    float ph = acos((*Z).z / distance) + alpha;

    if (fabs(ph) > 0.5f * M_PI_F) ph = sign(ph) * M_PI_F - ph;

    (*Z).x = rp * cos(th * power) * sin(ph * power);
    (*Z).y = rp * sin(th * power) * sin(ph * power);
    (*Z).z = rp * cos(ph * power);

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Sierpinski3D. made from Darkbeam's Sierpinski code from M3D - Sierpinski3dIteration
static void sierpinski3dIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float scale, const float3 offset, const float3 rot)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    float3 temp = *Z;

    if ((*Z).x - (*Z).y < 0.0f) (*Z).xy = (*Z).yx;
    if ((*Z).x - (*Z).z < 0.0f) (*Z).xz = (*Z).zx;
    if ((*Z).y - (*Z).z < 0.0f) (*Z).yz = (*Z).zy;
    if ((*Z).x + (*Z).y < 0.0f)
    {
        temp.x = -(*Z).y;
        (*Z).y = -(*Z).x;
        (*Z).x = temp.x;
    }
    if ((*Z).x + (*Z).z < 0.0f)
    {
        temp.x = -(*Z).z;
        (*Z).z = -(*Z).x;
        (*Z).x = temp.x;
    }
    if ((*Z).y + (*Z).z < 0.0f)
    {
        temp.y = -(*Z).z;
        (*Z).z = -(*Z).y;
        (*Z).y = temp.y;
    }

    *Z *= scale;
    *de *= scale;

    *Z -= offset;

    *Z = mtxPtMult( mtxRotate(rot) , *Z );

    if (julia.x == 1.0f)
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Menger Smooth - mengerSmoothIteration
static float mengerSmoothIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float scale, const float offset_s, const float3 offset_c, const float3 rot)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);
	
    float sc1 = scale - 1.0f;
	float sc2 = sc1 / scale;

    *Z = (float3)(sqrt((*Z).x * (*Z).x + offset_s), sqrt((*Z).y * (*Z).y + offset_s), sqrt((*Z).z * (*Z).z + offset_s));    

	float t;

	t = (*Z).x - (*Z).y;
	t = 0.5f * (t - sqrt(t * t + offset_s));
	(*Z).x = (*Z).x - t;
	(*Z).y = (*Z).y + t;

	t = (*Z).x - (*Z).z;
	t = 0.5f * (t - sqrt(t * t + offset_s));
	(*Z).x = (*Z).x - t;
	(*Z).z = (*Z).z + t;

	t = (*Z).y - (*Z).z;
	t = 0.5f * (t - sqrt(t * t + offset_s));
	(*Z).y = (*Z).y - t;
	(*Z).z = (*Z).z + t;

	(*Z).z = (*Z).z - offset_c.z * sc2;
	(*Z).z = -sqrt((*Z).z * (*Z).z + offset_s);
	(*Z).z = (*Z).z + offset_c.z * sc2;

	(*Z).x = scale * (*Z).x - offset_c.x * sc1;
	(*Z).y = scale * (*Z).y - offset_c.y * sc1;
	(*Z).z = scale * (*Z).z;

	*de *= scale;

    *Z = mtxPtMult( mtxRotate(rot) , *Z );

    if (julia.x == 1.0f)
    {
        *Z += *P_in;
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}


// testing new formulas

// amazing surf from M3D
// does nothing
static float amazingSurf(float3 P, float size)
{
//void AmazingSurfIteration(CVector4 &z, const sFractal *fractal, sExtendedAux &aux)
//{
//    // update aux.actualScale
//    aux.actualScale =
//        fractal).mandelbox.scale + fractal).mandelboxVary4D.scaleVary * (fabs(aux.actualScale) - 1.0);

//    CVector4 c = aux.const_c;
//    z.x = fabs(z.x + fractal).transformCommon.additionConstant111.x)
//                - fabs(z.x - fractal).transformCommon.additionConstant111.x) - z.x;
//    z.y = fabs(z.y + fractal).transformCommon.additionConstant111.y)
//                - fabs(z.y - fractal).transformCommon.additionConstant111.y) - z.y;
//    // no z fold

//    double rr = z.Dot(z);
//    if (fractal).transformCommon.functionEnabledFalse) // force cylinder fold
//        rr -= z.z * z.z;

//    double sqrtMinR = sqrt(fractal).transformCommon.minR05);
//    double dividend = rr < sqrtMinR ? sqrtMinR : min(rr, 1.0);

//    // use aux.actualScale
//    double m = aux.actualScale / dividend;

//    z *= (m - 1.0) * fractal).transformCommon.scale1 + 1.0;
//    // z *= m * fractal).transformCommon.scale1 + 1.0 * (1.0 - fractal).transformCommon.scale1);
//    aux.DE = aux.DE * fabs(m) + 1.0;

//    if (fractal).transformCommon.addCpixelEnabledFalse)
//        z += CVector4(c.y, c.x, c.z, c.w) * fractal).transformCommon.constantMultiplier111;

//    z = fractal).transformCommon.rotationMatrix.RotateVector(z);
//}
    //P /= size;
    float3 z = P;
    float dr = 1.0;
    int Iterations = 10;
    //int Bailout = 6;

    float actualScale = 1.39;
    float scaleVary = 0;
    float2 fold = (float2)(1.076562, 1.05);
    float minRad = 0.18;
    float scaleInf = 1;
    float auxScale = 1;

    for (int i = 0; i < Iterations ; i++)
    {
        //update aux.actualScale
        actualScale = actualScale + scaleVary * (fabs(actualScale) - 1.0f);
    
        //CVector4 c = aux.const_c;
        z.x = fabs(z.x + fold.x) - fabs(z.x - fold.x) - z.x;
        z.y = fabs(z.y + fold.y) - fabs(z.y - fold.y) - z.y;
        // no z fold
    
        float rr = z.x*z.x + z.y*z.y + z.z*z.z;
        //if (fractal).transformCommon.functionEnabledFalse) // force cylinder fold
        //    rr -= z.z * z.z;
    
        float sqrtMinR = sqrt(minRad);
        float dividend = rr < sqrtMinR ? sqrtMinR : min(rr, 1.0f);
    
        // use aux.actualScale
        //float m = auxScale / dividend;
        float m = actualScale / dividend;
    
        z *= (m - 1.0f) * scaleInf + 1.0f;
        dr = dr * fabs(m) + 1.0;
    
        //if (fractal).transformCommon.addCpixelEnabledFalse)
        //    z += CVector4(c.y, c.x, c.z, c.w) * fractal).transformCommon.constantMultiplier111;
    
        //z = fractal).transformCommon.rotationMatrix.RotateVector(z);
        float16 rotMtx = mtxRotate( (float3)(8.414060, 3.340000, 18.125000) );
        //z = mtxPtMult(rotMtx, z);

        z += P;
    }

    //float out = 0.5 * log(r) * r/dr;
    //return out * size;
    return dr;
}

// quaternion fractals, kind of works, but not sure what to do with z.w component :)
static float quaternion(float3 P, float size)
{
//void QuaternionIteration(CVector4 &z, const sFractal *fractal, sExtendedAux &aux)
//{
//    Q_UNUSED(fractal);

//    aux.r_dz = aux.r_dz * 2.0 * aux.r;
//    double newx = z.x * z.x - z.y * z.y - z.z * z.z - z.w * z.w;
//    double newy = 2.0 * z.x * z.y;
//    double newz = 2.0 * z.x * z.z;
//    double neww = 2.0 * z.x * z.w;
//    z.x = newx;
//    z.y = newy;
//    z.z = newz;
//    z.w = neww;
//}
    P /= size;

    float4 z = (float4)(P, 1);
    float dr = 1.0;
    float r = 0.0;
    int Iterations = 50;
    int Bailout = 40;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;

        dr = dr * 2.0f * r;
        float newx = z.x * z.x - z.y * z.y - z.z * z.z - z.w * z.w;
        float newy = 2.0f * z.x * z.y;
        float newz = 2.0f * z.x * z.z;
        float neww = 2.0f * z.x * z.w;
        z.x = newx;
        z.y = newy;
        z.z = newz;
        //z.w = neww;
    }

    float out = 0.5f * log(r) * r/dr;
    return out * size;
}

// quaternion3d
// kind of works, but does not exactly match M2 visual, but parameters deform it in a similar manner, I hardcoded some of the parameters into values, there are some buggy areas, missing parts, noisy normals etc.
static float quaternion3d(float3 P, float size)
{
//void Quaternion3dIteration(CVector4 &z, const sFractal *fractal, sExtendedAux &aux)
//{
//
//    aux.r_dz = aux.r_dz * 2.0 * aux.r;
//    z = CVector4(z.x * z.x - z.y * z.y - z.z * z.z, z.x * z.y, z.x * z.z, z.w);
//
//    double tempL = z.Length();
//    z *= fractal).transformCommon.constantMultiplier122;
//    // if (tempL < 1e-21) tempL = 1e-21;
//    CVector4 tempAvgScale = CVector4(z.x, z.y / 2.0, z.z / 2.0, z.w);
//    double avgScale = tempAvgScale.Length() / tempL;
//    double tempAux = aux.r_dz * avgScale;
//    aux.r_dz = aux.r_dz + (tempAux - aux.r_dz) * fractal).transformCommon.scaleA1;
//
//    if (fractal).transformCommon.rotationEnabled)
//        z = fractal).transformCommon.rotationMatrix.RotateVector(z);
//
//    z += fractal).transformCommon.additionConstant000;
//}
    P /= size;

    float4 z = (float4)(P, 1);
    float dr = 1.0;
    float r = 0.0;
    int Iterations = 250;
    int Bailout = 20;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;
        
        dr = dr * 2.0 * r;
        z = (float4)(z.x * z.x - z.y * z.y - z.z * z.z, z.x * z.y, z.x * z.z, z.w);

        float tempL = r;
        z *= (float4)(1,1,1,1);

        float4 tempAvgScale = (float4)(z.x, z.y / 2.0, z.z / 2.0, z.w);
        float avgScale = length(tempAvgScale) / tempL;
        float tempAux = dr * avgScale;
        dr = dr + (tempAux - dr) * 1.0f;

        //if (fractal).transformCommon.rotationEnabled)
        //    z = fractal).transformCommon.rotationMatrix.RotateVector(z);

        z += (float4)(0,0,0,0);
    }

    //float out = 0.5f * log(r) * r/dr;
    float out = r / dr;
    return out * size;
}

#endif