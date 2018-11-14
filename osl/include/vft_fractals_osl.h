#ifndef VFT_FRACTALS_OSL
#define VFT_FRACTALS_OSL

// Porting OpenCL fractal functions from vft_fractals.h, vft_shading.h ...
// it means removing 
//                  " keyword, 
//                  "f" literal, 
//                  "const" keyword, 
//                  fixing expressions with vector swizzling, 
//                  removing pointers, 
//                  removing casts
//                  removing de, log_lin variables

// sphere: position, radius
float sphere( float3 P, float rad )
{
    float dist = LENGTH(P) - rad;
    return dist;
}

// box: position, size
float box( float3 P, float3 b )
{
  float3 d = fabs(P) - b;
  return min( max( d.x, max(d.y, d.z) ), 0.0) + LENGTH( max(d, 0.0) );
}

// round box: position, size, roundness
float roundBox( float3 P, float3 b, float r )
{
  return LENGTH( max( fabs(P) - b, 0.0 ) )-r;
}

// torus: position, size (radius, thickness)
float torus( float3 P, float2 t )
{
  float2 q = float2(LENGTH( float2(P.x, P.z) )-t.x,P.y);
  return LENGTH(q)-t.y;
}

// infinite cone
float cone( float3 P, float2 c )
{
    c = NORMALIZE(c);
    float q = LENGTH(float2(P.x, P.y));
    return dot( c, float2(q, P.z) );
}

// [WEB] - http://blog.hvidtfeldts.net/index.php/2011/09/
void mandelbulbIter(float3 Z, float3 P_in, float weight, float4 julia, float power)
{
    float3 Z_orig = Z;
    
    float distance = LENGTH(Z);

    // convert to polar coordinates
    float theta = acos( DIV(Z.z, distance));
    float phi = atan2(Z.y, Z.x);
    
    // scale and rotate the point
    float zr = POWR(distance, power);
    theta *= power;
    phi *= power;
    
    // convert back to cartesian coordinates
    float3 newZ = zr * float3( SIN(theta)*COS(phi), SIN(phi)*SIN(theta), COS(theta) );

    if (julia.x == 0.0)
    {
        Z = newZ + P_in;
    }
    else 
    {
        Z = newZ + float3(julia.y, julia.z, julia.w);
    }
}

#endif
