#ifndef VFT_OSL
#define VFT_OSL

// Trying to make OSL more compatible with OpenCL syntax
// Some features like vector swizzling are not supported in OSL, so they won't be fully compatible with the current code,
// but I am at least creating some helper macros and structs to speed up porting my OpenCL code
// Note that list of types and functions here is not complete, I implemented only those which were required by OpenCL functions

// functions

#define SIN(x) sin(x)
#define COS(x) cos(x)
#define POWR(x, y) pow((x), (y))
#define SQRT(x) sqrt(x)
#define LOG(x) log(x)
#define EXP(x) exp(x)
#define DIV(x, y) ((x) / (y))
#define NORMALIZE(x) normalize(x)
#define LENGTH(x) length(x)

#define M_PI_F          3.14159265358979323846
#define LARGE_NUMBER            10e10

// float3

struct float3
{
    float x, y, z;
};

float length(float3 vec)
{
    return length( vector(vec.x, vec.y, vec.z) );
}

float3 fabs(float3 a)
{
    vector vec = fabs( vector(a.x, a.y, a.z) );
    return float3(vec[0], vec[1], vec[2]);
}

float3 max(float3 a, float3 b)
{
    vector vec = max( vector(a.x, a.y, a.z), vector(b.x, b.y, b.z) );
    return float3(vec[0], vec[1], vec[2]);
}

float3 max(float3 a, float b)
{
    vector vec = max( vector(a.x, a.y, a.z), vector(b) );
    return float3(vec[0], vec[1], vec[2]);
}

float3 mix(float3 a, float3 b, float t)
{
    vector vec = mix( vector(a.x, a.y, a.z), vector(b.x, b.y, b.z), t );
    return float3(vec[0], vec[1], vec[2]);
}

float3 __operator__neg__ (float3 a)
{
    return float3(-a.x, -a.y, -a.z);
}

float3 __operator__add__ (float3 a, float3 b)
{
    return float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

float3 __operator__sub__ (float3 a, float3 b)
{
    return float3(a.x-b.x, a.y-b.y, a.z-b.z);
}

float3 __operator__mul__ (float3 a, float3 b)
{
    return float3(a.x*b.x, a.y*b.y, a.z*b.z);
}

float3 __operator__mul__ (float a, float3 b)
{
    return float3(a*b.x, a*b.y, a*b.z);
}

float3 __operator__div__ (float3 a, float3 b)
{
    return float3(a.x/b.x, a.y/b.y, a.z/b.z);
}

float3 __operator__mod__ (float3 a, float3 b)
{
    return float3(fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z));
}

int __operator__eq__ (float3 a, float3 b)
{
    if (a.x == b.x && a.y == b.y && a.z == b.z)
        return 1;
    else
        return 0;
}

int __operator__ne__ (float3 a, float3 b)
{
    return 1 - (a == b);
}

// float2

struct float2
{
    float x, y;
};

float length(float2 vec)
{
    return length( vector(vec.x, vec.y, 0) );
}

float2 normalize(float2 a)
{
    vector vec = normalize( vector(a.x, a.y, 0) );
    return float2(vec[0], vec[1]);
}

float dot(float2 a, float2 b)
{
    return dot( vector(a.x, a.y, 0), vector(b.x, b.y, 0) );
}

// float4

struct float4
{
    float x, y, z, w;
};


// functions

float length2(vector vec)
{
    return dot(vec, vec);
}

// point to plane distance
float distPointPlane(vector pt, vector plane_n, vector plane_point)
{
    float sb, sn, sd;
    vector point_proj;

    sn = -dot( plane_n, (pt - plane_point));
    sd = dot(plane_n, plane_n);
    sb = DIV(sn, sd);

    point_proj = pt + sb * plane_n;
    
    return LENGTH(pt - point_proj);
}

// returns rotation matrix specified by euler rotations in rot vector in degrees
matrix euler_rotation(vector rot)
{
    vector rot_rad = radians(-rot);

    float cosx = COS(rot_rad[0]);
    float sinx = SIN(rot_rad[0]);
    matrix x = matrix(1,            0,              0,              0,
                      0,            cosx,           -sinx,          0,
                      0,            sinx,          cosx,           0,
                      0,            0,              0,              1 );

    float cosy = COS(rot_rad[1]);
    float siny = SIN(rot_rad[1]);
    matrix y = matrix(cosy,         0,              siny,           0,
                      0,            1,              0,              0,
                      -siny,        0,              cosy,           0,
                      0,            0,              0,              1 );

    float cosz = COS(rot_rad[2]);
    float sinz = SIN(rot_rad[2]);
    matrix z = matrix(cosz,         -sinz,          0,              0,
                      sinz,         cosz,           0,              0,
                      0,            0,              1,              0,
                      0,            0,              0,              1 );

    matrix xform = 1;
    xform *= x;
    xform *= y;
    xform *= z;

    return xform;
}

#endif
