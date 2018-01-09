#ifndef _VFT_MATH
#define _VFT_MATH



// identity 4x4 matrix
static float16 ident()
{
    return (float16)(1,0,0,0,
                     0,1,0,0,
                     0,0,1,0,
                     0,0,0,1);
}

// transpose a 4x4 matrix
static float16 trans(float16 m)
{
    float16 x;
    x = (float16)( m.s0,m.s4,m.s8,m.sc,
                   m.s1,m.s5,m.s9,m.sd,
                   m.s2,m.s6,m.sa,m.se,
                   m.s3,m.s7,m.sb,m.sf);
    return x;
}

// multiplication of two 4x4 matrices
static float16 mtxMult(float16 a, float16 b)
{
    // float16 to float4 array
    float4 ma[4];
    ma[0] = a.s0123;
    ma[1] = a.s4567;
    ma[2] = a.s89ab;
    ma[3] = a.scdef;

    float4 mb[4];
    mb[0] = b.s048c;
    mb[1] = b.s159d;
    mb[2] = b.s26ae;
    mb[3] = b.s37bf;
    
    // compute matrix multiplication as a table of dot products
    float4 x[4];
    x[0] = (float4)( dot(ma[0], mb[0]), dot(ma[0], mb[1]), dot(ma[0], mb[2]), dot(ma[0], mb[3]) );
    x[1] = (float4)( dot(ma[1], mb[0]), dot(ma[1], mb[1]), dot(ma[1], mb[2]), dot(ma[1], mb[3]) );
    x[2] = (float4)( dot(ma[2], mb[0]), dot(ma[2], mb[1]), dot(ma[2], mb[2]), dot(ma[2], mb[3]) );
    x[3] = (float4)( dot(ma[3], mb[0]), dot(ma[3], mb[1]), dot(ma[3], mb[2]), dot(ma[3], mb[3]) );
    
    // float4 array to float16
    float16 xOut = (float16)(x[0].x,x[0].y,x[0].z,x[0].w,
                             x[1].x,x[1].y,x[1].z,x[1].w,
                             x[2].x,x[2].y,x[2].z,x[2].w,
                             x[3].x,x[3].y,x[3].z,x[3].w );    

    return xOut;
}

// generates a 4x4 scaling matrix
static float16 mtxScale(float3 s)
{
    float16 x;
    x = (float16)(s.x,   0,   0,   0,
                  0  , s.y,   0,   0,
                  0  ,   0, s.z,   0,
                  0  ,   0,   0,   1);
    return x;
}

// generates a 4x4 rotation matrix in XYZ order, in degrees
static float16 mtxRotate(float3 rot)
{
    rot = radians(-rot);

    float16 x;
    float cosx = cos(rot.x);
    float sinx = sin(rot.x);
    x = (float16)(1,            0,              0,              0,
                  0,            cosx,           -sinx,          0,
                  0,            sinx,          cosx,           0,
                  0,            0,              0,              1 );

    float16 y;
    float cosy = cos(rot.y);
    float siny = sin(rot.y);
    y = (float16)(cosy,         0,              siny,           0,
                  0,            1,              0,              0,
                  -siny,        0,              cosy,           0,
                  0,            0,              0,              1 );

    float16 z;
    float cosz = cos(rot.z);
    float sinz = sin(rot.z);
    z = (float16)(cosz,         -sinz,          0,              0,
                  sinz,         cosz,           0,              0,
                  0,            0,              1,              0,
                  0,            0,              0,              1 );

    float16 xform = ident();
    xform = mtxMult(xform, x);
    xform = mtxMult(xform, y);
    xform = mtxMult(xform, z);

    return xform;
}


// multiplciation of a 4x4 matrix and a point (homogeneous)
static float3 mtxPtMult(float16 mtx, float3 vec)
{
    float4 m[4];
    m[0] = mtx.s048c;
    m[1] = mtx.s159d;
    m[2] = mtx.s26ae;
    m[3] = mtx.s37bf;

    float4 v = (float4)(vec.x, vec.y, vec.z, 1);
    float4 x = (float4)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]), dot(v, m[3]) );
    
    if (x.w != 1 && x.w != 0) {
        x.x /= x.w;
        x.y /= x.w;
        x.z /= x.w;
    }
    
    return (float3)(x.x, x.y, x.z);
}

// multiplciation of a 4x4 matrix and a vector
static float3 mtxVecMult(float16 mtx, float3 vec)
{
    float4 m[3];
    m[0] = mtx.s048c;
    m[1] = mtx.s159d;
    m[2] = mtx.s26ae;

    float4 v = (float4)(vec.x, vec.y, vec.z, 1);
    float3 x = (float3)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]) );

    return x;
}




#endif