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
    const float4 ma[4] = { (float4)(a.s0123),
                           (float4)(a.s4567),
                           (float4)(a.s89ab),
                           (float4)(a.scdef) };
    
    const float4 mb[4] = { (float4)(b.s048c),
                           (float4)(b.s159d),
                           (float4)(b.s26ae),
                           (float4)(b.s37bf) };
    
    // compute matrix multiplication as a table of dot products
    float4 x[4] = {
            (float4)( dot(ma[0], mb[0]), dot(ma[0], mb[1]), dot(ma[0], mb[2]), dot(ma[0], mb[3]) ),
            (float4)( dot(ma[1], mb[0]), dot(ma[1], mb[1]), dot(ma[1], mb[2]), dot(ma[1], mb[3]) ),
            (float4)( dot(ma[2], mb[0]), dot(ma[2], mb[1]), dot(ma[2], mb[2]), dot(ma[2], mb[3]) ),
            (float4)( dot(ma[3], mb[0]), dot(ma[3], mb[1]), dot(ma[3], mb[2]), dot(ma[3], mb[3]) )
    };
    
    // float4 array to float16
    return (float16)(x[0].x,x[0].y,x[0].z,x[0].w,
                     x[1].x,x[1].y,x[1].z,x[1].w,
                     x[2].x,x[2].y,x[2].z,x[2].w,
                     x[3].x,x[3].y,x[3].z,x[3].w );
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

    const float cosx = cos(rot.x);
    const float sinx = sin(rot.x);
    const float16 x = (float16)(1,            0,              0,              0,
                                0,            cosx,           -sinx,          0,
                                0,            sinx,          cosx,           0,
                                0,            0,              0,              1 );

    const float cosy = cos(rot.y);
    const float siny = sin(rot.y);
    const float16 y = (float16)(cosy,         0,              siny,           0,
                                0,            1,              0,              0,
                                -siny,        0,              cosy,           0,
                                0,            0,              0,              1 );

    const float cosz = cos(rot.z);
    const float sinz = sin(rot.z);
    const float16 z = (float16)(cosz,         -sinz,          0,              0,
                                sinz,         cosz,           0,              0,
                                0,            0,              1,              0,
                                0,            0,              0,              1 );

    float16 xform = ident();
    xform = mtxMult(xform, x);
    xform = mtxMult(xform, y);
    xform = mtxMult(xform, z);

    return xform;
}


// multiplciation of a 4x4 matrix and a position vector (homogeneous)
static float3 mtxPtMult(float16 mtx, float3 vec)
{
    const float4 m[4] = { (float4)(mtx.s048c),
                          (float4)(mtx.s159d),
                          (float4)(mtx.s26ae),
                          (float4)(mtx.s37bf) };

    const float4 v = (float4)(vec.x, vec.y, vec.z, 1);
    const float w = 1 / ( dot(v, m[3]) );
    float3 x = (float3)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]) );    
    x *= w;

    return x;
}

// multiplciation of a 4x4 matrix and a directional vector, assumes affine matrix
static float3 mtxDirMult(float16 mtx, float3 vec)
{
    const float4 m[3] = { (float4)(mtx.s048c),
                          (float4)(mtx.s159d),
                          (float4)(mtx.s26ae) };

    const float4 v = (float4)(vec.x, vec.y, vec.z, 1);
    float3 x = (float3)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]) );

    return x;
}




#endif