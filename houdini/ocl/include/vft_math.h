#ifndef VFT_MATH
#define VFT_MATH

// some of the functions here are ported/derived from Three.js library for WebGL:
// https://github.com/mrdoob/three.js/tree/master/src/math

// length2
static float length2(float3 vec)
{
    return (float)(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

// point to plane distance
static float distPointPlane(float3 point, float3 plane_n, float3 plane_point)
{
    float sb, sn, sd;
    float3 point_proj;

    sn = -dot( plane_n, (point - plane_point));
    sd = dot(plane_n, plane_n);
    sb = DIV(sn, sd);

    point_proj = point + sb * plane_n;
    
    return LENGTH(point - point_proj);
}

// identity 4x4 matrix
static float16 mtxIdent()
{
    return (float16)(1,0,0,0,
                     0,1,0,0,
                     0,0,1,0,
                     0,0,0,1);
}

// transpose a 4x4 matrix
static float16 mtxTranspose(float16 m)
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

// generates a 4x4 translation matrix
static float16 mtxTranslate(float3 t)
{
    float16 x;
    x = (float16)(1,   0,   0,   0,
                  0,   1,   0,   0,
                  0,   0,   1,   0,
                  t.x, t.y, t.z, 1);
    return x;
}

// generates a 4x4 rotation matrix in XYZ order, in degrees
static float16 mtxRotate(float3 rot)
{
    rot = radians(-rot);

    const float cosx = COS(rot.x);
    const float sinx = SIN(rot.x);
    const float16 x = (float16)(1,            0,              0,              0,
                                0,            cosx,           -sinx,          0,
                                0,            sinx,          cosx,           0,
                                0,            0,              0,              1 );

    const float cosy = COS(rot.y);
    const float siny = SIN(rot.y);
    const float16 y = (float16)(cosy,         0,              siny,           0,
                                0,            1,              0,              0,
                                -siny,        0,              cosy,           0,
                                0,            0,              0,              1 );

    const float cosz = COS(rot.z);
    const float sinz = SIN(rot.z);
    const float16 z = (float16)(cosz,         -sinz,          0,              0,
                                sinz,         cosz,           0,              0,
                                0,            0,              1,              0,
                                0,            0,              0,              1 );

    float16 xform = mtxIdent();
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
    float4 x = (float4)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]), dot(v, m[3]) );
    x = DIV(x, x.w);

    return (float3)(x.xyz);
}

// multiplciation of a 4x4 matrix and a directional vector, assumes affine matrix
static float3 mtxDirMult(float16 mtx, float3 vec)
{
    const float4 m[3] = { (float4)(mtx.s048c),
                          (float4)(mtx.s159d),
                          (float4)(mtx.s26ae) };

    const float4 v = (float4)(vec.x, vec.y, vec.z, 1.0f);
    float3 x = (float3)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]) );

    return x;
}

// inverts a 4x4 matrix
static float16 mtxInvert(float16 me)
{
    float16 te;

    float   
    n11 = me.s0, n21 = me.s1, n31 = me.s2, n41 = me.s3,
    n12 = me.s4, n22 = me.s5, n32 = me.s6, n42 = me.s7,
    n13 = me.s8, n23 = me.s9, n33 = me.sa, n43 = me.sb,
    n14 = me.sc, n24 = me.sd, n34 = me.se, n44 = me.sf,

    t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
    t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
    t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
    t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;
    
    float det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

    //if ( det == 0 ) return mtxIdent(); // assuming invertible matrices only

    float detInv = DIV(1.0f, det);

    te.s0 = t11 * detInv;
    te.s1 = ( n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44 ) * detInv;
    te.s2 = ( n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44 ) * detInv;
    te.s3 = ( n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43 ) * detInv;

    te.s4 = t12 * detInv;
    te.s5 = ( n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44 ) * detInv;
    te.s6 = ( n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44 ) * detInv;
    te.s7 = ( n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43 ) * detInv;

    te.s8 = t13 * detInv;
    te.s9 = ( n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44 ) * detInv;
    te.sa = ( n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44 ) * detInv;
    te.sb = ( n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43 ) * detInv;

    te.sc = t14 * detInv;
    te.sd = ( n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34 ) * detInv;
    te.se = ( n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34 ) * detInv;
    te.sf = ( n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33 ) * detInv;

    return te;

}

#endif
