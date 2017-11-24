#ifndef __RAYMARCHING_FUNCS_H__
#define __RAYMARCHING_FUNCS_H__


//////////////////////////////////////////// shape functions

// sphere: position, radius, center
static float sphere( float3 P, float rad, float3 center ) {
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

// mandelbulb, with size control: s
static float mandelbulb( float3 P, float s ) {
    P /= s;

    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    
    int Iterations = 11;
    int Bailout = 2;
    int Power = 5;
    
    for (int i = 0; i < Iterations ; i++) {
            r = length(z);
            if (r>Bailout) break;
            
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
    return out * s;
}

//////////////////////////////////////////// shape operations

// union
static float sdfUnion( float a, float b ) {
    return min(a, b);
}

// substraction
static float sdfSubstract( float a, float b ) {
    return max(-a, b);
}

// intersection
static float sdfIntersect( float a, float b ) {
    return max(a, b);
}

// infinitely repeat by a distance (c)
static float3 sdfRep( float3 p, float3 c ){
    p = fmod(p,c) - .5f*c  ;
    return p;

}

//////////////////////////////////////////// data, math functions

// export float attrib
static void vstore1(float dataIn, int i, global float* dataOut)
{
    dataOut[i] = dataIn;
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

// generates a 4x4 matrix that in mtx multiplication will scale the other matrix by float3
static float16 mtxScale(float3 s)
{
    float16 x;
    x = (float16)(s.x,   0,   0,   0,
                  0  , s.y,   0,   0,
                  0  ,   0, s.z,   0,
                  0  ,   0,   0,   1);
    return x;
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

// identity 4x4 matrix
static float16 ident()
{
    return (float16)(1,0,0,0,
                     0,1,0,0,
                     0,0,1,0,
                     0,0,0,1);
}


//////////////////////////////////////////// debug functions

// print a mtx
static void printMtx(float16 m) {
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "\n" );
        printf( "%2.2v4hlf\n", m.s0123 );
        printf( "%2.2v4hlf\n", m.s4567 );
        printf( "%2.2v4hlf\n", m.s89ab );
        printf( "%2.2v4hlf\n", m.scdef );
        printf( "\n" );
    }
}

// print a vector
static void printVec(float3 a) {
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "%2.2v3hlf\n", a );
    }
}

// print a float
static void printFloat(float a) {
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "%2.8f\n", a );
    }
}


#endif