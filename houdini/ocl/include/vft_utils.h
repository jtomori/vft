#ifndef _VFT_UTILS
#define _VFT_UTILS



// export float attrib
static void vstore1(float dataIn, int i, global float* dataOut)
{
    dataOut[i] = dataIn;
}


////////////// sdf operations


// polynomial smooth min
static float sminPoly( float a, float b, float k )
{
    float h = clamp( 0.5 + 0.5 * (b - a) / k, 0.0, 1.0 );
    return mix( b, a, h ) - k * h * (1.0 - h);
}

// exponential smooth min
static float sminExp( float a, float b, float k )
{
    float res = exp( -k * a ) + exp( -k * b );
    return -log( res ) / k;
}

// union
static float sdfUnion( float a, float b )
{
    return min(a, b);
}

// smooth union
static float sdfUnionSmooth( float a, float b, float k )
{
    return sminPoly(a, b, k);
}

// subtraction
static float sdfSubtract( float b, float a )
{
    return max(-a, b);
}

// intersection
static float sdfIntersect( float a, float b )
{
    return max(a, b);
}

// blend
static float sdfBlend(float d1, float d2, float a)
{
    return (1 - a) * d1 + a * d2;
}

// infinitely repeat by a distance (c)
static float3 sdfRep( float3 p, float3 c )
{
    p = fmod(p,c) - .5f*c  ;
    return p;

}


////////////// debug


// print a mtx
static void printMtx(float16 m)
{
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
static void printVec(float3 a)
{
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "%2.2v3hlf\n", a );
    }
}

// print a float
static void printFloat(float a)
{
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "%2.8f\n", a );
    }
}



#endif