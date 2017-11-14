//// implicit shape functions

// sphere
float sphere( float3 P, float rad ) {
    float dist = length(P) - rad;
    return dist;
}

// mandelbulb
float mandelbulb( float3 P ) {
    float3 z = P;
    float dr = 1.0;
    float r = 0.0;
    
    int Iterations = 7;
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
    return 0.5 * log(r) * r/dr;
}

// scene - construction of the scene for raymarching
float scene( float3 P ) {
    float dist;
    //dist = sphere(P, 1.5);
    dist = mandelbulb(P);
    return dist;
}

//// math functions

// import vector attrib
float3 getData3(global float* data)
{
    int i = get_global_id(0);
    return (float3)(data[i*3], data[i*3+1], data[i*3+2]);
}

// export vector attrib
void setData3(global float* dataOut, float3 dataIn)
{
    int i = get_global_id(0);
    dataOut[i * 3]     = dataIn.x;
    dataOut[i * 3 + 1] = dataIn.y;
    dataOut[i * 3 + 2] = dataIn.z;
}

// export float attrib
void setData1(global float* dataOut, float dataIn)
{
    int i = get_global_id(0);
    dataOut[i] = dataIn;
}

// transpose a 4x4 matrix
float16 trans(float16 m)
{
    float16 x;
    x = (float16)( m.s0,m.s4,m.s8,m.sc,
                   m.s1,m.s5,m.s9,m.sd,
                   m.s2,m.s6,m.sa,m.se,
                   m.s3,m.s7,m.sb,m.sf);
    return x;
}

// multiplication of two 4x4 matrices
float16 mtxMult(float16 a, float16 b)
{
    // float16 to float4 array
    float4 ma[4];
    ma[0] = (float4)(a.s0,a.s1,a.s2,a.s3);
    ma[1] = (float4)(a.s4,a.s5,a.s6,a.s7);
    ma[2] = (float4)(a.s8,a.s9,a.sa,a.sb);
    ma[3] = (float4)(a.sc,a.sd,a.se,a.sf);

    float4 mb[4];
    mb[0] = (float4)(b.s0,b.s4,b.s8,b.sc);
    mb[1] = (float4)(b.s1,b.s5,b.s9,b.sd);
    mb[2] = (float4)(b.s2,b.s6,b.sa,b.se);
    mb[3] = (float4)(b.s3,b.s7,b.sb,b.sf);
    
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
float16 mtxScale(float3 s)
{
    float16 x;
    x = (float16)(s.x,   0,   0,   0,
                  0  , s.y,   0,   0,
                  0  ,   0, s.z,   0,
                  0  ,   0,   0,   1);
    return x;
}

// multiplciation of a 4x4 matrix and a point (homogeneous)
float3 mtxPtMult(float16 mtx, float3 vec)
{
    float4 m[4];
    m[0] = (float4)(mtx.s0,mtx.s4,mtx.s8,mtx.sc);
    m[1] = (float4)(mtx.s1,mtx.s5,mtx.s9,mtx.sd);
    m[2] = (float4)(mtx.s2,mtx.s6,mtx.sa,mtx.se);
    m[3] = (float4)(mtx.s3,mtx.s7,mtx.sb,mtx.sf);

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
float3 mtxVecMult(float16 mtx, float3 vec)
{
    float4 m[3];
    m[0] = (float4)(mtx.s0,mtx.s4,mtx.s8,mtx.sc);
    m[1] = (float4)(mtx.s1,mtx.s5,mtx.s9,mtx.sd);
    m[2] = (float4)(mtx.s2,mtx.s6,mtx.sa,mtx.se);

    float4 v = (float4)(vec.x, vec.y, vec.z, 1);
    float3 x = (float3)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]) );

    return x;
}

// identity 4x4 matrix
float16 ident()
{
    return (float16)(1,0,0,0,
                     0,1,0,0,
                     0,0,1,0,
                     0,0,0,1);
}

//// debug functions

// print a mtx
void printMtx(float16 m) {
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "\n" );
        printf( "%2.2v4hlf\n", (float4)(m.s0, m.s1, m.s2, m.s3 ) );
        printf( "%2.2v4hlf\n", (float4)(m.s4, m.s5, m.s6, m.s7 ) );
        printf( "%2.2v4hlf\n", (float4)(m.s8, m.s9, m.sa, m.sb ) );
        printf( "%2.2v4hlf\n", (float4)(m.sc, m.sd, m.se, m.sf ) );
        printf( "\n" );
    }
}

// print a vector
void printVec(float3 a) {
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "%2.2v3hlf\n", a );
    }
}

// print a float
void printFloat(float a) {
    int idx = get_global_id(0);
    if (idx == 0) {
        printf( "%2.8f\n", a );
    }
}

//// main function

kernel void march(
        int P_length, global float* P,
        int planeZ_length, global float* planeZ,
        int width_length, global float* width,
        int height_length, global float* height,
        int px_length, global float* px,
        int camXform_length, global float* camXform,
        int camPos_length, global float* camPos,
        int N_length, global float* N,
        int iRel_length, global float* iRel,
        int Cd_length, global float* Cd
)
{
    // get current point id
    int idx = get_global_id(0);
    
    // if current point is not valid, then end
    if ( idx >= P_length )
        return;

    // read in P attrib
    float3 P_in = getData3(P);
    float3 P_out = P_in;

    //// transforming to near img plane
    
    // move to near img plane
    P_out.z = planeZ[0];
    
    // compute scale of near img plane
    float16 scale = mtxScale( (float3)(width[0]-px[0], height[0]-px[0], 1) );
    
    // read in cam world matrix
    float16 cam = (float16)(camXform[0],camXform[1],camXform[2],camXform[3],
                            camXform[4],camXform[5],camXform[6],camXform[7],
                            camXform[8],camXform[9],camXform[10],camXform[11],
                            camXform[12],camXform[13],camXform[14],camXform[15] );

    // create a mtx to hold transformations
    float16 xform = ident();
    
    // apply transformations
    xform = mtxMult(xform, scale);
    xform = mtxMult(xform, cam);

    // transform points into near img plane
    P_out = mtxPtMult(xform, P_out);    
    
    // get camera world space position and compute ray direction vector
    float3 camP = (float3)(camPos[0], camPos[1], camPos[2]);
    float3 rayDir = normalize(P_out - camP);
    
    //// raymarching
    
    // raymarch settings
    float dist;
    int i = 0;
    int max = 300;
    float stepSize = 0.5;
    float iso = 0.001;
    float t = planeZ[0];
    float maxDist = 20;
    
    float sphereR = 1.5;
    
    // raymarch
    for (i; i<max; i++) {
        dist = scene(P_out);
        //if ( dist <= iso * (t/300) || t >= maxDist ) break;
        if ( dist <= iso || t >= maxDist ) break;
        dist *= stepSize;
        t += dist;
        P_out += rayDir * dist;
    }
    
    // compute N
    float e = iso;
    float3 N_out = (float3)(0);
    float3 ePos[6] = { P_out + (float3)(e,0,0),
                       P_out - (float3)(e,0,0),
                       P_out + (float3)(0,e,0),
                       P_out - (float3)(0,e,0),
                       P_out + (float3)(0,0,e),
                       P_out - (float3)(0,0,e) };
                
    N_out = (float3)( scene( ePos[0] ) - scene( ePos[1] ),
                      scene( ePos[2] ) - scene( ePos[3] ),
                      scene( ePos[4] ) - scene( ePos[5] ) );
    
    N_out = normalize(N_out);
    
    // relative amount of steps
    float iRel_out = (float)(i)/(float)(max);
    
    // remove missed
    if ( dist > iso ) {
        //P_out = (float3)(0);
        iRel_out = -1;
    }
    
    // Cd for viz
    float3 Cd_out = fabs(N_out);
    Cd_out *= (1 - iRel_out);
    
    // export attribs
    setData3(P, P_out);
    setData3(N, N_out);
    setData3(Cd, Cd_out);
    setData1(iRel, iRel_out);
}