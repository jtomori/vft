#define NULL                    0

#define ORBITS_ARRAY_LENGTH     4
#define ENABLE_DELTA_DE         0

// because of weird bugs sometimes I have to use this custom mix function
float3 mix3(float3 x, float3 y, const float t)
{
    return x + (y - x) * t;
}

// export float attrib
static void vstore1(float dataIn, int i, global float* dataOut)
{
    dataOut[i] = dataIn;
}

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
    sb = sn / sd;

    point_proj = point + sb * plane_n;
    
    return length(point - point_proj);
}

// identity 4x4 matrix
static float16 mtxIdent()
{
    return (float16)(1,0,0,0,
                     0,1,0,0,
                     0,0,1,0,
                     0,0,0,1);
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

// multiplciation of a 4x4 matrix and a position vector (homogeneous)
static float3 mtxPtMult(float16 mtx, float3 vec)
{
    const float4 m[4] = { (float4)(mtx.s048c),
                          (float4)(mtx.s159d),
                          (float4)(mtx.s26ae),
                          (float4)(mtx.s37bf) };

    const float4 v = (float4)(vec.x, vec.y, vec.z, 1);
    float4 x = (float4)( dot(v, m[0]), dot(v, m[1]), dot(v, m[2]), dot(v, m[3]) );
    x = x/x.w;

    return (float3)(x.xyz);
}

// [M2] - Menger Sponge formula created by Knighty - MengerSpongeIteration
static void mengerSpongeIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
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

    *de *= 3.0f;

    if (julia.x == 1.0f)
    {
        *Z += julia.yzw;
    }

    *Z = mix3(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Bristorbrot formula - BristorbrotIteration
void bristorbrotIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    float3 new;
    new.x = Z->x * Z->x - Z->y * Z->y - Z->z * Z->z;
    new.y = Z->y * (2.0f * Z->x - Z->z);
    new.z = Z->z * (2.0f * Z->x + Z->y);
    
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

static float hybrid(float3 P_in, const int max_iterations, const float max_distance, const float size, const int calculate_orbits, float* orbit_colors)
{
    P_in /= size;
    float3 Z = P_in;
    float de = 1.0f;
    float distance;
    float out_de;
    int log_lin = 0; // positive is log, negative is lin

    // init orbit traps variables
    float3 orbit_pt = (float3)(0.0f,0.0f,0.0f);
    float orbit_pt_dist = 1e20f;
    float2 orbit_plane = (float2)(1.0f,0.0f);
    float3 orbit_plane_origin = (float3)(0.0f);
    float3 orbit_plane_dist = (float3)(1e20f);
    float orbit_coord_dist = 1e20f;
    float orbit_sphere_rad = 1.0f;
    float orbit_sphere_dist = 1e20f;
    float3 orbit_axis_dist = 1e20f;

    // fractal loop
    for (int i = 0; i < max_iterations; i++)
    {
        distance = length(Z);
        if (distance > max_distance) break;
        
        // fractals stack
        //bristorbrotIter(&Z, &de, &P_in, &log_lin, 1.0f, (float4)(0.0f, 1.3f, 3.3f, 0.0f)); // log
        mengerSpongeIter(&Z, &de, &P_in, &log_lin, 1.0f, (float4)(0.0f, 0.0f, 1.0f, 0.0f)); // lin

        // orbit traps calculations
        if (calculate_orbits == 1)
        {
            orbit_pt_dist = min(orbit_pt_dist, length2(Z - orbit_pt));
            orbit_plane_dist.x = min( orbit_plane_dist.x, distPointPlane(Z, orbit_plane.xyy, orbit_plane_origin) );
            orbit_plane_dist.y = min( orbit_plane_dist.y, distPointPlane(Z, orbit_plane.yxy, orbit_plane_origin) );
            orbit_plane_dist.z = min( orbit_plane_dist.z, distPointPlane(Z, orbit_plane.yyx, orbit_plane_origin) );
        }
    }

    // outputting orbit traps
    if (calculate_orbits == 1)
    {
        orbit_colors[0] = sqrt(orbit_pt_dist); // distance to point at specified coordinates
        orbit_colors[1] = orbit_plane_dist.x; // distance to YZ plane
        orbit_colors[2] = orbit_plane_dist.y; // distance to XZ plane
        orbit_colors[3] = orbit_plane_dist.z; // distance to XY plane
    }

    // automatic determining DE mode based on log_lin value
    if (log_lin >= 0) out_de = 0.5f * log(distance) * distance/de;
    else out_de = (distance) / fabs(de);

    return out_de * size;
}

// scene setup
static float scene( float3 P, const int final, float* orbit_colors, float3* N ) {
    float dist_out;

    float shape1 = hybrid(P, 10, 10.0f, 1.0f, final, orbit_colors);     

    dist_out = shape1;

    return dist_out;
}

//// main function
kernel void marchPerspCam(
        int P_length, global float* P,
        int planeZ_length, global float* planeZ,
        int width_length, global float* width,
        int height_length, global float* height,
        int px_length, global float* px,
        int camXform_length, global float* camXform,
        int camPos_length, global float* camPos,
        int N_length, global float* N,
        int iRel_length, global float* iRel,
        int Cd_length, global float* Cd,
        int orbits_length, global int* orbits_index, global float* orbits
        )
{
    // get current point id
    const int idx = get_global_id(0);

    // if current point is not valid, then end
    if ( idx >= P_length ) return;

    // read in P attrib
    const float3 pixel_P_origin = vload3(idx, P);
    float3 pixel_P_world = pixel_P_origin;

    //// transforming to near img plane

    // move to near img plane
    pixel_P_world.z = planeZ[0];

    // compute scale of near img plane
    const float16 near_plane_scale = mtxScale( (float3)(width[0]-px[0], height[0]-px[0], 1.0f) );

    // read in cam world matrix
    const float16 cam_xform_world = (float16)(camXform[0],camXform[1],camXform[2],camXform[3],
                                              camXform[4],camXform[5],camXform[6],camXform[7],
                                              camXform[8],camXform[9],camXform[10],camXform[11],
                                              camXform[12],camXform[13],camXform[14],camXform[15] );

    // create a mtx to hold transformations
    float16 near_plane_xform = mtxIdent();

    // apply transformations, also produce alternative matrix with scaled near plane
    near_plane_xform = mtxMult(near_plane_xform, near_plane_scale);
    float16 near_plane_xform_scaled = mtxMult(near_plane_xform, mtxScale( (float3)(100000.0f) ) );
    near_plane_xform = mtxMult(near_plane_xform, cam_xform_world);
    near_plane_xform_scaled = mtxMult(near_plane_xform_scaled, cam_xform_world);

    // create a scaled near plane position for more accurate ray_dir calculation
    float3 pixel_P_world_scaled = mtxPtMult(near_plane_xform_scaled, pixel_P_world);

    // transform pixels into near img plane
    pixel_P_world = mtxPtMult(near_plane_xform, pixel_P_world);

    // get camera world space position and compute ray direction vector
    const float3 cam_P_world = (float3)(camPos[0], camPos[1], camPos[2]);
    const float3 ray_dir = normalize(pixel_P_world_scaled - cam_P_world);

    //// raymarching

    // raymarch settings, initialize variables
    float3 color = (float3)(0.0f);    
    float orbit_colors[ORBITS_ARRAY_LENGTH];
    float3 Cd_out = (float3)(1.0f);
    float3 N_grad;

    float3 ray_P_world = pixel_P_world;
    float cam_dist = scene(cam_P_world, 0, NULL, NULL);
    float de = 0.0f;
    int i = 0;

    // quality settings
    float step_size = 0.6f;
    float iso_limit_mult = 1.0f;
    float ray_dist = planeZ[0];
    const int max_steps = 700;
    const float max_dist = 1000.0f;

    float iso_limit = cam_dist * 0.0001f * iso_limit_mult;  

    // coloring settings
    float Cd_mix_N = 0.9f;
    float Cd_mix_orbit = 0.9f;
    float Cd_mix_AO = 0.9f;

    // raymarching loop
    #pragma unroll
    for (i=0; i<max_steps; i++)
    {
        de = scene(ray_P_world, 0, NULL, NULL) * step_size;

        if ( de <= iso_limit || ray_dist >= max_dist )
        {
            de = scene(ray_P_world, 1, orbit_colors, &N_grad) * step_size;
            break;
        }

        ray_dist += de;
        ray_P_world += ray_dir * de;
    }

    // relative amount of steps
    float i_rel = (float)(i)/(float)(max_steps);
    i_rel = 1.0f-pow(i_rel, 1.0f/3.0f);

    // remove missed
    if ( de > iso_limit )
    {
        i_rel = -1.0f;
    }
    else
    {
        // AO
        float AO = 1.0f;
        float AO_occ = 0.0f;
        float AO_sca = 1.0f;

        color.x = orbit_colors[1];
        color.y = orbit_colors[2];
        color.z = orbit_colors[3];
        color *= orbit_colors[0];

        Cd_out = mix(Cd_out, Cd_out * fabs(N_grad), Cd_mix_N);
        Cd_out = mix(Cd_out, color, Cd_mix_orbit);
        Cd_out = mix(Cd_out, Cd_out * AO, Cd_mix_AO);
    }

    // export attribs
    vstore3(ray_P_world, idx, P);
    vstore3(N_grad, idx, N);
    vstore3(Cd_out, idx, Cd);
    vstore1(i_rel, idx, iRel);

    int orbits_idx_start = orbits_index[idx];
    int orbits_idx_end = orbits_idx_start + ORBITS_ARRAY_LENGTH;
    #pragma unroll
    for (int j=orbits_idx_start; j<orbits_idx_end; j++)
    {
        orbits[j] = orbit_colors[j-orbits_idx_start];
    }
}