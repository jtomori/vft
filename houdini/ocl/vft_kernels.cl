#include "vft_utils.h"
#include "vft_math.h"
#include "vft_fractals.h"

#define NULL                    0

#define LARGE_NUMBER            10e10

#define ORBITS_ARRAY_LENGTH     9
#define ENABLE_DELTA_DE         0

// forward func declarations
float3 compute_N(const float*, const float3*);
float compute_AO(const float3*, const float3*);
float hybrid(float3, const int, const float, const float, const int, float*, float3*, const int);

// contains fractal combinations
void fractal_stack(float3* Z, float* de, const float3* P_in, int* log_lin, const int stack)
{
    if (stack == 0)
    {
        //mandelbulbPower2Iter(Z, de, P_in, log_lin, 1.0f, (float4)(1.0f, 0.3f, 0.5f, 0.2f)); // log
        //bristorbrotIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 1.3f, 3.3f, 0.0f)); // log
        //xenodreambuieIter(Z, de, P_in, log_lin, 1.0f, (float4)(1.0f, 1.0f, 0.0f, 0.0f), 9.0f, 0.0f, 0.0f); // log
        //mandelboxIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 1.0f, 3.0f, 4.0f), 3.0f); // lin
        //mandelbulbIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 1.0f, 0.0f, 0.0f), 8.0f); // log
        mengerSpongeIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 1.0f, 0.0f)); // lin
        //sierpinski3dIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.5f), 2.0f, (float3)(1.0f, 1.0f, 1.0f), (float3)(0.0f, 0.0f, 0.0f) ); // lin
    }
    if (stack == 1)
    {
        sierpinski3dIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.5f), 2.0f, (float3)(1.0f, 1.0f, 1.0f), (float3)(0.0f, 0.0f, 0.0f) ); // lin        
    }
    if (stack == 2)
    {

    }
}

// scene setup - setting of coordinates and shapes in them
float scene( float3 P, const int final, float* orbit_colors, float3* N ) {
    float dist_out;

    float shape1 = hybrid(P, 10, 10.0f, 1.0f, final, orbit_colors, N, 0);
    float shape2 = hybrid(P - (float3)(2.0f), 10, 10.0f, 1.0f, final, orbit_colors, N, 1);

    dist_out = sdfUnion(shape1, shape2);

    return dist_out;
}

// main function
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
    float Cd_mix_AO = 0.6f;
    float AO = 1.0f;

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
        // compute N and AO only when not using DELTA DE mode
        #if !ENABLE_DELTA_DE
            N_grad = compute_N(&iso_limit, &ray_P_world);
            AO = compute_AO(&N_grad, &ray_P_world);
        #endif

        // default color assignment
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


//// shading, shape functions


// compute Normals
// based on "Modeling with distance functions" article from Inigo Quilez
float3 compute_N(const float* iso_limit, const float3* ray_P_world)
{
    float3 N_grad;
    float2 e = (float2)(1.0f, -1.0f) * (*iso_limit) * 0.01f;

    N_grad = normalize( e.xyy * scene( (*ray_P_world) + e.xyy, 0, NULL, NULL) + 
                        e.yyx * scene( (*ray_P_world) + e.yyx, 0, NULL, NULL) + 
                        e.yxy * scene( (*ray_P_world) + e.yxy, 0, NULL, NULL) + 
                        e.xxx * scene( (*ray_P_world) + e.xxx, 0, NULL, NULL) );
    return N_grad;
}

// compute Ambient Occlusion
// based on "Modeling with distance functions" article from Inigo Quilez
float compute_AO(const float3* N, const float3* ray_P_world)
{
    float AO = 1.0f;
    float AO_occ = 0.0f;
    float AO_sca = 1.0f;

    #pragma unroll
    for(int j=0; j<5; j++)
    {
        float AO_hr = 0.01f + 0.12f * (float)(j)/4.0f;
        float3 AO_pos =  (*N) * AO_hr + (*ray_P_world);
        float AO_dd = scene(AO_pos, 0, NULL, NULL);
        AO_occ += -(AO_dd-AO_hr)*AO_sca;
        AO_sca *= 0.95f;
    }
    
    AO = clamp( 1.0f - 3.4f * AO_occ, 0.0f, 1.0f );
    AO = pow(AO, 0.8f);

    return AO;
}

// primitive shape function
/*float primitive(float3 P, const int final, float* orbit_colors, float3* N)
{
    float out_distance;

    out_distance = sphere(P, 1, (float3)(1.5f, 0.0f, 0.0f));

    if (final == 1)
    {
        #pragma unroll
        for (int i=0; i<ORBITS_ARRAY_LENGTH; i++)
        {
            orbit_colors[i] = 1.0f;
        }

        // delta DE mode does not work well with multiple fractal objects
        #if ENABLE_DELTA_DE
            float iso_limit = 0.0001f;
            *N = compute_N(&iso_limit, &P);
        #endif
    }

    return out_distance;
}*/

// hybrid shape function - contains fractal loop, fractals combination is defined in fractals_stack(), in case of DELTA DE mode this also outputs N
float hybrid(float3 P_in, const int max_iterations, const float max_distance, const float size, const int final, float* orbit_colors, float3* N, const int stack)
{
    P_in /= size;
    float3 Z = P_in;
    float de = 1.0f;
    float distance;
    float out_distance;
    int log_lin = 0;

    // init orbit traps variables
    float3 orbit_pt = (float3)(0.0f,0.0f,0.0f);
    float orbit_pt_dist = LARGE_NUMBER;
    float2 orbit_plane = (float2)(1.0f,0.0f);
    float3 orbit_plane_origin = (float3)(0.0f);
    float3 orbit_plane_dist = (float3)(LARGE_NUMBER);
    float orbit_coord_dist = LARGE_NUMBER;
    float orbit_sphere_rad = 1.0f;
    float orbit_sphere_dist = LARGE_NUMBER;
    float3 orbit_axis_dist = LARGE_NUMBER;

    // fractal loop
    int i = 0;
    #pragma unroll
    for (i = 0; i < max_iterations; i++)
    {
        distance = length(Z);
        if (distance > max_distance) break;
        
        fractal_stack(&Z, &de, &P_in, &log_lin, stack);

        // orbit traps calculations
        if (final == 1)
        {
            orbit_pt_dist = min(orbit_pt_dist, length2(Z - orbit_pt));
            orbit_plane_dist.x = min( orbit_plane_dist.x, distPointPlane(Z, orbit_plane.xyy, orbit_plane_origin) );
            orbit_plane_dist.y = min( orbit_plane_dist.y, distPointPlane(Z, orbit_plane.yxy, orbit_plane_origin) );
            orbit_plane_dist.z = min( orbit_plane_dist.z, distPointPlane(Z, orbit_plane.yyx, orbit_plane_origin) );
            orbit_coord_dist = min( orbit_coord_dist, fabs(dot(Z, P_in)) );
            orbit_sphere_dist = min( orbit_sphere_dist, length2(Z - normalize(Z)*orbit_sphere_rad) );
            orbit_axis_dist.x = min(orbit_axis_dist.x, Z.y*Z.y + Z.z*Z.z);
            orbit_axis_dist.y = min(orbit_axis_dist.y, Z.x*Z.x + Z.z*Z.z);
            orbit_axis_dist.z = min(orbit_axis_dist.z, Z.x*Z.x + Z.y*Z.y);
        }
    }
    distance = length(Z);

    // delta DE method
    // based on Makin/Buddhi 4-point Delta-DE formula
    #if ENABLE_DELTA_DE
        float delta = 0.000005f;

        Z = P_in + (float3)(delta, 0.0f, 0.0f);
        #pragma unroll
        for (int j=0; j<i; j++)
        {
            fractal_stack(&Z, &de, &P_in, &log_lin, stack);
        }
        float rx = length(Z);
        float drx = (distance - rx) / delta;

        Z = P_in + (float3)(0.0f, delta, 0.0f);

        #pragma unroll
        for (int j=0; j<i; j++)
        {
            fractal_stack(&Z, &de, &P_in, &log_lin, stack);
        }
        float ry = length(Z);
        float dry = (distance - ry) / delta;

        Z = P_in + (float3)(0.0f, 0.0f, delta);

        #pragma unroll
        for (int j=0; j<i; j++)
        {
            fractal_stack(&Z, &de, &P_in, &log_lin, stack);
        }
        float rz = length(Z);
        float drz = (distance - rz) / delta;

        float3 dist_grad = (float3)(drx, dry, drz);
        de = length(dist_grad);

        dist_grad = normalize(dist_grad);
    #endif

    // automatic determining DE mode based on log_lin value
    if (log_lin >= 0) out_distance = 0.5f * log(distance) * distance/de;
    else out_distance = distance / fabs(de);

    // outputting orbit traps and N in case of DELTA DE mode
    if (final == 1)
    {
        orbit_colors[0] = sqrt(orbit_pt_dist); // distance to point at specified coordinates
        orbit_colors[1] = orbit_plane_dist.x; // distance to YZ plane
        orbit_colors[2] = orbit_plane_dist.y; // distance to XZ plane
        orbit_colors[3] = orbit_plane_dist.z; // distance to XY plane
        orbit_colors[4] = orbit_coord_dist; // dot(Z, world coords)
        orbit_colors[5] = sqrt(orbit_sphere_dist); // distance to sphere
        orbit_colors[6] = sqrt(orbit_axis_dist.x); // distance to X axis
        orbit_colors[7] = sqrt(orbit_axis_dist.y); // distance to Y axis
        orbit_colors[8] = sqrt(orbit_axis_dist.z); // distance to Z axis

        // delta DE mode does not work well with multiple fractal objects    
        #if ENABLE_DELTA_DE
            *N = dist_grad;
        #endif
    }

    return out_distance * size;
}