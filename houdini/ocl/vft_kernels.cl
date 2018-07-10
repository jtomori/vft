#include "vft_utils.h"
#include "vft_math.h"
#include "vft_fractals.h"
#include "vft_defines.h"
#include "vft_shading.h"

// contains primitives
float primitive_stack(float3 P, const int stack)
{
    float out_distance;
    switch (stack)
    {
        case 0:
        {
            out_distance = sphere(P, 1.0f);

            break;
        }
        case 1:
        {
            out_distance = torus(P, (float2)(1.0f, 0.5f));

            break;
        }
    }

    return out_distance;
}

// contains fractal combinations for all shapes
void fractal_stack(float3* Z, float* de, const float3* P_in, int* log_lin, const int stack)
{
    switch (stack)
    {
        case 0:
        {
            iqBulbIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.0f), 8.0f, 8.0f);
            //idesIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.0f), (float3)(1.0f, 2.0f, 1.0f), (float2)(0.5f, 0.5f));
            //benesiIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.0f));
            //hypercomplexIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.0f));            
            //josKleinianIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 1.0f, 1.0f, 3.0f), 2.0f, 0.0f, (float3)(1.0f, 1.0f, 1.0f));            
            //amazingSurfIter(Z, de, P_in, log_lin, 1.0f, (float4)(1.0f, 0.0f, 0.0f, 0.0f), 1.0f, 1.0f, 0, 0.5f, 2.0f, 1.0f, (float3)(0.0f, 0.0f, 0.0f), 1, (float3)(1.0f, 1.0f, 1.0f));
            //mandelbulbPower2Iter(Z, de, P_in, log_lin, 1.0f, (float4)(1.0f, 0.3f, 0.5f, 0.2f)); // log
            //bristorbrotIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 1.3f, 3.3f, 0.0f)); // log
            //xenodreambuieIter(Z, de, P_in, log_lin, 1.0f, (float4)(1.0f, 1.0f, 0.0f, 0.0f), 9.0f, 0.0f, 0.0f); // log
            //mandelboxIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 1.0f, 3.0f, 4.0f), 3.1f); // lin
            //mandelbulbIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 1.0f, 0.0f, 0.0f), 8.0f); // log
            //mengerSpongeIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 1.0f, 0.0f)); // lin
            //sierpinski3dIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.5f), 2.0f, (float3)(1.0f, 1.0f, 1.0f), (float3)(0.0f, 0.0f, 0.0f) ); // lin

            break;
        }

        case 1:
        {
            sierpinski3dIter(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f, 0.0f, 0.0f, 0.5f), 2.0f, (float3)(1.0f, 1.0f, 1.0f), (float3)(0.0f, 0.0f, 0.0f) ); // lin

            break;
        }

        case 2:
        {
            mandelbulbPower2Iter(Z, de, P_in, log_lin, 1.0f, (float4)(1.0f, 0.3f, 0.5f, 0.2f)); // log

            break;
        }
        case 3:
        {
            #define PY_FRACTAL_STACK

            break;
        }
    }
}

// scene setup - setting of coordinates and shapes in them

// scene with multiple unions of prims and fractals
/*float scene( float3 P, const int final, float* orbit_colors, float3* N ) {
    float dist_out;
    float orbit_closest = LARGE_NUMBER;

    float shape1 = hybrid(P,                  10, 10.0f, 1.0f, final, &orbit_closest, orbit_colors, N, 0);
    float shape2 = hybrid(P - (float3)(2.0f), 10, 10.0f, 1.0f, final, &orbit_closest, orbit_colors, N, 1);
    float shape3 = hybrid(P + (float3)(2.0f), 10, 10.0f, 1.0f, final, &orbit_closest, orbit_colors, N, 2);

    float shape4 = primitive(P - (float3)(2.0f, 0.0f, 0.0f), 0.3f,      final, &orbit_closest, orbit_colors, (float3)(1.0f,0.0f,1.0f),  N, 0);
    float shape5 = primitive(P - (float3)(-2.0f, 0.0f, 0.0f), 0.6f,     final, &orbit_closest, orbit_colors, (float3)(1.0f,1.0f,0.0f), N, 1);

    dist_out = sdfUnion( sdfUnion( sdfUnion( sdfUnion(shape1, shape2) , shape3 ) , shape4 ) , shape5 );

    return dist_out;
}*/

// scene showing unions and subtractions
/*float scene( float3 P, const int final, float* orbit_colors, float3* N ) {
    float dist_out;
    float orbit_closest = LARGE_NUMBER;

    float shape1 = hybrid(P,                  10, 10.0f, 1.0f, final, &orbit_closest, orbit_colors, N, 0);
    float shape2 = hybrid(P - (float3)(0.4f), 10, 10.0f, 1.0f, final, &orbit_closest, orbit_colors, N, 1);
    //float shape3 = hybrid(P - (float3)(-0.5f, 0.0f, 0.5f), 10, 10.0f, 1.0f, final, &orbit_closest, orbit_colors, N, 2);

    float shape4 = primitive(P - (float3)(0.7f, 0.0f, 0.0f), 1.0f,      final, &orbit_closest, orbit_colors, (float3)(1.0f,0.0f,1.0f),  N, 0);

    dist_out = sdfSubtract( sdfUnion(shape1, shape2) , shape4 );

    return dist_out;
}*/

float scene( float3 P, const int final, float* orbit_colors, float3* N ) {
    float dist_out;
    float orbit_closest = LARGE_NUMBER;

    float shape1 = hybrid(P, 25, 10.0f, 1.0f, final, &orbit_closest, orbit_colors, N, 3);

    dist_out = shape1;

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
    float16 near_plane_xform_scaled = mtxMult(near_plane_xform, mtxScale( (float3)(10000000.0f) ) );
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
    float AO = 1.0f;
    float orbit_colors[ORBITS_ARRAY_LENGTH];
    float3 Cd_out = (float3)(1.0f);
    float3 N_grad;

    float3 ray_P_world = pixel_P_world;
    float cam_dist = scene(cam_P_world, 0, NULL, NULL);
    float de = 0.0f;
    int i = 0;

    // quality settings
    float step_size = 0.4f;
    float iso_limit_mult = 0.5f;
    float ray_dist = planeZ[0];
    const int max_steps = 300;
    const float max_dist = 1000.0f;

    float iso_limit = cam_dist * 0.0001f * iso_limit_mult;  

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

        // output shading for viewport preview
        color.x = AO;
        color.y = orbit_colors[0];
        color.z = 1.0f;

        Cd_out = color;
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

kernel void computeSdf( 
    int surface_stride_x, 
    int surface_stride_y, 
    int surface_stride_z, 
    int surface_stride_offset, 
    float16 surface_xformtoworld, 
    global float * surface
    )
{
    int gidx = get_global_id(0);
    int gidy = get_global_id(1);
    int gidz = get_global_id(2);
    int idx = surface_stride_offset + surface_stride_x * gidx + surface_stride_y * gidy + surface_stride_z * gidz;

    float3 P_vol = (float3)(gidx, gidy, gidz);
    float3 P_world = mtxPtMult(surface_xformtoworld, P_vol);

    float de = 0.0f;
    de = scene(P_world, 0, NULL, NULL);
    vstore1(de, idx, surface);
}


kernel void computeSdfColors( 
    int color_0_stride_x, 
    int color_0_stride_y, 
    int color_0_stride_z, 
    int color_0_stride_offset, 
    float16 color_0_xformtoworld, 
    global float * color_0,
    global float * color_1,
    global float * color_2,
    global float * color_3,
    global float * color_4,
    global float * color_5,
    global float * color_6,
    global float * color_7,
    global float * color_8
    )
{
    int gidx = get_global_id(0);
    int gidy = get_global_id(1);
    int gidz = get_global_id(2);
    int idx = color_0_stride_offset + color_0_stride_x * gidx + color_0_stride_y * gidy + color_0_stride_z * gidz;

    float3 P_vol = (float3)(gidx, gidy, gidz);
    float3 P_world = mtxPtMult(color_0_xformtoworld, P_vol);

    float orbit_colors[ORBITS_ARRAY_LENGTH];    
    float de = 0.0f;

    de = scene(P_world, 1, orbit_colors, NULL);
    
    vstore1(orbit_colors[0], idx, color_0);
    vstore1(orbit_colors[1], idx, color_1);
    vstore1(orbit_colors[2], idx, color_2);
    vstore1(orbit_colors[3], idx, color_3);
    vstore1(orbit_colors[4], idx, color_4);
    vstore1(orbit_colors[5], idx, color_5);
    vstore1(orbit_colors[6], idx, color_6);
    vstore1(orbit_colors[7], idx, color_7);
    vstore1(orbit_colors[8], idx, color_8);
}
