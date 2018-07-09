#ifndef VFT_FRACTALS
#define VFT_FRACTALS

// mapping of variables
// aux.r_dz -> de
// aux.r -> distance
// aux.de -> de
// dr -> de
// r -> distance
// Bailout -> max_distance
// Iterations -> max_iterations
// positive log_lin -> log, negative -> lin

// fractal formula sources:
// [M2] - Mandelbulber v2
// [M3D] - Mandelbulb 3D
// [WEB] - From the web

////////////// primitives
// [WEB] - http://iquilezles.org/www/articles/distfunctions/distfunctions.htm

// sphere: position, radius
static float sphere( float3 P, float rad )
{
    float dist = length(P) - rad;
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


////////////// fractals

// [WEB] - http://blog.hvidtfeldts.net/index.php/2011/09/
static void mandelbulbIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float power)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    // convert to polar coordinates
    float theta = acos((*Z).z/distance);
    float phi = atan2((*Z).y, (*Z).x);

    *de =  pow(distance, power-1) * power * (*de) + 1.0f;
    
    // scale and rotate the point
    float zr = pow(distance, power);
    theta *= power;
    phi *= power;
    
    // convert back to cartesian coordinates
    float3 newZ = zr * (float3)( sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta) );

    if (julia.x == 0.0f)
    {
        *Z = newZ + *P_in;
    }
    else {
        *Z = newZ + julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [WEB] - http://www.fractalforums.com/index.php?topic=2785.msg14893#msg14893
static void mandelboxIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float scale)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float fixedRadius = 1.0f;
    float fR2 = fixedRadius * fixedRadius;
    float minRadius = 0.5f;
    float mR2 = minRadius * minRadius;

    if ((*Z).x > 1.0f) (*Z).x = 2.0f - (*Z).x;
    else if ((*Z).x < -1.0f) (*Z).x = -2.0f - (*Z).x;

    if ((*Z).y > 1.0f) (*Z).y = 2.0f - (*Z).y;
    else if ((*Z).y < -1.0f) (*Z).y = -2.0f - (*Z).y;

    if ((*Z).z > 1.0f) (*Z).z = 2.0f - (*Z).z;
    else if ((*Z).z < -1.0f) (*Z).z = -2.0f - (*Z).z;

    float r2 = (*Z).x*(*Z).x + (*Z).y*(*Z).y + (*Z).z*(*Z).z;

    if (r2 < mR2)
    {
        *Z = *Z * fR2 / mR2;
        *de = *de * fR2 / mR2;
    }
    else if (r2 < fR2)
    {
        *Z = *Z * fR2 / r2;
        *de *= fR2 / r2;
    }

    *de *= scale;

    if (julia.x == 0.0f)
    {
        *Z = *Z * scale + *P_in;
    }
    else
    {
        *Z = *Z * scale + julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Classic Mandelbulb Power 2 fractal - MandelbulbPower2Iteration
static void mandelbulbPower2Iter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia) 
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    *de = *de * 2.0f * distance;
    float x2 = (*Z).x * (*Z).x;
    float y2 = (*Z).y * (*Z).y;
    float z2 = (*Z).z * (*Z).z;
    float temp = 1.0f - z2 / (x2 + y2);
    float3 new;
    new.x = (x2 - y2) * temp;
    new.y = 2.0f * (*Z).x * (*Z).y * temp;
    new.z = -2.0f * (*Z).z * sqrt(x2 + y2);

    if (julia.x == 0.0f)
    {
        *Z = new + *P_in;
    }
    else
    {
        *Z = new + julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Menger Sponge formula created by Knighty, modulus modification by mancoast - MengerSpongeIteration
static void mengerSpongeIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const int modulus)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    (*Z).x = fabs((*Z).x);
    (*Z).y = fabs((*Z).y);
    (*Z).z = fabs((*Z).z);

    if ((*Z).x - (*Z).y < 0.0f) (*Z).xy = (*Z).yx;
    if ((*Z).x - (*Z).z < 0.0f) (*Z).xz = (*Z).zx;
    if ((*Z).y - (*Z).z < 0.0f) (*Z).yz = (*Z).zy;

    *Z *= 3.0f;

    (*Z).x -= 2.0f;
    (*Z).y -= 2.0f;
    
    if (modulus == 0)
    {
        if ((*Z).z > 1.0f) (*Z).z -= 2.0f;
    }
    else
    {
        if (fmod((*Z).z, M_PI_F) > 2.0f) (*Z).z -= 2.0f;
    }

    *de *= 3.0f;

    if (julia.x == 1.0f)
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Bristorbrot formula - BristorbrotIteration
static void bristorbrotIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    float3 new;
    new.x = (*Z).x * (*Z).x - (*Z).y * (*Z).y - (*Z).z * (*Z).z;
    new.y = (*Z).y * (2.0f * (*Z).x - (*Z).z);
    new.z = (*Z).z * (2.0f * (*Z).x + (*Z).y);
    
    *de = *de * 2.0f * distance;
    *Z = new;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Xenodreambuie - XenodreambuieIteration
static void xenodreambuieIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float power, float alpha, float beta)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    alpha = radians(alpha);
    beta = radians(beta);

    float rp = pow(distance, power - 1.0f);
    *de = rp * (*de) * power + 1.0f;
    rp *= distance;

    float th = atan2((*Z).y, (*Z).x) + beta;
    float ph = acos((*Z).z / distance) + alpha;

    if (fabs(ph) > 0.5f * M_PI_F) ph = sign(ph) * M_PI_F - ph;

    (*Z).x = rp * cos(th * power) * sin(ph * power);
    (*Z).y = rp * sin(th * power) * sin(ph * power);
    (*Z).z = rp * cos(ph * power);

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Sierpinski3D. made from Darkbeam's Sierpinski code from M3D - Sierpinski3dIteration
static void sierpinski3dIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float scale, const float3 offset, const float3 rot)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float3 temp = *Z;

    if ((*Z).x - (*Z).y < 0.0f) (*Z).xy = (*Z).yx;
    if ((*Z).x - (*Z).z < 0.0f) (*Z).xz = (*Z).zx;
    if ((*Z).y - (*Z).z < 0.0f) (*Z).yz = (*Z).zy;
    if ((*Z).x + (*Z).y < 0.0f)
    {
        temp.x = -(*Z).y;
        (*Z).y = -(*Z).x;
        (*Z).x = temp.x;
    }
    if ((*Z).x + (*Z).z < 0.0f)
    {
        temp.x = -(*Z).z;
        (*Z).z = -(*Z).x;
        (*Z).x = temp.x;
    }
    if ((*Z).y + (*Z).z < 0.0f)
    {
        temp.y = -(*Z).z;
        (*Z).z = -(*Z).y;
        (*Z).y = temp.y;
    }

    *Z *= scale;
    *de *= scale;

    *Z -= offset;

    *Z = mtxPtMult( mtxRotate(rot) , *Z );

    if (julia.x == 1.0f)
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Menger Smooth - mengerSmoothIteration
static void mengerSmoothIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float scale, const float offset_s, const float3 offset_c, const float3 rot)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    	
    float sc1 = scale - 1.0f;
	float sc2 = sc1 / scale;

    *Z = (float3)(sqrt((*Z).x * (*Z).x + offset_s), sqrt((*Z).y * (*Z).y + offset_s), sqrt((*Z).z * (*Z).z + offset_s));    

	float t;

	t = (*Z).x - (*Z).y;
	t = 0.5f * (t - sqrt(t * t + offset_s));
	(*Z).x = (*Z).x - t;
	(*Z).y = (*Z).y + t;

	t = (*Z).x - (*Z).z;
	t = 0.5f * (t - sqrt(t * t + offset_s));
	(*Z).x = (*Z).x - t;
	(*Z).z = (*Z).z + t;

	t = (*Z).y - (*Z).z;
	t = 0.5f * (t - sqrt(t * t + offset_s));
	(*Z).y = (*Z).y - t;
	(*Z).z = (*Z).z + t;

	(*Z).z = (*Z).z - offset_c.z * sc2;
	(*Z).z = -sqrt((*Z).z * (*Z).z + offset_s);
	(*Z).z = (*Z).z + offset_c.z * sc2;

	(*Z).x = scale * (*Z).x - offset_c.x * sc1;
	(*Z).y = scale * (*Z).y - offset_c.y * sc1;
	(*Z).z = scale * (*Z).z;

	*de *= scale;

    *Z = mtxPtMult( mtxRotate(rot) , *Z );

    if (julia.x == 1.0f)
    {
        *Z += *P_in;
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Amazing Surf from Mandelbulber3D, formula proposed by Kali, with features added by Darkbeam - AmazingSurfIteration
static void amazingSurfIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float fold_x, const float fold_y, const int force_cylindrical_fold, const float min_radius, float scale, const float scale_fold_influence, const float3 rot, const int multiply_c, const float3 c_multiplier)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float actual_scale = scale;

    (*Z).x = fabs((*Z).x + fold_x) - fabs((*Z).x - fold_x) - (*Z).x;
    (*Z).y = fabs((*Z).y + fold_y) - fabs((*Z).y - fold_y) - (*Z).y;

    float rr = dot(*Z, *Z);

    // get rid of if afterwards
    if (force_cylindrical_fold == 1)
    {
        rr -= (*Z).z * (*Z).z;
    }

	float sqrtMinR = sqrt(min_radius);
	float dividend = rr < sqrtMinR ? sqrtMinR : min(rr, 1.0f);

    float m = actual_scale / dividend;

	*Z *= (m - 1.0f) * scale_fold_influence + 1.0f;
	*de = *de * fabs(m) + 1.0f;

    if (multiply_c == 1)
    {
        *Z += (float3)((*P_in).y, (*P_in).x, (*P_in).z) * c_multiplier;
    }
    
    *Z = mtxPtMult( mtxRotate(rot) , *Z );

    if (julia.x == 0.0f)
    {
        *Z += (float3)((*P_in).y, (*P_in).x, (*P_in).z);
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Benesi formula invented by Benesi - BenesiIteration
static void benesiIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    *de = *de * 2.0f * distance;
    float r1 = (*Z).y * (*Z).y + (*Z).z * (*Z).z;
    float newx;
	if ((*P_in).x < 0.0f || (*Z).x < sqrt(r1))
	{
		newx = (*Z).x * (*Z).x - r1;
	}
	else
	{
		newx = -(*Z).x * (*Z).x + r1;
	}
	r1 = -1.0f / sqrt(r1) * 2.0f * fabs((*Z).x);
	float newy = r1 * ((*Z).y * (*Z).y - (*Z).z * (*Z).z);
	float newz = r1 * 2.0f * (*Z).y * (*Z).z;

	*Z = (float3)(newx, newy, newz);

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Mandelbulb 2 fractal formula created by Buddhi - Mandelbulb2Iteration
static void mandelbulb2Iter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

	*de = *de * 2.0f * distance;

	float tempR = sqrt((*Z).x * (*Z).x + (*Z).y * (*Z).y); //+ 1e-061
	*Z *= 1.0f / tempR;
	float temp = (*Z).x * (*Z).x - (*Z).y * (*Z).y;
	(*Z).y = 2.0f * (*Z).x * (*Z).y;
	(*Z).x = temp;
	*Z *= tempR;

	tempR = sqrt((*Z).y * (*Z).y + (*Z).z * (*Z).z); //+ 1e-061
	*Z *= 1.0f / tempR;
	temp = (*Z).y * (*Z).y - (*Z).z * (*Z).z;
	(*Z).z = 2.0f * (*Z).y * (*Z).z;
	(*Z).y = temp;
	*Z *= tempR;

	tempR = sqrt((*Z).x * (*Z).x + (*Z).z * (*Z).z); //+ 1e-061
	*Z *= 1.0f / tempR;
	temp = (*Z).x * (*Z).x - (*Z).z * (*Z).z;
	(*Z).z = 2.0f * (*Z).x * (*Z).z;
	(*Z).x = temp;
	*Z *= tempR;

	(*Z) *= distance;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)--;
}

// [M2] - Mandelbulb 3 fractal formula created by Buddhi - Mandelbulb3Iteration
static void mandelbulb3Iter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

	
    *de = *de * 2.0f * distance;

	float temp, tempR;

	float sign = 1.0f;
	float sign2 = 1.0f;

	if ((*Z).x < 0.0f) sign2 = -1.0f;
	tempR = sqrt((*Z).x * (*Z).x + (*Z).y * (*Z).y); //+ 1e-061
	*Z *= 1.0f / tempR;
	temp = (*Z).x * (*Z).x - (*Z).y * (*Z).y;
	(*Z).y = 2.0f * (*Z).x * (*Z).y;
	(*Z).x = temp;
	*Z *= tempR;

	if ((*Z).x < 0.0f) sign = -1.0f;
	tempR = sqrt((*Z).x * (*Z).x + (*Z).z * (*Z).z); //+ 1e-061
	*Z *= 1.0f / tempR;
	temp = (*Z).x * (*Z).x - (*Z).z * (*Z).z;
	(*Z).z = 2.0f * (*Z).x * (*Z).z * sign2;
	(*Z).x = temp * sign;
	*Z *= tempR;

	*Z *= distance;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Mandelbulb 4 fractal formula created by Buddhi - Mandelbulb4Iteration
static void mandelbulb4Iter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float power, const float3 angles)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

    float rp = pow(distance, power - 1.0f);
	*de = rp * (*de) * power + 1.0f;

	float angZ = degrees(atan2((*Z).y, (*Z).x)) + angles.x;
	float angY = degrees(atan2((*Z).z, (*Z).x)) + angles.y;
	float angX = degrees(atan2((*Z).z, (*Z).y)) + angles.z;

    float16 rotM = mtxRotate( (float3)((angX * (power - 1.0f)), (angY * (power - 1.0f)), (angZ * (power - 1.0f))) );

	*Z = mtxPtMult( rotM , *Z ) * rp;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Ides formula made by Trafassel, the original Ide's Formula thread - IdesIteration
static void idesIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float3 multiplier, const float2 sub_multiplier)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
	if (fabs((*Z).x) < 2.5f) (*Z).x = (*Z).x * 0.9f;
	if (fabs((*Z).z) < 2.5f) (*Z).z = (*Z).z * 0.9f;

	float3 z2 = (*Z) * (*Z);
	float3 newZ;
	newZ.x = multiplier.x * z2.x - sub_multiplier.x * (z2.y + z2.z);
	newZ.y = multiplier.y * (*Z).x * (*Z).y * (*Z).z;
	newZ.z = multiplier.z * z2.z - sub_multiplier.y * (z2.x + z2.y);
	*Z = newZ;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
}


// [M2] - Ides 2 formula made by Trafassel, the original Ide's Formula thread - Ides2Iteration
static void ides2Iter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float3 multiplier, const float2 sub_multiplier)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
	float3 z2 = (*Z) * (*Z);
	float3 newZ;
	newZ.x = multiplier.x * z2.x - sub_multiplier.x * (z2.y + z2.z);
	newZ.y = multiplier.y * (*Z).x * (*Z).y * (*Z).z;
	newZ.z = multiplier.z * z2.z - sub_multiplier.y * (z2.x + z2.y);

	*Z = newZ + (*Z);

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
}

// [M2] - IQ Bulb from Mandelbulb 3D and Inigo Quilez - IqBulbIteration
static void iqBulbIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float power, const float zpower)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

	// extract polar coordinates
	float wr = distance;
	float wo = acos((*Z).y / wr);
	float wi = atan2((*Z).x, (*Z).z);

	// scale and rotate the point
	wr = pow(wr, power - 1.0f);
	*de = wr * *de * power + 1.0f;
	wr *= distance;
	wo *= power;
	wi *= zpower;

	// convert back to cartesian coordinates
	(*Z).x = sin(wo) * sin(wi);
	(*Z).y = cos(wo);
	(*Z).z = sin(wo) * cos(wi);

	*Z *= wr; // then add Cpixel constant

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - Quaternion3DE fractal with extended controls - Quaternion3dIteration
static void quaternion3dIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float3 scale, const float3 offset, const float3 rot, const float de_influence)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);

	*de = (*de) * 2.0f * distance;
	*Z = (float3)((*Z).x * (*Z).x - (*Z).y * (*Z).y - (*Z).z * (*Z).z, (*Z).x * (*Z).y, (*Z).x * (*Z).z);

	float tempL = length(*Z);
	*Z *= scale;
	float3 tempAvgScale = (float3)((*Z).x, (*Z).y / 2.0f, (*Z).z / 2.0f);
	float avgScale = length(tempAvgScale) / tempL;
	float tempAux = *de * avgScale;
	*de = *de + (tempAux - *de) * de_influence;

    *Z = mtxPtMult( mtxRotate(rot) , *Z );

	*Z += offset;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - JosLeys-Kleinian - JosKleinianIteration
// will not work in this case, requires different DE computation
static void josKleinianIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia, const float r, const float l, const float3 box_size)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
	float a = r;
	float b = l;
	float f = sign(b);

	float3 box1 = (float3)(2.0f * box_size.x, a * box_size.y, 2.0 * box_size.z);
	float3 box2 = (float3)(-box_size.x, -box_size.y + 1.0f, -box_size.z);
	float3 wrapped = wrap(*Z, box1, box2);

    *Z = (float3)(wrapped.x, wrapped.y, wrapped.z);

    if ((*Z).y >= a * (0.5f + 0.2f * sin(f * M_PI_F * ((*Z).x + b * 0.5f) / box_size.x)))
		*Z = (float3)(-b, a, 0.0f) - *Z;

    float z2 = dot(*Z, *Z);

    float iR = 1.0f / z2;
	*Z *= -iR;
	(*Z).x = -b - (*Z).x;
	(*Z).y = a + (*Z).y;
	*de *= iR;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}

// [M2] - T>rotation - TransfRotationIteration
static void rotationIter(float3* Z, const float3 rot)
{
    *Z = mtxPtMult( mtxRotate(rot) , *Z );
}

// [M2] - T>Box Fold - TransfBoxFoldIteration
static void boxFoldIter(float3* Z, const int3 axis_enable, const float folding_limit, const float folding_value, const float z_scale)
{
    if (axis_enable.x == 1)
    {
        if (fabs((*Z).x) > folding_limit)
        {
            (*Z).x = sign((*Z).x) * folding_value - (*Z).x;
        }
    }
    if (axis_enable.y == 1)
    {
        if (fabs((*Z).y) > folding_limit)
        {
            (*Z).y = sign((*Z).y) * folding_value - (*Z).y;
        }
    }
    if (axis_enable.z == 1)
    {
        float zLimit = folding_limit * z_scale;
        float zValue = folding_value * z_scale;
        if (fabs((*Z).z) > zLimit)
        {
            (*Z).z = sign((*Z).z) * zValue - (*Z).z;
        }
    }
}

// [M2] - Rotation Folding: rotatedAbs & Rotated Folding transform from M3D - TransfRotationFoldingIteration
static void fabsFoldIter(float3* Z, const int3 axis_enable, const float3 offset)
{
    if (axis_enable.x == 1)
    {
        (*Z).x = fabs((*Z).x + offset.x) - offset.x;
    }
    if (axis_enable.y == 1)
    {
        (*Z).y = fabs((*Z).y + offset.y) - offset.y;
    }
    if (axis_enable.z == 1)
    {
        (*Z).z = fabs((*Z).z + offset.z) - offset.z;
    }
}


// [M2] - Rotation Folding: rotatedAbs & Rotated Folding transform from M3D - TransfRotationFoldingIteration
static void tgladFoldIter(float3* Z, const int3 axis_enable, const float3 offset)
{
    if (axis_enable.x == 1)
    {
        (*Z).x = fabs((*Z).x + offset.x) - fabs((*Z).x - offset.x) - (*Z).x;
    }
    if (axis_enable.y == 1)
    {
        (*Z).y = fabs((*Z).y + offset.y) - fabs((*Z).y - offset.y) - (*Z).y;
    }
    if (axis_enable.z == 1)
    {
        (*Z).z = fabs((*Z).z + offset.z) - fabs((*Z).z - offset.z) - (*Z).z;
    }
}

// [M2] - T>Scale - TransfScaleIteration
static void scaleIter(float3* Z, float* de, const float3 scale)
{
    float distance_init = length(*Z);

    *Z = mtxPtMult( mtxScale(scale) , *Z );
    
    float distance_scaled = length(*Z);
    *de *= distance_scaled / distance_init;
}

static void translateIter(float3* Z, const float3 translate)
{
    *Z += translate;
}

static void addCOffsetIter(float3* Z, const float3* P_in, const float3 offset)
{
    *Z += *P_in + offset;
}

/*
// testing new formulas


// [M2] - Hypercomplex 3D Mandelbrot formula invented by David Makin - HypercomplexIteration
static void hypercomplexIter(float3* Z, float* de, const float3* P_in, int* log_lin, const float weight, const float4 julia)
{
    float3 Z_orig = *Z;
    float de_orig = *de;
    
    float distance = length(*Z);
    float4 Z_complex = (float4)((*Z).x, (*Z).y, (*Z).z, 1.0f);
    float4* Z_ = &Z_complex;

	*de = *de * 2.0f * distance;
	float newx = (*Z_).x * (*Z_).x - (*Z_).y * (*Z_).y - (*Z_).z * (*Z_).z - (*Z_).w * (*Z_).w;
	float newy = 2.0f * (*Z_).x * (*Z_).y - 2.0f * (*Z_).w * (*Z_).z;
	float newz = 2.0f * (*Z_).x * (*Z_).z - 2.0f * (*Z_).y * (*Z_).w;
	float neww = 2.0f * (*Z_).x * (*Z_).w - 2.0f * (*Z_).y * (*Z_).z;
	//(*Z).x = newx;
	//(*Z).y = newy;
	//(*Z).z = newz;
	//(*Z).w = neww;
    
    (*Z).x = newx * newx - newy * newy - newz * newz - neww * neww;
    (*Z).y = newx * newy;
    (*Z).z = newx * newz;

    if (julia.x == 0.0f)
    {
        *Z += *P_in;
    }
    else 
    {
        *Z += julia.yzw;
    }

    *Z = mix(Z_orig, *Z, weight);
    *de = mix(de_orig, *de, weight);
    (*log_lin)++;
}



// quaternion fractals, kind of works, but not sure what to do with z.w component :)
static float quaternion(float3 P, float size)
{
//void QuaternionIteration(CVector4 &z, const sFractal *fractal, sExtendedAux &aux)
//{
//    Q_UNUSED(fractal);

//    aux.r_dz = aux.r_dz * 2.0 * aux.r;
//    double newx = z.x * z.x - z.y * z.y - z.z * z.z - z.w * z.w;
//    double newy = 2.0 * z.x * z.y;
//    double newz = 2.0 * z.x * z.z;
//    double neww = 2.0 * z.x * z.w;
//    z.x = newx;
//    z.y = newy;
//    z.z = newz;
//    z.w = neww;
//}
    P /= size;

    float4 z = (float4)(P, 1);
    float dr = 1.0;
    float r = 0.0;
    int Iterations = 50;
    int Bailout = 40;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;

        dr = dr * 2.0f * r;
        float newx = z.x * z.x - z.y * z.y - z.z * z.z - z.w * z.w;
        float newy = 2.0f * z.x * z.y;
        float newz = 2.0f * z.x * z.z;
        float neww = 2.0f * z.x * z.w;
        z.x = newx;
        z.y = newy;
        z.z = newz;
        //z.w = neww;
    }

    float out = 0.5f * log(r) * r/dr;
    return out * size;
}

// quaternion3d
// kind of works, but does not exactly match M2 visual, but parameters deform it in a similar manner, I hardcoded some of the parameters into values, there are some buggy areas, missing parts, noisy normals etc.
static float quaternion3d(float3 P, float size)
{
//void Quaternion3dIteration(CVector4 &z, const sFractal *fractal, sExtendedAux &aux)
//{
//
//    aux.r_dz = aux.r_dz * 2.0 * aux.r;
//    z = CVector4(z.x * z.x - z.y * z.y - z.z * z.z, z.x * z.y, z.x * z.z, z.w);
//
//    double tempL = z.Length();
//    z *= fractal).transformCommon.constantMultiplier122;
//    // if (tempL < 1e-21) tempL = 1e-21;
//    CVector4 tempAvgScale = CVector4(z.x, z.y / 2.0, z.z / 2.0, z.w);
//    double avgScale = tempAvgScale.Length() / tempL;
//    double tempAux = aux.r_dz * avgScale;
//    aux.r_dz = aux.r_dz + (tempAux - aux.r_dz) * fractal).transformCommon.scaleA1;
//
//    if (fractal).transformCommon.rotationEnabled)
//        z = fractal).transformCommon.rotationMatrix.RotateVector(z);
//
//    z += fractal).transformCommon.additionConstant000;
//}
    P /= size;

    float4 z = (float4)(P, 1);
    float dr = 1.0;
    float r = 0.0;
    int Iterations = 250;
    int Bailout = 20;

    for (int i = 0; i < Iterations ; i++)
    {
        r = length(z);
        if (r > Bailout) break;
        
        dr = dr * 2.0 * r;
        z = (float4)(z.x * z.x - z.y * z.y - z.z * z.z, z.x * z.y, z.x * z.z, z.w);

        float tempL = r;
        z *= (float4)(1,1,1,1);

        float4 tempAvgScale = (float4)(z.x, z.y / 2.0, z.z / 2.0, z.w);
        float avgScale = length(tempAvgScale) / tempL;
        float tempAux = dr * avgScale;
        dr = dr + (tempAux - dr) * 1.0f;

        //if (fractal).transformCommon.rotationEnabled)
        //    z = fractal).transformCommon.rotationMatrix.RotateVector(z);

        z += (float4)(0,0,0,0);
    }

    //float out = 0.5f * log(r) * r/dr;
    float out = r / dr;
    return out * size;
}

*/
#endif