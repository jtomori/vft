#ifndef VFT_FRACTALS_OSL
#define VFT_FRACTALS_OSL

// Porting OpenCL fractal functions from vft_fractals.h
// it means removing 
//                  " keyword, 
//                  "f" literal, 
//                  "const" keyword, 
//                  fixing expressions with vector swizzling, 
//                  removing pointers, 
//                  removing casts
//                  removing de, log_lin variables

// sphere: position, radius
float sphere( float3 P, float rad )
{
    float dist = LENGTH(P) - rad;
    return dist;
}

// box: position, size
float box( float3 P, float3 b )
{
  float3 d = fabs(P) - b;
  return min( max( d.x, max(d.y, d.z) ), 0.0) + LENGTH( max(d, 0.0) );
}

// round box: position, size, roundness
float roundBox( float3 P, float3 b, float r )
{
  return LENGTH( max( fabs(P) - b, 0.0 ) )-r;
}

// torus: position, size (radius, thickness)
float torus( float3 P, float2 t )
{
  float2 q = float2(LENGTH( float2(P.x, P.z) )-t.x,P.y);
  return LENGTH(q)-t.y;
}

// infinite cone
float cone( float3 P, float2 c )
{
    c = NORMALIZE(c);
    float q = LENGTH(float2(P.x, P.y));
    return dot( c, float2(q, P.z) );
}

// [WEB] - http://blog.hvidtfeldts.net/index.php/2011/09/
void mandelbulbIter(point Z, point P_in, float weight, float4 julia, float power)
{
    point Z_orig = Z;
    
    float distance = LENGTH(Z);

    // convert to polar coordinates
    float theta = acos( DIV(Z[2], distance));
    float phi = atan2(Z[1], Z[0]);
    
    // scale and rotate the point
    float zr = POWR(distance, power);
    theta *= power;
    phi *= power;
    
    // convert back to cartesian coordinates
    point new_pZ = zr * point( SIN(theta)*COS(phi), SIN(phi)*SIN(theta), COS(theta) );

    if (julia.x == 0.0)
    {
        Z = new_pZ + P_in;
    }
    else 
    {
        Z = new_pZ + point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

// [WEB] - http://www.fractalforums.com/index.php?topic=2785.msg14893#msg14893
void mandelboxIter(point Z, point P_in, float weight, float4 julia, float scale)
{
    point Z_orig = Z;
    
    float fixedRadius = 1.0;
    float fR2 = fixedRadius * fixedRadius;
    float minRadius = 0.5;
    float mR2 = minRadius * minRadius;

    if (Z[0] > 1.0) Z[0] = 2.0 - Z[0];
    else if (Z[0] < -1.0) Z[0] = -2.0 - Z[0];

    if (Z[1] > 1.0) Z[1] = 2.0 - Z[1];
    else if (Z[1] < -1.0) Z[1] = -2.0 - Z[1];

    if (Z[2] > 1.0) Z[2] = 2.0 - Z[2];
    else if (Z[2] < -1.0) Z[2] = -2.0 - Z[2];

    float r2 = Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2];

    if (r2 < mR2)
    {
        Z = Z * DIV(fR2, mR2);
    }
    else if (r2 < fR2)
    {
        Z = Z * DIV(fR2, r2);
    }

    if (julia.x == 0.0)
    {
        Z = Z * scale + P_in;
    }
    else
    {
        Z = Z * scale + point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

// [M2] - Classic Mandelbulb Power 2 fractal - MandelbulbPower2Iteration
void mandelbulbPower2Iter(point Z, point P_in, float weight, float4 julia) 
{
    point Z_orig = Z;
    
    float x2 = Z[0] * Z[0];
    float y2 = Z[1] * Z[1];
    float z2 = Z[2] * Z[2];
    float temp = DIV(1.0 - z2, (x2 + y2));
    point new_p;
    new_p[0] = (x2 - y2) * temp;
    new_p[0] = 2.0 * Z[0] * Z[1] * temp;
    new_p[0] = -2.0 * Z[2] * SQRT(x2 + y2);

    if (julia.x == 0.0)
    {
        Z = new_p + P_in;
    }
    else
    {
        Z = new_p +  point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

// [M2] - Bristorbrot formula - BristorbrotIteration
void bristorbrotIter(point Z, point P_in, float weight, float4 julia)
{
    point Z_orig = Z;
    
    point new_p;
    new_p[0] = Z[0] * Z[0] - Z[1] * Z[1] - Z[2] * Z[2];
    new_p[1] = Z[1] * (2.0 * Z[0] - Z[2]);
    new_p[2] = Z[2] * (2.0 * Z[0] + Z[1]);
    
    Z = new_p;

    if (julia.x == 0.0)
    {
        Z += P_in;
    }
    else 
    {
        Z += point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

// [M2] - Xenodreambuie - XenodreambuieIteration
void xenodreambuieIter(point Z, point P_in, float weight, float4 julia, float power, float alpha_parm, float beta_parm)
{
    point Z_orig = Z; 

    float distance = LENGTH(Z);

    float alpha = radians(alpha_parm);
    float beta = radians(beta_parm);

    float rp = POWR(distance, power - 1.0);
    rp *= distance;

    float th = atan2(Z[1], Z[0]) + beta;
    float ph = acos( DIV(Z[2], distance)) + alpha;

    if (fabs(ph) > 0.5 * M_PI_F) ph = sign(ph) * M_PI_F - ph;

    Z[0] = rp * COS(th * power) * SIN(ph * power);
    Z[1] = rp * SIN(th * power) * SIN(ph * power);
    Z[2] = rp * COS(ph * power);

    if (julia.x == 0.0)
    {
        Z += P_in;
    }
    else 
    {
        Z += point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

// [M2] - Benesi formula invented by Benesi - BenesiIteration
void benesiIter(point Z, point P_in, float weight, float4 julia)
{
    point Z_orig = Z;
    
    float distance = LENGTH(Z);

    float r1 = Z[1] * Z[1] + Z[2] * Z[2];
    float newx;
	if (P_in[0] < 0.0 || Z[0] < SQRT(r1))
	{
		newx = Z[0] * Z[0] - r1;
	}
	else
	{
		newx = -Z[0] * Z[0] + r1;
	}
	r1 = DIV(-1.0, SQRT(r1)) * 2.0 * fabs(Z[0]);
	float newy = r1 * (Z[1] * Z[1] - Z[2] * Z[2]);
	float newz = r1 * 2.0 * Z[1] * Z[2];

	Z = point(newx, newy, newz);

    if (julia.x == 0.0)
    {
        Z += P_in;
    }
    else 
    {
        Z += point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

// [M2] - Ides formula made by Trafassel, the original Ide's Formula thread - IdesIteration
void idesIter(point Z, point P_in, float weight, float4 julia, point multiplier, float2 sub_multiplier)
{
    point Z_orig = Z;
    
	if (fabs(Z[0]) < 2.5) Z[0] = Z[0] * 0.9;
	if (fabs(Z[2]) < 2.5) Z[2] = Z[2] * 0.9;

	point z2 = Z * Z;
	point newZ;
	newZ[0] = multiplier[0] * z2[0] - sub_multiplier.x * (z2[1] + z2[2]);
	newZ[1] = multiplier[1] * Z[0] * Z[1] * Z[2];
	newZ[2] = multiplier[2] * z2[2] - sub_multiplier.y * (z2[0] + z2[1]);
	Z = newZ;

    if (julia.x == 0.0)
    {
        Z += P_in;
    }
    else 
    {
        Z += point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

// [M2] - IQ Bulb from Mandelbulb 3D and Inigo Quilez - IqBulbIteration
void iqBulbIter(point Z, point P_in, float weight, float4 julia, float power, float zpower)
{
    point Z_orig = Z;
    
    float distance = LENGTH(Z);

	// extract polar coordinates
	float wr = distance;
	float wo = acos( DIV(Z[1], wr));
	float wi = atan2(Z[0], Z[2]);

	// scale and rotate the point
	wr = POWR(wr, power - 1.0);
	wr *= distance;
	wo *= power;
	wi *= zpower;

	// convert back to cartesian coordinates
	Z[0] = SIN(wo) * SIN(wi);
	Z[1] = COS(wo);
	Z[2] = SIN(wo) * COS(wi);

	Z *= wr; // then add Cpixel constant

    if (julia.x == 0.0)
    {
        Z += P_in;
    }
    else 
    {
        Z += point(julia.y, julia.z, julia.w);
    }

    Z = mix(Z_orig, Z, weight);
}

#endif
