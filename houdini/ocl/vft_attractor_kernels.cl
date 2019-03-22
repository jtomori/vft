kernel void lorenzAttractor( 
    int P_length, 
    global float* P,
    int k_length,
    global int* k_index,
    global float* k
)
{
    int idx = get_global_id(0);
    if (idx >= P_length)
        return;

    float3 p = vload3(idx, P);
    float a = k[0];
    float b = k[1];
    float c = k[2];
    float dt = 0.0005;

    float3 d;
    d.x = a * (p.y - p.x);
    d.y = p.x * (b - p.z) - p.y;
    d.z = p.x * p.y - c * p.z;

    p += d * dt;
    vstore3(p, idx, P);
}
