#ifndef VFT_BLINK_FUNCS
#define VFT_BLINK_FUNCS

static float fit(float value, float src_min, float src_max, float dst_min, float dst_max) {
    return dst_min + (dst_max - dst_min) * (value - src_min) / (src_max - src_min);
}

static float distPtLine(float2 pt, float2 line_pt_1, float2 line_pt_2)
{
    float numerator = fabs( pt.x*( line_pt_2.y - line_pt_1.y ) - pt.y*( line_pt_2.x - line_pt_1.x ) + line_pt_2.x*line_pt_1.y - line_pt_2.y * line_pt_1.x );
    float denominator = sqrt( pow( line_pt_2.y - line_pt_1.y , 2.0f) + pow( line_pt_2.x - line_pt_1.x , 2.0f) );
    return numerator / denominator;
}

#endif
