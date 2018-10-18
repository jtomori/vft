#ifndef VFT_BLINK_FUNCS
#define VFT_BLINK_FUNCS

static float fit(float value, float src_min, float src_max, float dst_min, float dst_max) {
    return dst_min + (dst_max - dst_min) * (value - src_min) / (src_max - src_min);
}

#endif
