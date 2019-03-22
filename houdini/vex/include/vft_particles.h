#ifndef VFT_PARTICLES
#define VFT_PARTICLES

function vector transform_right(vector vx; int i) {
    int i_mod = i % 3;
    
    if (i_mod == 0) return vx;
    if (i_mod == 1) return set(vx.z,vx.x,vx.y);
    if (i_mod == 2) return set(vx.y,vx.z,vx.x);
    else return vx;
}

function vector transform_left(vector vx; int i) {
    int i_mod = i % 3;

    if (i_mod == 0) return vx;
    if (i_mod == 1) return set(vx.y,vx.z,vx.x);
    if (i_mod == 2) return set(vx.z,vx.x,vx.y);
    else return vx;
}

function vector transform_full(vector vx; int i) {
    int i_mod = i % 6;

    if (i_mod == 0) return vx;
    if (i_mod == 1) return set(vx.y,vx.z,vx.x);
    if (i_mod == 2) return set(vx.z,vx.x,vx.y);
    if (i_mod == 3) return set(vx.x,vx.z,vx.y);
    if (i_mod == 4) return set(vx.z,vx.y,vx.x);
    if (i_mod == 5) return set(vx.y,vx.x,vx.z);

    else return vx;
}

function vector transform_straight(vector vx; int i) {
    return vx;
}

#endif
