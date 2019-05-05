Renderers and techniques comparison
==========================

## Techniques

### A) Polygons, converted from sdf
<img src="img/comparison/m_poly_10min.jpg" width="300px"> <img src="img/comparison/a_poly_2min.jpg" width="300px"> <img src="img/comparison/rs_poly_1min.jpg" width="300px">

Mantra ~ 10 min, Arnold ~ 2 min, Redshift ~ 1 min

<br>

### B) Points (sphere trace), ~ 50.6 M
<img src="img/comparison/m_pts_11min.jpg" width="300px"> <img src="img/comparison/a_pts_4min.jpg" width="300px"> <img src="img/comparison/rs_pts_4min.jpg" width="300px">

Mantra ~ 11 min, Arnold ~ 4 min, Redshift ~ 4 min

<br>

### C) Points (perspective camera trace), ~ 10.3 M (4 samples/pixel)
<img src="img/comparison/m_pts_cam_8min.jpg" width="300px"> <img src="img/comparison/a_pts_cam_1min.jpg" width="300px"> <img src="img/comparison/rs_pts_cam_2min.jpg" width="300px">

Mantra ~ 8 min, Arnold ~ 1 min, Redshift ~ 2 min

<br>

### D) Volume, ~ 90 M voxels
<img src="img/comparison/m_fog_spec_13min.jpg" width="300px"> <img src="img/comparison/m_fog_8min.jpg" width="300px"> <img src="img/comparison/a_fog_9min.jpg" width="300px"> <img src="img/comparison/rs_fog_1min.jpg" width="300px">

Mantra (with GGX) ~ 13 min, Mantra ~ 8 min, Arnold ~ 9 min, Redshift ~ 1 min

<br>

### E) SDF, ~ 70 M voxels
<img src="img/comparison/m_sdf_16min.jpg" width="300px"> <img src="img/comparison/a_sdf_3min.jpg" width="300px">

Mantra ~ 16 min, Arnold ~ 3 min

<br>

## Mantra sdf vs polygons vs points (sphere traced) vs points (camera traced) vs fog vs fog (with GGX)
<img src="img/comparison/m_sdf_16min.jpg" width="300px"> <img src="img/comparison/m_poly_10min.jpg" width="300px"> <img src="img/comparison/m_pts_11min.jpg" width="300px"> <img src="img/comparison/m_pts_cam_8min.jpg" width="300px"> <img src="img/comparison/m_fog_8min.jpg" width="300px"> <img src="img/comparison/m_fog_spec_13min.jpg" width="300px">

<br>

## Arnold points (camera traced), ~ 1 min
<img src="img/comparison/a_pts_cam_1min.jpg" width="300px"> <img src="img/comparison/a_pts_cam_3spp_1min.jpg" width="300px"> <img src="img/comparison/a_pts_cam_2spp_1min.jpg" width="300px"> <img src="img/comparison/a_pts_cam_1_5spp_1min.jpg" width="300px">

4 spp (~ 10.3 M points), 3 spp (~ 5.8 M points), 2 spp (~ 2.5 M points), 1.5 spp (~ 1.5 M points)

<br>

## Points (sphere traced) with GGX speculars, ~ 87 M points
<img src="img/comparison/m_pts_spec_3min.jpg" width="300px"> <img src="img/comparison/a_pts_spec_1min.jpg" width="300px"> <img src="img/comparison/a_pts_spec_too_long.jpg" width="300px"> <img src="img/comparison/rs_pts_spec_spheres_incorrect.jpg" width="300px">

Mantra ~ 3 min, Arnold ~ 1 min (particles as discs with correct normals), Arnold (using instancing, too slow), Redshift ~ 3 min (as spheres - incorrect speculars, instancing not possible because of VRAM)

<br>

## Grid vs Frustum volumes comparison of volume size and render time
<img src="img/comparison/grid_hd.jpg" width="300px"> <img src="img/comparison/frustum_hd.jpg" width="300px"> <img src="img/comparison/frustum_hd_75.jpg" width="300px"> <img src="img/comparison/frustum_hd_5.jpg" width="300px"> <img src="img/comparison/frustum_hd_25.jpg" width="300px">

Grid volume (~ 800MB, 272M voxels) - 3:30s, frustum volume 1.0 z (~ 800MB, 272M voxels) - 5:20s, frustum volume 0.75 z (~ 600MB, 204M voxels) - 5:21s, frustum volume 0.5 z (~ 400MB, 136M voxels) - 5:27s, frustum volume 0.25 z (~ 200MB, 68M voxels) - 5:42s

**1.0 z** means scaling of resolution along Z axis (inside of the image)