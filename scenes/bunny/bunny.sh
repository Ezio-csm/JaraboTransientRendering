#!/bin/bash

bunnykiller=../../build/BunnyKiller

camera_samples=4
bdpt_samples=1
max_bounces=4

height=600
width=600

time_res=900
toffset=0
exp=0.025

camera_pos='1.6 0.5 -1.6'
camera_at='-0.5 0.0 1.5'
camera_fov=60

u_a='0.1'
u_s='0.1'

name='bunny_image'

$bunnykiller \
    -camera-spp $camera_samples -camera-fov $camera_fov \
    -bidirectional-path-tracing $bdpt_samples \
    -max-nb-bounces $max_bounces -steady-state \
    -film-name $name -film-size-x $width -film-size-y $height \
    -film-size-t $time_res -film-exposure $exp -film-offset $toffset \
    -film-store-depth -film-store-normals -film-store-positions \
    -log-name $name.txt -camera-position $camera_pos -camera-focus $camera_at \
    -point-light-source 0.0 1.0 -0.5 10.0 10.0 10.0 \
    -name-mesh 'geometry/bunny.obj' -lambertian 0.5 0.0 0.5 \
    -name-mesh 'geometry/wall.obj'  -lambertian 1.0 0.0 0.0 \
    -name-mesh 'geometry/wall2.obj' -lambertian 0.7 0.7 0.7 \
    -name-mesh 'geometry/floor.obj' -lambertian 0.7 0.7 0.7

