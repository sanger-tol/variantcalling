#!/bin/sh

## render static maps
nf-metro render variantcalling.mmd \
    -o sanger-tol-variantcalling_metro_map_dark.svg \
    --logo sanger-tol-variantcalling_logo_dark.png
nf-metro render variantcalling.mmd \
    -o sanger-tol-variantcalling_metro_map_light.svg \
    --theme light \
    --logo sanger-tol-variantcalling_logo_light.png

## render animated maps
nf-metro render variantcalling.mmd \
    --animate \
    -o sanger-tol-variantcalling_metro_map_dark_animated.svg \
    --logo sanger-tol-variantcalling_logo_dark.png
nf-metro render variantcalling.mmd \
    --animate \
    -o sanger-tol-variantcalling_metro_map_light_animated.svg \
    --theme light \
    --logo sanger-tol-variantcalling_logo_light.png