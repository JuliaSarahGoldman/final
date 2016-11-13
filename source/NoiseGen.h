#pragma once
/**
  \file Noise.h

  Creates noise for the heightfields and colors used in terrain creation.
 */
#pragma once
#include <G3D/G3DAll.h>

/** \brief Application framework. */
class NoiseGen : public G3D::Noise {
protected:

    Random rng1;
    Random rng2;

public:

    NoiseGen(float seed1, float seed2);

    float sampleFloat(int x, int y, int numOctaves = 1);

    float noise2(float nx, float ny);

};