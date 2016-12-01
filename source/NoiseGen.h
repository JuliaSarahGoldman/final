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

    NoiseGen();

    float sampleFloat(int x, int y, int numOctaves = 1);

    void generateMountainImage(shared_ptr<Image> image, float frequency, float multiplier);

    void generateLandImage(shared_ptr<Image> image, float frequency);

    void generateSeaImage(shared_ptr<Image> image, float frequency);

    void colorMountainImage(shared_ptr<Image> noise, shared_ptr<Image> colorMap);

    void colorLandImage(shared_ptr<Image> noise, shared_ptr<Image> colorMap, float oceanLevel);
};