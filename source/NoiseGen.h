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

    void generateMountainImage(shared_ptr<Image> image, int frequency);
    void generateLandImage(shared_ptr<Image> image, int frequency);
    void generateSeaImage(shared_ptr<Image> image, int frequency);

public:

    NoiseGen();

    float sampleFloat(int x, int y, int numOctaves = 1);

    float noise2(float nx, float ny);
    
    void generateNoisyImage(shared_ptr<Image> image, int frequency);

    void generateNoisyImage(shared_ptr<Image> image, int type, int frequency);
};