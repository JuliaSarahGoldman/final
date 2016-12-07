/**
  \file Noise.h

  Creates noise for the heightfields and colors used in terrain creation.
 */
#pragma once
#include <G3D/G3DAll.h>

/** \brief Application framework. */
class NoiseGen : public G3D::Noise {
public:

    NoiseGen();

    void generateMountainImage(shared_ptr<Image> image, float frequency, float multiplier);

    void generateLandImage(shared_ptr<Image> image, float frequency);

    void generateSeaImage(shared_ptr<Image> image, float frequency);

    void colorMountainImage(shared_ptr<Image> noise, shared_ptr<Image> colorMap);

    void landMapImage(shared_ptr<Image> land, shared_ptr<Image> mountain, shared_ptr<Image> landMap, float oceanLevel, float power, float multiplier);
};