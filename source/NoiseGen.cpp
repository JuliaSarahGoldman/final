/** \file NoiseGen.cpp */
#include "NoiseGen.h"

using namespace G3D;

float NoiseGen::sampleFloat(int x, int y, int numOctaves) {
    return 1.0;
}


void NoiseGen::generateMountainImage(shared_ptr<Image> image, float frequency, float multiplier){
    for(int x(0); x < image->width(); ++x){
        for(int y(0); y < image->height(); ++y) {
            Color1 col;
            image->get(Point2int32(x, y), col);
            Color1 col2 = multiplier * Color1(unorm8::fromBits(sampleUint8( (int) (x * frequency) << 12, (int) (y * frequency) <<12, 0)));
            image->set(x, y, col + col2);
        }
    }
}

void NoiseGen::generateLandImage(shared_ptr<Image> image, float frequency){
    for(int x(0); x < image->width(); ++x){
        for(int y(0); y < image->height(); ++y) {
            image->set(x, y, Color1unorm8(unorm8::fromBits(sampleUint8( (int) (x * frequency) << 12, (int) (y * frequency) <<12, 0))));
        }
    }
}

void NoiseGen::generateSeaImage(shared_ptr<Image> image, float frequency){
    for(int x(0); x < image->width(); ++x){
        for(int y(0); y < image->height(); ++y) {
            image->set(x, y, Color1unorm8(unorm8::fromBits(sampleUint8( (int) (x * frequency) << 12, (int) (y * frequency) <<12, 0))));
        }
    }
}

NoiseGen::NoiseGen() {
}