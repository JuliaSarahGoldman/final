/** \file NoiseGen.cpp */
#include "NoiseGen.h"

using namespace G3D;

float NoiseGen::sampleFloat(int x, int y, int numOctaves) {
    return 1.0;
}


void NoiseGen::generateMountainImage(shared_ptr<Image> image, int frequency){
    for(int x(0); x < image->width(); ++x){
        for(int y(0); y < image->height(); ++y) {
            image->set(x, y, Color1unorm8(unorm8::fromBits(sampleUint8(x << 12, y <<12, 0))));
            //image->set(x, y, Color3(Noise::sampleFloat(x, y, 0, 4)));
        }
    }
}

void NoiseGen::generateLandImage(shared_ptr<Image> image, int frequency){
    for(int x(0); x < image->width(); ++x){
        for(int y(0); y < image->height(); ++y) {
            image->set(x, y, Color1unorm8(unorm8::fromBits(sampleUint8((x * frequency) << 12, (y * frequency) <<12, 0))));
        }
    }
}

void NoiseGen::generateSeaImage(shared_ptr<Image> image, int frequency){
    for(int x(0); x < image->width(); ++x){
        for(int y(0); y < image->height(); ++y) {
            image->set(x, y, Color1unorm8(unorm8::fromBits(sampleUint8((x * frequency) << 12, (y * frequency) <<12, 0))));
        }
    }
}

void NoiseGen::generateNoisyImage(shared_ptr<Image> image, int frequency){
    for(int x(0); x < image->width(); ++x){
        for(int y(0); y < image->height(); ++y) {
            image->set(x, y, Color1unorm8(unorm8::fromBits(sampleUint8((x * frequency) << 12, (y * frequency) <<12, 0))));
        }
    }
}

void NoiseGen::generateNoisyImage(shared_ptr<Image> image, int type, int frequency) {
    if (type == 0) {
        generateSeaImage(image, frequency);
    }
    else if (type == 1) {
        generateLandImage(image, frequency);
    }
    else if (type == 2) {
        generateMountainImage(image, frequency);
    }
    else {
        generateNoisyImage(image, frequency);
    }
}

float NoiseGen::noise2(float nx, float ny) {
    return 1.0;
}

NoiseGen::NoiseGen() {
}