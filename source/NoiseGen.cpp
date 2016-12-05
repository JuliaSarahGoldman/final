/** \file NoiseGen.cpp */
#include "NoiseGen.h"

using namespace G3D;

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

void NoiseGen::colorMountainImage(shared_ptr<Image> noise, shared_ptr<Image> colorMap) {
    for(int x(0); x < noise->width(); ++x){
        for(int y(0); y < noise->height(); ++y) {
            Color1unorm8 color;
            noise->get(Point2int32(x, y), color);
            if((float) color.value < 0.65f) {
                colorMap->set(x, y, Color3(1.0f));
            } else {
                colorMap->set(x, y, Color3(0.5f));
            }
        }
    }
}

void NoiseGen::colorLandImage(shared_ptr<Image> noise, shared_ptr<Image> colorMap, float oceanLevel) {
    for(int x(0); x < noise->width(); ++x){
        for(int y(0); y < noise->height(); ++y) {
            Color1unorm8 color;
            noise->get(Point2int32(x, y), color);
            /*if((float) color.value > (oceanLevel - 0.1f)) {
                colorMap->set(x, y, Color3(204, 204, 0));
            } else {*/
                float red = G3D::Random::threadCommon().uniform(0.0f, 0.5f);
                float green = G3D::Random::threadCommon().uniform(0.5f, 1.0f);
                float blue = G3D::Random::threadCommon().uniform(0.0f, 0.5f);
                colorMap->set(x, y, Color3(red, green, blue));
           // }
        }
    }
}

NoiseGen::NoiseGen() {
}