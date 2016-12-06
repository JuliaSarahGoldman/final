#include "Planet.h"
#include "Mesh.h"
#include "NoiseGen.h"

Planet::Planet() {

}

Planet::Planet(const String& name, const Any& planetSpec) {
    m_planetSpec = planetSpec;
    m_planetName = name;
    readSpec(planetSpec);
}

bool Planet::readSpec(const Any& planetSpec) {
    try {
        AnyTableReader x(planetSpec);

        x.getIfPresent("recursions", m_recursionLevel);
        x.getIfPresent("landBevel", m_landBevel);
        x.getIfPresent("mountainBevel", m_mountainBevel);
        x.getIfPresent("mountainHeight", m_mountainHeight);
        x.getIfPresent("mountainDiversity", m_mountianDiversity);
        x.getIfPresent("mountainNoise1", m_mountianNoise1);
        x.getIfPresent("mountainNoise2", m_mountianNoise2);
        x.getIfPresent("mountainNoise3", m_mountianNoise3);
        x.getIfPresent("oceanLevel", m_oceanLevel);
        x.getIfPresent("landNoise", m_landNoise);
        x.getIfPresent("oceanNoise", m_oceanNoise);


        x.getIfPresent("scale", m_scale);
        x.getIfPresent("planetName", m_planetName);
        x.getIfPresent("orbitPlanet", m_objectToOrbit);
        x.getIfPresent("orbitDistance", m_orbitDistance);

        x.getIfPresent("collapsingEnabled", m_collapsingEnabled);
        x.getIfPresent("oceanCollapsing", m_oceanEdgesToCollapse);
        x.getIfPresent("landCollapsing", m_landEdgesToCollapse);
        x.getIfPresent("mountainCollapsing", m_mountainEdgesToCollapse);

        x.getIfPresent("waterMount", m_waterMount);

        x.getIfPresent("useMountainTexture", m_useMTexture);
        x.getIfPresent("useWaterTexture", m_useWTexture);
        x.getIfPresent("useLandTexture", m_useLTexture);
        x.getIfPresent("mountainTexture", m_mountainTextureFile);
        x.getIfPresent("landTexture", m_landTextureFile);
        x.getIfPresent("waterTexture", m_waterTextureFile);

        float xPos, yPos, zPos;
        x.getIfPresent("xPos", xPos);
        x.getIfPresent("yPos", yPos);
        x.getIfPresent("zPos", zPos);
        m_position = Point3(xPos, yPos, zPos);

        float red, green, blue;
        x.getIfPresent("mRed", red);
        x.getIfPresent("mGreen", green);
        x.getIfPresent("mBlue", blue);
        m_mountainColor = Color3(red, green, blue);

        x.getIfPresent("lRed", red);
        x.getIfPresent("lGreen", green);
        x.getIfPresent("lBlue", blue);
        m_landColor = Color3(red, green, blue);

        x.getIfPresent("wRed", red);
        x.getIfPresent("wGreen", green);
        x.getIfPresent("wBlue", blue);
        m_waterColor = Color3(red, green, blue);

        float glossBase, glossPower;
        x.getIfPresent("mBase", glossBase);
        x.getIfPresent("mPow", glossPower);
        m_mountainGloss = Color4(Color3(glossBase), glossPower);

        x.getIfPresent("lBase", glossBase);
        x.getIfPresent("lPow", glossPower);
        m_landGloss = Color4(Color3(glossBase), glossPower);

        x.getIfPresent("wBase", glossBase);
        x.getIfPresent("wPow", glossPower);
        m_waterGloss = Color4(Color3(glossBase), glossPower);


        return true;
    }
    catch (...) {
        return false;
    }
}

bool Planet::generatePlanet() {
    try {
        Array<Vector3> vertices = Array<Vector3>();
        Array<Vector3int32> faces = Array<Vector3int32>();
        NoiseGen noise;
        AnyTableReader planetReader(m_planetSpec);

        shared_ptr<Image> noiseImage = Image::create(1024, 1024, ImageFormat::RGBA8());
        shared_ptr<Image> colorImage = Image::create(1024, 1024, ImageFormat::RGBA8());
        shared_ptr<Image> testImage = Image::create(1024, 1024, ImageFormat::RGBA8());

        //noise.generateSeaImage(noiseImage, m_oceanNoise);
        writeSphere(m_planetName, 12.5f, 5, vertices, faces);
        //applyNoiseWater(vertices, noiseImage);
        //noiseImage->save(m_planetName + "water.png");

        Mesh mesh(vertices, faces);
        //mesh.bevelEdges2(0.1f);
        m_waterObjFile = m_planetName + "water";

        int width, height;
        if (m_useWTexture){
            shared_ptr<Image> image = Image::fromFile(m_waterTextureFile);
            width = image->width();
            height = image->height();
        }
        else{
            width = 1;
            height = 1;
        }

        mesh.toObj(m_waterObjFile, width, height);


        noiseImage = Image::create(1024, 1024, ImageFormat::RGBA8());
        colorImage = Image::create(1024, 1024, ImageFormat::RGBA8());
        vertices.clear();
        faces.clear();

        noise.generateLandImage(noiseImage, m_landNoise);
        noiseImage->save(m_planetName + "land.png");

        writeSphere(m_planetName + "land", 12.0f, m_recursionLevel, vertices, faces);

        noise.colorLandImage(noiseImage, colorImage, m_oceanLevel);
        applyNoiseLand(vertices, noiseImage, testImage, m_oceanLevel);

        colorImage->save(m_planetName + "landColor.png");
        applyNoiseLand(vertices, noiseImage, testImage, m_oceanLevel);
        testImage->save(m_planetName + "landTest.png");

        Mesh mesh2(vertices, faces);
        mesh2.bevelEdges2(m_landBevel);
        m_landObjFile = m_planetName + "land";

        if (m_useLTexture){
            shared_ptr<Image> image = Image::fromFile(m_landTextureFile);
            width = image->width();
            height = image->height();
        }
        else{
            width = 1;
            height = 1;
        }

        mesh2.toObj(m_landObjFile, width, height);


        noiseImage = Image::create(1024, 1024, ImageFormat::RGBA8());
        colorImage = Image::create(1024, 1024, ImageFormat::RGBA8());
        vertices.clear();
        faces.clear();

        writeSphere(m_planetName + "mountain", 11.5f, m_recursionLevel, vertices, faces);

        noise.generateMountainImage(noiseImage, m_mountianNoise1, 1.0f);
        noise.generateMountainImage(noiseImage, m_mountianNoise2, 0.5f);
        noise.generateMountainImage(noiseImage, m_mountianNoise1, 0.25f);

        noise.colorMountainImage(noiseImage, colorImage);
        noiseImage->save(m_planetName + "mountain.png");
        colorImage->save(m_planetName + "mountainColor.png");
        applyNoiseMountain(vertices, noiseImage, testImage, m_waterMount, m_mountianDiversity, m_mountainHeight);

        Mesh mesh3(vertices, faces);
        mesh3.bevelEdges2(m_mountainBevel);
        m_mountainObjFile = m_planetName + "mountain";

        if (m_useMTexture){
            shared_ptr<Image> image = Image::fromFile(m_mountainTextureFile);
            width = image->width();
            height = image->height();
        }
        else{
            width = 1;
            height = 1;
        }

        mesh3.toObj(m_mountainObjFile, width, height);

    }
    catch (...) {
        return false;
    }
    return true;
}

void Planet::createWaterAnyFile(Any& waterModel, Any& waterEntity) {

    String anyStr = (m_useWTexture && !m_waterTextureFile.empty()) ? (String) "UniversalMaterial::Specification { glossy = Color4" + m_waterGloss.toString() + "; "
        "lambertian = Texture::Specification{ filename = \"" +
        m_waterTextureFile + "\"; encoding = Texture::Encoding { readMultiplyFirst = Color3" + m_waterColor.toString() + "}; }; }" :
        "UniversalMaterial::Specification { lambertian = Color3" + m_waterColor.toString() + "; glossy = Color4" + m_waterGloss.toString() + "; }";

    /*String anyStr = (m_useWTexture && !m_waterTextureFile.empty()) ? (String) "UniversalMaterial::Specification { glossy = Color4" + m_waterGloss.toString() + "; "
        "lambertian = Texture::Specification{ filename = \"" +
        m_waterTextureFile + "\"; encoding = Texture::Encoding { readMultiplyFirst = " + m_waterColor.toString() + "}; }; }" :
        "UniversalMaterial::Specification { lambertian = Color3" + m_waterColor.toString() + "; glossy = Color4" + m_waterGloss.toString() + "; }";*/

    String preprocess = "{ setMaterial(all()," + anyStr + ");" +
        "transformGeometry(all(), Matrix4::scale(" +
        (String)std::to_string(m_scale) + ", " + (String)std::to_string(m_scale) + "," + (String)std::to_string(m_scale) + ")); }";

    waterModel["filename"] = m_waterObjFile + ".obj";

    waterModel["preprocess"] = Any::parse(preprocess);

    Point3 position = getPosition();
    waterEntity["frame"] = CFrame::fromXYZYPRDegrees(position.x, position.y, position.z);

    String planetToOrbit;
    float orbitDistance;
    getPlanetOrbit(planetToOrbit, orbitDistance);

    if (planetToOrbit.empty()) {
        waterEntity["track"] = Any::parse("combine(orbit(0.0f, 10.0f), Point3" + position.toString() + " )");
    }
    else {
        waterEntity["track"] = Any::parse("transform(transform(entity(\"" + planetToOrbit + "\"), orbit(" + (String)std::to_string(orbitDistance) + ", 10.0f)), orbit(0, 7))");
    }
}

void Planet::createMountainAnyFile(Any& mountainModel, Any& mountainEntity, const String& waterEntity) {
    String anyStr = (m_useMTexture && !m_mountainTextureFile.empty()) ? (String) "UniversalMaterial::Specification { glossy = Color4" + m_mountainGloss.toString() + "; "
        "lambertian = Texture::Specification{ filename = \"" +
        m_mountainTextureFile + "\"; encoding = Texture::Encoding { readMultiplyFirst = Color3" + m_mountainColor.toString() + "}; }; }" :
        "UniversalMaterial::Specification { lambertian = Color3" + m_mountainColor.toString() + "; glossy = Color4" + m_mountainGloss.toString() + "; }";

    String preprocess = "{ setMaterial(all()," + anyStr + ");" +
        "transformGeometry(all(), Matrix4::scale(" +
        (String)std::to_string(m_scale) + ", " + (String)std::to_string(m_scale) + "," + (String)std::to_string(m_scale) + ")); }";

    mountainModel["filename"] = m_mountainObjFile + ".obj";

    mountainModel["preprocess"] = Any::parse(preprocess);
    Point3 position = getPosition();
    mountainEntity["frame"] = CFrame::fromXYZYPRDegrees(position.x, position.y, position.z);
    mountainEntity["track"] = Any::parse("entity(" + waterEntity + ")");
}

void Planet::createLandAnyFile(Any& landModel, Any& landEntity, const String& waterEntity) {
    String anyStr = (m_useLTexture && !m_landTextureFile.empty()) ? (String) "UniversalMaterial::Specification { glossy = Color4" + m_landGloss.toString() + "; "
        "lambertian = Texture::Specification{ filename = \"" +
        m_landTextureFile + "\"; encoding = Texture::Encoding { readMultiplyFirst = Color3" + m_landColor.toString() + "}; }; }" :
        "UniversalMaterial::Specification { lambertian = Color3" + m_landColor.toString() + "; glossy = Color4" + m_landGloss.toString() + "; }";

    String preprocess = "{ setMaterial(all()," + anyStr + ");" +
        "transformGeometry(all(), Matrix4::scale(" +
        (String)std::to_string(m_scale) + ", " + (String)std::to_string(m_scale) + "," + (String)std::to_string(m_scale) + ")); }";

    landModel["filename"] = m_landObjFile + ".obj";

    landModel["preprocess"] = Any::parse(preprocess);
    Point3 position = getPosition();
    landEntity["frame"] = CFrame::fromXYZYPRDegrees(position.x, position.y, position.z);
    landEntity["track"] = Any::parse("entity(" + waterEntity + ")");

}

void Planet::addCloudToPlanet(Any& cloudEntity, const String& name, const String& planetName, const Point3& position, const float scale) {

    int t[5] = { 2, 3, 5, 7, 11 };
    cloudEntity["model"] = name;
    cloudEntity["track"] = Any::parse((String)
        "transform(" +
        "Matrix4::rollDegrees(" + (String)std::to_string(Random::threadCommon().integer(-89, 89)) + "), " +
        "transform("
        "orbit(" +
        (String)std::to_string(Random::threadCommon().integer(40, 70)*getScale()) + ", " + (String)std::to_string(t[Random::threadCommon().integer(0, 4)]) +
        "), " +
        "combine(" +
        "Matrix4::pitchDegrees(90), entity(" + planetName + ")" +
        ")"
        "), " +
        ");");
}

void Planet::writeSphere(String filename, float radius, int depths, Array<Vector3>& vertices, Array<Vector3int32>& faces) {

    makeIcohedron(radius, vertices, faces);

    for (int i(0); i < depths; ++i) {
        subdivideIcoHedron(radius, vertices, faces);
    }

    Array<Vector2> texture;
    Array<Vector3> normals;
    Array<Vector3> verts = vertices;

    Array<int> indices;
    for (int i(0); i < faces.length(); ++i) {
        Vector3int32 face = faces[i];
        indices.append(face.x, face.y, face.z);
    }

    // # WELDING
    Welder::weld(verts, texture, normals, indices, G3D::Welder::Settings());

    faces = Array<Vector3int32>();

    for (int i(0); i < indices.size() - 2; i += 3) {
        faces.append(Vector3int32(indices[i], indices[i + 1], indices[i + 2]));
    }

    vertices = verts;
}

void Planet::applyNoiseWater(Array<Vector3>& vertices, shared_ptr<Image> image) {
    for (int i(0); i < vertices.size(); ++i) {
        Vector3 vertex = vertices[i];

        Vector3 d = (vertex - Vector3(0, 0, 0)).unit();

        float nx = image->width() * (0.5f + atanf(d.z / d.x) / (2.0f*pif()));
        float ny = image->height() * (0.5f - asinf(d.y) * 1 / pif());

        int ix = (int)abs((int)nx % image->width());
        int iy = (int)abs((int)ny % image->height());

        Color3 color = Color3();

        image->get(Point2int32(ix, iy), color);

        float bump = 1.0f - color.average();
        /*if (bump > 0.3f && bump < 0.6f) bump = 0.5f;
        else if (bump < 0.3f) bump = 0.1f;
        else bump = 1.0f;*/

        if (bump < 0.4f) {
            bump = -0.2f;
        }
        else if (bump > 0.6f) {
            bump = 0.2f;
        }
        else {
            bump = 0.0f;
        }

        vertices[i] += vertex.unit() * bump;
    }
}

void Planet::applyNoiseLand(Array<Vector3>& vertices, shared_ptr<Image> noise, shared_ptr<Image> test, float oceanLevel) {
    for (int i(0); i < vertices.size(); ++i) {
        Vector3 vertex = vertices[i];

        Vector3 d = (vertex - Vector3(0, 0, 0)).unit();

        float nx = noise->width() * (0.5f + atanf(d.z / d.x) / (2.0f*pif()));
        float ny = noise->height() * (0.5f - asinf(d.y) * 1 / pif());

        int ix = (int)abs((int)nx % noise->width());
        int iy = (int)abs((int)ny % noise->height());

        Color3 color = Color3();

        noise->get(Point2int32(ix, iy), color);

        float bump = 1.0f - color.average();
        /*if (bump > 0.3f && bump < 0.6f) bump = 0.5f;
        else if (bump < 0.3f) bump = 0.1f;
        else bump = 1.0f;*/

        if (bump < oceanLevel) {
            bump = 0.0f;
            test->set(Point2int32(ix, iy), Color1(1.0f));
            /*} else if(bump < 0.7f) {
                bump *= 2.857f;
                test->set(Point2int32(ix, iy), Color1(0.0f));*/
        }
        else {
            bump *= 4.0f;
            test->set(Point2int32(ix, iy), Color1(0.0f));
        }

        vertices[i] += vertex.unit() * bump;
    }
}

void Planet::applyNoiseMountain(Array<Vector3>& vertices, shared_ptr<Image> noise, shared_ptr<Image> test, bool waterMount, float power, float multiplier) {
    for (int i(0); i < vertices.size(); ++i) {
        Vector3 vertex = vertices[i];

        Vector3 d = (vertex - Vector3(0, 0, 0)).unit();

        float nx = noise->width() * (0.5f + atanf(d.z / d.x) / (2.0f*pif()));
        float ny = noise->height() * (0.5f - asinf(d.y) * 1 / pif());

        int ix = (int)abs((int)nx % noise->width());
        int iy = (int)abs((int)ny % noise->height());

        Color3 color = Color3();

        noise->get(Point2int32(ix, iy), color);

        float bump = 1.0f - color.average();
        /*if (bump > 0.3f && bump < 0.6f) bump = 0.5f;
        else if (bump < 0.3f) bump = 0.1f;
        else bump = 1.0f;*/

        /*e =    1 * noise(1 * nx, 1 * ny);
                +  0.5 * noise(2 * nx, 2 * ny);
                 + 0.25 * noise(4 * nx, 4 * ny);
                elevation[y][x] = Math.pow(e, 3.00);*/
        test->get(Point2int32(ix, iy), color);

        if (!waterMount && (bump > 0.75f || color.average() != 0.0f)) {
            bump = 0.0f;
        }
        //else bump *= 1.25f;

        vertices[i] += vertex.unit() * pow(bump, power) * multiplier;
    }
}

//Creates an initial icohedron with the given radius to be tessellated to create a sphere
void Planet::makeIcohedron(float radius, Array<Vector3>& vertices, Array<Vector3int32>& faces) {
    float t = (1.0f + sqrt(5.0f)) / 2.0f;

    vertices.append(radius * Vector3(t, 1, 0));
    vertices.append(radius * Vector3(-t, 1, 0));
    vertices.append(radius * Vector3(t, -1, 0));
    vertices.append(radius * Vector3(-t, -1, 0));
    vertices.append(radius * Vector3(1, 0, t));
    vertices.append(radius * Vector3(1, 0, -t));
    vertices.append(radius * Vector3(-1, 0, t));
    vertices.append(radius * Vector3(-1, 0, -t));
    vertices.append(radius * Vector3(0, t, 1));
    vertices.append(radius * Vector3(0, -t, 1));
    vertices.append(radius * Vector3(0, t, -1));
    vertices.append(radius * Vector3(0, -t, -1));

    faces.append(Vector3int32(0, 8, 4));
    faces.append(Vector3int32(1, 10, 7));
    faces.append(Vector3int32(2, 9, 11));
    faces.append(Vector3int32(7, 3, 1));
    faces.append(Vector3int32(0, 5, 10));
    faces.append(Vector3int32(3, 9, 6));
    faces.append(Vector3int32(3, 11, 9));
    faces.append(Vector3int32(8, 6, 4));
    faces.append(Vector3int32(2, 4, 9));
    faces.append(Vector3int32(3, 7, 11));
    faces.append(Vector3int32(4, 2, 0));
    faces.append(Vector3int32(9, 4, 6));
    faces.append(Vector3int32(2, 11, 5));
    faces.append(Vector3int32(0, 10, 8));
    faces.append(Vector3int32(5, 0, 2));
    faces.append(Vector3int32(10, 5, 7));
    faces.append(Vector3int32(1, 6, 8));
    faces.append(Vector3int32(1, 8, 10));
    faces.append(Vector3int32(6, 1, 3));
    faces.append(Vector3int32(11, 7, 5));
}

void Planet::getMiddle(float radius, Vector3& v1, Vector3& v2, Vector3& newVector) {
    float t = (1.0f + sqrt(5.0f)) / 2.0f;
    newVector = (v2 + v1) / 2.0f;
    newVector = newVector.unit();
    newVector *= sqrt(t*t + 1)*radius;
}

void Planet::subdivideIcoHedron(float radius, Array<Vector3>& vertices, Array<Vector3int32>& faces) {
    float t = (1.0f + sqrt(5.0f)) / 2.0f;

    Vector3 newVec1;
    Vector3 newVec2;
    Vector3 newVec3;

    int numFaces = faces.length();

    Array<Vector3int32> newFaces = Array<Vector3int32>();

    for (int i(0); i < numFaces; ++i) {

        int numVertices = vertices.length();
        Vector3int32 face = faces[i];

        Vector3 vert1 = vertices[face.x];
        Vector3 vert2 = vertices[face.y];
        Vector3 vert3 = vertices[face.z];

        //Find the midpoints
        getMiddle(radius, vert1, vert2, newVec1);
        getMiddle(radius, vert2, vert3, newVec2);
        getMiddle(radius, vert3, vert1, newVec3);

        vertices.append(newVec1, newVec2, newVec3);

        int posVec1 = numVertices;
        int posVec2 = numVertices + 1;
        int posVec3 = numVertices + 2;

        //add the new faces
        newFaces.append(Vector3int32(posVec3, face.x, posVec1));
        newFaces.append(Vector3int32(posVec1, face.y, posVec2));
        newFaces.append(Vector3int32(posVec2, face.z, posVec3));
        newFaces.append(Vector3int32(posVec1, posVec2, posVec3));
    }
    faces = newFaces;
}

Point3 Planet::getPosition() {
    return m_position;
}

float Planet::getScale() {
    return m_scale;
}

void Planet::getPlanetOrbit(String& objectToOrbit, float& orbitDistance) {
    objectToOrbit = m_objectToOrbit;
    orbitDistance = m_orbitDistance;
}

String Planet::getName() {
    return m_planetName;
}


    /*cloud["canChange"] = false;
    cloud["particlesAreInWorldSpace"] = true;

    int t[5] = { 2, 3, 5, 7, 11 };

    cloud["model"] = (String) "cloud" + (String)std::to_string(Random::threadCommon().integer(1, 3));
    cloud["scale"] = scale;
    cloud["track"] = Any::parse((String)
        "transform(" +
        "Matrix4::rollDegrees(90), " +
        "transform("
        "orbit(" +
        (String)std::to_string(Random::threadCommon().integer(10, 20)) + ", " + (String)std::to_string(t[Random::threadCommon().integer(0, 4)]) +
        "), " +
        "combine(" +
        "Matrix4::pitchDegrees(" + (String)std::to_string(Random::threadCommon().integer(-89, 89)) + "), entity(" + name + ")" +
        ")"
        "), " +
        ");");
        */

/*
    Table<Vector3, int> vertexPositions;
    for(int i(0); i < vertices->length(); ++i) {
        vertexPositions.set(vertices->operator[](i), i);
    }

            if(!vertexPositions.get(newVec1, posVec1)) {
            vertices->append(newVec1);
            vertexPositions.set(newVec1, numVertices);
            posVec1 = numVertices;
        }
        if(!vertexPositions.get(newVec2, posVec2)) {
            vertices->append(newVec2);
            vertexPositions.set(newVec2, numVertices + 1);
            posVec2 = numVertices + 1;
        }
        if(!vertexPositions.get(newVec3, posVec3)) {
            vertices->append(newVec3);
            vertexPositions.set(newVec3, numVertices + 2);
            posVec3 = numVertices + 2;
        }
*/


/*int numVert = vertices->size();
int numFace = faces->size();

TextOutput output("model/" + filename + ".obj");

for (int i = 0; i < numVert; ++i) {
    Vector3 vertex = vertices->operator[](i);
    output.printf("v %f %f %f\n", vertex[0], vertex[1], vertex[2]);
}

for (int i = 0; i < numFace; ++i) {
    Vector3int32 face = faces->operator[](i);
    output.printf("f %d %d %d\n", face[0]+1, face[1]+1, face[2]+1);
}
output.commit(true);
G3D::ArticulatedModel::clearCache();*/


/*for(int i(0); i < 6; ++i){
    subdivideIcoHedron(radius, vertices, faces);
}
//const G3D::Welder::Settings settings;
Array<Vector3> normals;
Array<Vector2> texture;
Array<int> indices;
Array<Vector3> verts;

for (int i(0); i < vertices->length(); ++i) {
    verts.append(vertices->operator[](i));
}

for(int i(0); i < faces->length(); ++i){
    Vector3int32 face = faces->operator[](i);
    indices.append(face.x, face.y, face.z);
}

Welder::weld(verts, texture, normals, indices, G3D::Welder::Settings());

faces = std::make_shared<Array<Vector3int32>>();
for(int i(0); i < indices.size()-3; i += 3) {
    faces->append(Vector3int32(indices[i], indices[i+1], indices[i+2]));
}

for(int i(0); i < verts.size(); ++i) {
    verts[i] += normals[i] * Random::threadCommon().uniform(1.0, 1.2);
}
vertices = std::make_shared<Array<Vector3>>(verts);*/

/*
       float d = noise.sampleFloat(nx, ny, nz, Random::threadCommon().integer(0,4))*Random::threadCommon().uniform(50.0f,100.0f);
       float o = noise.sampleFloat(nz, ny, nx, Random::threadCommon().integer(0,4))*Random::threadCommon().uniform(50.0f,100.0f);
       float p = noise.sampleFloat(ny, nx, nz, Random::threadCommon().integer(0,4))*Random::threadCommon().uniform(50.0f,100.0f);
       float e = noise.sampleFloat(nx, nz, ny, Random::threadCommon().integer(0,4))*Random::threadCommon().uniform(50.0f,100.0f);
       */


       /* float lat = (float)acosf(vertex.z / radius);
        float lon =  (float)atanf(vertex.x / vertex.y);

        int nx = (int)(radius * lon * image->width()) % (image->width());//freq * vertex.x;
        int ny = (int)(image->height() * radius * log(tanf((lat + pif()/2.0f) / 2.0f))) % (image->height());//freq * vertex.y;


        */



        /*    float frequency = 1.0f;
            for(int x(0); x < image->width(); x++){
                for(int y(0); y < image->width(); y++){
                    image->set(x, y, Color1unorm8(unorm8::fromBits(noise.sampleUint8( frequency * (x << 12), frequency * (y << 12), 0))));
                }
            }
            */


            /*for (int i(0); i < vertices.size(); ++i) {
                Vector3 vertex = vertices[i];

                Vector3 d = (vertex - Vector3(0, 0, 0)).unit();

                float nx = image->width() * (0.5f + atanf(d.z / d.x) / (2.0f*pif()));
                float ny = image->height() * (0.5f - asinf(d.y) * 1 / pif());

                int ix = (int)abs((int)nx % image->width());
                int iy = (int)abs((int)ny % image->height());

                Color3 color = Color3();

                image->get(Point2int32(ix, iy), color);

                float bump = color.average();
                if (bump > 0.3f && bump < 0.6f) bump = 0.5f;
                else if (bump < 0.3f) bump = 0.1f;
                else bump = 1.0f;

                vertices[i] += vertex.unit() * bump * radius;
            }*/