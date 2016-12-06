#include "SolarSystem.h"
#include "Planet.h"
#include "Mesh.h"
#include "NoiseGen.h"

SolarSystem::SolarSystem() {

}

void SolarSystem::onInit() {
    m_planetTable = Table<String, Planet>();
    m_scene = Any(Any::TABLE);
    m_entities = Any(Any::TABLE);
    m_models = Any(Any::TABLE);

    initializeSceneTable(m_scene);
    initializeEntityTable(m_entities);
    initializeModelsTable(m_models);
}

bool SolarSystem::addPlanet(const String& name, Planet& planet) {
    if (m_planetTable.containsKey(name)) return false;
    m_planetTable.set(name, planet);
    String objectToOrbit;
    float orbitDistance;
    planet.getPlanetOrbit(objectToOrbit, orbitDistance);

    addPlanetToScene(m_entities, m_models, name, planet);

    return true;
}

void SolarSystem::addPlanetToScene(Any& entities, Any& models, const String& name, Planet& planet) {
    Any waterModelDescription(Any::TABLE, "ArticulatedModel::Specification");
    Any waterEntityDescription(Any::TABLE, "VisibleEntity");

    planet.createWaterAnyFile(waterModelDescription, waterEntityDescription);
    waterEntityDescription["model"] = name;
    entities[name] = waterEntityDescription;
    models[name] = waterModelDescription;

    Any landModelDescription(Any::TABLE, "ArticulatedModel::Specification");
    Any landEntityDescription(Any::TABLE, "VisibleEntity");

    String levelName = name + "land";
    planet.createLandAnyFile(landModelDescription, landEntityDescription, name);
    landEntityDescription["model"] = levelName;
    entities[levelName] = landEntityDescription;
    models[levelName] = landModelDescription;

    Any mountainModelDescription(Any::TABLE, "ArticulatedModel::Specification");
    Any mountainEntityDescription(Any::TABLE, "VisibleEntity");

    levelName = name + "mountain";
    planet.createMountainAnyFile(mountainModelDescription, mountainEntityDescription, name);
    mountainEntityDescription["model"] = levelName;
    entities[levelName] = mountainEntityDescription;
    models[levelName] = mountainModelDescription;

    String preprocess = "{setMaterial(all(), UniversalMaterial::Specification{ lambertian = Color4(Color3(0.8), 1.0); emissive = Color3(0.1); } ); }";
    Any cloudModel(Any::TABLE, "ArticulatedModel::Specification");
    cloudModel["scale"] = 0.1 * planet.getScale();
    cloudModel["filename"] = "model/cloud/cloud.zip/altostratus00.obj";
    cloudModel["preprocess"] = Any::parse(preprocess);

    models["cloud1"] = cloudModel;

    for (int i(0); i < 10; ++i) {
        Any cloudEntityDescription(Any::TABLE, "VisibleEntity");
        planet.addCloudToPlanet(cloudEntityDescription, "cloud1", name, planet.getPosition(), planet.getScale());
        entities["cloud" + (String) std::to_string(i)] = cloudEntityDescription;
    }

    makeSceneTable(m_scene, models, entities, name);
}

void SolarSystem::makeSceneTable(Any& scene, const Any& models, const Any& entities, const String& name) {
    scene["entities"] = entities;
    scene["models"] = models;
}

void SolarSystem::initializeSceneTable(Any& scene) {

    scene["name"] = "Solar System";

    Any lightingEnvironment(Any::TABLE, "LightingEnvironment");

    String occlusionSettings = (String) "AmbientOcclusionSettings { bias = 0.12; blurRadius = 4; blurStepSize = 2;" +
        "depthPeelSeparationHint = 0.2; edgeSharpness = 1; enabled = false; intensity = 3; " +
        "monotonicallyDecreasingBilateralWeights = false; numSamples = 19; radius = 2; " +
        "temporalFilterSettings = TemporalFilter::Settings { hysteresis = 0.9; falloffEndDistance = 0.07; " +
        "falloffStartDistance = 0.05; }; temporallyVarySamples = true; useDepthPeelBuffer = true; useNormalBuffer = true; " +
        "useNormalsInBlur = true;}";

    lightingEnvironment["ambientOcclusionSettings"] = Any::parse(occlusionSettings);

    String environmentMapSettings = (String) "Texture::Specification { encoding = Texture::Encoding { readMultiplyFirst = 1.9; " +
        "}; filename = \"cubemap/majestic/majestic512_*.jpg\"; }";
    //\"cubemap/hipshot_m9_sky/16_*.png\"
    lightingEnvironment["environmentMap"] = Any::parse(environmentMapSettings);

    scene["lightingEnvironment"] = lightingEnvironment;
}

void SolarSystem::initializeEntityTable(Any& entities) {
    //Create the light source
    Any light(Any::TABLE, "Light");
    light["attenuation"] = Vector3(0, 0, 1);
    light["bulbPower"] = Power3(1e+006);
    light["castsShadows"] = true;
    light["shadowMapBias"] = 0.05f;
    light["track"] = Any::parse("lookAt(Point3(-15, 200, 40), Point3(0, 0, 0));");
    light["shadowMapSize"] = Vector2int16(2048, 2048);
    light["spotHalfAngleDegrees"] = 8;
    light["spotSquare"] = true;
    light["type"] = "SPOT";
    entities["sun"] = light;

    Any camera(Any::TABLE, "Camera");
    camera["frame"] = CFrame::fromXYZYPRDegrees(50, 50, 50);
    camera["projection"] = Any::parse("Projection { farPlaneZ = -inf; fovDegrees = 50; fovDirection = \"VERTICAL\"; nearPlaneZ = -0.1; }");

    String filmSettings = (String) "FilmSettings{" +
        "antialiasingEnabled = true; antialiasingFilterRadius = 0; antialiasingHighQuality = true; bloomRadiusFraction = 0.009;"
        + "bloomStrength = 0.2; debugZoom = 1; gamma = 2.2; sensitivity = 1; toneCurve = \"CELLULOID\";" +
        "vignetteBottomStrength = 0.05; vignetteSizeFraction = 0.17; vignetteTopStrength = 0.5; }";
    camera["filmSettings"] = Any::parse(filmSettings);
    entities["camera"] = camera;

    Any skybox(Any::TABLE, "Skybox");
    skybox["texture"] = "cubemap/hipshot_m9_sky/16_*.png";
    entities["skybox"] = skybox;
}

void SolarSystem::initializeModelsTable(Any& models) {
    /*
    Any cloud(Any::TABLE, "ParticleSystemModel::Emitter::Specification");
    cloud["material"] = "material/smoke/smoke.png";
    cloud["initialDensity"] = 8;
    cloud["radiusMean"] = 3;
    cloud["radiusVariance"] = 0.5;
    cloud["noisePower"] = 0;
    cloud["angularVelocityMean"] = 0.5;
    cloud["angularVelocityVariance"] = 0.25;

    Any cloud2 = cloud;
    Any cloud3 = cloud;

    Any cloudShape(Any::TABLE, "ArticulatedModel::Specification");
    cloudShape["scale"] = 0.05;

    Any cloudShape2 = cloudShape;
    Any cloudShape3 = cloudShape;

    cloudShape["filename"] = "model/cloud/cloud.zip/cumulus02.obj";
    cloudShape2["filename"] = "model/cloud/cloud.zip/cumulus01.obj";
    cloudShape3["filename"] = "model/cloud/cloud.zip/cumulus00.obj";

    cloud["shape"] = cloudShape;
    cloud2["shape"] = cloudShape2;
    cloud3["shape"] = cloudShape3;
    */
}

bool SolarSystem::containsPlanet(const String& name){
    return m_planetTable.containsKey(name);
}

bool SolarSystem::removePlanet(const String &name){
    return m_planetTable.remove(name);
}

bool SolarSystem::printSolarSystemToScene(const String& save){
    try {
        m_scene.save(save + ".Scene.Any");
        return true;
    } catch(...) {
        return false;
    }
}