#pragma once
/**
  \file SolarSystem.h

  Container Class for the Solar System being created.
 */

#include <G3D/G3DAll.h>
#include "Mesh.h"
#include "Planet.h"

/** \brief Application framework. */
class SolarSystem{
protected:
    Table<String, Planet> m_planetTable;
    Any m_scene;
    String m_systemName;
    Any m_models;
    Any m_entities;

    void addPlanetToScene(Any& entities, Any& models, const String& name, Planet& planet);

    void initializeEntityTable(Any& entities);
    void initializeModelsTable(Any& models);
    void initializeSceneTable(Any& scene);
    void makeSceneTable(Any& scene, const Any& models, const Any& entities, const String& name);

public:

    SolarSystem();

    bool printSolarSystemToScene(const String& save);
    bool addPlanet(const String& name, Planet& planet);
    bool containsPlanet(const String& name);
    bool removePlanet(const String &name);
    void onInit();
};