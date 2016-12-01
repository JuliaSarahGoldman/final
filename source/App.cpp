/** \file App.cpp */
#include "App.h"
#include "NoiseGen.h"
#include "Planet.h"
#include "Mesh.h"
#include "SimpleMesh.h"

// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();

int main(int argc, const char* argv[]) {
    {
        G3DSpecification g3dSpec;
        g3dSpec.audio = false;
        initGLG3D(g3dSpec);
    }

    GApp::Settings settings(argc, argv);

    // Change the window and other startup parameters by modifying the
    // settings class.  For example:
    settings.window.caption = argv[0];

    // Set enable to catch more OpenGL errors
    // settings.window.debugContext     = true;

    // Some common resolutions:
    // settings.window.width            =  854; settings.window.height       = 480;
    // settings.window.width            = 1024; settings.window.height       = 768;
    settings.window.width = 1280; settings.window.height = 720;
    //settings.window.width             = 1920; settings.window.height       = 1080;
    // settings.window.width            = OSWindow::primaryDisplayWindowSize().x; settings.window.height = OSWindow::primaryDisplayWindowSize().y;
    settings.window.fullScreen = false;
    settings.window.resizable = !settings.window.fullScreen;
    settings.window.framed = !settings.window.fullScreen;

    // Set to true for a significant performance boost if your app can't render at 60fps, or if
    // you *want* to render faster than the display.
    settings.window.asynchronous = false;

    settings.hdrFramebuffer.depthGuardBandThickness = Vector2int16(64, 64);
    settings.hdrFramebuffer.colorGuardBandThickness = Vector2int16(0, 0);
    settings.dataDir = FileSystem::currentDirectory();
    settings.screenshotDirectory = "../journal/";

    settings.renderer.deferredShading = true;
    settings.renderer.orderIndependentTransparency = false;

    return App(settings).run();
}


App::App(const GApp::Settings& settings) : GApp(settings) {
}

// Called before the application loop begins.  Load data here and
// not in the constructor so that common exceptions will be
// automatically caught.
void App::onInit() {
    debugPrintf("Target frame rate = %f Hz\n", realTimeTargetDuration());
    GApp::onInit();
    setFrameDuration(1.0f / 120.0f);

    // Call setScene(shared_ptr<Scene>()) or setScene(MyScene::create()) to replace
    // the default scene here.

    showRenderingStats = false;

    makeGUI();

    // For higher-quality screenshots:
    // developerWindow->videoRecordDialog->setScreenShotFormat("PNG");
    // developerWindow->videoRecordDialog->setCaptureGui(false);
    developerWindow->cameraControlWindow->moveTo(Point2(developerWindow->cameraControlWindow->rect().x0(), 0));
    loadScene(
        //"G3D Sponza"
        "Ground" // Load something simple
        //developerWindow->sceneEditorWindow->selectedSceneName()  // Load the first scene encountered 
    );
}

void App::makePentagon() {
    Array<Vector3> vertices = Array<Vector3>();
    Array<Vector3int32> faces = Array<Vector3int32>();
    Array<Vector3> vertices1(Vector3(0, 0, -2), Vector3(0, 2, -2), Vector3(-2, 0, -2), Vector3(-1, -2, -2), Vector3(1, -2, -2), Vector3(2, 0, -2));
    Array<Vector3int32> indices1(Vector3int32(0, 1, 2), Vector3int32(0, 2, 3), Vector3int32(0, 3, 4), Vector3int32(0, 4, 5), Vector3int32(0, 5, 1));
    //Planet planet;
    //planet.writeSphere("pentagon", 12.0f, 5, vertices, faces);
    //m_myMesh = Mesh::create(vertices, faces);

    m_myMesh = Mesh::create(vertices1, indices1);
    m_myMesh->toObj("pentagon");

    loadScene("Pentagon");
    GuiPane* pentPane = debugPane->addPane("Info", GuiTheme::ORNATE_PANE_STYLE);

    // Example of how to add debugging controls
    pentPane->addLabel("Pentagon - Collapse Edges");
    pentPane->addNumberBox("# edges", &m_edgesToCollapse, "",
        GuiTheme::LINEAR_SLIDER, 0, 20000)->setUnitsSize(1);

    pentPane->addButton("Collapse!", [this]() {
        m_myMesh->collapseEdges(m_edgesToCollapse);
        m_myMesh->toObj("pentagon");
        G3D::ArticulatedModel::clearCache();
        loadScene("Pentagon");

    });

    pentPane->addButton("Reset", [this]() {
        Array<Vector3> vertices1(Vector3(0, 0, -2), Vector3(0, 2, -2), Vector3(-2, 0, -2), Vector3(-1, -2, -2), Vector3(1, -2, -2), Vector3(2, 0, -2));
        Array<Vector3int32> indices1(Vector3int32(0, 1, 2), Vector3int32(0, 2, 3), Vector3int32(0, 3, 4), Vector3int32(0, 4, 5), Vector3int32(0, 5, 1));
        m_myMesh = Mesh::create(vertices1, indices1);
        m_myMesh->toObj("pentagon");
        G3D::ArticulatedModel::clearCache();
        loadScene("Little Heightfield");

    });

    pentPane->pack();
}

void App::makeBunny() {
    m_myMesh = Mesh::create("bunny.ifs");
    m_myMesh->toObj("bunny");

    loadScene("Bunny");
    GuiPane* bunnyPane = debugPane->addPane("Collapse Edges", GuiTheme::ORNATE_PANE_STYLE);

    // Example of how to add debugging controls
    bunnyPane->addLabel("Collapse Edges");

    bunnyPane->addNumberBox("# edges", &m_edgesToCollapse, "",
        GuiTheme::LINEAR_SLIDER, 0, 20000)->setUnitsSize(1);

    bunnyPane->addButton("Collapse!", [this]() {
        m_myMesh->collapseEdges(m_edgesToCollapse);
        m_myMesh->toObj("bunny");
        G3D::ArticulatedModel::clearCache();
        loadScene("Bunny");

    });

    bunnyPane->addButton("Reset", [this]() {
        m_myMesh = Mesh::create("bunny.ifs");
        m_myMesh->toObj("bunny");
        G3D::ArticulatedModel::clearCache();
        loadScene("Little Heightfield");

    });

    bunnyPane->pack();
}

void App::makeLittleHeightfield() {
    Array<Vector3> vertices;
    Array<Vector3int32> indices;

    shared_ptr<Image>base(Image::fromFile("heightmap.png"));
    for (int x = 0; x < base->width(); ++x) {
        int i(x*base->height());

        for (int z = 0; z < base->height(); ++z) {
            Color3 color;
            base->get(Point2int32(x, z), color);
            float y(color.average());
            vertices.append(Vector3(x, y * 5, z));

            if (x < base->width() - 1 && z < base->height() - 1) {
                int c(i + z);
                int h(base->height());
                indices.append(Vector3int32(c, c + 1, c + h), Vector3int32(c + 1, c + h + 1, c + h));
            }
        }
    }

    m_myMesh = Mesh::create(vertices, indices);
    m_myMesh->toObj("littleHf");

    loadScene("Little Heightfield");
    GuiPane* hfPane = debugPane->addPane("Collapse Edges", GuiTheme::ORNATE_PANE_STYLE);

    // Example of how to add debugging controls
    hfPane->addLabel("Collapse Edges");

    hfPane->addNumberBox("# edges", &m_edgesToCollapse, "",
        GuiTheme::LINEAR_SLIDER, 0, 1000)->setUnitsSize(1);

    hfPane->addButton("Collapse!", [this]() {
        m_myMesh->collapseEdges(m_edgesToCollapse);
        m_myMesh->toObj("littleHf");
        G3D::ArticulatedModel::clearCache();
        loadScene("Little Heightfield");

    });

    hfPane->addButton("Reset", [this]() {
        Array<Vector3> vertices;
        Array<Vector3int32> indices;

        shared_ptr<Image>base(Image::fromFile("heightmap.png"));
        for (int x = 0; x < base->width(); ++x) {
            int i(x*base->height());

            for (int z = 0; z < base->height(); ++z) {
                Color3 color;
                base->get(Point2int32(x, z), color);
                float y(color.average());
                vertices.append(Vector3(x, y * 5, z));

                if (x < base->width() - 1 && z < base->height() - 1) {
                    int c(i + z);
                    int h(base->height());
                    indices.append(Vector3int32(c, c + 1, c + h), Vector3int32(c + 1, c + h + 1, c + h));
                }
            }
        }
        m_myMesh = Mesh::create(vertices, indices);
        m_myMesh->toObj("littleHf");
        G3D::ArticulatedModel::clearCache();
        loadScene("Little Heightfield");

    });

    hfPane->pack();
}

void App::addPlanetToScene(Mesh& mesh, String name, Point3& position, String filename, Matrix3& rotation) {
    String anyStr = "UniversalMaterial::Specification { lambertian = \"" + filename + "\"; }";
    shared_ptr<Image> im(Image::fromFile(filename));
    addPlanetToScene(mesh, name, position, anyStr, im->width(), im->height(), rotation);
}

void App::addPlanetToScene(Mesh& mesh, String name, Point3& position, Color3& color, Matrix3& rotation) {
    String anyStr = "UniversalMaterial::Specification { lambertian = Color3(" + (String)std::to_string(color.r) + ", " + (String)std::to_string(color.g) + ", " + (String)std::to_string(color.b) + "); }";
    addPlanetToScene(mesh, name, position, anyStr, 1, 1, rotation);
}

void App::addPlanetToScene(Mesh& mesh, String name, Point3& position, String anyStr, int width, int height, Matrix3& rotation) {
    // Replace any existing torus model. Models don't 
    // have to be added to the model table to use them 
    // with a VisibleEntity.
    const shared_ptr<Model>& planetModel = mesh.toArticulatedModel(name + "Model", anyStr, width, height);
    if (scene()->modelTable().containsKey(planetModel->name())) {
        scene()->removeModel(planetModel->name());
    }
    scene()->insert(planetModel);

    // Replace any existing planet entity that has the wrong type
    shared_ptr<Entity> planet = scene()->entity(name);
    if (notNull(planet) && isNull(dynamic_pointer_cast<VisibleEntity>(planet))) {
        logPrintf("The scene contained an Entity named %s that was not a VisibleEntity\n", name);
        scene()->remove(planet);
        planet.reset();
    }

    if (isNull(planet)) {
        // There is no torus entity in this scene, so make one.
        //
        // We could either explicitly instantiate a VisibleEntity or simply
        // allow the Scene parser to construct one. The second approach
        // has more consise syntax for this case, since we are using all constant
        // values in the specification.
        //planet = scene()->createEntity();
        String anyStr("VisibleEntity { model = \"" + name + "Model\"; };");
        //frame = CFrame::fromXYZYPRDegrees(" + (String) std::to_string(position.x) +", " + (String) std::to_string(position.y) + ", " + (String) std::to_string(position.z) + ", 0,0,0);
        Any any = Any::parse(anyStr);
        planet = scene()->createEntity(name, any);
        Matrix3 startMat = Matrix3::identity();
        PhysicsFrame frame(startMat);
        PhysicsFrameSpline spline(frame.toAny());
        spline.extrapolationMode = SplineExtrapolationMode::CYCLIC;
        spline.interpolationMode = SplineInterpolationMode::LINEAR;

        for(int i(0); i < 128; i++){
            startMat *= rotation;
            spline.append(PhysicsFrame(startMat));
        }
        planet->setFrameSpline(spline);
    }
    else {
        // Change the model on the existing planet entity
        dynamic_pointer_cast<VisibleEntity>(planet)->setModel(planetModel);
    }

    planet->setFrame(CFrame::fromXYZYPRDegrees(position.x, position.y, position.z, 45.0f, 45.0f));
}

void App::makePlanetGUI() {

    GuiPane* planetPane = debugPane->addPane("Planet");

    planetPane->setNewChildSize(240);
    planetPane->addNumberBox("Recursion Level", &m_recursionLevel, "",
        GuiTheme::LINEAR_SLIDER, 1, 8)->setUnitsSize(1);

    planetPane->addNumberBox("Land Bevel", &m_landBevel, "",
        GuiTheme::LOG_SLIDER, 0.0001f, 1.0f)->setUnitsSize(1);

    planetPane->addNumberBox("Mount Bevel", &m_mountainBevel, "",
        GuiTheme::LOG_SLIDER, 0.0001f, 1.0f)->setUnitsSize(1);

    planetPane->addNumberBox("Mount Height", &m_mountainHeight, "",
        GuiTheme::LOG_SLIDER, 10.0f, 200.0f)->setUnitsSize(1);

    planetPane->addNumberBox("Mount Diversity", &m_mountianDiversity, "",
        GuiTheme::LOG_SLIDER, 1.0f, 10.0f)->setUnitsSize(1);

    planetPane->addNumberBox("Ocean Level", &m_oceanLevel, "",
        GuiTheme::LOG_SLIDER, 0.0f, 1.0f)->setUnitsSize(1);

    planetPane->addNumberBox("Ocean Noise", &m_oceanNoise, "",
        GuiTheme::LOG_SLIDER, 0.0f, 10.0f)->setUnitsSize(1);

    planetPane->addCheckBox("Water Mount", &m_waterMount);

    planetPane->addNumberBox("Land Noise", &m_landNoise, "",
        GuiTheme::LOG_SLIDER, -1.0f, 10.0f)->setUnitsSize(1);

    planetPane->addTextBox("Save to:", &m_planetSave);
    
   /* heightfieldPane->beginRow(); {


    /*planetPane->addNumberBox("Frequency", &m_frequency, "",
        GuiTheme::LOG_SLIDER, 0.0001f, 1.0f)->setUnitsSize(1);
    /*
    heightfieldPane->beginRow(); {
        heightfieldPane->addTextBox("Input Image", &m_heightfieldSource)->setWidth(210);
        heightfieldPane->addButton("...", [this]() {
            FileDialog::getFilename(m_heightfieldSource, "png", false);
        })->setWidth(30);
    } heightfieldPane->endRow();*/

    planetPane->addButton("Get Constants from Any File", [this]() {
        FileDialog::getFilename(m_planetSource, "", false);
        try{
            Any any(Any::TABLE, "Planet");
            any.load(m_planetSource);
            AnyTableReader x(any);
            x.getIfPresent("recursions", m_recursionLevel);
            x.getIfPresent("landBevel", m_landBevel);
            x.getIfPresent("mountainBevel", m_mountainBevel);
            x.getIfPresent("mountainHeight", m_mountainHeight);
            x.getIfPresent("mountainDiversity", m_mountianDiversity);
            x.getIfPresent("oceanLevel", m_oceanLevel);
            x.getIfPresent("landNoise", m_landNoise);
            x.getIfPresent("oceanNoise", m_oceanNoise);
            x.getIfPresent("waterMount", m_waterMount);
        }
        catch (...) {
            msgBox("Unable to load the image.");
        }
    });

    planetPane->addButton("Save Constants to Any File", [this]() {
        try{
            Any x(Any::TABLE, "Planet");
            x["recursions"] =  m_recursionLevel;
            x["landBevel"] = m_landBevel;
            x["mountainBevel"] = m_mountainBevel;
            x["mountainHeight"] = m_mountainHeight;
            x["mountainDiversity"] = m_mountianDiversity;
            x["oceanLevel"] = m_oceanLevel;
            x["landNoise"] = m_landNoise;
            x["oceanNoise"] = m_oceanNoise;
            x["waterMount"] = m_waterMount;
            x.save(m_planetSave);
        }
        catch (...) {
            msgBox("Unable to save the image.");
        }
    });

    planetPane->addButton("Generate", [this]() {
        shared_ptr<G3D::Image> image;
        //Noise noise = G3D::Noise::common();
        try {
            Planet planet;
            NoiseGen noise;

            Array<Vector3> vertices = Array<Vector3>();
            Array<Vector3int32> faces = Array<Vector3int32>();
            float freq = 3.0f;
            shared_ptr<Image> image = Image::create(1024, 1024, ImageFormat::RGBA8());
            noise.generateSeaImage(image, m_oceanNoise);
            //noise.generateSeaImage(image, freq);
            image->save("ocean.png");
            planet.writeSphere("ocean", 12.5f, m_recursionLevel, vertices, faces);
            planet.applyNoiseWater(vertices, image);
            Mesh mesh(vertices, faces);
            mesh.bevelEdges2(0.1f);

            vertices = Array<Vector3>();
            faces = Array<Vector3int32>();
            freq = 0.25f;
            image = Image::create(1024, 1024, ImageFormat::RGBA8());
            shared_ptr<Image> test = Image::create(1024, 1024, ImageFormat::RGBA8());
            noise.generateLandImage(image, m_landNoise);
            //noise.generateLandImage(image, freq);
            image->save("land.png");
            planet.writeSphere("land", 12.0f, m_recursionLevel, vertices, faces);
            planet.applyNoiseLand(vertices, image, test, m_oceanLevel);
            test->save("test.png");
            Mesh mesh2(vertices, faces);
            mesh2.bevelEdges2(m_landBevel);

            vertices = Array<Vector3>();
            faces = Array<Vector3int32>();
            planet.writeSphere("mountain", 11.5f, m_recursionLevel, vertices, faces);
            image = Image::create(1024, 1024, ImageFormat::RGBA8());
            noise.generateMountainImage(image, 0.125f, 1.0f);
            noise.generateMountainImage(image, 0.25f, 0.5f);
            noise.generateMountainImage(image, 0.5f, 0.25f);
            planet.applyNoiseMountain(vertices, image, test, m_waterMount, m_mountianDiversity, m_mountainHeight);
            Mesh mesh3(vertices, faces);
            mesh3.bevelEdges2(m_mountainBevel);

            image->save("image.png");

            loadScene("Ground");


            /*String material1 = "UniversalMaterial::Specification { lambertian = \"space.png\"; }";
            addPlanetToScene(mesh, "ocean", Point3(0, 0, 0), material1, 1000, 1000);
            String material2 = "UniversalMaterial::Specification { lambertian = \"texture2.jpg\"; }";
            addPlanetToScene(mesh2, "land", Point3(0, 0, 0), material2, 1000, 1000);
            String material3 = "UniversalMaterial::Specification { lambertian = \"texture.jpg\"; }";
            addPlanetToScene(mesh3, "mountain", Point3(0, 0, 0), material3, 1000, 1000);*/
                        
            float cosrot = cos(pif()/4.0f);
            float sinrot = sin(pif()/4.0f);

            Matrix3 rotation(cosrot, 0, sinrot, 0, 1, 0, -sinrot, 0, cosrot);

            cosrot = cos(pif()/5.0f);
            sinrot = sin(pif()/5.0f);

            Matrix3 waterRotation(cosrot, 0, sinrot, 0, 1, 0, -sinrot, 0, cosrot); 

            addPlanetToScene(mesh, "ocean", Point3(100, 0, 0), Color3(0, 0, 1), waterRotation);
            addPlanetToScene(mesh2, "land", Point3(100, 0, 0), Color3(0, 1, 0), rotation);
            addPlanetToScene(mesh3, "mountain", Point3(100, 0, 0), Color3(.5, .5, .5), rotation);

            addPlanetToScene(mesh, "lava", Point3(0, 0, 0), "orange.jpg", waterRotation);
            addPlanetToScene(mesh2, "twin", Point3(0, 0, 0), "deep_red.jpg", rotation);
            addPlanetToScene(mesh3, "mount", Point3(0, 0, 0), "texture2.jpg", rotation);

            /*addPlanetToScene(mesh, "ocean", Point3(0, 0, 0), "orange.jpg");
            addPlanetToScene(mesh2, "land", Point3(0, 0, 0), "deep_red.jpg");
            addPlanetToScene(mesh3, "mountain", Point3(0, 0, 0), "texture2.jpg");*/

            //loadScene("Planet");

            //Julia's Generate Code
            /*Planet planet;
            Array<Vector3> vertices1 = Array<Vector3>();
            Array<Vector3int32> faces1 = Array<Vector3int32>();
            planet.writeSphere("mountains", 10.1f, 5, vertices1, faces1);
            Mesh mesh(vertices1, faces1);
            mesh.bevelEdges2(.1);
            loadScene("Ground");

            addPlanetToScene(mesh, "sphere", Point3(-4,-4,-4), Color3(1,0,0));
            addPlanetToScene(mesh, "sphere2", Point3(10,12,10), Color3(0,0,1));*/
        }
        catch (...) {
            msgBox("Unable to load the image.", m_heightfieldSource);
        }
    });
}

//Creates a GUI that allows a user to generate a heightfield with a given xz and y scale based on a given image
void App::makeHeightfield() {

    GuiPane* heightfieldPane = debugPane->addPane("Heightfield");

    heightfieldPane->setNewChildSize(240);
    heightfieldPane->addNumberBox("Max Y", &m_heightfieldYScale, "m",
        GuiTheme::LOG_SLIDER, 0.0f, 100.0f)->setUnitsSize(30);

    heightfieldPane->addNumberBox("XZ Scale", &m_heightfieldXZScale, "m/px",
        GuiTheme::LOG_SLIDER, 0.001f, 10.0f)->setUnitsSize(30);

    heightfieldPane->beginRow(); {
        heightfieldPane->addTextBox("Input Image", &m_heightfieldSource)->setWidth(210);
        heightfieldPane->addButton("...", [this]() {
            FileDialog::getFilename(m_heightfieldSource, "png", false);
        })->setWidth(30);
    } heightfieldPane->endRow();

    heightfieldPane->addButton("Generate", [this]() {
        shared_ptr<G3D::Image> image;
        Noise noise = G3D::Noise::common();
        try {
            drawMessage("Generating Heightfield.");
            Noise n;
            image = Image::create(640, 640, ImageFormat::RGBA8());
            for (int y = 0; y < image->height(); ++y) {
                for (int x = 0; x < image->width(); ++x) {
                    image->set(x, y, lerp(Color3(0.2f, 0.3f, 0.7f), Color3(1.0f), noise.sampleFloat(x, y, x + y, 3)));
                }
            }
            image->save("../data-files/noise.png");

            TextOutput output("model/heightfield.off");
            output.writeSymbol("OFF\n");
            output.printf("%d %d 0\n", image->width() * image->height(), (image->width() - 1) * (image->height() - 1));

            for (int x = 0; x < image->width(); ++x) {
                for (int z = 0; z < image->height(); ++z) {
                    Color3 color;
                    image->get(Point2int32(x, z), color);
                    float y = color.average();
                    output.printf("%f %f %f\n", ((float)x)*m_heightfieldXZScale, y*m_heightfieldYScale, ((float)z)*m_heightfieldXZScale);
                }
            }

            for (int i = 1; i < image->height(); ++i) {
                for (int j = 1; j < image->width(); ++j) {
                    output.printf("4 %d %d %d %d\n", i + ((image->height())*j), i + ((image->height())*j) - 1, i + ((image->height())*(j - 1)) - 1, i + ((image->height())*(j - 1)));
                }
            }

            output.commit(true);
            G3D::ArticulatedModel::clearCache();
        }
        catch (...) {
            msgBox("Unable to load the image.", m_heightfieldSource);
        }
    });
}

void App::makeGUI() {
    // Initialize the developer HUD
    createDeveloperHUD();

    debugWindow->setVisible(true);
    developerWindow->videoRecordDialog->setEnabled(true);

    /*GuiPane* infoPane = debugPane->addPane("Info", GuiTheme::ORNATE_PANE_STYLE);

    // Example of how to add debugging controls
    infoPane->addLabel("You can add GUI controls");
    infoPane->addLabel("in App::onInit().");
    infoPane->addButton("Exit", [this]() { m_endProgram = true; });
    infoPane->pack();*/

    /*Planet planet;
    shared_ptr<Array<Vector3>> vertices = std::make_shared<Array<Vector3>>();
    shared_ptr<Array<Vector3int32>> faces = std::make_shared<Array<Vector3int32>>();
    planet.writeSphere("water.obj", 9.9f, 3, vertices, faces);

    vertices = std::make_shared<Array<Vector3>>();
    faces = std::make_shared<Array<Vector3int32>>();
    planet.writeSphere("land.obj", 10.0f, 5, vertices, faces);


    vertices = std::make_shared<Array<Vector3>>();
    faces = std::make_shared<Array<Vector3int32>>();
    planet.writeSphere("mountains.obj", 10.1f, 5, vertices, faces);*/

    //makeHeightfield();
    makePlanetGUI();

    //makeBunny();
    //makeLittleHeightfield();
    //makePentagon();

    Array<Vector3> verticeArray(Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(.5, 0, 1), Vector3(.5, 1, .5));
    Array<Vector3int32> triangles(Vector3int32(3, 1, 0), Vector3int32(1, 2, 0), Vector3int32(3, 2, 1), Vector3int32(3, 0, 2));

    //Mesh mesh(*vertices, *faces);
    //Mesh mesh(verticeArray, triangles);
    //mesh.bevelEdges(.3);
    //mesh.toObj("wtf.obj");

    //loadScene("Ground");
    //mesh.toArticulatedModel("test");
    //addPlanetToScene(mesh);

    debugWindow->pack();
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}


// This default implementation is a direct copy of GApp::onGraphics3D to make it easy
// for you to modify. If you aren't changing the hardware rendering strategy, you can
// delete this override entirely.
void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& allSurfaces) {
    if (!scene()) {
        if ((submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) && (!rd->swapBuffersAutomatically())) {
            swapBuffers();
        }
        rd->clear();
        rd->pushState(); {
            rd->setProjectionAndCameraMatrix(activeCamera()->projection(), activeCamera()->frame());
            drawDebugShapes();
        } rd->popState();
        return;
    }

    GBuffer::Specification gbufferSpec = m_gbufferSpecification;
    extendGBufferSpecification(gbufferSpec);
    m_gbuffer->setSpecification(gbufferSpec);
    m_gbuffer->resize(m_framebuffer->width(), m_framebuffer->height());
    m_gbuffer->prepare(rd, activeCamera(), 0, -(float)previousSimTimeStep(), m_settings.hdrFramebuffer.depthGuardBandThickness, m_settings.hdrFramebuffer.colorGuardBandThickness);

    m_renderer->render(rd, m_framebuffer, scene()->lightingEnvironment().ambientOcclusionSettings.enabled ? m_depthPeelFramebuffer : shared_ptr<Framebuffer>(),
        scene()->lightingEnvironment(), m_gbuffer, allSurfaces);

    // Debug visualizations and post-process effects
    rd->pushState(m_framebuffer); {
        // Call to make the App show the output of debugDraw(...)
        rd->setProjectionAndCameraMatrix(activeCamera()->projection(), activeCamera()->frame());
        drawDebugShapes();
        const shared_ptr<Entity>& selectedEntity = (notNull(developerWindow) && notNull(developerWindow->sceneEditorWindow)) ? developerWindow->sceneEditorWindow->selectedEntity() : shared_ptr<Entity>();
        scene()->visualize(rd, selectedEntity, allSurfaces, sceneVisualizationSettings(), activeCamera());

        // Post-process special effects
        m_depthOfField->apply(rd, m_framebuffer->texture(0), m_framebuffer->texture(Framebuffer::DEPTH), activeCamera(), m_settings.hdrFramebuffer.depthGuardBandThickness - m_settings.hdrFramebuffer.colorGuardBandThickness);

        m_motionBlur->apply(rd, m_framebuffer->texture(0), m_gbuffer->texture(GBuffer::Field::SS_EXPRESSIVE_MOTION),
            m_framebuffer->texture(Framebuffer::DEPTH), activeCamera(),
            m_settings.hdrFramebuffer.depthGuardBandThickness - m_settings.hdrFramebuffer.colorGuardBandThickness);
    } rd->popState();

    // We're about to render to the actual back buffer, so swap the buffers now.
    // This call also allows the screenshot and video recording to capture the
    // previous frame just before it is displayed.
    if (submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) {
        swapBuffers();
    }

    // Clear the entire screen (needed even though we'll render over it, since
    // AFR uses clear() to detect that the buffer is not re-used.)
    rd->clear();

    // Perform gamma correction, bloom, and SSAA, and write to the native window frame buffer
    m_film->exposeAndRender(rd, activeCamera()->filmSettings(), m_framebuffer->texture(0), settings().hdrFramebuffer.colorGuardBandThickness.x + settings().hdrFramebuffer.depthGuardBandThickness.x, settings().hdrFramebuffer.depthGuardBandThickness.x);
}


void App::onAI() {
    GApp::onAI();
    // Add non-simulation game logic and AI code here
}


void App::onNetwork() {
    GApp::onNetwork();
    // Poll net messages here
}


void App::onSimulation(RealTime rdt, SimTime sdt, SimTime idt) {
    GApp::onSimulation(rdt, sdt, idt);

    // Example GUI dynamic layout code.  Resize the debugWindow to fill
    // the screen horizontally.
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}


bool App::onEvent(const GEvent& event) {
    // Handle super-class events
    if (GApp::onEvent(event)) { return true; }

    // If you need to track individual UI events, manage them here.
    // Return true if you want to prevent other parts of the system
    // from observing this specific event.
    //
    // For example,
    // if ((event.type == GEventType::GUI_ACTION) && (event.gui.control == m_button)) { ... return true; }
    // if ((event.type == GEventType::KEY_DOWN) && (event.key.keysym.sym == GKey::TAB)) { ... return true; }
    // if ((event.type == GEventType::KEY_DOWN) && (event.key.keysym.sym == 'p')) { ... return true; }

    return false;
}


void App::onUserInput(UserInput* ui) {
    GApp::onUserInput(ui);
    (void)ui;
    // Add key handling here based on the keys currently held or
    // ones that changed in the last frame.
}


void App::onPose(Array<shared_ptr<Surface> >& surface, Array<shared_ptr<Surface2D> >& surface2D) {
    GApp::onPose(surface, surface2D);

    // Append any models to the arrays that you want to later be rendered by onGraphics()
}


void App::onGraphics2D(RenderDevice* rd, Array<shared_ptr<Surface2D> >& posed2D) {
    // Render 2D objects like Widgets.  These do not receive tone mapping or gamma correction.
    Surface2D::sortAndRender(rd, posed2D);
}


void App::onCleanup() {
    // Called after the application loop ends.  Place a majority of cleanup code
    // here instead of in the constructor so that exceptions can be caught.
}

        /*planet = scene()->createEntity("planet",
            PARSE_ANY(
                VisibleEntity {
                    model = "planetModel";
                };
            ));*/