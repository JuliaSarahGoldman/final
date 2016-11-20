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


void App::addPlanetToScene(Mesh& mesh, String name, Point3& position, Color3& color) {
    // Replace any existing torus model. Models don't 
    // have to be added to the model table to use them 
    // with a VisibleEntity.
    const shared_ptr<Model>& planetModel = mesh.toArticulatedModel(name + "Model", color);
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
        Any any = Any::parse(anyStr);
        planet = scene()->createEntity(name, any);
        /*planet = scene()->createEntity("planet",
            PARSE_ANY(
                VisibleEntity {
                    model = "planetModel";
                };
            ));*/
    } else {
        // Change the model on the existing planet entity
        dynamic_pointer_cast<VisibleEntity>(planet)->setModel(planetModel);
    }

    //planet->setFrame(CFrame::fromXYZYPRDegrees(0.0f, 1.8f, 0.0f, 45.0f, 45.0f));
    planet->setFrame(CFrame::fromXYZYPRDegrees(position.x, position.y, position.z, 45.0f, 45.0f));
}

void App::makePentagon() {
    Array<Vector3> vertices1(Vector3(0, 0, -2), Vector3(0, 2, -2), Vector3(-2, 0, -2), Vector3(-1, -2, -2), Vector3(1, -2, -2), Vector3(2, 0, -2));
    Array<Vector3int32> indices1(Vector3int32(0, 1, 2), Vector3int32(0, 2, 3), Vector3int32(0, 3, 4), Vector3int32(0, 4, 5), Vector3int32(0, 5, 1));
    m_myMesh = Mesh::create(vertices1, indices1);

    m_myMesh->toObj("pentagon.obj");

    GuiPane* pentPane = debugPane->addPane("Info", GuiTheme::ORNATE_PANE_STYLE);

    // Example of how to add debugging controls
    pentPane->addLabel("Collapse Edges");

    pentPane->addButton("Collapse one edge", [this]() {
        m_myMesh->collapseEdges(1);
        m_myMesh->toObj("pentagon.obj");
        G3D::ArticulatedModel::clearCache();
        loadScene("Planet");

    });

    pentPane->pack();
}

void App::makePlanetGUI() {

    GuiPane* planetPane = debugPane->addPane("Planet");

    planetPane->setNewChildSize(240);
    planetPane->addNumberBox("Recursion Level", &m_recursionLevel, "",
        GuiTheme::LINEAR_SLIDER, 0, 8)->setUnitsSize(1);

    planetPane->addNumberBox("Frequency", &m_frequency, "",
        GuiTheme::LOG_SLIDER, 0.001f, 8.0f)->setUnitsSize(1);
    /*
    heightfieldPane->beginRow(); {
        heightfieldPane->addTextBox("Input Image", &m_heightfieldSource)->setWidth(210);
        heightfieldPane->addButton("...", [this]() {
            FileDialog::getFilename(m_heightfieldSource, "png", false);
        })->setWidth(30);
    } heightfieldPane->endRow();*/

    planetPane->addButton("Generate", [this]() {
        shared_ptr<G3D::Image> image;
        //Noise noise = G3D::Noise::common();
        try {
            Planet planet;
            NoiseGen noise;

            Array<Vector3> vertices = Array<Vector3>();
            Array<Vector3int32> faces = Array<Vector3int32>();
            float freq = 3.0f;
            shared_ptr<Image> image = Image::create(1024, 1024, ImageFormat::RGBA8()); //Image::fromFile("noise.jpg"); 
            //noise.generateNoisyImage(image, 0, freq);
            image->save("ocean.png");
            planet.writeSphere("ocean", 12.5f, 3, vertices, faces);
            Mesh mesh(vertices, faces);
            //mesh.toObj("ocean");

            vertices = Array<Vector3>();
            faces = Array<Vector3int32>();
            freq = 0.25f;
            image = Image::create(1024, 1024, ImageFormat::RGBA8()); //Image::fromFile("noise.jpg"); 
            noise.generateNoisyImage(image, 1, freq);
            image->save("land.png");
            planet.writeSphere("land", 12.0f, 8, vertices, faces);
            planet.applyNoiseLand(vertices, image);
            Mesh mesh2(vertices, faces);
            //mesh2.toObj("land");

            loadScene("Ground");
            addPlanetToScene(mesh, "ocean", Point3(0,0,0), Color3(0,0,1));
            addPlanetToScene(mesh2, "land", Point3(0,0,0), Color3(0,1,0));

            //loadScene("Planet");

            //Julia's Generate Code
            /*Planet planet;
            Array<Vector3> vertices1 = Array<Vector3>();
            Array<Vector3int32> faces1 = Array<Vector3int32>();
            planet.writeSphere("mountains", 10.1f, 5, vertices1, faces1);
            Mesh mesh(vertices1, faces1);
            mesh.bevelEdges(1);
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
                    image->set(x, y, lerp(Color3(0.2f, 0.3f, 0.7f), Color3(1.0f), noise.sampleFloat(x, y, x+y, 3)));
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

    GuiPane* infoPane = debugPane->addPane("Info", GuiTheme::ORNATE_PANE_STYLE);

    // Example of how to add debugging controls
    infoPane->addLabel("You can add GUI controls");
    infoPane->addLabel("in App::onInit().");
    infoPane->addButton("Exit", [this]() { m_endProgram = true; });
    infoPane->pack();
    
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
    
    Array<Vector3> verticeArray(Vector3(0,0,0), Vector3(1,0,0), Vector3(.5, 0, 1), Vector3(.5, 1, .5));
    Array<Vector3int32> triangles(Vector3int32(3,1,0), Vector3int32(1,2,0), Vector3int32(3,2,1), Vector3int32(3,0,2));
    
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

