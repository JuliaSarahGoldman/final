/** \file App.cpp */
#include "App.h"

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
    settings.window.caption             = argv[0];

    // Set enable to catch more OpenGL errors
    // settings.window.debugContext     = true;

    // Some common resolutions:
    // settings.window.width            =  854; settings.window.height       = 480;
    // settings.window.width            = 1024; settings.window.height       = 768;
    settings.window.width               = 1280; settings.window.height       = 720;
    //settings.window.width             = 1920; settings.window.height       = 1080;
    // settings.window.width            = OSWindow::primaryDisplayWindowSize().x; settings.window.height = OSWindow::primaryDisplayWindowSize().y;
    settings.window.fullScreen          = false;
    settings.window.resizable           = ! settings.window.fullScreen;
    settings.window.framed              = ! settings.window.fullScreen;

    // Set to true for a significant performance boost if your app can't render at 60fps, or if
    // you *want* to render faster than the display.
    settings.window.asynchronous        = false;

    settings.hdrFramebuffer.depthGuardBandThickness = Vector2int16(64, 64);
    settings.hdrFramebuffer.colorGuardBandThickness = Vector2int16(0, 0);
    settings.dataDir                    = FileSystem::currentDirectory();
    settings.screenshotDirectory        = "../journal/";

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
    
    showRenderingStats      = false;

    makeGUI();

    // For higher-quality screenshots:
    // developerWindow->videoRecordDialog->setScreenShotFormat("PNG");
    // developerWindow->videoRecordDialog->setCaptureGui(false);
    developerWindow->cameraControlWindow->moveTo(Point2(developerWindow->cameraControlWindow->rect().x0(), 0));
    loadScene(
        //"G3D Sponza"
        "G3D Cornell Box" // Load something simple
        //developerWindow->sceneEditorWindow->selectedSceneName()  // Load the first scene encountered 
        );
}

void App::writeSphere(String filename, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3>>& faces) {
    makeIcohedron(5.0f, vertices, faces);

    int numVert = vertices->size();
    int numFace = faces->size();

    String sphere = "OFF\n";
    sphere += (String) std::to_string(numVert) + " " + (String) std::to_string(numFace) + " 0\n\n";

    for(int i = 0; i < numVert; ++i) {
        Vector3 vertex = vertices->operator[](i);
        sphere += (String) std::to_string(vertex[0]) + " " + (String) std::to_string(vertex[1]) + " " + (String) std::to_string(vertex[2]) + "\n";
    }

    for(int i = 0; i < numFace; ++i) {
        Vector3 face = faces->operator[](i);
        sphere += "3 " + (String) std::to_string(face[0]) + " " + (String) std::to_string(face[1]) + " " + (String) std::to_string(face[2]) + "\n";
    }

    TextOutput output ("model/" + filename + ".off");
    output.writeSymbol(sphere);
    output.commit(true);
    G3D::ArticulatedModel::clearCache();
}

void App::makeIcohedron(float radius, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3>>& faces) {
    float t = (1.0 + sqrt(5.0)) / 2.0;
    float s = sqrt(1 + square(t));

    vertices->append(radius * Vector3(t, 1, 0) / s);
    vertices->append(radius * Vector3(-t, 1, 0) / s);
    vertices->append(radius * Vector3(t, -1, 0) / s);
    vertices->append(radius * Vector3(-t, -1, 0) / s);
    vertices->append(radius * Vector3(1, 0, t) / s);
    vertices->append(radius * Vector3(1, 0, -t) / s);
    vertices->append(radius * Vector3(-1, 0, t) / s);
    vertices->append(radius * Vector3(-1, 0, -t) / s);
    vertices->append(radius * Vector3(0, t, 1) / s);
    vertices->append(radius * Vector3(0, -t, 1) / s);
    vertices->append(radius * Vector3(0, t, -1) / s);
    vertices->append(radius * Vector3(0, -t, -1) / s);

    faces->append(Vector3(0, 8, 4));
    faces->append(Vector3(1, 10, 7));
    faces->append(Vector3(2, 9, 11));
    faces->append(Vector3(7, 3, 1));
    faces->append(Vector3(0, 5, 10));
    faces->append(Vector3(3, 9, 6));
    faces->append(Vector3(3, 11, 9));
    faces->append(Vector3(8, 6, 4));
    faces->append(Vector3(2, 4, 9));
    faces->append(Vector3(3, 7, 11));
    faces->append(Vector3(4, 2, 0));
    faces->append(Vector3(9, 4, 6));
    faces->append(Vector3(2, 11, 5));
    faces->append(Vector3(0, 10, 8));
    faces->append(Vector3(5, 0, 2));
    faces->append(Vector3(10, 5, 7));
    faces->append(Vector3(1, 6, 8));
    faces->append(Vector3(1, 8, 10));
    faces->append(Vector3(6, 1, 3));
    faces->append(Vector3(11, 7, 5));
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
    
    heightfieldPane->addButton("Generate", [this](){
        shared_ptr<G3D::Image> image;
        try {
            drawMessage("Generating Heightfield.");
            image = Image::fromFile(m_heightfieldSource);

            TextOutput output ("model/heightfield.off");
            output.writeSymbol("OFF\n");
            output.printf("%d %d 0\n", image->width() * image->height(), (image->width() - 1) * (image->height() - 1));

            for(int x = 0; x < image->width(); ++x){
                for(int z = 0; z < image->height(); ++z){
                    Color3 color;
                    image->get(Point2int32(x, z), color);
                    float y = color.average();
                    output.printf("%f %f %f\n", ((float)x)*m_heightfieldXZScale, y*m_heightfieldYScale, ((float)z)*m_heightfieldXZScale);
                }
            }

            for(int i = 1; i < image->height(); ++i){
                for(int j = 1; j < image->width(); ++j){
                    output.printf("4 %d %d %d %d\n", i + ((image->height())*j), i + ((image->height())*j) - 1, i + ((image->height())*(j-1)) - 1, i + ((image->height())*(j-1)));
                }
            }

            output.commit(true);
            G3D::ArticulatedModel::clearCache();
        } catch (...) {
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

    shared_ptr<Array<Vector3>> vertices = std::make_shared<Array<Vector3>>();
    shared_ptr<Array<Vector3>> faces = std::make_shared<Array<Vector3>>();
    writeSphere("test", vertices, faces);
    makeHeightfield();
    // More examples of debugging GUI controls:
    // debugPane->addCheckBox("Use explicit checking", &explicitCheck);
    // debugPane->addTextBox("Name", &myName);
    // debugPane->addNumberBox("height", &height, "m", GuiTheme::LINEAR_SLIDER, 1.0f, 2.5f);
    // button = debugPane->addButton("Run Simulator");
    // debugPane->addButton("Generate Heightfield", [this](){ generateHeightfield(); });
    // debugPane->addButton("Generate Heightfield", [this](){ makeHeightfield(imageName, scale, "model/heightfield.off"); });

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
